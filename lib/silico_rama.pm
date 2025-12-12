#!/usr/bin/perl -w

#  COPYRIGHT NOTICE
#  
#  Silico - A Perl Molecular Toolkit
#  Copyright (C) 2008-25 David K. Chalmers and Benjamin P. Roberts,
#  Department of Medicinal Chemistry, Monash University
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#<
#! silico_rama.pm
#? Silico routines to calculate Ramachandran plots and peptide
#  dihedral angles
#. $Revision: 1.12.2.4.2.7 $
#>

use strict;
package Silico;

sub mol_rama {

	#<
	#? Find backbone dihedrals and phi, psi and omega angles
	#; Requires: Molecule, parent_ramadata(optional)
	#; Returns: Ramadata record
	#>
	
	my $mol = $_[0];
	my $parent_ramadata = $_[1];
	
	if ($parent_ramadata) {
	
		$mol->{RAMADATA} = deep_copy($parent_ramadata);
		
	} else {
	
		molecule_check_and_fix_connectivity($mol);
		label_aa_backbone($mol);
		rama_find_backbone_dihedrals($mol);
	}
	
	mol_rama_angles($mol);
	
	return $mol->{RAMADATA};
}

sub rama_find_backbone_dihedrals {

	#<
	#? Extract backbone dihedral information from a molecule
	#; Requires: Molecule
	#; Returns: RAMADATA record
	#>

	my $mol = $_[0];
	
	#molecule_printout($mol);
	
	foreach my $atom (atoms($mol)) {
	
		# Next atom if we are not looking at an alpha carbon
		next if !$atom->{FG}{AA_CA};
		
		my $c;
		my $ca = $atom;
		my $n;
		my $nextn;
		my $nextca;
		my $prevc;
		my $record;
	
		foreach my $con (connected($ca, $mol)) {	
			$c = $con if $con->{FG}{AA_C};
			$n = $con if $con->{FG}{AA_N};
		}
		
		next if (!defined $c || !defined $n);
		
		foreach my $con (connected($n, $mol)) {
			$prevc = $con if $con->{FG}{AA_C};
		}
		
		if (!defined $prevc) {
		
			$atom->{FG}{C_TERMINAL} = 1;
			$record->{C_TERMINAL} = 1;
		} else {
			
			foreach my $con (connected($c, $mol)) {
				$nextn = $con if $con->{FG}{AA_N};
			}
		}
		
		if (!defined $nextn) {
			$atom->{FG}{N_TERMINAL} = 1;
			$record->{N_TERMINAL} = 1;
			
		} else {
		
			foreach my $con (connected($nextn, $mol)) {
				$nextca = $con if $con->{FG}{AA_CA};
			}
		}
		
		my $res = $ca->{SUBNAME}.$ca->{SUBID};
		$res =~ s/\s+//g;
		
		$record->{SUBNAME} = $ca->{SUBNAME};
		$record->{SUBID} = $ca->{SUBID};
		$record->{CHAIN} = $ca->{CHAIN} || '';
		$record->{LABEL} = $res;
		
		my $prevc_idx;
		my $nextn_idx;
		my $nextca_idx;
		
		$prevc_idx = $prevc->{NUM}-1 if defined $prevc;
		$nextn_idx = $nextn->{NUM}-1 if defined $nextn;
		$nextca_idx = $nextca->{NUM}-1 if defined $nextca;

		@{$record->{ATOMNUMS}} = ($prevc_idx, $n->{NUM}-1, $ca->{NUM}-1, $c->{NUM}-1, $nextn_idx, $nextca_idx);
 		push @{$mol->{RAMADATA}}, $record;
	}
	
	if (!defined $mol->{RAMADATA}[0]) {
		silico_msg ('w', "No Ramachandran dihedrals found found in molecule: $mol->{NAME}\n");
		return undef;
	}
	
	return $mol->{RAMADATA};
}

sub mol_rama_angles {

	#<
	#? Caculate Ramachandran angles for RAMADATA
	#. Omega is undefined if the residue is C- or N-terminal
	#. Also stores mol->{ENERGY} in ramadata->[0]{ENERGY};
	#; Requires: molecule
	#; Returns: Nothing
	#; Sets: $mol->{RAMADATA}[$rescount]{PHI} and $mol->{RAMADATA}[$rescount]{PSI} = $psi
	#>
	
	my $mol = $_[0];
	
	my $atoms = $mol->{ATOMS};
	
	$mol->{RAMADATA}[0]{ENERGY} = $mol->{ENERGY};
	
	foreach my $res (@{$mol->{RAMADATA}}) {

		my ($prevc, $n, $ca, $c, $nextn, $nextca);

		$prevc = $atoms->[$res->{ATOMNUMS}[0]]  if defined $res->{ATOMNUMS}[0];
		$n =     $atoms->[$res->{ATOMNUMS}[1]]  if defined $res->{ATOMNUMS}[1];
		$ca =    $atoms->[$res->{ATOMNUMS}[2]]  if defined $res->{ATOMNUMS}[2];
		$c =     $atoms->[$res->{ATOMNUMS}[3]]  if defined $res->{ATOMNUMS}[3];
		$nextn = $atoms->[$res->{ATOMNUMS}[4]]  if defined $res->{ATOMNUMS}[4];
		$nextca = $atoms->[$res->{ATOMNUMS}[5]] if defined $res->{ATOMNUMS}[5];

		my $phi;
		if (defined $prevc && defined $n && defined $ca && defined $c) {
			$phi = dihedral($prevc, $n, $ca, $c);
			$phi = rad_to_deg($phi);
		}
		
		my $psi;
		if (defined $n && defined $ca && defined $c && defined $nextn) {
			$psi = dihedral($n, $ca, $c, $nextn);
			$psi = rad_to_deg($psi);
		}
		
		my $omega;
		my $cis = 0;
		if (defined $nextca && defined $nextn) {
			$omega = dihedral($ca, $c, $nextn, $nextca);
			$omega = rad_to_deg($omega);
		
			if ($omega > -90 && $omega < 90) {
				$cis = 1;
			}
		}
		
		#print "omega ", ($omega || '_'), "\n";
		
		$res->{PHI} = $phi;
		$res->{PSI} = $psi;
		$res->{OMEGA} = $omega;
		$res->{CIS} = $cis;
	}
	
	return undef;
}

sub ramaplot_open {

	#<
	#? Open .rama files and write headers
	#; Requires:output filename
	#; Returns: file record 
	#>
	
	my $filebase = $_[0];
	my $options = $_[1] || '';

	my $cart = 1 if $options =~ /\bCART\b/i;
	my $noheader = 1  if $options =~ /\bNOHEADER\b/i;
	my $outfile = $filebase;
	
	$outfile .= "_cart" if $cart;
	$outfile .= "_rama.txt";
	
	#$outfile = overwrite_check($outfile) unless $main::force;
	
	my $FH = new FileHandle;
	my $fr;
	$fr->{FH} = $FH;
	$fr->{FILENAME} = $outfile;
	
	open($FH, ">$outfile") || file_write_error($outfile, 1);
	
	if ($cart) {
		print $FH "# Ramachandran Cartesian Output\n" if !$noheader;
		print $FH "# Generated on ".`date` if !$noheader;
		print $FH "#Residue\tsin(Phi)\tcos(Phi)\tsin(Phi)\tcos(Phi)\tcos(Omega)\n";
	} else {
		print $FH "# Ramachandran Plot Output\n" if !$noheader;
		print $FH "# Generated on ".`date` if !$noheader;
		print $FH "#Residue\tPhi\tPsi\tOmega\n";
	}
	print $FH "#\n";
	
	return $fr;
}

sub ramaplot_by_mol_open {

	#<
	#? Open .ramol files and write headers
	#; Requires: ramadata record, output filename base, options
	#; Returns: file record
	#; Options: 'CART' will output sin(phi), cos(phi), sin(psi), cos(psi) instead of phi, phsi. 
	#  NOHEADER - don't write comment at top of file - required for R
	#>
	
	my $ramadata = $_[0];
	my $filebase = $_[1];
	my $options = $_[2];
	
	my $en = 1 if get_flag('energy', 's');
	my $FH;
	my $fr;
	my $outfile = "$filebase\_rama.txt";
	my $res;
	
	my $cart = 1 if $options =~ /\bCART\b/i;
	
	if ($cart) {
		$outfile = $filebase."_ramol_cart.txt";
	} else {
		$outfile = $filebase."_ramol.txt";
	}
	
	$fr->{FH} = $FH = new FileHandle;
	$fr->{FILENAME} = $outfile;
	
	open($FH, ">$outfile") || file_write_error($outfile, 1);
	
	if ($options !~ /\bNOHEADER\b/) {
		print $FH "# Ramachandran Plot Output\n";
		print $FH "# Generated on ".`date`;
		print $FH "#\n";
	}
	
	print $FH "Energy\t" if $en;
	
	my $i = 0;
	foreach (@$ramadata) {

		$res = $_->{LABEL};
	
		print $FH  "\t" if $i > 0;

		if (!$cart) {
			print $FH "$res\_phi\t$res\_psi";
		} else {
			print $FH "$res\_sin_phi\t$res\_cos_phi\t$res\_sin_psi\t$res\_cos_psi";
		}
		++$i;
	};

	print $FH "\n";
	
	return $fr;
}


sub write_rama_structure {

	#<
	#? Write out a _rama.txt file of a single structure suitable for the subroutine 'ramaplot'
	#; Requires: molecule, filehandle, optional residue number (indexed from 0).
	#; Returns: nothing
	#>
	
	my $ramadata = $_[0];
	my $fr = $_[1];
	my $residue = $_[2]; # Residue number (optional);
	my $options = $_[3] || '';

	my $cart = 1 if $options =~ /\bCART\b/i;
	my $en = 1 if get_flag('energy', 's');
	my $FH = $fr->{FH};
	
	my $list;

	if ($residue) {
		@$list = ($ramadata->[$residue]);
	} else {
		@$list = @{$ramadata}
	}
	
	foreach my $rec (@$list) {

		my $res = $rec->{LABEL};
		
		my $phi = $rec->{PHI};
		my $psi = $rec->{PSI};

		if (!defined $phi || !defined $psi) {
			
			print "Residue $res has missing phi or psi values. Skipping\n";;
			next;
		}

		if ($cart ) {

			my $a = deg_to_rad($phi);
			my $b = deg_to_rad($psi);
			
			printf $FH "$res\t%f\t%f\t%f\t%f", sin($a), cos($a), sin($b), cos($b);
			
			if (defined $rec->{OMEGA}) {
				my $c =  deg_to_rad($rec->{OMEGA});
				printf $FH "\t%f\t%f", sin($c), cos($c);
				
			} else {
				printf $FH "\t\t"
			}
			printf $FH "\n";

		} else {
	
			$phi = sprintf "%8.3f", $rec->{PHI};
			$psi = sprintf "%8.3f", $rec->{PSI};
			
			print $FH "$res\t$phi\t$psi";
			
			if (defined $rec->{OMEGA}) {
				my $omega = sprintf "%8.3f", $rec->{OMEGA};
				print $FH "\t$omega";
				
			} else {
				print $FH "\t";
			}
			print $FH "\n";
		}
	}
}



sub write_rama_structure_by_mol {

	#<
	#? Write out a .ramol file of a single structure with all phi and psi angles written on a single line
	#; Requires: molecule, file record, optional residue number (indexed from 0), options
	#; Options: 'CART' will output sin(phi), cos(phi), sin(psi), cos(psi) instead of phi, phsi
	#; Returns: nothing
	#>
	
	my $ramadata = $_[0];
	my $fr = $_[1];
	my $options = $_[2];
	
	my $cart = 1 if $options =~ /\bCART\b/i;
	my $en = 1 if get_flag('energy', 's');
	my $FH = $fr->{FH};
	my @list;
	
	my $energy = $ramadata->[0]{ENERGY} || 99999;
	printf $FH "%8.3f\t", $energy if $en;
	
	my $i = 0;
	foreach my $rec (@$ramadata) {

		my $res = $rec->{LABEL};
		my $phi = sprintf "%8.3f", $rec->{PHI};
		my $psi = sprintf "%8.3f", $rec->{PSI};
		
		my $omega;
		if (defined $rec->{OMEGA}) {
			$omega = sprintf "%8.3f", $rec->{OMEGA};
		} else {
			$omega ='';
		}
		
		print $FH "\t" if $i > 0;
		
		if (!$cart) {
		
			print $FH "$phi\t$psi\t$omega";
			
		} else {
		
			my $a = deg_to_rad($phi);
			my $b = deg_to_rad($psi);
			
			printf $FH "%f\t%f\t%f\t%f", sin($a), cos($a), sin($b), cos($b);
			
			if ($omega) {
				my $c = my $b = deg_to_rad($omega);
				printf $FH "\t%f\t%f", sin($c), cos($c)
			
			} else {
				printf $FH "\t\t";
			}
		}
		++$i;
	}
	
	print $FH "\n";
}



sub ramaplot {

	#<
	#? Writes a grace parameter file for the
	#  Ramachandran plot, and calls grace
	#; Requires: Ramachandran data file, format, hardcopy flag, residue label flag
	#; Returns: Nothing
	#>
	
	my $filebase = $_[0];    # File containing data to be plotted
	my $oformat = $_[1]; # Grace output format - see check_grace_output_format()
	my $hardcopy= $_[2]; # Send output directly to printer
	my $label = $_[3];   # Optionally label points with residue name

	$label = 0 if !defined $label;
	
	#
	# Check that xmgrace exists and that the output format is sane
	# 
	
	my $grace_status = system("which xmgrace");
	
	if ($grace_status > 0) {
		silico_msg('n', "Grace is not installed on this computer.\n",
				"Only data files will be generated, and no graphs.\n");
	}
	
	my ($format, $ext) = check_grace_output_format ($oformat);
	
	if (!defined $format) {

		silico_msg('d', "Unsupported Grace output format '$oformat'!\n",
			"Please choose another output format using the -o flag.\n");
	}
	
	if ($hardcopy && ($oformat ne 'ps')) {
		silico_msg('w', "Print option and output format $oformat are not compatible!\n",
				"Switching to PostScript format.\n");
	
		$oformat = 'ps';
	}
	
	$| = 0;
	print "Creating graph...\n";
	
	my $diff = length ($filebase) - length ("Ramachandran plot of molecule");
	my $offset = $diff * -0.3;
	
	my $time = `date`;
	chomp $time;
	
	open(PARAM, ">$filebase.gpf") || file_write_error("$filebase.gpf", 1);
	
	print PARAM "# Grace project file
#
version 50118
hardcopy device \"$format\"
page size 595,843
page scroll 5%
page inout 5%
link page off
map font 0 to \"Times-Roman\", \"Times-Roman\"
map font 1 to \"Times-Italic\", \"Times-Italic\"
map font 2 to \"Times-Bold\", \"Times-Bold\"
map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\"
map font 4 to \"Helvetica\", \"Helvetica\"
map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\"
map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\"
map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\"
map font 8 to \"Courier\", \"Courier\"
map font 9 to \"Courier-Oblique\", \"Courier-Oblique\"
map font 10 to \"Courier-Bold\", \"Courier-Bold\"
map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\"
map font 12 to \"Symbol\", \"Symbol\"
map font 13 to \"ZapfDingbats\", \"ZapfDingbats\"
map color 0 to (255, 255, 255), \"white\"
map color 1 to (0, 0, 0), \"black\"
map color 2 to (255, 0, 0), \"red\"
map color 3 to (0, 255, 0), \"green\"
map color 4 to (0, 0, 255), \"blue\"
map color 5 to (255, 255, 0), \"yellow\"
map color 6 to (188, 143, 143), \"brown\"
map color 7 to (220, 220, 220), \"grey\"
map color 8 to (148, 0, 211), \"violet\"
map color 9 to (0, 255, 255), \"cyan\"
map color 10 to (255, 0, 255), \"magenta\"
map color 11 to (255, 165, 0), \"orange\"
map color 12 to (114, 33, 188), \"indigo\"
map color 13 to (103, 7, 72), \"maroon\"
map color 14 to (64, 224, 208), \"turquoise\"
map color 15 to (0, 139, 0), \"green4\"
reference date 0
date wrap off
date wrap year 1950
default linewidth 1.0
default linestyle 1
default color 1
default pattern 1
default font 8
default char size 1.000000
default symbol size 1.000000
default sformat \"%.8g\"
background color 0
page background fill on
timestamp off
timestamp 0.03, 0.03
timestamp color 1
timestamp rot 0
timestamp font 8
timestamp char size 1.000000
timestamp def \"$time\"
r0 off
link r0 to g0
r0 type above
r0 linestyle 1
r0 linewidth 1.0
r0 color 1
r0 line 0, 0, 0, 0
r1 off
link r1 to g0
r1 type above
r1 linestyle 1
r1 linewidth 1.0
r1 color 1
r1 line 0, 0, 0, 0
r2 off
link r2 to g0
r2 type above
r2 linestyle 1
r2 linewidth 1.0
r2 color 1
r2 line 0, 0, 0, 0
r3 off
link r3 to g0
r3 type above
r3 linestyle 1
r3 linewidth 1.0
r3 color 1
r3 line 0, 0, 0, 0
r4 off
link r4 to g0
r4 type above
r4 linestyle 1
r4 linewidth 1.0
r4 color 1
r4 line 0, 0, 0, 0
g0 on
g0 hidden false
g0 type XY
g0 stacked false
g0 bar hgap 0.000000
g0 fixedpoint off
g0 fixedpoint type 0
g0 fixedpoint xy 0.000000, 0.000000
g0 fixedpoint format general general
g0 fixedpoint prec 6, 6
with g0
    world xmin -180
    world xmax 180
    world ymin -180
    world ymax 180
    stack world 0, 0, 0, 0
    znorm 1
    view xmin 0.150000
    view xmax 0.850000
    view ymin 0.150000
    view ymax 0.850000
    title \"Ramachandran Plot of molecule\\n\\h{$offset}$filebase\"
    title font 8
    title size 1.000000
    title color 1
    subtitle \"\"
    subtitle font 8
    subtitle size 1.000000
    subtitle color 1
    xaxes scale Normal
    yaxes scale Normal
    xaxes invert off
    yaxes invert off
    xaxis  on
    xaxis  type zero false
    xaxis  offset 0.000000 , 0.000000
    xaxis  bar on
    xaxis  bar color 1
    xaxis  bar linestyle 1
    xaxis  bar linewidth 1.0
    xaxis  label \"\\f{Symbol}f\"
    xaxis  label layout para
    xaxis  label place auto
    xaxis  label char size 1.000000
    xaxis  label font 8
    xaxis  label color 1
    xaxis  label place normal
    xaxis  tick on
    xaxis  tick major 90
    xaxis  tick minor ticks 1
    xaxis  tick default 6
    xaxis  tick place rounded true
    xaxis  tick in
    xaxis  tick major size 1.000000
    xaxis  tick major color 1
    xaxis  tick major linewidth 1.0
    xaxis  tick major linestyle 1
    xaxis  tick major grid off
    xaxis  tick minor color 1
    xaxis  tick minor linewidth 1.0
    xaxis  tick minor linestyle 1
    xaxis  tick minor grid off
    xaxis  tick minor size 0.500000
    xaxis  ticklabel on
    xaxis  ticklabel format general
    xaxis  ticklabel prec 5
    xaxis  ticklabel formula \"\"
    xaxis  ticklabel append \"\"
    xaxis  ticklabel prepend \"\"
    xaxis  ticklabel angle 0
    xaxis  ticklabel skip 0
    xaxis  ticklabel stagger 0
    xaxis  ticklabel place normal
    xaxis  ticklabel offset auto
    xaxis  ticklabel offset 0.000000 , 0.010000
    xaxis  ticklabel start type auto
    xaxis  ticklabel start 0.000000
    xaxis  ticklabel stop type auto
    xaxis  ticklabel stop 0.000000
    xaxis  ticklabel char size 1.000000
    xaxis  ticklabel font 8
    xaxis  ticklabel color 1
    xaxis  tick place both
    xaxis  tick spec type none
    yaxis  on
    yaxis  type zero false
    yaxis  offset 0.000000 , 0.000000
    yaxis  bar on
    yaxis  bar color 1
    yaxis  bar linestyle 1
    yaxis  bar linewidth 1.0
    yaxis  label \"\\f{Symbol}y\"
    yaxis  label layout para
    yaxis  label place auto
    yaxis  label char size 1.000000
    yaxis  label font 8
    yaxis  label color 1
    yaxis  label place normal
    yaxis  tick on
    yaxis  tick major 90
    yaxis  tick minor ticks 1
    yaxis  tick default 6
    yaxis  tick place rounded true
    yaxis  tick in
    yaxis  tick major size 1.000000
    yaxis  tick major color 1
    yaxis  tick major linewidth 1.0
    yaxis  tick major linestyle 1
    yaxis  tick major grid off
    yaxis  tick minor color 1
    yaxis  tick minor linewidth 1.0
    yaxis  tick minor linestyle 1
    yaxis  tick minor grid off
    yaxis  tick minor size 0.500000
    yaxis  ticklabel on
    yaxis  ticklabel format general
    yaxis  ticklabel prec 5
    yaxis  ticklabel formula \"\"
    yaxis  ticklabel append \"\"
    yaxis  ticklabel prepend \"\"
    yaxis  ticklabel angle 0
    yaxis  ticklabel skip 0
    yaxis  ticklabel stagger 0
    yaxis  ticklabel place normal
    yaxis  ticklabel offset auto
    yaxis  ticklabel offset 0.000000 , 0.010000
    yaxis  ticklabel start type auto
    yaxis  ticklabel start 0.000000
    yaxis  ticklabel stop type auto
    yaxis  ticklabel stop 0.000000
    yaxis  ticklabel char size 1.000000
    yaxis  ticklabel font 8
    yaxis  ticklabel color 1
    yaxis  tick place both
    yaxis  tick spec type none
    altxaxis  off
    altyaxis  off
    legend off
    legend loctype view
    legend 0.85, 0.8
    legend box color 1
    legend box pattern 1
    legend box linewidth 1.0
    legend box linestyle 1
    legend box fill color 0
    legend box fill pattern 1
    legend font 8
    legend char size 1.000000
    legend color 1
    legend length 4
    legend vgap 1
    legend hgap 1
    legend invert false
    frame type 0
    frame linestyle 1
    frame linewidth 1.0
    frame color 1
    frame pattern 1
    frame background color 0
    frame background pattern 0
    s0 hidden false
    s0 type xy
    s0 symbol 1
    s0 symbol size 0.240000
    s0 symbol color 2
    s0 symbol pattern 1
    s0 symbol fill color 2
    s0 symbol fill pattern 1
    s0 symbol linewidth 1.0
    s0 symbol linestyle 1
    s0 symbol char 67
    s0 symbol char font 8
    s0 symbol skip 0
    s0 line type 0
    s0 line linestyle 1
    s0 line linewidth 1.0
    s0 line color 1
    s0 line pattern 1
    s0 baseline type 0
    s0 baseline off
    s0 dropline off
    s0 fill type 0
    s0 fill rule 0
    s0 fill color 1
    s0 fill pattern 1
    s0 avalue on
    s0 avalue type 4
    s0 avalue char size 0.400000
    s0 avalue font 8
    s0 avalue color 1
    s0 avalue rot 0
    s0 avalue format general
    s0 avalue prec 3
    s0 avalue prepend \"\"
    s0 avalue append \"\"
    s0 avalue offset 0.020000 , -0.010000
    s0 errorbar on
    s0 errorbar place both
    s0 errorbar color 2
    s0 errorbar pattern 1
    s0 errorbar size 1.000000
    s0 errorbar linewidth 1.0
    s0 errorbar linestyle 1
    s0 errorbar riser linewidth 1.0
    s0 errorbar riser linestyle 1
    s0 errorbar riser clip off
    s0 errorbar riser clip length 0.100000
    s0 comment \"$filebase\"
    s0 legend  \"\"
";

	my $command = "xmgrace ";
	$command .= "-block $filebase\_rama.txt ";
	if ($label) {
		$command .= "-bxy 2:3:'{1}' ";
	} else {
		$command .= "-bxy 2:3 ";
	}
	$command .= "-param $filebase.gpf ";
	$command .= "-printfile $filebase\_rama.$ext " if !$hardcopy;
	$command .= "-hardcopy";  #  No interactive session, just print and quit
	
	print "$command\n";
	
	silico_msg('g', "\nCommand: $command\n");
	
	system ("$command");
	unlink "$filebase.gpf" unless $Silico::debug;
	
	print "done\n";
	$| = 1;
}

return 1;
