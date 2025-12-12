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
#! silico_mopac.pm
#? Silico Mopac and Divcon routines.
#. $Revision: 1.56.2.1.2.11 $
#>

use strict;

package Silico;

##################################################################
#
#	Mopac read routines
#
##################################################################

sub read_mopac_cart {

	#<
	#. Read in a Cartesian mopac .dat, .arc or .out file
	#; Requires: filename, (optional) unused, unused, options string
	#; Returns: ensemble or zero if read failed.
	#>

	my $infile  = $_[0];
	my $options = $_[3] || '';

	my $arc_flag  = 0;
	my $cart_flag = 0;

	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
	if (!defined $fr) {
		silico_msg('e', "Can not open file $infile for reading\n");
		return undef;
	}
	my $FH = $fr->{FH};

	silico_msg('c', "Reading Mopac Cartesian file: $infile\n") if $options !~ /\bQUIET\b/;

	# Decide if the file contains Cartesian coordinates
	while (<$FH>) {
		next if !/CARTESIAN COORDINATES/;
		$cart_flag = 1;
		last;
	}

	# Close and reopen file
	seek $FH, 0, 0;

	if (!$cart_flag) {
		while (<$FH>) {

			# Last line of header in arc file
			next if !/FINAL GEOMETRY OBTAINED/;
			$arc_flag = 1;
			last;
		}
	}

	# Close and reopen file
	seek $FH, 0, 0;
	
	my $mol;
	if ($arc_flag) {
		$mol = read_mopac_arc_single($fr, $options);
		#print "Reading Mopac .arc file\n";
	} elsif ($cart_flag) {
		$mol = read_mopac_out_single($fr, $options);
		#print "Reading Mopac .out file\n";
	} else {
		$mol = read_mopac_dat_single($fr, $options);
		#print "Reading Mopac .dat file\n";
	}
	
	my $mols;
	$mols->[0] = $mol;
	return $mols;
}

sub read_mopac_arc_single {

	#
	# Read a single structure from a mopac arc file and return molecule
	#

	my $fr  = $_[0];
	my $options = $_[1];

	my $mols;
	my $mol;
	my $head;
	my $comments;
	my $keywords;
	my %element_count;

	my $FH = $fr->{FH};

	while (<$FH>) {

		if (/HEAT OF FORMATION/) {

			my @f = split " ", $_;
			$mol->{ENERGY}                 = $f[7];
			$mol->{ENERGY_UNITS}           = "kJ/mol";
			$mol->{SDF_DATA}{ENERGY}       = $f[7];
			$mol->{SDF_DATA}{ENERGY_UNITS} = "kJ/mol";
		}

		if (/GRADIENT NORM/) {

			my @f = split " ", $_;
			$mol->{SDF_DATA}{GRADIENT_NORM} = $f[3];
		}

		# Last line of header in arc file
		next if !(/FINAL GEOMETRY OBTAINED/);

		# Read in keywords and next two comment lines
		$keywords = <$FH>;
		$comments = <$FH>;
		$comments .= <$FH>;

		if (!(defined $keywords && defined $comments)) {
			silico_msg('e', "Error reading file $fr->{FILENAME}!\n");
			return undef;
		}
		
		chomp $keywords;
		$mol->{SDF_DATA}{KEYWORDS} = $keywords;
		
		# Read in z-matrix type Cartesian coordinates
		# Read in body of input file
		my $i = 0;
		while (<$FH>) {

			last if !defined $_;

			chomp;

			#Stop if we don't match a non-white space character
			last if !/\S/;

			my @l = split;

			my $atom;

			$atom->{NUM}     = $i + 1;            #Atom number
			$atom->{ELEMENT} = $l[0];             #Atom element
			$atom->{SUBNAME} = "MOL";             #Residue name
			$atom->{SUBID}   = 1;                 #Residue number
			$atom->{X}       = $l[1];             #x
			$atom->{Y}       = $l[3];             #y
			$atom->{Z}       = $l[5];             #z
			$atom->{OPT_X}   = $l[2];             #flag to optimise x
			$atom->{OPT_Y}   = $l[4];             #flag to optimise y
			$atom->{OPT_Z}   = $l[6];             #flag to optimise z
			$atom->{CHARGE}  = $l[7] if $l[7];    #Charge
			
			push @{$mol->{ATOMS}}, $atom;
		}

		mopac_finalise($mol, $fr, $options);
	}

	return $mol;
}

sub read_mopac_out_single {

	#
	# Read a single structure from a mopac out file and return ensemble
	#

	my $fr  = $_[0];
	my $options = $_[1];

	my $mol;
	my $head;
	my $comments;
	my $keywords;
	my %element_count;

	my $FH = $fr->{FH};

	# Find start
	while (<$FH>) {

		next if !/CURRENT VALUE OF GEOMETRY/;

		# Skip 4 lines
		foreach (1..4) {
			<$FH>;
		}

		# Read in body of input file
		my $i = 0;
		while (<$FH>) {

			last if !defined $_;

			chomp;
			my @l = split;
			
			++$i;

			#Stop if there were no fields
			last if $#l < 0;

			my $atom;

			$atom->{NUM}     = $i;        #Atom number
			$atom->{ELEMENT} = $l[0];     #Atom element
			$atom->{SUBNAME} = "MOL";     #Residue name
			$atom->{SUBID}   = 1;         #Residue number
			$atom->{X}       = $l[1];     #x
			$atom->{Y}       = $l[3];     #y
			$atom->{Z}       = $l[5];     #z

			push @{ $mol->{ATOMS} }, $atom;
			
			atom_guess_element($atom);
		}

		last;
	}

	return $mol;
}

sub read_mopac_dat_single {

	#<
	#? Read a mopac .dat file
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $fr  = $_[0];
	my $options = $_[1];

	my $mol;
	my $head;
	my $comments;
	my $keywords;
	my %element_count;

	my $FH = $fr->{FH};

	# Read in z-matrix type Cartesian coordinates
	$keywords = <$FH>;
	$comments = <$FH>;
	$comments .= <$FH>;

	# Read in keywords and next two comment lines
	if (!(defined $keywords && defined $comments)) {
		silico_msg('e', "Error reading file $fr->{FILENAME}!\n");
		return undef;
	}

	chomp $keywords;

	# Read in body of input file
	my $i = 0;
	while (<$FH>) {

		last if !defined $_;

		chomp;

		#Stop if we don't match a non-white space character
		last if !/\S/;

		my @l = split;

		silico_msg('g', "@l\n");

		my $atom;

		$atom->{NUM}     = $i + 1;            #Atom number
		$atom->{ELEMENT} = $l[0];             #Atom element
		$atom->{SUBNAME} = "MOL";             #Residue name
		$atom->{SUBID}   = 1;                 #Residue number
		$atom->{X}       = $l[1];             #x
		$atom->{Y}       = $l[3];             #y
		$atom->{Z}       = $l[5];             #z
		$atom->{OPT_X}   = $l[2];             #flag to optimise x
		$atom->{OPT_Y}   = $l[4];             #flag to optimise y
		$atom->{OPT_Z}   = $l[6];             #flag to optimise z
		$atom->{CHARGE}  = $l[7] if $l[7];    #Charge

		# Assign atom record to molecule
		$mol->{ATOMS}[$i] = $atom;

		++$i;
	}

	mopac_finlalise($mol, $fr, $options);

	return $mol;
}

sub mopac_finalise {

	#<
	#? Fix up atomic elements from mopac file read and perform checks
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $fr  = $_[1];
	my $options = $_[2] || '';
	
	#fr_printout($fr);
	
	if ($#{$mol->{ATOMS}} < 0) {
		molecule_printout($mol);
		silico_msg('e', "The file $fr->{FILENAME} does not contain any atom records! Perhaps MOPAC failed to run.\n", "Skipping.\n");
		return undef;
	}

	# Check to see if we read in a z-matrix
	if (        defined $mol->{ATOMS}[0]
		&&  defined $mol->{ATOMS}[1]
		&&  defined $mol->{ATOMS}[2]
		&& $mol->{ATOMS}[0]{X} == 0
		&& $mol->{ATOMS}[0]{Y} == 0
		&& $mol->{ATOMS}[0]{Z} == 0
		&& $mol->{ATOMS}[1]{Y} == 0
		&& $mol->{ATOMS}[1]{Z} == 0
		&& $mol->{ATOMS}[2]{Z} == 0) {

		silico_msg('w', "You appear to be reading in a mopac z-matrix file as Cartesian coordinates!\n");
		molecule_printout($mol);
	}

	# Check and fix atomic elements and set atom names
	my $i = 0;
	my %element_count;
	foreach my $atom (atoms($mol)) {

		# Checks and calculated fields
		$atom->{ELEMENT} =~ /^\d/ && silico_msg('d', "Element starts with a digit! Can't read this yet.\n");
		$atom->{NAME} = $atom->{ELEMENT} . (++$element_count{ $atom->{ELEMENT} });    #Atom name

		# Set atom element number (a hash would be better)
		my $aname = ucfirst lc $atom->{ELEMENT};
		for (my $j = 1 ; $j <= $#Atomic_Elements ; ++$j) {
			if ($aname =~ /^$Atomic_Elements[$j]/) {
				$atom->{ELEMENT_NUM} = $j;
				last;
			}
		}

		if (!defined $atom->{ELEMENT_NUM}) {
			$atom->{ELEMENT}     = "Du";
			$atom->{ELEMENT_NUM} = 0;
		}
		++$i;
	}

	$mol->{NUMATOMS} = $i;
	$mol->{NUMBONDS} = 0;
	mol_set_source_data($mol, $fr, 'mopac');
	molecule_generate_bonds($mol) if ($options !~ /\bNOBOND\b/);
	$mol->{NAME} = $mol->{SOURCE_FILE_NAME};

}

##################################################################
#
#	Mopac write routines
#
##################################################################

sub write_mopac_cart {

	#<
	#? Print out a mopac Cartesian format file
	#. Currently prints out only the first molecule
	#  of the ensemble
	#; Requires: Ensemble or molecule, filename.
	#; Returns: undef if file open fails otherwise returns 1.
	#>

	my $molecules = ens($_[0]);
	my $outfile   = $_[1];
	my $keywords  = $_[2] || $molecules->[0]{MOPAC_KEYWORDS} || "AM1 NOINTER LET XYZ";

	my $atom;
	my $i;
	my $mol;
	my $name;
	my $success;

	$outfile = get_ofilebase($molecules) . ".dat" if !$outfile;

	$mol  = $molecules->[0];
	$name = $mol->{NAME} || "Mol";
	chomp($name);

	my $fr;
	$fr = open_molfile($fr, $outfile, undef, 'write');
	if (!defined $fr) {
		silico_msg('e', "Can not open file $outfile for writing\n");
		return undef;
	}
	++$fr->{MOL_WRITE_COUNT};
	my $FH = $fr->{FH};

	silico_msg('c', "Writing Mopac Cartesian file: $outfile\n");

	print $FH "$keywords\n";
	print $FH "$name\n";
	print $FH "\n";

	for ($i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

		$atom = $mol->{ATOMS}[$i];

		$atom->{OPT_X} || ($atom->{OPT_X} = 1);
		$atom->{OPT_Y} || ($atom->{OPT_Y} = 1);
		$atom->{OPT_Z} || ($atom->{OPT_Z} = 1);

		printf $FH " %-3s %10.6f %2d %10.6f %2d %10.6f %2d",
		  $atom->{ELEMENT}, $atom->{X}, $atom->{OPT_X},
		  $atom->{Y},       $atom->{OPT_Y},
		  $atom->{Z},       $atom->{OPT_Z};

		printf $FH " %10.6f", $atom->{CHARGE} if defined $atom->{CHARGE};
		print $FH "\n";
	}

	close($FH);

	return 1;
}

##################################################################
#
#	Run Mopac
#
##################################################################

sub setup_mopac_flags {

	#<
	#? Setup mopac flags
	#; Requires: Nothing
	#; Returns: Mopac keywords
	#>

	my $mopac_exe ||= `which mopac`;
	chomp $mopac_exe;

	my $ham    = make_flag('ham', 'hamiltonian',   "Method [AM1, PM7, etc]", 1, "PM7");
	my $charge = make_flag('c',   'formal-charge', "Formal charge",          1);
	my $radical_ok = make_flag('r',      'radicals-ok',           "Run calculations on radical species");
	my $geo_ok     = make_flag('geo',    'mopac-geo-ok',          "Use Mopac GEO-OK keyword");
	my $onescf     = make_flag('1scf',   'mopac-1scf',            "Calculate charge only (1SCF)");
	my $ef         = make_flag('ef',     'eigenvector-following', "Geometry optimisation using eigenvector following");
	my $esp        = make_flag('esp',    'esp-charges',           "Use ESP charges");
	my $mozyme     = make_flag('mozyme', 'mopac-mozyme',          "Use mozyme keyword");
	my $gnorm      = make_flag('gnorm',  'mopac-gnorm',           "Value for GNORM keyword", 1);
	my $cycles     = make_flag('cycles', 'mopac-cycles',          "Use CYCLES keyword", 1, 1000);
	my $key        = make_flag('k',      'mopac-keywords',        "Additional keywords", 1);
	
	
	my $noclean = make_flag('noclean',  'noclean',       "Do not delete output files");
	
	$mopac_exe = make_flag('exe', 'mopac-exe', "Location of MOPAC executable", 1, $mopac_exe);

	if (!-x $mopac_exe) {
		silico_msg('d', "$mopac_exe is not a valid executable file!\n");
	}

	chomp $key if $key;
	my $keywords = $ham . " MMOK XYZ";

	$keywords .= " ESP"            if $esp;
	$keywords .= " 1SCF"           if $onescf;
	$keywords .= " EF"             if !$onescf && $ef;
	$keywords .= " MOZYME"         if $mozyme;
	$keywords .= " GNORM=$gnorm"   if $gnorm;
	$keywords .= " CHARGE=$charge" if $charge;
	$keywords .= " GEO-OK"         if $geo_ok;
	$keywords .= " CYCLES=$cycles" if $cycles;
	$keywords .= " $key"           if $key;

	$Silico::mopac_keywords = $keywords;

	return $keywords;
}

sub run_mopac {

	#<
	#? Run mopac minimisation or mopac ESP charge calculation on molecule
	#. MOPAC keywords should be in $mol->{KEYWORDS}
	#; Requires: Molecule,  filebase, location of mopac_exe, flag to not clean up mopac files after running
	#; Returns: molecule
	#>

	my $mol       = $_[0];
	my $mopac_exe = $_[1] || get_lflag('mopac-exe');
	my $noclean   = $_[2] || get_lflag('noclean');

	#$noclean = 1;
	
	my @list = qw( ENERGY ENERGY_UNITS GRADIENT_NORM KEYWORDS);
	

	$Silico::temp_dir = "/tmp";
	$Silico::temp_dir = '.' if $noclean;

	my $mfilebase = get_tempname();
	my $minfile   = $mfilebase . "_mop.dat";
	my $marcfile  = $mfilebase . "_mop.arc";
	my $moutfile  = $mfilebase . "_mop.out";

	my $run = "rm -f $marcfile; rm -f $moutfile";
	
	
	if (system($run)) {
		silico_msg('d', "Execution of command \"$run\" failed!\n");
	}

	my $ens;
	$ens->[0] = $mol;
	if (!write_mopac_cart($ens, $minfile)) {
		silico_msg('d', "Could not write $minfile!\n");
	}

	$run = "$mopac_exe $minfile";
	$|   = 1;
	#print "$run\n" if (!get_flag('quiet', 'l'));
	print "$run\n";

	if (system($run)) {
		silico_msg('e', "Mopac failed to run!\n", "Skipping.\n\n");
		return;
	}

	# If we use the esp flag only read in the ESP charges
	# which are written out at the end of the file
	# and apply to the original molecule
	# otherwise read in whole molecule with read_mopac_cart
	if (get_flag('mopac-esp', 'l')) {
	
		if (!-e $marcfile) {

			silico_msg('e', "Mopac failed to produce an archive file: $marcfile\nSkipping.\n\n");
			return;
		}

		silico_msg('c', "Reading ESP charges only.  Coordinates not read in\n");

		open(INFILE, $marcfile) || file_read_error($marcfile, 1);

		while (<INFILE>) {
			if (/ELECTROSTATIC POTENTIAL CHARGES/) {
				<INFILE>;    # Blank line
				<INFILE>;    # HEADER

				my $atom;

				foreach $atom (atoms($mol)) {

					my @f;
					my $line;

					$line = <INFILE>;

					#print $atom->{NAME},$line;
					chomp $line;
					@f = split(' ', $line);
					$atom->{CHARGE} = $f[2];
				}
				last;
			} else {

				silico_msg('w', "No charges found in output file\n");
			}
		}
		close INFILE;
	} else {

		my $newens;
		
		if (-e $marcfile) {
			
			$newens = read_mopac_cart($marcfile, undef, undef, "NOBOND")
		  	|| silico_msg('d', "Failed to read molecule data from file $marcfile\n");
		
		} elsif (-e $moutfile) {
		
			$newens = read_mopac_cart($moutfile, undef, undef, "NOBOND")
		  	|| silico_msg('d', "Failed to read molecule data from file $moutfile\n");
		
		} else {
		
			silico_msg('w', "Failed to read molecule data from file $moutfile or $marcfile\n");
			return;
		}
		  
		my $newmol = $newens->[0];
		
		# Copy charges and geometry back to original molecule
		my $i = 0;
		foreach my $atom (atoms($mol)) {

			$atom->{CHARGE} = $newmol->{ATOMS}[$i]{CHARGE};
			$atom->{X}      = $newmol->{ATOMS}[$i]{X};
			$atom->{Y}      = $newmol->{ATOMS}[$i]{Y};
			$atom->{Z}      = $newmol->{ATOMS}[$i]{Z};
			++$i;
		}
		
		foreach (@list) {
			$mol->{SDF_DATA}{$_} = $newmol->{SDF_DATA}{$_};
		}

	}
	
	$mol->{SYBCHARGETYPE} = "AMPAC_CHARGES";
	
	my $string = get_mopac_sdf_data_name();
	$mol->{SDF_DATA}{$string} = $mol->{SDF_DATA}{ENERGY};

	if (-e $marcfile) {
		open(MINFILE, $marcfile) || file_read_error($marcfile, 1);
	} elsif (-e $moutfile) {
		open(MINFILE, $moutfile) || file_read_error($moutfile, 1);
	}

	# Check to see if the molecule was a radical
	while (<MINFILE>) {
		if (/OPEN LEVELS/ || /SINGLY OCCUPIED/) {
			silico_msg('w', "Uneven number of electrons in calculation!\n");
			$mol->{WARN_RADICAL} = 1;
			last;
		}
	}
	close MINFILE;
	
	
	# Clean up
	if (!$noclean) {
	
		$run = "rm -f $minfile $marcfile $moutfile ";

		if (system($run)) {
			silico_msg('d', "Execution of cleanup (\"rm\") command failed!\n", "Command: $run\n");
		}
	}

	return $mol;
}

sub get_mopac_sdf_data_name {

	my $ham = get_lflag('hamiltonian');
	my $string = uc("MOPAC_ENERGY_$ham");
	
	return $string;
}


sub get_mopac_outfile_name {

	my $mol = $_[0];

	my $filebase = get_ofilebase($mol);
	my $fileext = get_oformat($mol);
	my $exist;
		
	$filebase =~ s/mop.*mop/mop/;
	
	if (-e "$filebase.$fileext") {
		my $n = 1;
		$exist = "$filebase.$fileext";
		while(-e "$filebase\_$n.$fileext") {
			$exist = "$filebase\_$n.$fileext";
			++$n;
		}
		$filebase = "$filebase\_$n";
	}
	
	print "\nOutfile: $filebase.$fileext\n";		
	
	return $filebase, $fileext, $exist;
}
##################################################################
#
#	misc
#
##################################################################



sub mol_count_valence_electrons {

	#<
	#? Simplistic routine to calculate the total number of valence electrons and work out if molecule
	#  is a radical or not
	#; Requires: Molecule, charge
	#; Returns: total number of valence electrons, flag denoting odd number of valence electrons
	#>

	my $mol    = $_[0];
	my $charge = $_[1] || 0;

	my $atom;
	my $tot = 0;
	my $unpaired;
	my $val;
	my @Valence_electrons;

	@Valence_electrons = qw(

	  0

	  1  2
	  1  2                                                                          3  4  5  6  7  8
	  1  2                                                                          3  4  5  6  7  8
	  1  2  Sc                                           Ti V  Cr Mn Fe Co Ni Cu Zn 3  4  5  6  7  8
	  1  2  Y                                            Zr Nb Mo Tc Ru Rh Pd Ag Cd 3  4  5  6  7  8
	  1  2  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W  Re Os Ir Pt Au Hg 3  4  5  6  7  8
	  1  2  Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr

	);

	foreach $atom (atoms($mol)) {

		$val = $Valence_electrons[ $atom->{ELEMENT_NUM} ];

		if ($val =~ /A-Z/) {

			silico_msg('w', "Unable to calculate valence electrons for element $val. Setting to 0.\n");
			$val = 0;
		}

		$tot += $val;
	}

	# Subtract charge from valence electrons
	$tot -= $charge;

	# Determine whether the number of electrons is odd
	# (in which case we have at least one unpaired electron)

	# This uses a bitwise AND to determine oddness
	# If the last binary digit of $tot is 0 (i.e., an even number)
	# then this, when ANDed to the number 1, will give a zero result
	if (($tot & 1) != 0) {
		$unpaired = 1;
	}

	return $tot, $unpaired;
}

##################################################################
#
#	DivCon
#
##################################################################

sub read_divcon_cart {

	#<
	#? Read a DivCon Cartesian format file
	#. Reads only a single structure
	#; Requires: infile, undef, undef, options (optional)
	#; Returns: ensemble or undef on read error
	#>

	my $infile  = $_[0];
	my $options = uc($_[3] || '');

	my $anisotropy;
	my $atom;
	my $count;
	my @f;
	my $mol;
	my $nmrflag;
	my $shielding;
	my $shift;

	# Options
	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
	if (!defined $fr) {
		silico_msg('e', "Can not open file $infile for reading\n");
		return undef;
	}
	my $FH = $fr->{FH};

	silico_msg('c', "Reading DivCon file: $infile\n") if $options !~ /\bQUIET\b/;

	mol_set_source_data($mol, $fr, 'divcon');
	$mol->{NUMBONDS} = 0;

	while (<$FH>) {

		if (/^ HEAT OF FORMATION/) {

			@f             = split;
			$mol->{ENERGY} = $f[4];
			$mol->{NAME}   = $infile;
			next;
		}

		#
		# Coordinates
		#

		if (/NUMBER     SYMBOL            X           Y           Z/) {

			<$FH>;    # Blank line

			$count = 0;
			while (<$FH>) {

				if (!defined $_) {
					silico_msg('e', "Premature end to file\n");
					return undef;
				}

				my $atom;
				@f = split;
				last if !$f[0];    # Exit loop on blank line

				$atom->{NUM}         = $count;
				$atom->{ELEMENT}     = ucfirst(lc($f[1]));
				$atom->{ELEMENT_NUM} = element_symbol2number($f[1]);
				$atom->{X}           = $f[2];
				$atom->{Y}           = $f[3];
				$atom->{Z}           = $f[4];
				$atom->{NAME}        = $f[1] . ($count + 1);
				$atom->{SUBID}       = 1;
				$atom->{SUBNAME}     = 'UNK';
				$atom->{CHAIN}       = 'A';

				push @{ $mol->{ATOMS} }, $atom;
				++$count;
			}

			$mol->{NUMATOMS} = $count;
			next;
		}

		#
		# Charges
		#

		if (/CHARGE     CHARGE     CHARGE/) {

			<$FH>;    # Blank line

			$count = 0;
			while (<$FH>) {

				if (!defined $_) {
					silico_msg('e', "Premature end to file\n");
					return undef;
				}

				chomp;

				@f = split;
				last if !$f[0];    # Exit loop on blank line

				$mol->{ATOMS}[$count]{MULLIKEN_CHARGE} = $f[2];
				$mol->{ATOMS}[$count]{CHARGE}          = $f[2];    # Sets charge to Mulliken charge
				$mol->{ATOMS}[$count]{CM1_CHARGE}      = $f[3];
				$mol->{ATOMS}[$count]{CM2_CHARGE}      = $f[4];
				++$count;
			}
		}

		#
		# NMR shieldings
		#

		if (/NUCLEAR MAGNETIC SHIELDING TENSOR/) {

			$nmrflag = 1;
			$count   = 0;
			while (<$FH>) {

				if (!defined $_) {
					silico_msg('e', "Premature end to file\n");
					return undef;
				}

				chomp;

				@f = split;
				last if !$f[0];    # Exit loop on blank line

				if (/\*\*  ATOM(....)/) {

					#print "$_ $1\n";
					$count = $1;    # Indexed from 1
				}

				if (/\*\* AVERAGE/) {
					my $atom = $mol->{ATOMS}[ $count - 1 ];
					$atom->{SHIELDING}  = $f[3];
					$atom->{ANISOTROPY} = $f[6];
				}
			}
		}
	}

	if (!$mol->{NUMATOMS}) {
		silico_msg('w', "No atoms read from $infile\n");
		return undef;
	}

	molecule_make_residues($mol);

	if ($nmrflag) {

		# Add shift offset according to Patchkovskii and Thiel J. Comp Chem, 12, 1220, 1999
		# Divcon does not seem to entirely use shift offset from above paper so
		# modified some offsets using pddg-pm3 optimised structures

		foreach $atom (atoms($mol)) {

			next if !defined $atom->{SHIELDING};

			if ($atom->{ELEMENT} eq 'C') {

				# Referenced so that benzene C is 137.2
				$atom->{SHIFT} = 47.60 - $atom->{SHIELDING};
			}
			if ($atom->{ELEMENT} eq 'H') {

				# Note: separate OH and NH corrections

				my $con = $mol->{ATOMS}[ $atom->{CONNECT}[0] ];

				if ($con->{ELEMENT} eq 'N') {

					# Reference so that formamide NHs average 4.3
					$atom->{SHIFT} = 50.6 - $atom->{SHIELDING};

				} elsif ($con->{ELEMENT} eq 'O') {

					# Methanol OH referenced to 0
					$atom->{SHIFT} = 47.0 - $atom->{SHIELDING};

				} else {

					# MB3
					$atom->{SHIFT} = 50.64 - $atom->{SHIELDING};
				}
			}
			if ($atom->{ELEMENT} eq 'N') {

				# Note: magnetogyric ratio of N is -ve

				$atom->{SHIFT} = -172.20 - $atom->{SHIELDING} + 288;    # MB3
			}
			if ($atom->{ELEMENT} eq 'O') {

				# Referenced so that formamide O is 370
				$atom->{SHIFT} = 29 - $atom->{SHIELDING};
			}
		}
	}

	atom_data_to_sdf($mol, qw(SHIFT SHIELDING ANISOTROPY));

	my $ens;
	$ens->[0] = $mol;

	return $ens;
}

sub write_divcon_cart {

	#<
	#? Print out a DivCon Cartesian format file
	#. Currently prints out only the first molecule
	#  of the ensemble
	#; Options: CLUSTER
	#; Flags: nmr (peptide, all), opt (num optimization steps)
	#; Requires: Ensemble (or molecule), filename.
	#; Returns: undef if file open fails otherwise returns 1.
	#>

	require silico_split;

	my $molecules = ens($_[0]);
	my $outfile   = $_[1];
	my $options   = uc($_[3] || '');

	my $atom;
	my $buff1 = "4.0";
	my $buff2 = "2.0";
	my $chg;
	my $cluster;
	my $dumpfreq;
	my $ha;
	my $i;
	my @f;
	my $keywords;
	my $len;
	my $lig;
	my $mol;
	my $maxit;
	my $name;
	my $optmethod;
	my $residues;
	my $subatoms;
	my $success;
	my $total_charge;

	$outfile = get_ofilebase($molecules) . ".in" if !$outfile;

	my $fr;
	$fr = open_molfile($fr, $outfile, undef, undef, $options);
	if (!defined $fr) {
		silico_msg('e', "Can not open file $outfile for reading\n");
		return undef;
	}
	++$fr->{MOL_WRITE_COUNT};
	
	my $FH = $fr->{FH};

	silico_msg('c', "Writing DivCon Cartesian file: $outfile\n");

	# Kludge
	$mol = $molecules->[0];

	#
	# Options
	#
	$ha = uc get_sflag('ha') || 'MNDO';
	$ha = 'MNDO/d' if $ha eq 'MNDOD';
	$ha = 'PDDG-PM3' if $ha eq 'PDDG';

	$maxit = get_sflag('maxit');

	$dumpfreq  = 0;
	$dumpfreq  = 5 if get_sflag('opt');
	$optmethod = uc get_sflag('optm');

	# Clustering
	if (get_sflag('cluster')) {

		$cluster = get_sflag('cluster');

	} elsif (get_sflag('res')) {

	} else {
		$cluster = 1 if $mol->{NUMATOMS} > 300;
	}

	#
	# Label things
	#
	molecule_check_and_fix_connectivity($mol);
	molecule_fix_and_get_residues($mol);
	mol_renumber_substructures($mol);
	molecule_label_submols($mol);
	$residues = molecule_get_residues($mol, undef, undef, 1);

	$chg = get_sflag('c');
	if (defined $chg) {
		silico_msg('c', "\nSetting total charge on system to be $chg\n");
		$total_charge = $chg;
	} else {
		$total_charge = molecule_formal_charge($mol);
		silico_msg('c', "\nTotal charge on molecule is calculated to be $total_charge\n");
	}

	#
	# Keywords
	#
	$keywords = "CARTESIAN $ha DIRECT CHARGE=$total_charge CUTBOND=8 PDUMP=2";
	if ($cluster) {
		$keywords .= ' CLUSTER RESIDUE';
	} else {
		$keywords .= ' STANDARD';
	}

	$keywords .= ' NMR'           if get_sflag('nmr') || $mol->{DIVCON_NMR_ATOMS};
	$keywords .= ' DPSCF=1.0D-9 ' if get_sflag('nmr') || get_sflag('strict');

	if (get_sflag('opt')) {

		$maxit ||= 100;

		$keywords .= " OPT=$optmethod MAXOPT=" . get_sflag('opt') . " MAXIT=$maxit ETEST=0.1 GTEST=1.5 ";
		$keywords .= " FORCE-IT ECRIT=1.0E-3 DCRIT=1.0E-3" if !get_sflag('strict');
	}

	$keywords .= ' ADDMM' if $ha =~ /MNDO/i || $ha =~ /AM1/i || $ha =~ /PM3/i;

	# Not implemented in current divcon
	#$keywords .= ' SNAPGEOM=2' if get_sflag('opt') && $mol->{NUMATOMS} > 100;
	$keywords .= ' SCRF'                       if get_sflag('scrf');
	$keywords .= ' SHIFT=' . get_sflag('shift') if get_sflag('shift');
	$keywords .= ' RESTART'                    if get_sflag('restart');
	$keywords .= ' DOUBLE=0'                   if get_sflag('double') && $cluster;
	$keywords .= ' GUESS'                      if get_sflag('guess');
	$keywords .= " MAXIT=$maxit"               if $maxit;
	$keywords .= " DUMP=$dumpfreq";
	$keywords .= " SCREEN";
	$keywords .= " " . (get_sflag('key') || '');

	$keywords = divcon_wrap_keywords($keywords);

	print $FH "$keywords\n";

	silico_msg('c', heading("Keywords"));
	silico_msg('c', "$keywords\n");

	#
	# Comments and coordinates
	#
	$name = $mol->{NAME} || $mol->{SOURCE_FILE_NAME} || "Mol";

	if ($name) {
		chomp $name;
		print $FH "$name\n";
	} else {
		print $FH "Mol\n";
	}

	for ($i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

		$atom = $mol->{ATOMS}[$i];

		printf $FH "%5d  %4s %12.4f  %12.4f  %12.4f", ($i + 1), $atom->{ELEMENT}, $atom->{X}, $atom->{Y}, $atom->{Z};

		if ($atom->{STARTRES}) {
			print $FH "  RES " if $atom->{STARTRES};
		} else {
			print $FH "      ";
		}

		printf $FH "(%s;%s;%d;%d;%d)\n", ($atom->{SUBNAME} || 'UNK'), ($atom->{NAME} || 'X'), ($atom->{SUBID} || 1), 1, 1;
	}

	print $FH "END_COORD\n";

	#
	# Clusters
	#

	if ($cluster) {

		print $FH "CLUSTER\n";
		silico_msg('c', heading("Charged residues"));

		my $num_residues = $#{$residues} + 1;

		#print $FH "     NCORE=1 (1-$num_residues)\n";
		print $FH "     NCORE=1\n";

		print $FH "    DBUFF1=$buff1 DBUFF2=$buff2\n";
		print $FH "END_CLUSTER\n";
	}

	if (get_sflag('res')) {

		molecule_fix_and_get_residues($mol);

		my $count;
		my $scharge;

		print $FH "CLUSTER\n";
		silico_msg('c', heading("Charged residues"));

		$count = 1;
		foreach $subatoms (@{ $mol->{RESIDUES} }) {

			$scharge = 0;
			foreach $atom (@$subatoms) {
				$scharge += $atom->{FORMAL_CHARGE} || 0;
				(printf "Atom: %-4s %-4s %-4d %3d\n", $atom->{SUBNAME}, $atom->{NAME}, $atom->{SUBID}, $atom->{FORMAL_CHARGE})
				  if $atom->{FORMAL_CHARGE};
			}

			print $FH "    NCORE=1  ($count) [$scharge]\n";

			++$count;
		}
		print $FH "    DBUFF1=$buff1 DBUFF2=$buff2\n";
		print $FH "END_CLUSTER\n";
	}

	#
	# Guess
	#
	if (get_sflag('guess')) {

		print $FH "GUESS\n";
		print $FH "    " . (get_sflag('guess')) . "\n";
		print $FH "END_GUESS\n";
	}

	#
	# NMR
	#
	if (get_sflag('nmr') || $mol->{DIVCON_NMR_ATOMS}) {

		my $alist;
		my $count;
		my $last;

		silico_msg('c', heading("NMR"));

		# Select atoms for NMR cacluation using the option specified with the
		# 'nmr' command line flag.
		$alist = divcon_select_nmratoms($mol, get_sflag('nmr'));

		if ($#{$alist} >= 0) {

			print $FH "NMR\n";
			print $FH "    ATOMS ";

			$last  = $alist->[-1];
			$count = 0;
			foreach (@$alist) {
				++$count;
				print $FH "    " if $count == 1;
				print $FH "$_ ";
				if ($count == 10) {
					$count = 0;
					print $FH " &" if $_ != $last;
					print $FH "\n";
				}
			}
			print $FH "\n" if $count != 0;
			print $FH "END_NMR";

			silico_msg('c', "Calculating shift for " . ($#{$alist} + 1) . " atoms\n");
		} else {

			silico_msg('c', "No NMR atoms selected.  No NMR caculations will be run\n");
		}
	}

	close_molfile($fr);

	silico_msg('c', "\n");

	return 1;
}

sub divcon_select_nmratoms {

	# Can be run before running write_divcon_cart or within it.

	my $mol  = shift;
	my $flag = uc(shift) || 'ALL';    # String containing keywords: all, core

	my $alist;
	my $atom;
	my $count;

	#my $nmr_cutoff2 = get_sflag('cut') || 12; # Cutoff for determining nmr shifts if lig flag is used
	#my $nmr_cutoff1 = 7;			 # Close cutoff to determine these shifts first

	# List of elements for shift calculations
	my @pep = qw(AA_CA AA_CB AA_N AA_HN);

	# Return saved data if we have already run this subroutine
	return $mol->{DIVCON_NMR_ATOMS} if defined $mol->{DIVCON_NMR_ATOMS};

	silico_msg('c', "\nSelecting atoms for DIVCON NMR cacluation using flag: $flag\n");

	$count = 0;
      ATOM: foreach $atom (atoms($mol)) {
		++$count;

		# Select only CORE atoms
		if ($flag =~ /\bCORE\b/) {
			next ATOM if !$atom->{CORE};
		}

		# Select only peptide backbone atoms
		if ($flag =~ /\bPEPTIDE\b/) {
			foreach (@pep) {
				$flag = 0;
				next if !$atom->{FG}{$_};
				$flag = 1;
				last;
			}
			next ATOM if !$flag;
		}

		# Select only peptide backbone atoms
		if ($flag =~ /\bCHN\b/) {
			my $el = $atom->{ELEMENT};
			next ATOM if $el ne 'C' && $el ne 'H' && $el ne 'N';
		}

		push @$alist, $count;
	}

	print "Selected " . ($#{$alist} + 1) . " atoms\n";

	$mol->{DIVCON_NMR_ATOMS} = $alist;

	return $alist;
}

sub divcon_wrap_keywords {

	my @f;
	my $i;
	my $len;

	my $keywords = $_[0];

	@f        = split " ", $keywords;
	$keywords = '';
	$len      = 0;
	$i        = 0;
	foreach (@f) {

		if ($len + length($_) + 1 > 78) {
			$len = 0;
			$keywords .= "&";
			$keywords .= "\n";
		}

		$keywords .= $_ . " ";
		$len += length($_) + 1;
		++$i;
	}

	return $keywords;
}

return 1;
