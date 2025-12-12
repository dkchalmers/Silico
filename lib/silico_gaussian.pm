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
#! silico_gaussian.pm
#? Silico routines for reading and writing Gaussian Z-matrix and Cartesian files
#. $Revision: 1.8.2.1.2.2 $
#>

use strict;
package Silico;


##################################################################
#
#	Gaussian read routines
#
##################################################################

sub read_gaussian_cart {

	#<
	#? Simplistic routine to extract Cartesian coordinates from
	#  a Gaussian (98) output file using the standard orientation
	#. Note all structures are read by default.  To read only the 
	#  final structure, use the --last-structure flag
	#; Requires: filename.
	#; Returns: ensemble or zero if read failed.
	#>

	my $infile = $_[0];
	my $options = uc ($_[3] || '');
	
	my $aname;
	my $i;
	my @l;
	my %element_count;
	my $molecules;
	my $molcount = 0;
	
	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
        if (!defined $fr) {
                silico_msg('e', "Can not open file $infile for reading\n");
                return undef;
        }
	my $FH = $fr->{FH};
	
	silico_msg('c', "Reading gaussian Cartesian file: $infile\n") if $options !~ /\bQUIET\b/;

	while (<$FH>) {

		# Marker for start of Cartesian coordinates
		if (!/Standard orientation:/) {
			next;
		}
		<$FH>; # dashes
		<$FH>; # words
		<$FH>; # words
		<$FH>; # dashes

		my $mol;

		# Read in body of input file
		$i = 0;
		while (<$FH>) {

			chomp;

			#We have finished if we match a whole lot of dashes
			if (/-----/ ) {
				last;
			}

			@l = split;
			silico_msg('g', "@l\n");

			my $atom;

			$atom->{NUM} = $i+1;			#Atom number
			$atom->{ELEMENT_NUM} = $l[1];		#Atom element number
			$atom->{SUBNAME} = "MOL";		#Residue name
			$atom->{SUBID} = 1;			#Residue number
			$atom->{X} = $l[3];			#x
			$atom->{Y} = $l[4];			#y
			$atom->{Z} = $l[5];			#z

			# Calculated fields
			$atom->{ELEMENT} = element_number2symbol($atom->{ELEMENT_NUM});
			$atom->{NAME} = $atom->{ELEMENT}.(++$element_count{$atom->{ELEMENT}});
	
			# Set atom element number
			$aname = ucfirst lc $atom->{ELEMENT};

			# Assign atom record to molecule
			$mol->{ATOMS}[$i] = $atom;

			++$i;
		
			$mol->{NAME} = $infile;
			$mol->{NUMATOMS} = $i;
			$mol->{NUMBONDS} = 0;
		}
		
		mol_set_source_data($mol, $fr, 'gauss');

		$molecules->[$molcount] = $mol;
		++$molcount;
	}
	
	close ($FH);

	if (!$molcount) {
		silico_msg('e', "No strucures were read from file $infile!\n",
				"Skipping.\n");
		return undef;
	}
	
	if (get_flag('last-structure', 'l')) {
	
		silico_msg('w', "Retaining only last structure from input file\n");
		my $newmolecules;
		$newmolecules->[0] = $molecules->[-1];
		$molecules = $newmolecules;
	}
	
	my $mol;
	$i = 0;
	foreach $mol (@$molecules) {
		++$i;
		silico_msg('c', "$mol->{NUMATOMS} atoms read in molecule $i\n");
		molecule_check_and_fix_connectivity($mol) if ($options !~ /\bNOBOND\b/);
	}

	return $molecules;
}



##################################################################
#
#	gaussian write routines
#
##################################################################


sub write_gaussian_cart {

	#<
	#? Print out a gaussian Cartesian format file
	#; Requires: ensemble (or molecule), output filename
	#; Returns: undef if write failed
	#>
      
	my $molecules = ens($_[0]);
	my $outfile = $_[1];
	my $options = uc ($_[3] || '');
	
	my $atom;
	my $i;
	my $mol;
	my $name;
	my $success;

	$outfile = get_ofilebase($molecules).".dat" if !$outfile;

	$mol = $molecules->[0];
       
	$success = open (OUTFILE, ">$outfile");
	if (!defined $success) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing gaussian Cartesian file: $outfile\n");
	
	for ($i=0; $i < $mol->{NUMATOMS}; ++$i) {

		$atom = $mol->{ATOMS}[$i];

		printf OUTFILE " %-3s %12.8f %12.8f %12.8f",
			$atom->{ELEMENT_NUM},$atom->{X}, $atom->{Y},$atom->{Z};
			print OUTFILE "\n";
	}
	
	close (OUTFILE);

	return 1;
}

sub write_gaussian_zmatrix {

	#<
	#? Gaussian Z-matrix writer.
	#; Requires: ensemble or molecule, filename, forcefield file (optional) to
	#  type atoms.
	#; Returns: undef if file open failed otherwise returns 1.
	#>

	my $ensemble = ens($_[0]);
	my $outfile = $_[1];
	my $options = uc $_[2] || '';
	
	my $angle;
	my $atom;
	my $distance;
	my $ext = 'zmat';
	my $filebase;
	my $i;
	my $mol;
	my $molcount;
	my $outfile;
	my $string;
	my $success;
	my $torsion;
	
	$outfile = get_ofilebase($ensemble).".zmat" if !$outfile;
	$filebase = get_filebase($outfile);
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	ensemble_printout($ensemble, 'all') if $Silico::debug;
	
	$molcount = -1;
	
	foreach $mol (@$ensemble) {
		
		++$molcount;
		
		if (!defined $mol->{NUMATOMS}) {
			silico_msg('e', "\$mol->NUMATOMS is not defined!\n");
			return undef;
		}
		
		if ($molcount > 0) {
			$outfile = $filebase.sprintf "%04d", $molcount.".zmat";
		}
		
		$success = open (GAUOUT, ">$outfile");
		
		if (!$success) {
			silico_msg('e', "Can not create or open file $outfile for writing!\n");
			return undef;
		}
		
		if ($options !~ /\bQUIET\b/) {
			silico_msg('c', "Writing zmat file $outfile");
			silico_msg('c', " using $options") if $options;
			silico_msg('c',  "\n");
		}
		
		for ($i=0; $i <= $#{$mol->{ATOMS}}; ++$i) {
			
			# Reset variables
			$atom = $mol->{ATOMS}[$i];
			$distance = 999.999999;
			$angle = 999.999999;
			$torsion = 999.999999;
			$string = '';
			
			# Element
			$string .= sprintf "%3s", $atom->{ELEMENT};
			
			# Bond
			if ($i > 0) {
				$distance = distance($atom, $mol->{ATOMS}[$i-1]);
				$string .= sprintf "%5d %10.6f",  $i, $distance;
			}
			
			# Angle
			if ($i > 1) {
				$angle = bndangle($atom, $mol->{ATOMS}[$i-1], $mol->{ATOMS}[$i-2]);
				$string .= sprintf "%5d %10.6f", $i-1, $angle;
			}
			
			# Dihedral
			if ($i > 2) {
				$torsion = rad_to_deg(dihedral($atom, $mol->{ATOMS}[$i-1], $mol->{ATOMS}[$i-2], $mol->{ATOMS}[$i-3]));
				$string .= sprintf "%5d %10.6f", $i-2, $torsion;
			}
			
			$string .= "\n";
			
			print GAUOUT "$string";
		}
		close (GAUOUT);
	}

	return 1;
}


return 1;
