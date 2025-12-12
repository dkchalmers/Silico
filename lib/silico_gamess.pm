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
#! silico_gamess.pm
#? Silico routines for reading and writing Gamess files
#. $Revision: 1.7.2.1.2.1 $
#>

use strict;
package Silico;


##################################################################
#
#	Gamess read routines
#
##################################################################

sub read_gamess_cart {

	#<
	#? Simplistic routine to extract Cartesian coordinates from
	#  a Gamess output file 
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

	silico_msg('c', "Reading Gamess Cartesian file: $infile\n") if $options !~ /\bQUIET\b/;

	while (<$FH>) {

		# Marker for start of Cartesian coordinates
		next if !/^ COORDINATES OF ALL ATOMS/;
		
		<$FH>; # headers
		<$FH>; # dashes

		my $mol;

		# Read in body of input file
		$i = 0;
		while (<$FH>) {

			chomp;

			#We have finished if we match a blank line
			if ($_ eq '' ) {
				last;
			}

			@l = split;
			silico_msg('g', "@l\n");

			my $atom;

			$atom->{NUM} = $i+1;			#Atom number
			$atom->{ELEMENT} = ucfirst (lc $l[0]);	#Atom element 
			$atom->{SUBNAME} = "MOL";		#Residue name
			$atom->{SUBID} = 1;			#Residue number
			$atom->{X} = $l[2];			#x
			$atom->{Y} = $l[3];			#y
			$atom->{Z} = $l[4];			#z

			# Calculated fields
			$atom->{ELEMENT_NUM} = element_symbol2number($atom->{ELEMENT});
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
		
		mol_set_source_data($mol, $fr, 'gamess');

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




return 1;
