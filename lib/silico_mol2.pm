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
#! silico_mol2.pm
#? Silico routines to deal with Sybyl mol2 files.
#. $Revision: 1.15-46-gd5da50c $
#>

use strict;
package Silico;

##################################################################
#
#	Read Mol2 routines
#
##################################################################

sub read_mol2 {

	#<
	#? Read in a mol2 file and return as an ensemble
	#. Reads multiple structures from a single file. Information
	#  is currently read from the MOLECULE, ATOM, BOND,
	#  SUBSTRUCTURE and SET records.  Only static atom sets are supported.
	#  Comment lines preceeding the molecule data are preserved.
	#  Other records are ignored.
	#. DOCK binding energies are extracted from the Comment lines and stored
	#  as $mol->{DOCK_ENERGY}
	#; Requires: Filename, (optional), start record, maximum number of records, options string .
	#; Returns: ensemble or undef if no molecule is read.
	#. Options:
	#, Will respond to a -me (maximum energy) flag
	#, QUIET - do not print 'Reading' line
	#, NOSORT - do not sort atoms
	#>
	
	my $infile = $_[0];
	my $start = $_[1] || get_flag('ss', 's') || 0;  # Note, it is not clear that this option works!!
	my $max_molecules = $_[2] || get_flag('ms', 's') || 1000000;
	my $options = uc ($_[3] || '');
	
	my $energy_start; # Energy of first compound in file
	my $fr;
	my $mols;	  # Ensemble
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	# Maximum desired molecule energy
	my $energy_max = get_flag('me', 's');

	$fr = open_molfile($fr, $infile, undef, undef, $options);
        if (!defined $fr) {
		silico_msg('e', "Can not open file $infile for reading\n");
		return undef;
	}

	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Reading mol2 file: $infile");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}
	
	my $i = -1;
	my $molcount = 0;
        while (my $mol = read_mol2_single($fr)) {

                ++$i;
                next if $start and ($i < $start);
                ++$molcount;
                push @$mols, $mol;
                last if $molcount == $max_molecules;

		# Maxenergy option
		if (defined $energy_max) {
			if ($molcount == 0) {
				if (defined $mol->{ENERGY}) {
					$energy_start = $mol->{ENERGY};
				} else {
				
					silico_msg('w', "Could not obtain energy of initial structure.  Ignoring maxenergy option\n");
				} 
			}
			if (defined $energy_start && defined $mol->{ENERGY}) {
			
				print "en".($mol->{ENERGY} - $energy_start)."\n";
				if ($mol->{ENERGY} - $energy_start > $energy_max) {
					silico_msg('n', "Reached energy maximum: $energy_max max\n");
					last;
				}
			}
		}
        }

	# Check that we read in a molecule
	if ($molcount < 0) {
		silico_msg('e', "No molecule data was read in from file $infile!\n",
				"The file may be empty or is not a properly formatted mol2.\n");
		$fr->{ERROR} = 1;
		return undef;
	}

	#
	# Error messages
	#
	my $warn1 = 0;
	my $warn2 = 0;
	foreach my $mol (@$mols) {
		++$warn1 if $mol->{WARN}{MOLECULE_SORT};
		++$warn2 if $mol->{WARN}{FIXED_VMD_BONDORDERS};
	}
	silico_msg('n', "Substructures in $warn1 molecule".($warn1 == 1 ? '' : 's')." were out of order and have been sorted.\n") if $warn1;
	silico_msg('w', "Correcting bad Sybyl atomtypes and bond orders produced by VMD!\n") if $warn2;

	close_molfile($fr);

	return $mols;
}

sub read_mol2_single {

        #<
        #? Read in a single record from a mol2 file.
        #. See read_mol2 for general description
        #; Requires: file record,  options
        #; Returns: molecule or undef
        #>

        my $fr = $_[0];
        my $options = uc($_[1] || '');

	my $comments;
	my $flag = 0;
	my $line;
	my $stride;

	# Make a new molecule 
	my $mol = {};
	mol_set_source_data($mol, $fr, 'mol2');

	$stride = get_flag('stride', 's');
	if ($stride) {
		while (1) {
			last if $stride == 1;
	  		$line = silico_readline($fr);
               		if (!defined $line) {
        	   	    last;
               		}
			if ($line =~ /\@<TRIPOS>MOLECULE/ || $line =~ /^#/) {
				--$stride;
				++$fr->{MOL_READ_COUNT};
				silico_putline($fr, $line);
               		}
		}
	}

	# Read any comments that come before the @<TRIPOS>MOLECULE line
	while (1) {

		$line = silico_readline($fr);
		return undef if !defined $line;
		if ($line =~ /^#/) {
			push @$comments, $line;
			next;
		}

		if ($line =~ /\@<TRIPOS>MOLECULE/) {
			silico_putline ($line, $fr);

			# Extract any SDF data contained in comments
			$comments = extract_mol2_sdf_data($mol, $comments);
			# Process information in comments provided by Dock etc
			mol2_process_comments($mol, $comments);
			last;
		}
	}

	while (1) {

		$line = silico_readline($fr);
		last if !defined $line; # We have come to the end of the file

		if ($line =~ /^#/ || ($flag == 1 && $line =~ /\@<TRIPOS>MOLECULE/)) {
			silico_putline ($line, $fr);
			last;
		}
		if ($line =~ /\@<TRIPOS>MOLECULE/) {
			if (!defined mol2_molecule_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;
			}
			$flag = 1;	
			next;
		}
		if ($line =~ /\@<TRIPOS>ATOM/) {
			if (!defined mol2_atom_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;
			}
			next;
		}
		if ($line =~ /\@<TRIPOS>BOND/) {
			if (!defined mol2_bond_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;
			}
			next;
		}
		if ($line =~ /\@<TRIPOS>SUBSTRUCTURE/) {
			if (!defined mol2_substructure_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;	
			}
			next;
		}
		if ($line =~ /\@<TRIPOS>SET/) {
			if (!defined mol2_set_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;
			}
			next;
		}
		if ($line =~ /\@<TRIPOS>CRYSIN/) {
			if (!defined mol2_crystal_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;
			}
			next;
		}
		if ($line =~ /\@<TRIPOS>CENTROID/) {
			if (!defined mol2_centroid_record($mol, $fr)) {
				$fr->{ERROR} = 1;
				return undef;
			}
			next;
		}
	}

	# Fix bad VMD mol2 files which lack Sybyl atom types and bond orders
	if ($mol->{NAME} =~ /generated.*by.*VMD/) {
		mol_type_sybyl($mol);
		molecule_generate_bondorders($mol);
		$mol->{WARN}{FIXED_VMD_BONDORDERS} = 1;
	}
	$mol->{NUM} = $fr->{MOL_READ_COUNT}+1;
	
	# Test that we have not exeeded maximum molecule energy or countset by flags
	return undef if test_energy_molcount($fr, $mol);

	return $mol;
}
	
		
sub mol2_process_comments {

	#<
	#? Extract data from mol2 comments strings such as Dock energies and add it to a molecule
	#. Called from read_mol2.
	#; Requires: molecule, comments
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $comments = $_[1];

	my $dock_energy;

	foreach (@$comments) {
	
		my @f = split;

		# Dock energy score (Dock 4)
		# This value is also asignted to ENERGY to
		# retain compatablity with some programs
		# Use of this value this way is discouraged!
		# The vlaue is also added to SDF_DATA (use of this is preferred)
		if (/Energy score/) {
			
			$dock_energy = $f[4];
			$mol->{DOCK_ENERGY}= $dock_energy;
			$mol->{ENERGY}= $f[4];
			$mol->{ENERGY_TYPE} = "Dock4_energy";
			$mol->{ENERGY_UNITS} = "UNK";
			$mol->{SDF_DATA}{DOCK_ENERGY} = $f[4];
		}
		
		# Dock energy score (Dock 5)
		# This value is also asignted to ENERGY to
		# retain compatablity with some programs
		# Use of this value this way is discouraged!
		# The vlaue is also added to SDF_DATA (use of this is preferred)
		if (/Energy Score:/) {

			$dock_energy = $f[3];
			$mol->{DOCK_ENERGY}= $f[3];
			$mol->{ENERGY}= $f[3];
			$mol->{ENERGY_TYPE} = "Dock5_energy";
			$mol->{ENERGY_UNITS} = "UNK";
			$mol->{SDF_DATA}{DOCK_ENERGY} = $f[3];
		}
		if (/vdw:/) {

			$mol->{DOCK_VDW}= $f[2];
			$mol->{SDF_DATA}{DOCK_VDW} = $f[2];
		}
		if (/es:/) {

			$mol->{DOCK_ES} = $f[2];
			$mol->{SDF_DATA}{DOCK_ES} = $f[2];
		}
	}

	# Save comments in molecule
	foreach (@$comments) {
		$mol->{COMMENTS} .= ($_."\n");
	}

	return;
}

sub mol2_molecule_record {

	#<
	#? Read in a @<TRIPOS>MOLECULE record.
	#. Called from read_mol2.
	#; Requires: molecule
	#; Returns; 1 or undef on failure
	#>
	
	my $mol = $_[0];
	my $fr = $_[1];
	
	my $error;
	
	# Molecule name
	defined (my $line = silico_readline($fr)) or return undef;
	$mol->{NAME} = $line;
	
	# Numatoms, numbonds, numsubstruct, num features (unused), numsets
	defined ($line = silico_readline($fr)) or return undef;
	my @f = split(' ', $line);
	
	$mol->{NUMATOMS} = $f[0] || 0;
	$mol->{NUMBONDS} = $f[1] || 0;
	$mol->{NUMSUBSTRUCT} = $f[2] || 0;
	#$mol->{NUM_FEATURES} = $f[3] || 0;
	$mol->{NUMSETS} = $f[4] || 0;

	# Small molecule or protein
	defined ($line = silico_readline($fr)) or return undef;
	$mol->{SYBTYPE} = $line || 'SMALL';
	
	# Charge type
	defined ($line = silico_readline($fr)) or return undef;
	$mol->{SYBCHARGETYPE} = $line || 'USER_CHARGES';
	
	# Read in two lines following the charge type line
	# Take into account that these lines are optional by checking that
	# the line is not a TRIPOS<ATOM> line
	defined ($line = silico_readline($fr)) or return undef;
	if ($line =~ /^@<TRIPOS>ATOM/) {
		$error = &mol2_atom_record;
		return undef if !defined $error;
		return 1;
	}
	#$mol->{SYBSTATUS} = $line || '';
	
	defined ($line = silico_readline($fr)) or return undef;
	if ($line =~ /^@<TRIPOS>ATOM/) {
		$error = &mol2_atom_record;
		return undef if !defined $error;
		return 1;
	}
	$mol->{SYBCOMMENT} = $line || '';

	return 1;
}

sub mol2_atom_record {

	#<
	#? Read in a @<TRIPOS>ATOM record.
	#. Called from read_mol2. Checks to see if more than half the atoms have
	#  zero Z coordinates, and sets mol->{HAS_BAD_GEOM} if this is the case.
	#; Requires: molecule
	#; Returns: 1 if successful, undef otherwise
	#>
	
	my $mol = $_[0];
	my $fr = $_[1];

	my $z_zerocount = 0;	# Counter to find if molecule is 2D

	for (my $i = 0; $i < $mol->{NUMATOMS}; ++$i) {
	
		
		my $line = silico_readline($fr);
		if (!defined $line) {
			silico_msg('e', "Premature end of file!\n");
			return undef;
		}
		my @f = split(' ', $line);
		
		if ($#f < 5) {
			silico_msg('e', "Too few fields in atom record ".($i+1)."!\n",
					"Line: \"$line\"\n");
			return undef;
		}
		
		# Note atom numbers are indexed from zero
		my $atom;
		$atom->{NUM} = $f[0];
		$atom->{NAME} = $f[1];
		$atom->{X} = $f[2];
		$atom->{Y} = $f[3];
		$atom->{Z} = $f[4];
		$atom->{TYPE} = $f[5] || 'Du';
		$atom->{SUBCOUNT} = $f[6] || 1;
		$atom->{RESFIELD} = $f[7];
		$atom->{CHARGE} = $f[8] if defined $f[8];
		
		++$z_zerocount if $atom->{Z} == 0;
		
		# Set atom element
		($atom->{ELEMENT}, $atom->{ELEMENT_NUM}) = atom_guess_element_mol2($atom, $mol);
		
		$mol->{ATOMS}[$i] = $atom;
	}
	
	if ($mol->{WARN}{GUESSED_ATOM_ELEMENTS}) {
	
		silico_msg('w', "Guessed element type for $mol->{WARN}{GUESSED_ATOM_ELEMENTS} atoms\n");
	}
	
	# Set bad geometry flag if the z value was zero for the majority of atoms
	if ($z_zerocount > $mol->{NUMATOMS}/2 + 2) {
		$mol->{HAS_BAD_GEOM} = 1;
		silico_msg('n', "Molecule appears to be 2D\n");
	}
	
	return 1;
}

sub mol2_bond_record {

	#<
	#? Read in @<TRIPOS>BOND record.
	#. Called from read_mol2.
	#; Requires: molecule
	#; Returns; 1 or undef on failure
	#>
	
	my $mol = $_[0];
	my $fr = $_[1];

	for (my $i = 1; $i <= $mol->{NUMBONDS}; ++$i) {
	
		my $line = silico_readline($fr);
		if (!defined $line) {
			silico_msg('e', "Premature end of file!\n");
			return undef;
		}
		
		my @f = splitline($line);
		
		if ($#f < 3) {
			silico_msg('e', "Missing information in bond record $i!\n",
					"Molecule: '$mol->{NAME}'\n",
					"Line: '$line'\n");
			return undef;
		}
		
		# This error doesn't matter
		#if ($f[0] != $i) {
			#silico_msg('w', "Bond record number $i has incorrect ID (\"$f[0]\").\n");
		#}
		
		my $batom1 = $f[1];
		my $batom2 = $f[2];
		my $border = $f[3];
		
		$border = 4 if $border eq "ar";
		$border = 1 if $border eq "am";
		$border = 0 if $border eq "du";
		
		
		if ($border eq 'un') {
			silico_msg('w', "Found Sybyl 'un' (Unknown) bond between atoms $batom1 and $batom2.\n",
					"    Converting to a single bond.\n");
			$border = 1;
		}
		if ($border eq 'nc') {
			silico_msg('w', "Found Sybyl 'nc' (Not connected) bond between atoms $batom1 and $batom2.\n",
					"    Converting to a zero order bond.\n");
			$border = 1;
		}
		
				
		# Make connection table
		# Note that connected atoms are indexed from 0
		# This differs from the Sybyl input file which is indexed from 1
		# mol->NUMBONDS is not incremented because this is already set
		
		my $val = bond_create($mol, $batom1-1, $batom2-1, $border, 1);
		
		if (!$val) {
			silico_msg('w', "Could not create a bond based on bond record $i!\n",
					"Line: \"$line\"\n",
					"Skipping.\n");
		}
	}
	
	remove_aromatic_bondorders($mol);
	
	return 1;
}

sub remove_aromatic_bondorders {

	#<
	#? Convert aromatic bondorders in carboxylate and guanidinium groups
	#  to single and double bonds
	#; Requires: molecule
	#; Returns; 1 or undef on failure
	#>
	
	my $mol = $_[0];
	
	my $arO;
	my $arN;
	my $flag;
	my $order;
	my $total;
	
	# Ensure that atom numbering is correct for bond_modify_order
	my $i = 1;
	foreach my $atom (atoms($mol)) {
		$atom->{NUM} = $i;
		++$i;
	}

	foreach my $atom (atoms($mol)) {
	
		next if $atom->{ELEMENT} ne 'C';
		
		$arO = 0;
		$arN = 0;
		$i = -1;
		
		foreach my $bo (@{$atom->{BORDERS}}) {
		
			++$i;
			next if $bo != 4;
			
			my $con = $mol->{ATOMS}[$atom->{CONNECT}[$i]];
		
			++$arO if $con->{ELEMENT} eq 'O';
			++$arN if $con->{ELEMENT} eq 'N';
		}
		
		# Convert aromatic carboxylate to C=O(-O-)
		if ($arO == 2) {
		
			$order = 1;
			foreach my $con (connected($atom, $mol)) {
				next if ($con->{ELEMENT} ne 'O');
				
				bond_modify_order($mol, $atom, $con, $order);
				$con->{FORMAL_CHARGE} = -1 if $order == 1;
				$con->{FORMAL_CHARGE} = 0 if $order == 2;
				++$order;
			}
			next;
		}
		
		# Convert aromatic guanidine to double and single bonds
		# Note that converting them to a representation with a formal
		# charge on carbon causes problems (eg in the 'add hydrogens' routine)
		if ($arN == 3) {
		
			$total = 0;
			
			# Make all bonds single
			foreach my $con (connected($atom, $mol)) {
				bond_modify_order($mol, $atom, $con, 1);
				$con->{FORMAL_CHARGE} = 0;
				$total += valence ($con);
			}
					
			if ($total == 9) {
			
				# Make double bond to the first connected atom
				my $con = $mol->{ATOMS}[$atom->{CONNECT}[0]];
				bond_modify_order($mol, $atom, $con, 2);
				$con->{FORMAL_CHARGE} = 1;
			
			} elsif ($total <= 8) {
				foreach my $con (connected($atom, $mol)) {
				
					# Make double bond to first connected atom with valence of 2
					# or less
					if (valence($con) <= 2) {
						bond_modify_order($mol, $atom, $con, 2);
						last;
					}
				}
			
			} else {
				silico_msg('w', "Badly formed guanidine\n");
			}
		}
	}
}

sub mol2_substructure_record {

	#<
	#? Read in a @<TRIPOS>SUBSTRUCTURE record.
	#. Called from read_mol2.  Creates the following arrays: mol->{SUB_SUBNAME}, an array of
	#  substructure names; mol->{SUB_ROOTATOM}, an array of root atoms; mol->{SUB_CHAIN} (optional)
	#  an array of chains
	#; Requires: molecule
	#; Returns; 1 or undef on failure
	#>

	my $mol = $_[0];
	my $fr = $_[1];

	my @chain;
	my @dicttype;
	my @resnum;
	my $restype;
	my @restype;
	my @rootatom;
	my $subnum;
	my @subname;
	my @subtype;
	my $warn_improper_restype;
	
	for (my $i = 0; $i < $mol->{NUMSUBSTRUCT}; ++$i) {
		
		my $line = silico_readline($fr);
		if (!defined $line) {
			silico_msg('w', "Not enough substructures were found in SUBSTRUCTURE record!\n",
					"Number of substructures declared: $mol->{NUMSUBSTRUCT}, number actually found: $i\n");
			last;
		}
		
		my @f = split(' ', $line);

		$subname[$i] = $f[1];		# Residue name concatenated with residue number
		$rootatom[$i] = $f[2];		# Root atom of residue (CA?)
		$subtype[$i] = $f[3] || '';     # TEMP, PERM, RESIDUE, GROUP, DOMAIN
		#$dicttype[$i] = $f[4] || 0;    # INTEGER
		$chain[$i]  = $f[5] || ' ';     # Chain
		$restype[$i] = $f[6] || '';     # Residue type?
		#$interbonds[$i]  = $f[7] || ''; # Number of interresidue bonds
		#$status[$i]  = $f[8] || '';	# ROOT
		#$comment[$i] =  $f[9] || '';   #  All remaining text
		
		$chain[$i] =~ s/\*//g;
		$chain[$i] = ' ' if !$chain[$i];
		
		# Check to make sure the root atom is valid
		if ($rootatom[$i] > $mol->{NUMATOMS}) {
			silico_msg('w', "Substructure $i (\"$subname[$i]\") has as its root atom\n",
					"a nonexistent atom: $rootatom[$i]!\n");
		}
	}

	# Assign chain name to all atoms in residue
	my $j = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {

		++$j;
		
		$subnum = $atom->{SUBCOUNT}-1;
		$atom->{CHAIN} = $chain[$subnum];

		$restype = $restype[$subnum];
		if (!defined $restype) {
		} else {
			$restype =~ s/\*|//g;
		}
		
		# Attempt to disentangle residue name and residue number in a reliable manner
		if ($restype && $atom->{RESFIELD} =~ m/^$restype(.*)/) {
		
			$atom->{SUBID} = $1;
			$atom->{SUBNAME} = $restype;
			
		} else {
		
			# Obtain substructure number from trailing digits on substructure name
			# Remove angle brackets
			$atom->{SUBID} = $atom->{RESFIELD};

			# Strip trailing numbers off substructure names
			$atom->{SUBNAME} = $atom->{RESFIELD};
			$atom->{SUBNAME} =~ s/\d+$//;
			
			# Fix Sybyl undefined substructure names
			$atom->{SUBNAME} = 'UNK' if $atom->{SUBNAME} =~ '<';
		}
				
		# Remove non-numeric characters from the substructure ID.
		$atom->{SUBID} =~ s/\<|\>//g;
		# Remove leading nondigit characters
		$atom->{SUBID} =~ s/^\D+//g;
		# If there is a trailing alpha character this belongs in
		# $atom->{EXT}
		if ($atom->{SUBID} =~ /\D$/) {
			$atom->{EXT} = $atom->{SUBID};
			$atom->{EXT} =~ s/^\d+//;
			$atom->{SUBID} =~ s/\D$//;
		}
		# Remove any nondigit characters from SUBID
		$atom->{SUBID} =~ s/\D//g;
			
		# Fallback value for substructure ID: substructure count.
		# but beware of residues actually numbered '0' which can be 
		# generated by Sybyl
		$atom->{SUBID} = $atom->{SUBCOUNT} if !$atom->{SUBID} && $atom->{SUBID} ne '0';
		
		# Fallback value for substructure name: UNK.
		$atom->{SUBNAME} ||= 'UNK';
		
		# Throw away the RESFIELD.
		delete $atom->{RESFIELD};
	}

	if ($warn_improper_restype) {

		silico_msg('w', "$warn_improper_restype atom were found belonging to improperly defined residues\n");
	}

	return 1;
}

sub mol2_set_record {

	#<
	#? Read in a @<TRIPOS>SET record.
	#. Called from read_mol2.  Only reads static sets of atoms.
	#  Ignores dynamic sets or sets of bonds or substructures.
	#. Atoms are placed in a record: mol->{ATOMSET}{<setname>}[<atoms>]
	#  and also $mol->{ATOMSET_BY_NUMBER}[$i]
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $fr = $_[1];
	
	my $count;	# Number of static sets
	my $numatoms;	# Number of atoms
	my $string;	# String to accumulate multi-line records
	my $setname;	# Set name
	
	SET: for (my $i = 0 ; $i < $mol->{NUMSETS}; ++$i) {
	
		my $line = silico_readline($fr);
		my @f = split(' ', $line);
	
		# Read static atom sets
		if ($f[1] eq 'STATIC' && $f[2] eq 'ATOMS') {
	
			$setname = $f[0];
			$string = '';
		
			# Read line(s) until line does not end in a backslash
			while (1) {
		
				$line = silico_readline($fr);
				if (!defined $line) {
					silico_msg('w', "Not enough sets found in SET record!\n",
							"Number of sets declared: $mol->{NUMSETS}, number actually found: $i\n");
					last;
				}
			
				# If line ends in a backslash then read another line
				if ($line =~ /\\\s*$/) {
					$line =~ s/\\//;
					$string .= $line;
					next;
				} else {
					$string .= $line;
					++$count;
					last;
				}
			}
		
			my @f = split(' ', $string);
			$numatoms = shift @f;

			my $atomset;
			my $j = 0;
			foreach my $atnum (@f) {
		
				push @{$atomset}, $mol->{ATOMS}[$atnum-1];
				++$j;
			}

			$mol->{ATOMSET}{$setname} = $atomset;
			$mol->{ATOMSET_BY_NUMBER}[$i] = $atomset;
		
			if ($numatoms != $j) {
				silico_msg('w', "Number of atoms specified in set ($j) does not match number listed ($numatoms)!\n");
			}
		
		} else {
	
			# Read line(s) until line does not end in a backslash
			while (1) {
		
				$line = silico_readline($fr);
				last if $line !~ /\\\s*/;
			}
		}
	}
	
	# Correct number of sets read in
	$mol->{NUMSETS} = keys (%{$mol->{ATOMSET}}) +1;
	
	return 1;
}

sub mol2_set_formal_charges {

	#<
	#? Assign atomic formal charges based on FORMAL_CHARGE sets
	#. FORMAL_CHARGE atomsets are then deleted
	#; Requires: molecule
	#>

	my $mol = $_[0];

	my $set;
	my $string;
	
	foreach $set (keys %{$mol->{ATOMSET}}) {

		if ($set =~ /^FORMAL_CHARGE_MINUS_([1-9])/) {

			my $charge = -$1;

			my @atoms = @{$mol->{ATOMSET}{$set}};

			foreach my $atom (@atoms) {

				$atom->{FORMAL_CHARGE} = $charge;
			}

			delete $mol->{ATOMSETS}{$set};
			next;
		}
		
		if ($set =~ /^FORMAL_CHARGE_PLUS_([1-9])/) {

			my $charge = $1;

			my @atoms = @{$mol->{ATOMSET}{$set}};

			foreach my $atom (@atoms) {

				$atom->{FORMAL_CHARGE} = $charge;
			}

			delete $mol->{ATOMSETS}{$set};
			next;
		}
	}
}

sub mol2_crystal_record {
	
	#<
	#? Read in a @<TRIPOS>CRYSIN record
	#. Called from read_mol2().
	#. Unit cell data is placed in hash keys within $mol.
	#; Requires: molecule
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $fr = $_[1];
	
	space_groups();
	
	my $line = silico_readline($fr);
	my @f = splitline($line);
	
	# First three values are cell sides
	$mol->{CELL}[0] = shift @f;
	$mol->{CELL}[1] = shift @f;
	$mol->{CELL}[2] = shift @f;
	
	# Values 4 to 6 are cell angles
	$mol->{CELL_ALPHA} = shift @f;
	$mol->{CELL_BETA} = shift @f;
	$mol->{CELL_GAMMA} = shift @f;
	
	# Values 7 and 8 define the space group:
	# - value 7 is the space group number
	# - value 8 is the axial orientation
	# These are joined together to form an internal Silico code.
	my $sg = shift @f;
	my $set = shift @f;
	
	$mol->{SPACE_GROUP} = (defined $sg ? $sg : '000');
	$mol->{SPACE_GROUP} .= (defined $set ? $set : '0');
	
	# Add extra zeros to pad out the space group so it matches entries in %Silico::Space_Groups.
	while (length $mol->{SPACE_GROUP} < 4) {
		$mol->{SPACE_GROUP} = "0".$mol->{SPACE_GROUP};
	}
		
	# Check to make sure a space group is actually matched.
	my $flag;
	foreach my $key (keys %Silico::Space_Groups) {
		if ($mol->{SPACE_GROUP} eq $key) {
			$flag = 1;
			silico_msg('g', "Space group: $sg    Axis setting: $set    Hermann-Mauguin symbol: $Silico::Space_Groups{$key}\n");
		}
	}
	
	# If no space group is matched, print a warning and set the space group to 1 and axis setting to 1.
	# This is done by using the internal Silico notation: 0011 (see silico_definitions for more info).
	if (!$flag) {
		
		my $sgm = (defined $sg ? "\"$sg\"" : "Undefined");
		my $asm = (defined $set ? "\"$set\"" : "Undefined");		
		
		silico_msg('w', "Invalid space group entry in CRYSIN record!\n",
				"Space group: $sgm    Axis setting: $asm\n",
				"Setting space group to 001 and axis setting to 1 (Hermann-Mauguin: P1).\n");
		$mol->{SPACE_GROUP} = "0011";
	}
	
	# If there's anything we don't recognise on the crystal record line, skip it.
	if (defined $f[0]) {
		silico_msg('w', "Unrecognised entries on CRYSIN record line!\n",
				"These have been ignored.\n");
	}
	
	return 1;
}

sub mol2_centroid_record {

	#<
	#? Read a tripos centroid feature record
	#. The centroid is stored in $mol->{FEATURES}
	#. Note:  This code assumes that the CENTROID record is located
	#  after the SET record.  This is normally true but not strictly
	#  required.
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $fr = $_[1];

	while (1) {

		my $line = silico_readline($fr);
		last if !defined $line;

		# We have come to the end of the features
		if ($line =~ /^@/) {
			silico_putline($line, $fr);
			last;
		}

		my @f = split ' ', $line;

		my $feature;
		$feature->{TYPE} = 'centroid';
		$feature->{NAME} = shift @f;
		$feature->{COMMENT} = join ' ', @f;

		# Centroid atom number and associated set number
		$line = silico_readline($fr);
		last if !defined $line;

		silico_msg('c', "line $line\n");
		@f = splitline($line);
		my $atomnum  = shift @f;

		my $atom = $mol->{ATOMS}[$atomnum-1];

		push @{$feature->{ATOMS}}, $atom;

		my $setnum  = shift @f;
		--$setnum;
		my $atomset = $mol->{ATOMSET_BY_NUMBER}[$setnum];

		$feature->{ATOMNUMBERS} = $atomset;

		push @{$mol->{FEATURES}}, $feature;
	}

	return 1;
}

sub extract_mol2_sdf_data {

	#<
	#? Read the sdf type data from mol2 file comment lines and put into molecule
	#. Removes this data from the comments array
	#; Requires: molecule, comments array
	#; Returns: comments lines with sdf data removed or undef on failure
	#>

	my $mol = $_[0];
	my $comments = $_[1];
	
	my $newcomments;
	my $sdf_data;
	
	@$newcomments = ();
	
	my $key;
	my $flag = 0;
	foreach my $line (@$comments) {
	
		if ($line=~ /^#SDF_DATA/) {
			$flag = 1;
			next;
		}
		if ($line =~ /^#END_SDF_DATA/) {
			$flag = 0;
			next;
		}
		
		# Save non sdf data
		if (!$flag) {
		
			push @$newcomments, $line;
			next;
		}

		if ($line =~ /^#KEY/) {

			my @f = split "'", $line;
			$key = $f[1];
			next;
		}

		if ($line =~ /^#DATA/ && defined $key) {
			$line =~ /'(.*)'/;
			push @$sdf_data, $1;
			next;
		}
		
		if ($line =~ /^#END_KEY/) {
		
			# Save data into SDF_DATA as either a scalar or a
			# hash depending on how many elements were read
			# This is not an ideal data structure and should
			# be fixed
			if ($#{$sdf_data} == 0) {
				$mol->{SDF_DATA}{$key} = $sdf_data->[0];
			} elsif ($#{$sdf_data} >= 0) {
				#print "key $key @$sdf_data\n";
				@{$mol->{SDF_DATA}{$key}} = @{$sdf_data};
			}
			undef $sdf_data;
			undef $key;
			next;
		}

		chomp $line;
		silico_msg('e', "Error reading SDF data!\n",
				"Line: \"$line\"\n");
		
		return undef;
	}

	return $newcomments;
}

sub atom_guess_element_mol2 {

	#<
	#? Guess the element of an atom from a mol2 file.
	#. Uses the Sybyl atom type to generate the atom element. 
	#; Requires: atom.
	#; Returns: Atomic symbol, atomic number.
	#>
	
	my $atom = $_[0];
	my $mol = $_[1];
	
	my @f = split /\./, $atom->{TYPE};
	
	defined $f[0] || return ('Du', 0);
	$f[0] =~ ucfirst (lc $f[0]);
	
	# Shortcuts!
	return ('H', 1) if $f[0] eq 'H';
	return ('C', 6) if $f[0] eq 'C';
	return ('N', 7) if $f[0] eq 'N';
	return ('O', 8) if $f[0] eq 'O';
	
	$f[0] = 'Du' if $f[0] eq 'Lp';
		
	# Test to see if we have a valid element.  Otherwise
	# use available info to guess it.
	# Note that VMD provides bad atom types in Mol2 files
	if (!defined $Silico::Atomic_Elements{$f[0]}) {
	
		# Guess atomic element using quiet and CHONSP options
		my ($el, undef) = atom_guess_element($atom, $mol, 1, 1);
		
		silico_msg('g', "Unknown atomic element '$f[0]'.  Guessing $el.\n");
		
		# Count number of errors
		++$mol->{WARN}{GUESSED_ATOM_ELEMENTS};
		
		return ($el, element_symbol2number($el));
	}
		
	return ($f[0], element_symbol2number($f[0]));
}


##################################################################
#
#	Write Mol2 routines
#
##################################################################

sub write_mol2 {

	#<
	#? Write an ensemble as a mol2 file.
	#. Carries out the following steps:
	#, Sort atoms so that substructures are written out in numerical order
	#, Generate bonds if there seem to be too few (less than (NUMATOMS -1)/2)
	#, Convert kekule representation of aromatic rings to aromatic type bonds
	#, Convert guanidinium and carboxylates to aromatic type bonds
	#, Generate @<TRIPOS>SUBSTRUCTURE record
	#, Generate Sybyl atom types
	#.
	#; Requires: ensemble (or molecule), filename, optional options string.
	#; Returns: undef if file open failed or other error, otherwise 1.
	#.
	#. Options:
	#, BOND - Force bond generation
	#, DIVIDE - Write each molecule to a separate file
	#, FAST - Equivalent to NOSORT, NOBOND, NOTYPE. This option can be used
	#  if the file read was already a good quality mol2 file
	#, FLAT - atom type molecule using routines that do not require good geometry
	#, MOLMOL - produce nonstandard mol2 file that is compatible with molmol
	#, NOSORT - do not sort atom order before writing
	#, NOTYPE - do not perform Sybyl atom typing
	#, NOBOND - do not create any bonds
	#, NOSDF_DATA - do not write SDF_DATA into comments
	#, NO_AROM_DELOC_BONDS - do not convert carboxylate and guanidinium bonds to aromatic
	#, PROTEIN - write in Sybyl protein format
	#, QUIET - do not print 'Writing' line
	#, SMALL - write in Sybyl small molecule format
	#.
	#. Note that by default write_mol2 will change the order of the atoms in the
	#  molecule data structure.  This may upset subsequent routines.
	#>
	
	my $mols = ens($_[0]);
	my $outfile = $_[1];
	my $options = uc ($_[2] || '');
	
	silico_msg('d', "Argument 0 is not an ensemble of molecules\n") if ref($_[0] ne 'ARRAY');
	silico_msg('d', "Incorrect arguments\n") if $_[3];
	
	silico_msg('t', "Starting ".subname()."\n");

	$outfile = get_ofilebase($mols).".mol2" if !$outfile;
	
        # Options
	if ($options =~ /\bFAST\b/) {
		add_option("NOSORT", $options);
		add_option("NOBOND", $options);
		add_option("NOTYPE", $options);
	}
	
	# Add standard options
	foreach (qw(append bond divide nosort nobond notype protein)) {
		add_option(uc $_, $options) if get_lflag($_);
	}

        add_option("QUIET", $options) if (quiet());
        remove_option("QUIET", $options) if $Silico::debug;
        my $quiet = 1 if ($options =~ /\bQUIET\b/);

	my $fr = open_molfile(undef, $outfile, undef, 'write', $options);
        return undef if !$fr;
        return undef if $fr->{ERROR};
	my $FH = $fr->{FH};
	
	silico_msg('c', "Writing mol2 file $outfile".($options ? " using $options" : '')."\n") if !$quiet;

	my $i = 0;
	foreach	my $mol (@$mols) {

		++$i;

		if (!defined $mol->{NUMATOMS}) {
			silico_msg('e', "\$mol->{NUMATOMS} is not defined!\n", "Aborting write.\n");
			return undef;
		}
		my $ret = write_mol2_molecule($fr, $mol, $options);
		
		silco_msg('e', "Molecule write_failed\n") if !defined $ret;
	}
	
	close_molfile($fr);
	
	return 1;
}

sub write_mol2_molecule {

	#<
	#? Write a single molecule record to a mol2 file
	#; Requires: file record, molecule, options
	#. See read_mol2 for general description
	#; Returns: 1 or undef if failed
	#>

	my $fr = $_[0]  || silico_msg('d', "File record not defined\n");
	my $mol = $_[1] || silico_msg('d', "Molecule not defined\n");
	my $options = $_[2] || '';

	my $FH = $fr->{FH};
	
	++$fr->{MOL_WRITE_COUNT};
	#fr_printout($fr);

	# Set SYBTYPE to protein if we came from a pdb file
	$mol->{SYBTYPE} = 'PROTEIN' if (defined $mol->{SOURCE_FILE_TYPE} && ($mol->{SOURCE_FILE_TYPE} eq 'pdb'));
	
	# options line overides other assignments
	$mol->{SYBTYPE} = 'PROTEIN' if $options =~ /\bPROTEIN\b/;
	$mol->{SYBTYPE} = 'SMALL' if $options =~ /\bSMALL\b/;
	$mol->{SYBTYPE} = 'SMALL' if $options =~ /\bMOLMOL\b/;
	$mol->{HAS_BAD_GEOM} = 1 if $options =~ /\bFLAT\b/;

	 # Add standard options 
        foreach (qw(append bond divide nosort nobond notype protein)) {
                add_option(uc $_, $options) if get_lflag($_);
        }

	# Make sure the substructure records are OK and the atoms are sorted
	sybyl_sort($mol) if !($options =~ /\bNOSORT\b/);

	# Renumber atoms (this is necessary for substructure cleanup)
	molecule_renumber($mol);
	
	molecule_check_and_fix_connectivity($mol, $options);

        make_atom_names_unique($mol) if get_flag('unique-atom-names', 'l');

	# Ensure that aromatic rings have aromatic bond types
	# This is required by Sybyl to function correctly
	make_aromatic_bonds($mol) if ($options !~ /\bNOBOND\b/);
	
	# Type atoms
	mol_type_sybyl($mol) if ($options !~ /\bNOTYPE\b/);
		
	# Create substructure records
	molecule_fix_sybyl_substructures($mol);

	# Write formal charges to atom sets
	make_atomset_formal_charge($mol);

	# Molecule Record
	#----------------
	write_mol2_comments($mol, $FH);
	write_mol2_sdfdata($mol, $FH) if (defined $mol->{SDF_DATA} && $options !~ /\bNOSDF_DATA\b/);
	write_mol2_atoms($mol, $FH, $options);
	write_mol2_bonds($mol, $FH, $options);
	write_mol2_substructures($mol, $FH, $options);
	write_mol2_set($mol, $FH) if $mol->{NUMSETS};
		
	if (	   defined $mol->{CELL}[0]
		|| defined $mol->{CELL}[1]
		|| defined $mol->{CELL}[2]
		|| defined $mol->{CELL_ALPHA}
		|| defined $mol->{CELL_BETA}
		|| defined $mol->{CELL_GAMMA}
		|| defined $mol->{SPACE_GROUP}
		|| defined $mol->{SPACE_GROUP_SETTING}) {
		
		write_mol2_crystal($mol, $FH);
	}
	
	return 1;
}

sub write_mol2_comments {

	#<
	#? Write mol2 comments
	#; Requires: molecule, filehandle
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $FH  = $_[1];
	
	if (defined $mol->{COMMENTS} && $mol->{COMMENTS} =~ /\w/) {
		if ($mol->{COMMENTS} !~ /^#/) {
			$mol->{COMMENTS} = '# '.$mol->{COMMENTS};
		}
		print $FH $mol->{COMMENTS};
	}

	return;
}

sub write_mol2_sdfdata {

	#<
	#? Write SDF datainto mol2 comment lines
	#; Requires: molecule, filehandle
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $FH  = $_[1];
	
	my $key;

	return if !defined $mol->{SDF_DATA};
	
	print $FH "#SDF_DATA\n";
	foreach $key (keys %{$mol->{SDF_DATA}}) {
				
		my $type;
		next if !defined $mol->{SDF_DATA}{$key};
		$type = ref($mol->{SDF_DATA}{$key}) || "SCALAR";
		
		print $FH  "#KEY '$key'\n";
		
		if ($type eq 'SCALAR') {
		
			my @f = split "\n", $mol->{SDF_DATA}{$key};
			foreach (@f) {
				print $FH "#DATA '$_'\n";
			}

		} elsif ($type eq 'HASH') {
		
			foreach (sort keys %{$mol->{SDF_DATA}{$key}}) {
				print $FH "#DATA '$_ $mol->{SDF_DATA}{$key}{$_}'\n";
			}

		} elsif ($type eq 'ARRAY') {
		
			foreach (@{$mol->{SDF_DATA}{$key}}) {
				print $FH "#DATA '$_'\n";
			}
		}
		
		else {
			silico_msg('w', "Don't know how to handle SDF DATA field of type $type!\n");
		}
		
		print $FH "#END_KEY\n";
	}
	
	print $FH "#END_SDF_DATA\n";

	return;
}

sub make_atomset_formal_charge {

	#<
	#? Make an ATOMSET record from atoms carrying a formal charge
	#; Requires: molecule
	#>

	my $mol = shift;
	
	my $charge;
	my $key;
	my $string;
	
	# Remove any existing sets
	foreach $key (keys %{$mol->{ATOMSET}}) {
		delete $mol->{ATOMSET}{$key} if $key =~ /^FORMAL_CHARGE_PLUS_/;
		delete $mol->{ATOMSET}{$key} if $key =~ /^FORMAL_CHARGE_MINUS_/;
	}

	foreach my $atom (atoms($mol)) {
	
		$charge = $atom->{FORMAL_CHARGE};
	
		if (defined $charge) {

			$charge =~ s/ //g;
			next if !$charge; # Empty
			next if $charge == 0; # Zero
		
			# Note that the character '-' can not be used and hence we must
			# use the word 'minus' for negative charges
			$string = "FORMAL_CHARGE_PLUS_".$charge if $charge > 0;
			$string = "FORMAL_CHARGE_MINUS_".abs($charge) if $charge < 0;
			
			push @{$mol->{ATOMSET}{$string}}, $atom;
		}
	}
}


sub write_mol2_atoms {

	#<
	#? Write out mol2 atom records
	#; Requires: molecule, filehandle, options;
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $FH  = $_[1];
	my $options = $_[2];
	
	my $atom;
	my @f;			# Temp string array
	my $i;
	my $nonamecount = 0;
	my $warn_atype = 0;
	my $warn_prime = 0;
	my $quiet = 1 if $options =~ /QUIET/;
	
	# Comments line.  REREG from is obtained from rdf format files files for DOCK
	$mol->{REREG} ||= '';
	$mol->{SYBCOMMENT} ||= '';
	
	$mol->{NUMBONDS} ||= 0;
	$mol->{NUMSETS} = keys %{$mol->{ATOMSET}} || 0;
	$mol->{SYBTYPE} = 'SMALL' if !$mol->{SYBTYPE};
	
	# Tripos Header
	print $FH "\@<TRIPOS>MOLECULE\n";
	print $FH ($mol->{NAME}||'Mol')."\n";
	print $FH "$mol->{NUMATOMS}\t$mol->{NUMBONDS}\t$mol->{NUMSUBSTRUCT}\t0\t$mol->{NUMSETS}\n";
	print $FH "$mol->{SYBTYPE}\n";
	print $FH ($mol->{SYBCHARGETYPE} || 'USER_CHARGES')."\n";
	print $FH "\n";
	print $FH "$mol->{SYBCOMMENT}$mol->{REREG}\n" || "\n";
	
	# Atom Records
	print $FH "\@<TRIPOS>ATOM\n";

	$i = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {

		++$i;
		
		$f[0] = $i+1;
		
		$f[1] = $atom->{NAME} || 'X';
		
		# Move any leading digits to end of atom name (Sybyl does not like these)
		if ($f[1] =~ /^\d/) {
			$f[1] =~ s/(^\d*)(.*)/$2$1/;
			++$warn_atype if $1;
		}
		
		# Replace '*' with 'S' and "'" with 'P' (Doesn't like these either!)
		if ($f[1] =~ /\*|\'/) {
			$f[1] =~ s/\*/S/g;
			$f[1] =~ s/\'/P/g;
			++$warn_prime;
		}

		# and remove spaces
		$f[1] =~ s/ //g;
		if ($f[1] eq "") {

			++$nonamecount;
			$f[1] = uc ($atom->{ELEMENT}) . $nonamecount;
		}
		
		$f[2] = $atom->{X} || 0;
		$f[3] = $atom->{Y} || 0;
		$f[4] = $atom->{Z} || 0;
		
		$f[5] = $atom->{TYPE} || $atom->{ELEMENT};
		
		$f[6] = 1;
		$f[6] = $atom->{SUBCOUNT} if defined $atom->{SUBCOUNT};
		
		# Append atom SUBID and EXT to substructure name
		if (defined $atom->{SUBNAME}) {
			$f[7] = $atom->{SUBNAME};
			# Remove digits and all characters after them
			# ie SER27A becomes SER
			$f[7] =~ s/\d.*//;
			$f[7] = "UNK" if $f[7] eq "";
			$f[7] .= $atom->{SUBID} if defined $atom->{SUBID};
			$f[7] =~ s/ //g;
		} else {
			$f[7] = 'SUB';
		}
		if (defined $atom->{EXT}) {
			$f[7] .=$atom->{EXT};
		}
		
		$f[8] = $atom->{CHARGE} || 0;

		$f[9] = $atom->{SYBYL_STATUS} || '';
		
		printf $FH " %6d %4s %12.6f %12.6f %12.6f %5s %5s %5s %8.5f %s\n", @f;
	}
		
	if (!$quiet) {
		if ($warn_atype == 1) {
			silico_msg('n', "Renamed one atom with a name containing leading digits.\n");
		} elsif ($warn_atype > 1) {
			silico_msg('n', "Renamed $warn_atype atoms with names containing leading digits.\n");
		}
	
		if ($warn_prime == 1) {
			silico_msg('n', "Renamed one atom with a name containing \"'\" and/or \"*\" characters.\n");
		} elsif ($warn_prime > 1) {
			silico_msg('n', "Renamed $warn_prime atoms with names containing \"'\" and/or \"*\" characters.\n");
		}
	}

	return;
}

sub write_mol2_bonds {

	#<
	#? Write a mol2 bond records
	#. Note.  This routine converts carboxylate and guanidinium groups to
	#  aromatic bond order (in the output file)
	#; Requires: molecule, filehandle, options
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $FH  = $_[1];
	my $options = $_[2];
	
	my $a1;
	my $a2;
	my $atom1;
	my $atom2;
	my $bo;
	my $bondorder;
	my $connum;
	my $count;
	my $i;
	my $j;
	my $status;
	my @types;
	
	# Bond Records
	# 1 single, 2 double, 3 triple, 4 aromatic,
	# 5 single or double, 6 single or aromatic,
	# 7 double or aromatic,
	# 8 any (undefined), 0 zero order
	
	print $FH "\@<TRIPOS>BOND\n";

	$count = 0;
	$i = -1;
	foreach (@{$mol->{ATOMS}}) {
		
		++$i;
		$atom1 = $mol->{ATOMS}[$i];
		
		$j =-1;
		foreach $connum (@{$atom1->{CONNECT}}) {
		
			++$j;
			
			if (!defined $connum) {
				++$mol->{WARN}{"Attempting to create bond to missing atom"};				
				next;
			}

			$atom2 = $mol->{ATOMS}[$connum];
			$a1 = $connum+1;
			$a2 = $i+1;
			
			# Each bond is represented twice. Print each bond once by
			# only printing if atom1 < atom2
			next if ($a1 >= $a2);
			
			++$count;
				
			$bo = ${$atom1->{BORDERS}}[$j];
				
			@types = ($mol->{ATOMS}[$i]{TYPE} || '', $mol->{ATOMS}[$connum]{TYPE} || '');
			@types = sort @types;
				
			# No bond order defined
			if (!defined $bo) {
				silico_msg('w', "Undefined bond order between atoms $a1 and $a2!\n",
						"Setting bond order to Single.\n");
				$bo = 1;
			}
				
			# Sybilify the bond orders
			$bondorder = $bo;
			$bondorder = "ar" if ($bo == 4);
			$bondorder = "du" if ($bo == 0);
			$bondorder = 1 if ($bo == 5 || $bo == 6 || $bo == 8);
			$bondorder = 2 if ($bo == 7);
			
			
			if ($options !~ /\bNO_AROM_DELOC_BONDS\b/) {
				# Sybyl wants aromatic bonds for delocalised bonds in
				# carboxylate and guanidinium groups
				if (($types[0] eq 'C.2' && $types[1] eq 'O.co2') || ($types[0] eq  'C.cat' && $types[1] eq 'N.pl3')){
				
					$bondorder = 'ar';
				}
			}
			
			$status = '';
			$status = 'BACKBONE' if ($atom1->{SYBYL_STATUS} && $atom2->{SYBYL_STATUS}
				&& $atom1->{SYBYL_STATUS} =~ 'BACKBONE' && $atom1->{SYBYL_STATUS} =~ 'BACKBONE');
				
			printf $FH "%6d %5d %5d %1s %s\n",$count,$a1,$a2,$bondorder, $status;
		}
	}
}

sub write_mol2_substructures {

	#<
	#? Write a mol2 substructure record
	#; Requires: molecule, filehandle
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $FH  = $_[1];
	my $options = $_[2];
	
	my $chain;
	my $i;
	my $r;
	
	# Definition of SUBSTRUCTURE record
	#subst_id subst_name root_atom [subst_type [dict_type [chain [sub_type [inter_bonds [status [comment]]]]]]]
	
	print $FH "\@<TRIPOS>SUBSTRUCTURE\n";

	for ($i=0; $i < $mol->{NUMSUBSTRUCT}; ++$i) {
		
		print $FH $i+1 ."\t$mol->{SUB_SUBNAME}[$i]\t";
		
		if (defined $mol->{SUB_ROOTATOM}[$i]) {
			print $FH "$mol->{SUB_ROOTATOM}[$i]";
		} else {
			silico_msg('w', "List of substructure ROOT ATOMs contains no entry for this residue!\n",
					"Using 1 instead.\n");
			print $FH "1";
		}

		# PROTEIN format files
		
		if ($options =~ /\bMOLMOL\b/) {
			
			# Write specific format for MolMol which is broken
			# The following line is nessesary for MolMol - but breaks Sybyl(!)
			print $FH " **** 0 **** ****";
			
		} elsif ($mol->{SYBTYPE} eq "BIOPOLYMER" || $mol->{SYBTYPE} eq "PROTEIN") {
			
			# Protein format
			$r = $mol->{SUB_ROOTATOM}[$i] || 1; # r is the atom number of the root atom of the substructure
			$chain = $mol->{ATOMS}[$r-1]{CHAIN} || 'Z';
			$chain =~ s/ //g if defined $chain;
			
			print $FH "\t$mol->{SUB_SUBTYPE}[$i]";
			print $FH "\t".($mol->{SUB_DICT}[$i] || 0);
			print $FH "\t".($chain || 'Z');
			print $FH "\t".($mol->{SUB_INTERBONDS}[$i]);
			print $FH "\t".($mol->{SUB_REALSUBNAME}[$i]);
			print $FH "\t".($mol->{SUB_STATUS}[$i]);

			# Comment

		} else  {
		
			# Small molecule format
		}
		
		print $FH  "\n";
	}

	return;
}

sub write_mol2_set {

	#<
	#? Write mol2 static sets
	#; Requires: molecule, filehandle
	#>
	
	my $mol = $_[0];
	my $FH  = $_[1];
	
	my $res;
	my $setname;
	my $setatomcount;
	
	return if !$mol->{NUMSETS};
	
	print $FH "\@<TRIPOS>SET\n";
	
	# Fix atom->{NUM};
	my $i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
		++$i;
		$atom->{NUM} = $i;
	}
	
	foreach $setname (sort keys %{$mol->{ATOMSET}}) {
		
		# Take care of any deleted atoms
		$setatomcount = 0;
		foreach my $atom (@{$mol->{ATOMSET}{$setname}}) {
			next if !defined $atom;
			++$setatomcount;
		}
			
		print $FH "$setname STATIC ATOMS <user> **** \"\"\n";
		print $FH "$setatomcount  ";
			
		$i = 0;
		foreach my $atom (@{$mol->{ATOMSET}{$setname}}) {
			
			next if !defined $atom; # Skip if it has been deleted
			
			++$i;	
			if ($i == 10) {
				print $FH "\\\n";
				$i = 0;
			}
			
			print $FH $atom->{NUM}." ";
		}
			
		print $FH "\n";
	}
	
	#
	# Substructure sets
	#
	foreach $setname (sort keys %{$mol->{SUBSET}}) {
		
		print $FH "$setname STATIC SUBSTS <user> **** \"\"\n";
		print $FH ($#{$mol->{SUBSET}{$setname}}+1)."  ";
		
		$i = 0;	
		foreach $res (@{$mol->{SUBSET}{$setname}}) {
		
			++$i;	
			if ($i == 10) {
				print $FH "\\\n";
				$i = 0;
			}
			
			# Write atomnumber of first atom in set
			print $FH ($res->[0]{NUM} || 0). " ";
		}
			
		print $FH "\n";
	}
	return;
}

sub write_mol2_crystal {
	
	#<
	#? Write out the CRYSTAL record to a mol2 file
	#; Requires: Molecule, filehandle
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $FH  = $_[1];
	
	my @f;
	my $i;
	
	$f[0] = $mol->{CELL}[0];
	$f[1] = $mol->{CELL}[1];
	$f[2] = $mol->{CELL}[2];
	$f[3] = $mol->{CELL_ALPHA};
	$f[4] = $mol->{CELL_BETA};
	$f[5] = $mol->{CELL_GAMMA};
	# Make sure no values are undefined
	for ($i = 0; $i <= 5; ++$i) {
		$f[$i] = 0 if (!defined $f[$i]);
	}
	
	if (defined $mol->{SPACE_GROUP} && length $mol->{SPACE_GROUP} >= 4) {
		$f[6] = substr($mol->{SPACE_GROUP}, 0, 3);
		$f[7] = substr($mol->{SPACE_GROUP}, 3, 1);
		$f[6] =~ s/^0+//;
	} else {
		# If there is no space group, we set it to the default,
		# that is, 0011 (which is P1 by Hermann-Mauguin).
		$f[6] = '1';
		$f[7] = '1';
	}
	
	print $FH "\@<TRIPOS>CRYSIN\n";
	
	printf $FH "%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3s %1s\n", @f;
	
	return;
}


sub molecule_fix_sybyl_substructures {
	
	#<
	#? Creates records for sybyl substructures. Also generates
	#  BACKBONE status flags if amino acids are identified.
	#; Requires: molecule.
	#; Returns: nothing.
	#>
	
	my $mol = $_[0];
	
	my $aa;
	my $oldsegid = 'xxxxxx';
	my $oldchain = 'xxxxxx';
	my $residues;
	my $residue;
	my $subcount = 0;

	$residues = molecule_get_residues($mol);

	foreach $residue(@$residues) {
		
		my $realsubname;
		my $resatom;
		my $reference_atom;
		my $subname;
		
		# Set a reference atom: The first atom in the residue.
		$reference_atom = $residue->[0];
		
		# The text-only residue name: Whatever the reference atom's
		# residue is. Or "SUB".
		$realsubname = $reference_atom->{SUBNAME} || "SUB";
		$realsubname =~ s/ //g;
		$realsubname ||= 'SUB';
		$mol->{SUB_REALSUBNAME}[$subcount] = $realsubname;
		
		# The full residue name is a concatenation of the text-only name,
		# the substructure ID (or residue count if the ID is unknown) and the
		# extension. Spaces are removed.
		$subname = $realsubname;
		if (defined $reference_atom->{SUBID}) {
			$subname .= $reference_atom->{SUBID}
		} else {
			$subname .= $subcount;
		}
		$subname .= $reference_atom->{EXT} if defined $reference_atom->{EXT};
		$subname =~ s/ //g;
		$mol->{SUB_SUBNAME}[$subcount] = $subname;

		# If the substructure is an amino or nucleic acid,
		# apply various dictionary values.
		if ($Silico::Amino_Nucleic_Acids{$realsubname}) {
			
			# At least one amino acid has been found.
			$aa = 1;

			# This residue appears in a dictionary.
			$mol->{SUB_DICT}[$subcount] = 1;

			# If this residue contains an alpha-carbon (name "CA"),
			# use this atom as the "root atom".
			# Also set SYBYL_STATUS flag to BACKBONE
			foreach $resatom (@$residue) {
				
				my $aname;
				
				# Interrogate the atom name. Start by removing spaces.
				$aname = $resatom->{NAME};
				$aname =~ s/\s*//g;
				
				# Is the present atom the alpha carbon?
				# If so, this becomes the root atom.
				# We can use $atom->{NUM} here because ensemble_renumber was run
				# when write_mol2 commenced.
				$mol->{SUB_ROOTATOM}[$subcount] = $resatom->{NUM} if ($aname eq 'CA');
				
				# If the present atom is the carbonyl carbon, alpha carbon, amide nitrogen
				# or amide oxygen, it is a backbone atom.
				if ($aname eq 'C' || $aname eq 'CA' || $aname eq 'N' || $aname eq 'O') {
					$resatom->{SYBYL_STATUS} = 'BACKBONE';
				}
			}

			# If we have changed chain or segment...
			if ((defined $reference_atom->{CHAIN} && $reference_atom->{CHAIN} ne $oldchain) ||
				 (defined $reference_atom->{SEGID} && $reference_atom->{SEGID} ne $oldsegid)) {

				# This residue has "ROOT" status, it has one bond to it (instead of
				# the default two bonds), and if it is not the first residue in the
				# entire molecule, the previous residue (being therefore a chain end)
				# only has one bond to it as well.
				$mol->{SUB_STATUS}[$subcount] = 'ROOT';
				$mol->{SUB_INTERBONDS}[$subcount] = 1;
				$mol->{SUB_INTERBONDS}[$subcount-1] = 1 if $subcount > 0;
			
			# Otherwise...
			} else {
			
				# This residue has no particular status, and it has two bonds to it
				# (the default protein number, one each side).
				# Note that if it is an end of a chain or segment, this is no cause for
				# concern: the next residue will start a new chain or segment, and the
				# code above will reduce the SUB_INTERBONDS of this residue to 1.
				$mol->{SUB_INTERBONDS}[$subcount] = 2;
				$mol->{SUB_STATUS}[$subcount] = '';
			}
			
			# Whether the chain or segment has been changed or not,
			# this residue is of type RESIDUE (since it is in a protein).
			# Non-protein residues are of type GROUP (see also below).
			$mol->{SUB_SUBTYPE}[$subcount] = "RESIDUE";
			
			
		# In all other cases (not amino or nucleic acids)...
		} else {
			
			# This residue is of type GROUP (since it is not in a protein).
			# It also has no dictionary type (SUB_DICT = 0).
			# It also is said to only have one bond to neighbouring residues.
			# (This, however, is not particularly important - all actual bonds
			# are specified in the BONDS section).
			$mol->{SUB_SUBTYPE}[$subcount] = "GROUP";
			$mol->{SUB_DICT}[$subcount] = 0;
			$mol->{SUB_INTERBONDS}[$subcount] = 1;

			# If this residue is the first in the entire molecule, it has "ROOT"
			# status. Otherwise, it has no particular status, and its predecessor
			# has only one residue attached to it. (Presumably this is to address
			# ends of protein chains that aren't followed in the file by another
			# amino acid.)
			if ($subcount == 0) {
				$mol->{SUB_STATUS}[$subcount] = 'ROOT';
			} else {
				$mol->{SUB_STATUS}[$subcount] = '';
				$mol->{SUB_INTERBONDS}[$subcount-1] = 1;
			}
		}
		
		# The root atom may have been set to the alpha carbon, but this will
		# only happen if the residue was a correct amino acid. So we need to
		# ensure that the root atom is actually set. If it is not, we use the
		# reference atom as the root atom.
		# We can use $atom->{NUM} here because ensemble_renumber was run
		# when commencing write_mol2.
		if (!defined $mol->{SUB_ROOTATOM}[$subcount]) {
			$mol->{SUB_ROOTATOM}[$subcount] = $reference_atom->{NUM};
		}
		
		# Each atom in the residue belongs to substructure number ($subcount+1),
		# since $subcount is indexed from zero.
		# Note that $subcount + 1 is not necessarily the same as Sub ID.
		foreach $resatom (@$residue) {
			$resatom->{SUBCOUNT} = $subcount + 1;
		}
		
		# In any case, use the current reference atom's chain and SegID details
		# as the previous ones, so we know when we've found a new chain or segment.
		$oldchain = $reference_atom->{CHAIN} || '';
		$oldsegid = $reference_atom->{SEGID} || '';
		
		# Finally, increment the substructure count.
		++$subcount;
	}

	# If any amino acids have been found, we have a protein.
	$mol->{SYBTYPE} = 'PROTEIN' if $aa;
	
	# Set the total number of substructures.
	$mol->{NUMSUBSTRUCT} = $subcount;
	
	# Leave the subroutine.
	return;
}

sub sybyl_sort {

	#<
	#?
	#  Ensure that molecule is written out with substructures
	#  in numerical order.  This can cause problems with other
	#  software packages if it is not the case.  Sybyl doesn't
	#  mind this though
	#; Requires: ensemble.
	#; Returns: nothing.
	#>

	my $mol = $_[0];
	
	my $i = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {

		++$i;
		next if $i == 0;
		
		my $prevatom = $mol->{ATOMS}[$i-1];
		
		if (defined $atom->{SUBID} && defined $prevatom->{SUBID} && $atom->{SUBID} < $prevatom->{SUBID}) {

			# OK if new chain or segid
			next if (defined $atom->{SEGID} && defined $prevatom->{SEGID} && $atom->{SEGID} ne $prevatom->{SEGID});
			next if (defined $atom->{CHAIN} && defined $prevatom->{CHAIN} && $atom->{CHAIN} ne $prevatom->{CHAIN});
							
			pdb_molecule_sort_atoms($mol);
			$mol->{WARN}{MOLECULE_SORT} = 1;
			last;
		}
	}

	return;
}
	
sub mol_type_sybyl {

	#<
	#? Generate Sybyl atom types for a molecule
	#. Types are put in to $atom->{TYPE}.
	#. This routine does not assume that hydrogens are present.
	#. Calls molecule_find_rings to rings and planar rings
	#. Identifies carboxylate, amide and  guanidinium carbon atoms and
	#  marks them with the flags $atom->{FG}{CARBOXYLATE_C},  $atom->{FG}{AMIDE_C}
	#  and $atom->{FG}{GUANIDINIUM_C} respectively.
	#. If hydrogens are not present then this routine uses simple rules
	#  to try to put molecules in physiologialy relevant ionisation state.
	#  Makes primary amines N.4 and makes O-C=O
	#  into carboxylate groups.
	#. This routine makes Ti,Cr and Co octahedral.
	#; Requires: atom, molecule.
	#; Returns: -1 if there was an atom typing error.
	#>
	
	my $mol = $_[0];
	
	my $i;
	my $con;
	my $connum;
	my $error = 0;
	my $bo;
	my $numsingle;
	my $numdouble;
	my $numtriple;
	my $numaromatic;
	my $numbonds;
	my $aromatic_6;
	my $con_amide;
	my $con_carboxylate;
	my $con_guanidinium;
	my $con_planar;
	my $con_o;
	my $con_h;
	my $con_s;
	my $el;
	my $ring;
	
	return if $mol->{HAS_NO_COORDS};

	# Label amide, carboxylate and guanidinium carbons
	# This routine calls find_rings
	mol_label_functional_group($mol);
	
	ATOM: foreach my $atom (@{$mol->{ATOMS}}) {
	
		$numsingle = 0;
		$numdouble = 0;
		$numtriple = 0;
		$numaromatic = 0;
		$numbonds = 0;
		$aromatic_6 = 0;
		$con_amide = 0;
		$con_carboxylate = 0;
		$con_guanidinium = 0;
		$con_planar = 0;
		$con_o = 0;
		$con_h = 0;
		$con_s = 0;
	
		$el = $atom->{ELEMENT};
				
		if (!$el) {
			silico_msg('w', "Atom number $atom->{NUM} has no defined element!\n",
					"Setting type to DUMMY.\n");
			$atom->{TYPE} = 'Du';
			next;
		}

		# Make Co, Cr and Ti octahedral
		if ($el eq 'Co' || $el eq 'Cr' || $el eq 'Ti') {
			$atom->{TYPE} =  $el.".oh";
			next;
		}
		
		# Most elements
		if (!($el eq 'C' || $el eq 'O' || $el eq 'N'|| $el eq 'S' || $el eq 'P')) {
			$atom->{TYPE} = $el;
			next;
		}

		# Is the current atom part of a flat, 6-membered ring
		foreach $ring (@{$atom->{PRINGS}}) {
			++$aromatic_6 if $#{$ring} == 5;
		}
		
		# Get types of bonds and connected atoms
		$i = 0;
		foreach $connum (@{$atom->{CONNECT}}) {
	
			# Skip if zero order or no bonds defined
			next if (!defined $bo);
			next if !defined $connum;

			$bo = $atom->{BORDERS}[$i];
			next if ($bo == 0);
		
			$con = $mol->{ATOMS}[$connum];

			if ($con->{FG}{AMIDE_C}) {
				$con_amide = 1;
			}
			if ($con->{FG}{CARBOXYLATE_C}) {
				$con_carboxylate = 1;
			}
			if ($con->{ELEMENT} eq 'O') {
				++$con_o;
			}
			if ($con->{ELEMENT} eq 'H') {
				++$con_h;
			}
			if ($con->{ELEMENT} eq 'S') {
				++$con_s;
			}
			if ($con->{PLANAR_ATOM}) {
				++$con_planar;
			}
			if ($con->{FG}{GUANIDINIUM_C}) {
				++$con_guanidinium;
			}
		
			++$numsingle if ($bo == 1);
			++$numdouble if ($bo == 2);
			++$numtriple if ($bo == 3);
			++$numaromatic if ($bo == 4);
			++$numbonds;
			++$i;
		}
		
		# Carbon
		if ($el eq 'C') {
		
			if ($atom->{FG}{GUANIDINIUM_C}) {
				$atom->{TYPE} = 'C.cat';
				next;
			}
			if ($numaromatic ||$aromatic_6) {
				$atom->{TYPE} = 'C.ar';
				next;
			}
			if ($numdouble == 0 && $numtriple == 0) {
				$atom->{TYPE} = 'C.3';
				next;
			}
			if ($numdouble >= 1) {
				$atom->{TYPE} = 'C.2';
				next;
			}
			if ($numtriple == 1) {
				$atom->{TYPE} = 'C.1';
				next;
			}
			
			print_atom_type_error($atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{TYPE} = 'C.unk';
			$error = 1;
			next;
		}
		
		# Nitrogen
		if ($el eq 'N') {
			if ($numaromatic || $aromatic_6) {
				$atom->{TYPE} = 'N.ar';
				next;
			}
			if ($con_amide) {
				$atom->{TYPE} = 'N.am';
				next;
			}
			if ($numsingle == 4) {
				$atom->{TYPE} = 'N.4';
				next;
			}
			if ($con_guanidinium) {
				$atom->{TYPE} = 'N.pl3';
				next;
			}
			# Sulfonamides
			if ($con_s) {
				$atom->{TYPE} = 'N.pl3';
				next;
			}
			if ($numdouble == 0 && $numtriple == 0) {
				if ($atom->{PLANAR_ATOM}) {
					$atom->{TYPE} = 'N.pl3';
					next;
				}
				if ($atom->{PLANAR_RING}) {
					$atom->{TYPE} = 'N.pl3';
					next;
				}
				if ($con_planar) {
				$atom->{TYPE} = 'N.pl3';
					next;
				}
				# Promote nitrogens without attached hydrogens to N.4
				if (!$con_h) {
					$atom->{TYPE} = 'N.4';
					next;
				}
				$atom->{TYPE} = 'N.3';
				next;
			}
			
			if ($numdouble >= 1) {
				$atom->{TYPE} = 'N.2';
				next;
			}
			if ($numtriple == 1) {
				$atom->{TYPE} = 'N.1';
				next;
			}
			
			print_atom_type_error($atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{TYPE} = 'N.unk';
			$error = 1;
			next;
		}
		
		# Oxygen
		if ($el eq 'O') {
		
			# Return O.co2 for O attached to S or P
			foreach $connum (@{$atom->{CONNECT}}) {
				$con = $mol->{ATOMS}[$connum];
				if ($con->{ELEMENT} eq 'S' || $con->{ELEMENT} eq 'P') {
					$atom->{TYPE} = 'O.co2';
					next ATOM;
				}
			}
		
			if ($con_carboxylate) {
				$atom->{TYPE} = 'O.co2';
				next;
			}
				
			if ($numdouble == 1) {
				$atom->{TYPE} = 'O.2';
				next;
			}
			if ($numdouble == 0) {
				$atom->{TYPE} = 'O.3';
				next;
			}
			
			print_atom_type_error($atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{TYPE} = 'O.unk';
			$error = 1;
			next;
		}
		
		# Sulfur
		if ($el eq 'S') {
			if ($con_o == 1) {
				$atom->{TYPE} = 'S.O';
				next;
			}
			if ($con_o == 2) {
				$atom->{TYPE} = 'S.O2';
				next;
			}
			if ($con_o > 2) {
				$atom->{TYPE} = 'S.3';
				next;
			}
			if ($numdouble >= 1) {
				$atom->{TYPE} = 'S.2';
				next;
			}
			if ($numdouble == 0) {
				$atom->{TYPE} = 'S.3';
				next;
			}
			
			print_atom_type_error($atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{TYPE} = 'S.unk';
			$error = 1;
			next;
		}
		
		# Phosphorus
		if ($el eq 'P') {
			$atom->{TYPE} =  'P.3';
			next;
		}
	}
	
	return -1 if $error;
	return 0;
}


##################################################################
#
#	Sets
#
##################################################################

sub make_atomset {

	#<
	#? Make an ATOMSET record from atoms carrying a particular flag
	#. Flag must be defined and be set to true (ie not zero or undefined).  
	#  Atomsets are written out in mol2 format as static sets
	#; Requires: molecule, set name, flag1, flag2, ...
	#; Returns: count of atoms having that flag
	#>
	
	my $mol = shift;
	my $setname = shift;
	
	my $flag;
	my @flags = @_;
	my $count = 0;
	my @list;
	my %present;
	
	foreach my $atom (atoms($mol)) {
		foreach $flag (@flags) {
	
			next if !$atom->{$flag};
			++$count;
			push @list, $atom;
			next;
		}
	}
	
	# Remove duplicates
	foreach (@list) {
		next if ($present{$_});
		
		$present{$_} = 1;
		push @{$mol->{ATOMSET}{$setname}}, $_;
	}
	
	# Delete the atomset if at this point it is empty
	if ($count == 0) {
		undef $mol->{ATOMSET}{$setname};
	}
		
	return $count;
}

sub make_atomset_fg {

	#<
	#? Make an ATOMSET record from atoms carrying a particular functional group flag
	#. Flag must be defined and be set to true (ie not zero or undefined).  The set name
	# '*' will make an atomset of all functional groups 
	#  Atomsets are written out in mol2 format as static sets
	#; Requires: molecule, set name, flag1, flag2, ...
	#; Returns: count of atoms having that flag
	#>
	
	my $mol = shift;
	my $setname = shift;
	
	my $count = 0;
	my $flag;
	my @flags = @_;
	my $key;

	if ($setname eq '*') {
	
		foreach my $atom (atoms($mol)) {
	
			foreach $key (keys %{$atom->{FG}} ) {
				++$count;
				push @{$mol->{ATOMSET}{$key}}, $atom;
			}
		}
		
	} else {
	
		foreach my $atom (atoms($mol)) {
		
			next if !defined $atom->{FG};
		
			foreach $flag (@flags) {
				
				next if !$atom->{FG}{$flag};
				++$count;
				push @{$mol->{ATOMSET}{$setname}}, $atom;
			}
		}
	}
	
	return $count;
}

sub atom_colour {

	#<
	#? Colour an atom. 
	#. Compatible with Mol2 files
	#; Requres: atom, colour
	#>
	
	my $atom = $_[0];
	my $colour = uc $_[1];
	
	my $sybcolours = " WHITE RED MAGENTA YELLOW ORANGE BLUE VIOLET GREEN BLUEGREEN REDORANGE CYAN PURPLE ";
	
	silico_msg('w', "Colour $colour does not match sybyl colours\n") if $sybcolours !~ /\b$colour\b/;
	
	$atom->{COLOUR} = uc $colour;
}

sub make_atomset_colour {

	#<
	#? Make an ATOMSET record from atoms carrying an atom->{COLOUR} flag.  This produces
	#  atomsets compatible with Sybyl
	#. Atomsets have names of the type 'ATOM$RED'.  Sybyl colours are WHITE RED MAGENTA YELLOW ORANGE 
	#  BLUE VIOLET GREEN BLUEGREEN REDORANGE CYAN PURPLE (In no particular order).
	#; Requires: molecule
	#; Returns: count of atoms having that flag
	#>

	my $mol = shift;
	
	my $atomsets;
	my $count = 0;
	my $colour;
	my $i;
	my $sybcolours = " WHITE RED MAGENTA YELLOW ORANGE BLUE VIOLET GREEN BLUEGREEN REDORANGE CYAN PURPLE ";
	my $setname;
	my $warn = 0;
	my $warnhash;

	# Create atomset hash if required
	$mol->{ATOMSET} = {} if (!defined $mol->{ATOMSET});

	$atomsets = $mol->{ATOMSET};

	$i = -1;
	foreach my $atom (atoms($mol)) {

		++$i;
		
		$colour = $atom->{COLOUR} if $atom->{COLOUR};

		next if !$colour;

		if ($sybcolours !~ / $colour /) {
			++$warn;
			++$warnhash->{$colour};
		}

		$setname = 'ATOM$'.$colour;

		push @{$atomsets->{$setname}}, $atom;

		++$count;
	}

	if ($warn) {
		silico_msg('w', "Found $warn atoms with nonstandard colour names: ".(keys %$warnhash)."\n");
	}

	return $count;
}

sub make_subset {

	#<
	#? Make a SUBSET (substructure set) record from atoms carrying a particular flag
	#. Flag must be defined and be set to true (ie not zero or undefined).  
	#  Atomsets are written out in mol2 format as static sets
	#; Requires: molecule, set name, flag1, flag2, ...
	#; Returns: count of substructures having that flag
	#>
	
	my $mol = shift;
	my $setname = shift;
	
	my $count = 0;
	my @flags = @_;
	my @list;
	my %present;
	
	# Create $atom->{SUBCOUNT}
	my $residues = molecule_get_residues($mol);
	
	foreach my $atom (atoms($mol)) {
	
		my $subid = $atom->{SUBCOUNT}-1;
		
		foreach my $flag (@flags) {
	
			next if !$atom->{$flag};
			
			next if $present{$subid};
			$present{$subid} = 1 ;
			
			# Save all atoms in residue to substructure set
			push @{$mol->{SUBSET}{$setname}}, $residues->[$subid];
			++$count;
		}
	}
	
	# Delete the subset if at this point it is empty
	if ($count == 0) {
		delete $mol->{SUBSET}{$setname};
	}
		
	return $count;
}


sub delete_atomset {
	
	#<
	#? Delete an atomset on the basis of its label
	#; Requires: Molecule, atomset label
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $label = $_[1];
	
	# For each atom set that has been defined in the molecule...
	foreach my $atomset (keys %{$mol->{ATOMSET}}) {			
		
		my $satom;
		my $flag;
					
		# ... if the current group's label is the atom set label...
		if ($label eq $atomset) {
		
			# ...for each atom in the set...
			foreach $satom (@{$mol->{ATOMSET}{$atomset}}) {
			
				# skip this atom if it is part of the group
				next if (defined $satom->{"$label"}
					&& $satom->{"$label"} == 1);
				
				# otherwise, set a flag to print out the note,
				# and since we need to print a note, we don't need
				# to go through any more atoms
				$flag = 1;
				last;
			}
			
			# If the flag to print is set...
			if ($flag) {
				# print a note
				silico_msg('n', "Found atomset with label $atomset. Deleting.\n");
			}
			
			# In either case, delete the atom set
			delete $mol->{ATOMSET}{$atomset};
		}
	}
}

sub sybyl_spectrum {

	#< 
	#? Return a Sybyl colour based on a value and a maximum range
	#. Assumes that the value is between 0 and $range.  Returns
	#  one of a spectrum of 10 colours (RED REDORANGE ORANGE YELLOW GREEN GREENBLUE CYAN BLUE PURPLE VIOLET)
	#; Reqires: value, minimum value of colour range (default 0), maximum value of colour range (default 100)
	#; Returns: colour (string)
	#>

	my $value = $_[0];
	my $minval = $_[1] || 0;
	my $maxval = $_[2] || 100;  

	my $bin;
	my @spectrum;
	my $range = $maxval - $minval;

	# Avoid divide by zero
	return 'WHITE' if !$range;

	$bin = int(($value-$minval)/$range*10);
	
	# Squash the value '10' into the last bin
	return 'VIOLET' if $bin >= 10;

	# Out of range
	if ($bin < 0 || $bin > 10) {
		return 'WHITE';
	}

	@spectrum = qw(RED REDORANGE ORANGE YELLOW GREEN GREENBLUE CYAN BLUE PURPLE VIOLET);
	
	return $spectrum[$bin];
}
return 1;
