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
#! silico_check.pm
#? Silico routines to check the quality of molecules
#. $Revision: 1.34.2.1.2.3 $
#>

use strict;
package Silico;

##################################################################
#
#	Molecule checking routines
#
##################################################################

sub ensemble_check {
	
	#<
	#? Check all molecules in an ensemble.
	#. Calls 'mol_check'.
	#; Requires: ensemble.
	#>
	
	my $molecules = $_[0];
	
	my $i;
	my $nummols = $#{$molecules}+1;
	
	for ($i = 1; $i <= ($#{$molecules}+1); ++$i) {
		my $mol = $molecules->[$i];
		silico_msg('q', "Checking molecule $i of $nummols\n");
		mol_check($mol);
	}
}

sub mol_check_all {

	#<
	#? Call all the mol_check routines
	#; Requires: molecule.
	#. Returns; 0 if everything is OK. Otherwise returns the minor error
	#  count and returns -1 if there is a major problem
	#>

	my $mol = $_[0];

	my $error = 0;
	my $major_error = 0;
	my $val;

	$val = mol_check($mol);
	$major_error = 1 if $val == -1;
	$error += $val;
	
	$val= mol_check_atom_overlap($mol);
	$major_error = 1 if $val == -1;
	$error += $val;
	
	$val = mol_check_valences($mol);
	$major_error = 1 if $val == -1;
	$error += $val;

	$val = mol_check_bondorders($mol);
	$major_error = 1 if $val == -1;
	$error += $val;
	
	return -1 if $major_error;
	return $error;
}

sub mol_check {

	#<
	#? Check a molecule for self-consistency
	#. This routine checks for problems with the internal
	#  representation of the molecule.
	#. Checks for undefined fields - ELEMENT, ELEMENT_NUM, NUM, NAME,
	#  that all atoms connected to another atom actually
	#  exist and that atoms are not connected to themselves.
	#. Atoms with problems are marked with the flag $atom->{ERROR}
	#
	#; Requires: molecule.
	#; Returns: minor error count or -1 if there is a major problem
	#>

	my $mol = $_[0];
	
	my $atom_error;
	my $numbondcount;
	my $error = 0; # minor errors
	my $merror = 0; # major errors

	if (!defined $mol->{NAME}) {
		silico_msg('w', "Molecule has no defined name.\n",
				"Name will be set to \"<unnamed>\".\n");
		$mol->{NAME} = '<unnamed>';
	}
	
	silico_msg('q', "Checking integrity of molecule \"$mol->{NAME}\"\n",
			"\n");
	
	# First check: Source file type
	if (defined $mol->{SOURCE_FILE_TYPE}) {
		silico_msg('q', "Source file type: $mol->{SOURCE_FILE_TYPE}\n");
	} else {
		silico_msg('w', "Source file type is undefined.\n");
	}
	
	# Second check: number of atoms
	#	- is $mol->{NUMATOMS} defined?
	#	- is it an acceptable value?
	#	- is it the same as the number of atoms listed in @{$mol->{ATOMS}}?
	
	if (!defined $mol->{NUMATOMS}) {
		silico_msg('e', "\$mol->{NUMATOMS} is not defined.\n");
		return -1;
	} elsif ($mol->{NUMATOMS} < 0) {
		silico_msg('w', "Value of \$mol->{NUMATOMS} is negative ($mol->{NUMATOMS}).\n");
		++$error;
	} elsif ($mol->{NUMATOMS} == 0) {
		silico_msg('w', "Value of \$mol->{NUMATOMS} is zero.\n");
	} elsif ($mol->{NUMATOMS} != $#{$mol->{ATOMS}}+1) {
		silico_msg('w', "\$mol->{NUMATOMS} does not match actual atom count.\n",
				"\$mol->{NUMATOMS}: $mol->{NUMATOMS}    Number of atoms: ".($#{$mol->{ATOMS}}+1)."\n");
		$merror = -1;
	} else {
		silico_msg('q', "Number of atoms: $mol->{NUMATOMS}\n");
	}
	
	# Third check: Number of bonds
	#	- Is the number of bonds defined?
	#	- Is the number of bonds positive?
	#	- Does the value of $mol->{NUMBONDS} match the number of listed bonds?
	
	$numbondcount = molecule_count_numbonds($mol);
	
	if (!defined $mol->{NUMBONDS}) {
		silico_msg('w', "\$mol->{NUMBONDS} is not defined!\n");
		$merror = -1;
	} elsif ($mol->{NUMATOMS} < 0) {
		silico_msg('w', "Value of \$mol->{NUMBONDS} is negative ($mol->{NUMBONDS}).\n");
		++$error;
	} elsif ($mol->{NUMATOMS} == 0) {
		silico_msg('w', "Value of \$mol->{NUMBONDS} is zero.\n");
	} elsif ($numbondcount != $mol->{NUMBONDS}) {
		silico_msg('w', "\$mol->{NUMBONDS} does not match actual bond count.\n",
				"\$mol->{NUMBONDS}: $mol->{NUMBONDS}    Number of bonds: $numbondcount\n");
		$merror = -1;
	} else {
		silico_msg('q', "Number of bonds: $mol->{NUMBONDS}\n");
	}
	
	$atom_error = mol_check_atoms($mol);

	$error += $atom_error;
	
	return $error;
}

sub mol_check_atoms {

	my $mol = $_[0];
	
	my $atom;
	my $error;
	my $i = -1;
	my $con;
	my $conat;
	my $connum;
	my $j;
	my $flag;
	
	foreach $atom (atoms($mol)) {
	
		++$i;
		
		if (!defined $atom->{NUM}){
			silico_msg('w', "\$atom->{NUM} is not defined for atom $i.\n");
			$atom->{ERROR} = 1;
			++$error;
		} elsif ($atom->{NUM} != $i+1) {
			silico_msg('w', "\$atom->{NUM} has an incorrect value: $atom->{NUM}. Should be: ".($i+1)."\n");
			$atom->{ERROR} = 1;
			++$error;
		}

		if (!defined $atom->{NAME}){
			silico_msg('w', "\$atom->{NAME} is not defined for atom $i.\n");
			$atom->{ERROR} = 1;
			++$error;
		}
		
		if (!defined $atom->{SUBID}){
			silico_msg('w', "\$atom->{SUBID} is not defined for atom $i.\n");
			$atom->{ERROR} = 1;
			++$error;
		}
		
		# Each atom should have at least the following defined
		if (!defined $atom->{ELEMENT}) {
			silico_msg('w', "\$atom->{ELEMENT} is not defined for atom $i.\n");
			$atom->{ERROR} = 1;
			++$error;
		}
		
		if (!defined $atom->{ELEMENT_NUM}){
			silico_msg('w', "\$atom->{ELEMENT_NUM} is not defined for atom $i.\n");
			$atom->{ERROR} = 1;
			++$error;
		}

		# If $mol->{HAS_NO_COORDS} is not set,
		# every atom should have coordinates

		if (!$mol->{HAS_NO_COORDS}) {

			if (!defined $atom->{X}){
				silico_msg('w', "Atom $i has no defined X coordinate.\n");
				$atom->{ERROR} = 1;
				++$error;
			}
			if (!defined $atom->{Y}){
				silico_msg('w', "Atom $i has no defined Y coordinate.\n");
				$atom->{ERROR} = 1;
				++$error;
			}
			if (!defined $atom->{Z}){
				silico_msg('w', "Atom $i has no defined Z coordinate.\n");
				$atom->{ERROR} = 1;
				++$error;
			}
		}

		# Check connected atoms
		for ($j=0; $j<=$#{$atom->{CONNECT}}; ++$j) {
			
			$connum = $mol->{ATOMS}[$i]{CONNECT}[$j];
			$con = $mol->{ATOMS}[$connum];
			
			# Check each bond only once
			next if $connum > $i;

			# Check that all connection records are defined
			if (!defined $connum) {

				silico_msg('w', "A connection record is not properly defined.\n");
				print_con ($atom, $con);
				$atom->{ERROR} = 1;
				++$error;
			}
			
			# Check to see if this atom is connected to itself
			if ($connum == $i) {
				silico_msg('w', "Atom $i is connected to itself.\n");
				print_con ($atom, $con);
				$atom->{ERROR} = 1;
				++$error;
			}
			
			# Check that the connected atom has its own CONNECT record back to
			# this atom (the 'parent')
			$flag = 0;
			foreach $conat (@{$con->{CONNECT}}) {
				
				if ($conat == $i) {
					$flag = 1;
					last;
				}
			}
			
			if ($flag == 0) {
				silico_msg('w', "Atom $i has a connection record to atom $connum,\n",
						"but atom $connum has no connection record to atom $i.\n");
				$atom->{ERROR} = 1;
				++$error;
			}
		}
	
	}
}


sub mol_check_bondlengths {

	#<
	#? Check molecule bondlengths
	#. Sets $atom->{BONDLENGTH} on short and long bonds
	# (generally those more than 10% longer than standard bond lengths)
	#; Requires: molecule.
	#. Returns; 1 if OK or -1 if there is a major problem
	#>

	my $mol = $_[0];
 
	my $atom;
	my $bl;
	my $bond;
	my $con;
	my $connum;
	my $distance;
	my $error = 0;
	my $i;
	my $j;
	my $sbond;
	
	if ($mol->{HAS_NO_COORDS}) {
	
		silico_msg('w', "Molecule has no coordinates!\n",
				"Can not check bond lengths.\n");
		return 1;
	}

	#These are maximum allowable lengths and are generally 10% longer than
	# standard bond lengths
	$bond = bondlengths();
	
	$i = -1;
	foreach $atom (atoms($mol)) {
	
		++$i;
	
		# Connected atoms
		for ($j=0; $j<=$#{$atom->{CONNECT}}; ++$j) {
			
			$connum = $mol->{ATOMS}[$i]{CONNECT}[$j];
			$con = $mol->{ATOMS}[$connum];
			
			# Check each bond only once
			next if $connum > $i;

			$bl = $bond->{$atom->{ELEMENT}."_".$con->{ELEMENT}};
		
			# Check bond length
			$distance = distance($atom, $con);
			
			if (defined $bl) {
					
				if ($distance > $bl) {
					
					my $string = sprintf "Long bond: $atom->{ELEMENT} $con->{ELEMENT}. Distance: %f5.2 Max allowed length: %f5.2\n", $distance, $bl;
					
					silico_msg('w', "$string");
					print_con ($atom, $con);
					$atom->{BONDLENGTH} = 1;
					$con->{BONDLENGTH} = 1;
					++$error;
				}
				
				$sbond = $bl*0.75;
				if ($distance < $sbond) {
					
					my $string = sprintf "Short bond: $atom->{ELEMENT}. $con->{ELEMENT} Distance: %f5.2 Min allowed length: %f5.2\n", $distance, $sbond;
					
					silico_msg('w', "$string");
					print_con ($atom, $con);
					$atom->{BONDLENGTH} = 1;
					$con->{BONDLENGTH} = 1;
					++$error;
				}
				
			} else {
			
				if ($distance > 2 ) {
					
					my $string = sprintf "Long bond: $atom->{ELEMENT}. $con->{ELEMENT} Distance %f5.2 (longer than 2 Angstroms)\n", $distance;
					
					silico_msg('w', "$string");
					print_con ($atom, $con);
					$atom->{BONDLENGTH} = 1;
					$con->{BONDLENGTH} = 1;
					++$error;
				}
				
				if ($distance < 0.7) {
					
					my $string = sprintf "Short bond: $atom->{ELEMENT}. $con->{ELEMENT} Distance %f5.2 (shorter than 0.7 Angstroms)\n", $distance;
					
					silico_msg('w', "$string");
					print_con ($atom, $con);
					$atom->{BONDLENGTH} = 1;
					$con->{BONDLENGTH} = 1;
					++$error;
				}
			}
		}
	}
	
	return $error;
}


sub mol_check_atom_overlap {

	#<
	#? Pariwise check that no atoms in molecule overlap
	#. Checks that atoms are not closer than 90% of the min VDW distance.
	#; Requires: molecule, flag to supress warning messages and atom labelling
	#. Returns: Overlap count
	#>
	
	my $mol = $_[0];
	my $suppress = $_[1];
	
	my $atom1;
	my $atom2;
	my $d2;
	my $error = 0;
	my $i;
	my $j;
	my $thresh2 = 0.9**2;
	my ($x2,$y2,$z2);
	
	if ($mol->{HAS_NO_COORDS}) {
	
		silico_msg('w', "Molecule has no coordinates!\n",
				"Can not check for atom overlap.\n");
		return 0;
	}
	
	for ($i = 0; $i < $#{$mol->{ATOMS}}; ++$i) {
		
		$atom1 = $mol->{ATOMS}[$i];
		
		for ($j = $i+1; $j <= $#{$mol->{ATOMS}}; ++$j) {
		
			$atom2 = $mol->{ATOMS}[$j];
			
			$x2 = ($atom1->{X} - $atom2->{X})**2;
			next if $x2 > $thresh2;
			$y2 = ($atom1->{Y} - $atom2->{Y})**2;
			next if $y2 > $thresh2;
			$z2 = ($atom1->{Z} - $atom2->{Z})**2;
			next if $z2 > $thresh2;
			
			$d2 = $x2 + $y2 + $z2;
	
			if ($d2 < $thresh2) {
			
				if (!$suppress) {
					silico_msg('w', "Atoms overlap!\n",
							"Atom 1 -- Count: ".($i+1)."  Number: $atom1->{NUM}  Name:\"$atom1->{NAME}\"  Residue: $atom1->{SUBNAME}$atom1->{SUBID}\n",
							"Atom 2 -- Count: ".($j+1)."  Number: $atom2->{NUM}  Name:\"$atom2->{NAME}\"  Residue: $atom2->{SUBNAME}$atom2->{SUBID}\n",
							"Distance: ".distance($atom1, $atom2)." Angstroms\n");
					$atom1->{OVERLAP} = 1;
					$atom2->{OVERLAP} = 1;
				}
				++$error;
			}
		}
	}
	
	return $error;
}

sub mol_check_distorted_rings {

	#<
	#? Look for distorted aromatic rings
	#; Requires: molecule
	#. Returns: Number of bent aromatic rings
	#>
	
	my $mol = $_[0];

	my $atom;
	my $atomnum;
	my @bolist;
	my $error = 0;
	my $j;
	my $flag;
	my $ring;
	
	if ($mol->{HAS_NO_COORDS}) {
	
		silico_msg('w', "Molecule has no coordinates!\n",
				"Can not check for distorted rings.\n");
		return 0;
	}
	
	molecule_find_rings($mol, 6);
	
	foreach $ring (@{$mol->{RINGS}}) {
	
		# Only do 6-membered rings
		next if $#{$ring} != 5;
	
		$flag = 0;
		RINGATOM: foreach $atomnum (@$ring) {
		
			$atom = $mol->{ATOMS}[$atomnum];
			
			# Skip if we are in a planar ring
			next RINGATOM if $atom->{PRINGS}[0];
			
			@bolist = find_bondorders_ring($ring, $mol);
			
			$j = join '', @bolist;
						
			if ($j eq '444444') {
				$flag = 1;
				last RINGATOM;
			}
			
			if ($j eq '121212' || $j eq '212121') {
				$flag = 1;
				last RINGATOM;
			}
		}
		
		# Mark bad rings
		if ($flag) {
		
			++$error;
		
			foreach $atomnum (@$ring) {
			
				$atom = $mol->{ATOMS}[$atomnum];
			
				$atom->{ERROR} = 1;
				$atom->{BENT_RING} = 1;
			}
		}
	}
	
	silico_msg('q', "Found $error distorted rings\n") if $error;
	
	return $error;
}

sub mol_check_aromatic_bonds {

	#<
	#? Subroutine to look for atoms, not in aromatic rings, which
	#  have aromatic bonds
	#. This benefits from mol_check_distorted_rings being run first
	#  though that is not strictly necessary
	#; Requires: molecule
	#. Returns: Number of atoms with out-of-place aromatic bonds
	#>
	
	my $mol = $_[0];
	
	my $atom;
	my $bo;
	my $error = 0;
	my $flag;
	
	foreach $atom (atoms($mol)) {
	
		$flag = 0;
	
		# We have a bona fide aromatic atom if it is
		# in a planar ring or bent aromatic ring
		next if $atom->{PRINGS}[0];
		next if $atom->{BENT_RING};
		
		# Set a flag if this atom has an aromatic bond
		foreach $bo (@{$atom->{BORDERS}}) {
		
			next if $bo != 4;
			
			$flag = 1;
			last;
		}
		
		# If a flag has been set for this atom, mark it
		if ($flag) {
		
			++$error;
			
			$atom->{ERROR} = 1;
			$atom->{BAD_AROM_BOND} = 1;
		}
	}
	
	if ($error == 0) {
		silico_msg('q', "No inappropriate aromatic bonds were found.\n");
	} elsif ($error == 1) {
		silico_msg('q', "One atom was found with inappropriate aromatic bonds.\n");
	} else {
		silico_msg('q', "$error atoms were found with inappropriate aromatic bonds.\n");
	}
	
	return $error;
}

sub mol_check_amides {

	#<
	#? Subroutine that checks for non-trans amides
	#. Checks two dihedrals: C-C-N-C and O-C-N-H, as each of these
	#  should be 180 degrees
	#; Requires: molecule, name of output file
	#. Returns: Number of atoms in distorted amide groups
	#>
	
	my $mol = $_[0];
	
	my $ostring; # String containing data for output file

	my $cis = 0;
	my $dcis = 0;
	my $dtrans = 0;
	my $trans = 0;
	
	if ($mol->{HAS_NO_COORDS}) {
	
		silico_msg('w', "Molecule has no coordinates!\n",
				"Can not check amide geometry.\n");
		return 0;
	}
	
	# Before beginning, make sure all functional groups are labelled
	mol_label_functional_group ($mol) if !$mol->{FG};

	CARBON: foreach my $atom (atoms($mol)) {
	
		my $hcount = 0;
		my $ncount = 0;
		
		my $amide_c;
		my $amide_o;
		my $amide_n;
		my $amide_h;
		my $other_on_c;
		my $other_on_n;

		my $atomnum_on_c;
		my $atom_on_c;
		my $atomnum_on_n;
		my $atom_on_n;
		
		# Only look at amide carbons
		next CARBON if (!$atom->{FG}{S_AMIDE_C});
		
		# Set this atom as being the 'amide carbon'
		$amide_c = $atom;
		
		# For each atom connected to the amide carbon...
		ATOM_ON_C: foreach $atomnum_on_c (@{$amide_c->{CONNECT}}) {
		
			# Get the full information for this atom
			$atom_on_c = $mol->{ATOMS}[$atomnum_on_c];
			
			# If it is marked as an amide nitrogen then...
			if ($atom_on_c->{FG}{S_AMIDE_N}) {
				
				# Set it as being _the amide nitrogen
				$amide_n = $atom_on_c;
				
				# Count the number of amide nitrogens found so far
				++$ncount;
				
				# Move to the next amide carbon if there are 2 or more
				# (this is a copout - we aren't doing ureas now)
				if ($ncount == 2) {next CARBON}
				
				# Otherwise keep on getting info about this amide
				else {next ATOM_ON_C}
			}
			
			# If it is marked as an amide oxygen then...
			if ($atom_on_c->{FG}{S_AMIDE_O}) {
			
				# Set it as being _the amide oxygen
				$amide_o = $atom_on_c;
				next ATOM_ON_C;
			}
			
			# If it is neither an amide oxygen or an amide nitrogen,
			# set it as being the 'other' atom
			$other_on_c = $atom_on_c;
		}
			
		# Now we look at the amide nitrogen
		
		# For each atom on the amide nitrogen...
		ATOM_ON_N: foreach $atomnum_on_n (@{$amide_n->{CONNECT}}) {
			
			# Get the full information for this atom
			$atom_on_n = $mol->{ATOMS}[$atomnum_on_n];
				
			# Skip it if it is the amide carbon - we don't need to deal with these here
			next ATOM_ON_N if $atom_on_n == $amide_c;

			# If it is an amide hydrogen (ie, on the N), then...
			if ($atom_on_n->{ELEMENT} eq 'H') {
			
				#Mark it as an amide hydrogen
				$amide_h = $atom_on_n;
				
				# Count the number of amide hydrogens found so far
				++$hcount;
				
				# If there are two move on to the next amide carbon,
				# as we aren't considering primary amides here (there is no point)
				if ($hcount == 2) {next CARBON}
				
				# Otherwise, move on to the next atom
				else {next ATOM_ON_N}
			}
				
			# If an atom isn't the amide carbon or an amide nitrogen,
			# mark it as being the 'other' atom
			$other_on_n = $atom_on_n;
		}
		
		# Skip this amide now if we don't have an amide N (unlikely), or it is a tertiary amide (quite possible)
		next CARBON if ($ncount == 0 || $hcount == 0);
	
		# Calculate the dihedral angles
		# Dihedral1 is the "standard" dihedral, notionally C-C(O)-N(H)-C
		# Dihedral2 is the O-C-N-H dihedral which is used to check for distortion but not cis-trans
		my $omega = dihedral ($other_on_c, $amide_c, $amide_n, $other_on_n);
		my $omega1 = dihedral ($amide_o, $amide_c, $amide_n, $amide_h);
	
		my $omega_deg = rad_to_deg($omega);
		$omega_deg = 360 + $omega_deg if $omega_deg < 0;
		my $omega1_deg = rad_to_deg($omega1);
		$omega1_deg = 360 + $omega1_deg if $omega1_deg < 0;
		
		my  $cres;
		$cres = $amide_c->{SUBNAME}.($amide_c->{SUBID}||'undef');
		$cres =~ s/\s+//g;
		
		my  $nres;
		$nres = $amide_n->{SUBNAME}.($amide_n->{SUBID}||'undef');
		$nres =~ s/\s+//g;
		
		# Amides are said to be not distorted if the dihedrals are within 0.2 radians
		# (about 11 or 12 degrees) of being planar
		
		# Print to the output file
                $ostring .= "$amide_c->{NUM}($cres)\t$amide_n->{NUM}($nres)\t$omega_deg\t$omega1_deg\n";
		
		# An undistorted cis amide: Both dihedral angles are close to 0
		if ((-0.2 < $omega && $omega < 0.2) && (-0.2 < $omega1 && $omega1 < 0.2)) {
			
			$amide_c->{CIS_AMIDE} = 1;
			$amide_n->{CIS_AMIDE} = 1;
			$amide_h->{CIS_AMIDE} = 1;
			$amide_o->{CIS_AMIDE} = 1;
			$other_on_c->{CIS_AMIDE} = 1;
			$other_on_n->{CIS_AMIDE} = 1;
			
			++$cis; # Count the number of undistorted cis amides
			next CARBON;
		}
		
		# An undistorted trans amide (best case): Both dihedral angles are close to either Pi or -Pi
		if ((0.2 - $Pi > $omega || $Pi - 0.2 < $omega)
		&& (0.2 - $Pi > $omega1 || $Pi - 0.2 < $omega1)) {
		
			++$trans; # Count the number of undistorted trans amides
			next CARBON;
		}
		
		# A distorted cis amide: Dihedral 1 is between 0.2 and Pi/2 (about 12 degrees and 90 degrees)
		# This takes no account of Dihedral 2
		if ((-0.2 >= $omega && $omega > $Pi/-2) || (0.2 <= $omega && $omega < $Pi/2)) {
		
			$amide_c->{DISTORTED_CIS_AMIDE} = 1;
			$amide_n->{DISTORTED_CIS_AMIDE} = 1;
			$amide_h->{DISTORTED_CIS_AMIDE} = 1;
			$amide_o->{DISTORTED_CIS_AMIDE} = 1;
			$other_on_c->{DISTORTED_CIS_AMIDE} = 1;
			$other_on_n->{DISTORTED_CIS_AMIDE} = 1;
			
			++$dcis; # Count the number of distorted cis amides
			next CARBON;
		}
		
		# A distorted trans amide: Anything else (ie Dihedral1 is between Pi/2 and Pi - 0.2)
		# This acts as a catch all so it will also get if dihedral fails, maybe?
		else {
			
			$amide_c->{DISTORTED_TRANS_AMIDE} = 1;
			$amide_n->{DISTORTED_TRANS_AMIDE} = 1;
			$amide_h->{DISTORTED_TRANS_AMIDE} = 1;
			$amide_o->{DISTORTED_TRANS_AMIDE} = 1;
			$other_on_c->{DISTORTED_TRANS_AMIDE} = 1;
			$other_on_n->{DISTORTED_TRANS_AMIDE} = 1;
			
			++$dtrans; # Count the number of distorted trans amides
			next CARBON;
		}
	}
	
	# Print out a summary of amides (reports back on each type found)
	
	my $error = $dtrans + $cis + $dcis; # Calculate the overall number of bad amides
	my $total = $trans + $error; # Calculate the overall number of amides (good or bad)
	
	
	return $total, $error, $trans, $dtrans, $cis, $dcis, $ostring;
}

sub mol_check_valences {

	#<
	#? Simple check that atom valences are OK
	#. Calls atom_check_valence
	#; Requires: molecule, flag to mark that molecule is missing hydrogen atoms,
	#  flag to convert non-ring aromatic bonds to single/double bonds
	#; Returns: Negative of number of atoms with abberant valences
	#. An atom->{VALENCE} flag is set on each offending atom
	#>

	require silico_hydrogens;
	
	my $mol = $_[0];
	my $noh = $_[1]; # Flag to denote hydrogens are missing from molecule
	my $convert = $_[2]; # Flag to convert non-ring aromatic bonds to double/single
	
	my $atom;
	my $errorcount = 0;
	my $val;
	
	if ($noh) {
		silico_msg('q', "This molecule contains no hydrogen atoms.\n",
				"Atoms will only be checked for exceeding maximum allowed valences.\n");
	}
	
	convert_aromatic_bonds($mol) if $convert;
	
	foreach $atom (atoms($mol)) {
	
		$val = atom_check_valence($atom, $mol, $noh);
		if ($val == -1) {
			$atom->{VALENCE} = 1;
			
			--$errorcount;
		}
	}
	
	return $errorcount;
}

sub atom_check_valence {

	#<
	#? Simple check that atom valences are OK
	#; Requires: atom, molecule, (optional: flag to denote hydrogens not present)
	#; Returns: -1 if error, otherwise returns valence
	#>

	my $atom = $_[0];
	my $mol = $_[1];
	my $noh = $_[2] || 0; # Molecule is missing hydrogens
	
	my $arcount = 0;
	my $el;
	my $valence;
	
	# Calculate valence
	
	$el = $atom->{ELEMENT};
	$valence = 0; # No bonds to start with
	foreach (@{$atom->{BORDERS}}) {
	
		if ($_ == 1 || $_ == 2 || $_ == 3) {
			$valence += $_;
			next;
		}
		
		# Aromatic bonds. A total of 3 aromatic bonds to carbon is OK
		if ($_ == 4) {
			$valence += 1.5;	# Aromatic
			++$arcount;
			$valence = 4 if $arcount == 3 ;
			next;
		}
		
		silico_msg('w', "Atom $atom->{NUM} has poorly defined bond orders!\n");
	}
	
	# Take formal charge into account (Note: this has not always been set)
	if ($atom->{FORMAL_CHARGE}) {
		$valence -= $atom->{FORMAL_CHARGE};
	}
	
	if ($el eq "C") {
	
		return 4 if $valence == 4; # Neutral
		if ($valence == 0) {
			silico_msg('w', "Atom $atom->{NUM} appears to be an isolated carbon!\n");
			check_printerror($mol, $atom, $valence);
			return -1;
		}
		
		$noh && return $valence if $valence < 4;
		check_printerror($mol, $atom, $valence);

		return -1;
	}

	# Hydrogen
	if ($el eq 'H') {
	
		return 1 if $valence == 1; # Neutral
		
		if ($valence == 0) {
			silico_msg('w', "Atom $atom->{NUM} appears to be an isolated hydrogen!\n");
			check_printerror($mol, $atom, $valence);
			return -1;
		}
	}
	
	
	if ($el eq "N") {
	
		if ($valence == 2) {
			return 2; # -1
		} elsif  ($valence == 3) {
			return 3; # Neutral
		} elsif ($valence == 3.5) {
			return 3.5; # Guanidinium with aromatic bond orders
		} elsif ($valence == 4) {
			return 4; # +1
		} else {
			check_printerror($mol, $atom, $valence);
			return -1;
		}
	}

	if ($el eq "O") {
	
		return 2 if $valence == 2; # Neutral
		# Carboxylate with aromatic bonds
		if ($valence == 1.5) {
			return 1.5;
		}
		return 1 if $valence == 1; # -1
		$noh && return $valence if $valence < 1;
		check_printerror($mol, $atom, $valence);
		return -1;
	}

	if ($el eq "S") {
	
		
		return 6 if $valence == 6; # Neutral
		return 4 if $valence == 4; # Neutral
		return 3 if $valence == 3; # S+
		return 2 if $valence == 2; # Neutral
		return 1 if $valence == 1; # S-
		
		$noh && return $valence if $valence < 2;
		check_printerror($mol, $atom, $valence);
		return -1;
	}

	if ($el eq "P") {
	
		return 5 if $valence == 5; # Neutral
		return 4 if $valence == 4; # P+?
		return 3 if $valence == 3; # Neutral
		return 2 if $valence == 2; # P-
		check_printerror($mol, $atom, $valence);
		
		return -1;
	}
	
	# Halogen things
	if ($el eq 'Cl'||$el eq 'F'||$el eq 'Br'||$el eq 'I') {
	
		return 1 if $valence == 1; # Neutral
		return 0 if $valence == 0; # -1 ion
		$noh && return $valence if $valence < 1;
		check_printerror($mol, $atom, $valence);
		return -1;
	}
	
	if ($el eq "Na" || $el eq "K") {
	
		return 0 if $valence == 0; # +1 ion
		check_printerror($mol, $atom, $valence);
		return -1;
	}
	
	silico_msg('w', "Can not check valence on element: $atom->{ELEMENT} valence: $valence\n");
	return 0;
}

sub mol_check_bondorders {

	#<
	#? Simple check that bondorders are OK
	#. An atom->{ERROR} flag is set on each offending atom
	#; Requires: molecule, 
	#; Returns: 0 if OK,  otherwise number of unusual bonds found
	#>

	my $mol = $_[0];
	
	my $atom;
	my $bond;
	my $error_bond_du = 0;
	my $error_bond_ambiguous = 0;
	my $error_undefined_bond = 0;

	foreach $atom (atoms($mol)) {

		foreach $bond (@{$atom->{BORDERS}}) {

			if (!defined $bond) {

				++$error_undefined_bond;
				$atom->{ERROR} = 1;
				next;
			}
			
			next if $bond == int $bond && $bond >= 1 && $bond <= 4;

			++$error_bond_du if $bond == 0;
			++$error_bond_ambiguous if $bond >4;
			++$error_bond_ambiguous if $bond <0;
			$atom->{ERROR} = 1;
		}
	}
			
	$error_bond_du /= 2;
	$error_bond_ambiguous /= 2;
	$error_undefined_bond /= 2;

	silico_msg('w', "Found $error_undefined_bond bond(s) undefined bondorder\n") if $error_undefined_bond != 0;
	silico_msg('w', "Found $error_bond_ambiguous bond(s) with dummy (zero order) bonds\n") if $error_bond_ambiguous != 0;
	silico_msg('w', "Found $error_bond_du bond(s) with unusual bond orders (ie > 4 or < 0)\n") if $error_bond_du != 0;

	return ($error_bond_du +  $error_bond_ambiguous + $error_undefined_bond);
}

sub mol_check_elements {

	#<
	#? Check to make sure that all atom elements are sensible
	#; Requires: molecule
	#. An atom->{ERROR} flag is set on each offending atom
	#; Returns: Number of errors
	#>
	
	my $mol = $_[0];
	
	my $atom;
	my $error = 0;

	use vars qw (@Atomic_Elements %Atomic_Elements);
	
	# Set up hash of elements (but not Du)
	if (!defined $Atomic_Elements{H}) {
	
		# Dummy atom is element 0
		foreach (@Atomic_Elements[1..$#Atomic_Elements]) {
			$Atomic_Elements{$_} = 1;
		}
	}
	
	foreach $atom (atoms($mol)) {
	
		if ($atom->{ELEMENT} eq 'Du') {
		
			$atom->{ERROR} = 1;
			silico_msg('w', "A dummy atom was found!\n",
					"Number: ".($atom->{NUM}||'undefined').", Name: ".($atom->{NAME}||'undefined').", SubID: ".($atom->{SUBID}||'undefined').", SubName: ".($atom->{SUBNAME}||'undefined')."\n");
			++$error;
		}
		
		if (!$Atomic_Elements{$atom->{ELEMENT}}) {
		
			$atom->{ERROR} .= "Unknown_element";
			silico_msg('w', "An atom of an unknown element was found!\n",
					"Number: ".($atom->{NUM}||'undefined').", Name: ".($atom->{NAME}||'undefined').", SubID: ".($atom->{SUBID}||'undefined').", SubName: ".($atom->{SUBNAME}||'undefined').", Element: ".($atom->{ELEMENT}||'undefined')."\n");
			++$error;
		}
	}
	
	if ($error == 0) {
		silico_msg('q', "All atom elements are sensible.\n");
	} elsif ($error == 1) {
		silico_msg('q', "One atom was found with a bad element (dummy or unknown).\n");
	} else {
		silico_msg('q', "$error atoms were found with bad elements (dummy or unknown).\n");
	}

	return $error;
}

sub check_printerror {

	#<
	#? Print an error line for an atom with incorrect valences
	#; Requires: molecule, atom, valence
	#>

	my $mol = $_[0];
	my $atom = $_[1];
	my $valence = $_[2];
	
	my $conatom;
	my $i;
	
	silico_msg('w', "Atom $atom->{NUM} appears to have a bad valence ($valence)!\n",
			"Name: ".($atom->{NAME}||'undef').", SubID: ".($atom->{SUBID}||'undef').", SubName: ".($atom->{SUBNAME}||'undef').", Element: ".($atom->{ELEMENT}||'undef').", Fchg: ".($atom->{FORMAL_CHARGE}||'undef')."\n");
	
	$i = 0;
	
	foreach $conatom (@{$atom->{CONNECT}}) {
		
		my $string;
		my $con;
		
		$con = $mol->{ATOMS}[$conatom];
		
		$string = sprintf "    --> bo: %-2s num: %-5s name: %-4s subid: %-3s sname: %-4s el: %-2s ",
			($atom->{BORDERS}[$i]|| 'undef'), ($con->{NUM}||'undef'), ($con->{NAME}||'undef'), ($con->{SUBID}||'undef'), ($con->{SUBNAME}||'undef'), ($con->{ELEMENT}||'undef');
		
		silico_msg('q', $string);
		
		if (defined $atom->{X} && defined $atom->{Y} && defined $atom->{Z}
			&& defined $con->{X} && defined $con->{Y} && defined $con->{Z}) {
			
			silico_msg('q', "Distance: ".distance($atom, $con)."\n");
		} elsif (!$mol->{HAS_NO_COORDS}) {
			silico_msg('q', "Distance: undefined\n");
		} else {
			silico_msg('q', "\n");
		}
		
		++$i;
	}
}

return 1;
