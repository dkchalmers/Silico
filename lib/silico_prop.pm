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
#! silico_prop.pm
#? Silico routines to calculate molecular and atomic properties including
#  molecular interface interactions
#. $Revision: 1.15-72-gda9b5d8 $
#>

use strict;
package Silico;


##################################################################
#
#	Atom/bond/molecule properties
#
##################################################################

sub molecule_formal_charge {

	#<
	#? Calculate the formal charge on a molecule.
	#  Note. Formal charges are recalcuated for each atom which may REMOVE formal charges that are already present.
	#  This happens particularly if no hydrogens are present in the molecule
	#. Calls atom_formal_charge
	#; Requires: molecule
	#; Optional: atom_list, optional flag to denote missing hydrogens, flag to print atoms with nonzero formal charge
	#; Returns: total formal charge 
	#; Sets: $mol->{FORMAL_CHARGE}, $mol->{CHARGED_ATOMS};
	#>

	my $mol = $_[0];
	my $list = $_[1] || $mol->{ATOMS};
	my $noh = $_[2];
	my $print = $_[3];

	my $total = 0;
	
	#
	# Warning: subrotine has problems and can give the wrong answer
	# Needs work
	#
	
	silico_msg('w', "Recalculating charges for a molecule that already has formal charges\n") if $mol->{HAS_FORMAL_CHARGES};
	
	molecule_check_and_fix_connectivity($mol);
	
	mol_has_h($mol) if !defined $mol->{HAS_HYDROGENS};
	
	if ($mol->{HAS_HYDROGENS} < 0) {
		silico_msg ('e', "Molecule has unbonded hydrogens. Skipping\n");
		return;
	}
	if ($Silico::debug && $mol->{HAS_HYDROGENS} < 2) {
		silico_msg ('n', "Molecule is missing hydrogens\n");
	} 
	
	@{$mol->{CHARGED_ATOMS}} = ();
	
	my $i = 0;
	foreach my $atom (@$list) {
	
                my $val = atom_formal_charge($atom, $mol, $noh);
		
		my $c = $#{$atom->{CONNECT}}+1;
		#print "el $atom->{ELEMENT} chg $val ncon $c\n";
		
		$total += $val;
		if ($val != 0) {
			push @{$mol->{CHARGED_ATOMS}}, $atom;
			
			if ($print) {
				print "Num\tName\tSubid\tChain\tSubname\tCharge\n" if $i == 0;
				printf "%s\t%s\t%s\t%s\t%s\t%s\n", $atom->{NUM}, $atom->{NAME} || '-', $atom->{SUBID} || '-', $atom->{CHAIN} || '-', $atom->{SUBNAME} || '-', $atom->{FORMAL_CHARGE} || '-';
			}
			++$i;
		}
		
        }
	
	$mol->{FORMAL_CHARGE} = $total;
	$mol->{SDF_DATA}{FORMAL_CHARGE} = $total;
	$mol->{HAS_FORMAL_CHARGES} = 1;

	return $total;
}

sub atom_formal_charge {
	
	#<
	#? Calculate the formal charge on an atom.
	#. Uses the acceptable oxidation states of that atom,
	#  and whether that element is more intrinsically likely
	#  to be a cation or an anion.  Assumes that all valences are
	#  filled with hydrogens if 'noh' flag is set
	#; Requires: atom, molecule, optional flag to denote missing hydrogens, flag to print charged atoms
	#; Returns: formal charge 
	#>
	
	my $atom = $_[0];
	my $mol = $_[1];
	my $noh = $_[2] || 0;
		
	if (!defined $atom->{ELEMENT}) {
		
		silico_msg('d', "Element of atom is not defined!\n",
				"Could not determine formal charge.\n",
				"Skipping.\n");
		atom_printout($atom);
		return undef;
	}
	
	# Molecule is marked as being bad already
	if (defined $mol->{HAS_HYDROGENS}) {
		$noh = 1 if $atom->{ELEMENT} eq 'C' && $mol->{HAS_HYDROGENS} < 2;
		$noh = 1 if $mol->{HAS_HYDROGENS} < 1;
	}
	
	acceptable_neutral_valences();
	my @nvalences;
	if (defined $Silico::Valences->{$atom->{ELEMENT}}) {

		@nvalences = (@{$Silico::Valences->{$atom->{ELEMENT}}});

	} elsif ($atom->{ELEMENT} eq 'A' || $atom->{ELEMENT} eq 'Q' || $atom->{ELEMENT} eq 'Du') {

		# Allow query atom types 'A' and 'Q', and Du atoms
		$atom->{FORMAL_CHARGE} = 0;
		return 0;

	} else {

		silico_msg('w', "Formal charge of an atom of element $atom->{ELEMENT} is not known.  Setting charge to 0!\n");
		$atom->{FORMAL_CHARGE} = 0;
		return 0;
	}
	
	my $valence = valence($atom);

	# This looks at each of the Acceptable Neutral Valences
	# for the atom (corresponding to oxidation states, etc.).
	# If an Acceptable Neutral Valence is found that corresponds
	# to the atom's actual valence, the atom is formally neutral.
	
	my $lanv = $nvalences[0];
	my $hanv;
	
	foreach my $anv (@nvalences) {

		# If this a.n.v. is greater than the valence, we have
		# run off the end, and don't need to bother with the rest
		if ($anv > $valence) {
			last;
		} elsif ($anv == $valence) {
			$atom->{FORMAL_CHARGE} = 0;
			return 0;
		} elsif ($anv < $valence) {

			# This sets a property called the Highest
			# Acceptable Neutral Valence ($hanv). This
			# is the highest valence corresponding to a
			# neutral atom in a valid oxidation state
			# that is less than the actual valence.

			$hanv = $anv if (!defined $hanv || $anv > $hanv);
		}
	}
		
	# If atom is missing hydrogens then assume that it has
	# the lowest acceptable neutral valence and that the
	# remaining valences are filled with hydrogen
	if ($noh && $valence <= $lanv) {
		$atom->{FORMAL_CHARGE} = 0;
		return 0;
	}
	
	my $lp;
	if ($hanv) {
		
		# If we have a lone pair available given the current
		# environment, the atom is held to have given it away and thus
		# to be formally positive. Otherwise, it is held to have
		# received a dative bond and to be formally negative.

		$lp = atom_lone_pairs($atom);
		my $fc = ($lp >= 0 ? $valence-$hanv : $hanv-$valence);
		$atom->{FORMAL_CHARGE} = $fc;
		return $fc;
	}
	
	# In the last resort (atom has less connections to it than even
	# the lowest Acceptable Neutral Valence), the formal charge is
	# worked out by considering the oxidation state of the atom, which
	# is computed using electronegativities. If the oxidation state
	# is negative, the formal charge is the result when the Lowest
	# Acceptable Neutral Valence is deducted from the actual valence.
	# Otherwise, the formal charge is the reverse. 

	my $os = atom_oxidation_state($atom, $mol);
	my $fc = ($os < 0 ? $valence-$lanv : $lanv-$valence);
	
	if (!defined $fc) {
		silico_msg('w', "Could not determine an atom's formal charge!\n",
				"Atom: $atom->{NUM} (\"$atom->{NAME}\")\n",
				"Setting this atom's formal charge to zero.\n");
		$fc = 0;
	}
	
	$atom->{FORMAL_CHARGE} = $fc;
	return $fc;
}

sub count_formal_charges {

	#<
	#? Calculate the formal charge on a molecule. Also count total numbers of positive and negative atoms
	#; Requires: Molecule
	#; Returns: Total charge. Counts of positive and negative atoms (not counting nitros, azides etc with adjacent +ve and -ve charges)
	#; Sets: $mol->{FORMAL_CHARGE}, $mol->{CHARGED_ATOMS}
	#>
	
	my $mol = $_[0];
	
	my $pos = 0;
	my $neg = 0;
	my $total = 0;
	
	@{$mol->{CHARGED_ATOMS}} = ();
	
	ATOM: foreach my $atom (atoms ($mol)) {
	
		my $c = $atom->{FORMAL_CHARGE};
		
		next if !$c;
		
		$total += $c;
		push @{$mol->{CHARGED_ATOMS}}, $atom;
		
		# Skip groups such as nitro, azide, etc which have adjacent 
		# opposite charges
		foreach my $con (connected($atom, $mol)) {
			next ATOM if $con->{FORMAL_CHARGE} && $con->{FORMAL_CHARGE} == -$c;
		}
		
		++$pos if $c > 0;
		++$neg if $c < 0;	
	}
	
	$mol->{FORMAL_CHARGE} = $total;
	$mol->{HAS_FORMAL_CHARGES} = 1;
	
	return $total, $pos, $neg;
}


sub atom_lone_pairs {
	
	#<
	#? Calculate the number of lone pairs present on an atom.
	#; Requires: atom
	#; Returns: number of lone pairs
	#>
	
	my $atom = $_[0];
	
	my $lp;
	my $shell;
	my @shells = (2, 8, 8, 18, 18, 32, 32, 50, 50);
	my $sum = 0;
	
	my $en = element_symbol2number($atom->{ELEMENT});
	return 0 if $en == 0;
	
	my $i = -1;
	while ($en > $sum) {
		++$i;
		$sum += $shells[$i];
	}
	
	$shell = $shells[$i];
	
	if (($sum-$en)*2 < $shell) {
		my $maxv = $sum - $en;
		$lp = $shell/2 - $maxv;
		$lp -= (valence($atom)-$maxv) if ($maxv < valence($atom));
	} else {
		$lp = 0;
	}
	
	# Incidentally, if $lp is less than zero, this implies that the
	# atom has received one or more dative bonds.
	return $lp;	
}

sub atom_oxidation_state {
	
	#<
	#? Calculate the oxidation state of an atom, not including
	#  a formal-charge component.
	#; Requires: atom, molecule
	#; Returns: number
	#>
		
	my $atom = $_[0];
	my $mol = $_[1];
	
	my $os = 0;
	
	acceptable_neutral_valences();
	electronegativities();

	my $anum = element_symbol2number($atom->{ELEMENT});
	return 0 if $anum == 0;
	my $aen = $Silico::Electronegativities[$anum];
	return 0 if $aen == 0;
	
	if ($#{$atom->{CONNECT}} >= 0) {
		foreach my $con (connected($atom, $mol)) {
			my $cnum = element_symbol2number($con->{ELEMENT});
			
			next if $cnum == 0;
			
			my $cen = $Silico::Electronegativities[$cnum];
			# Fudge: count nothing if the electronegativity of the
			# connected atom is not established (e.g., noble gases)
			$cen = $aen if $cen == 0;
			my $border = find_bondorder2($atom, $con, $mol);
			
			if ($aen > $cen) {
				$os -= $border;
			} elsif ($aen < $cen) {
				$os += $border;
			}
		}
	# If nothing is connected, the oxidation state is derived from the lowest
	# acceptable neutral valence and the electronegativity compared to that of
	# carbon. If the electronegativity is less than or equal to that of carbon,
	# the oxidation state is the lowest acceptable neutral valence. If the
	# electronegativity is greater than that of carbon, the oxidation state
	# is the negative of the lowest acceptable neutral valence. This will work
	# for most things (but will fail in the case of high-oxidation-state
	# transition metals).
	} else {
		return (-1*$Silico::Valences->{$atom->{ELEMENT}}[0]) if ($aen > $Silico::Electronegativities[6]);
		return $Silico::Valences->{$atom->{ELEMENT}}[0] if ($aen <= $Silico::Electronegativities[6]);
	}
	
	return $os;
}

sub find_polar_atoms {

	#<
	#? Mark atoms as polar or nonpolar
	#. Polar atoms are defined as being non C or H atoms or C
	#  or H atoms that are bonded to an element
	#  other than C or H.
	#; Requires: Molecule
	#; Returns: 
	#>

	my $mol = $_[0];
	
	foreach my $atom (atoms($mol)) {
	
		my $flag= 0;
		
		if ($atom->{ELEMENT} ne 'C' && $atom->{ELEMENT} ne 'H') {
		
			$atom->{POLAR} = 1;
			next;
		}
	
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			
			next if $con->{ELEMENT} eq 'C' || $con->{ELEMENT} eq 'H';
			
			$flag = 1;
			last;
		}
			
		if ($flag == 1) {
		
			$atom->{POLAR} = 1;
		} else {
			$atom->{NONPOLAR} = 1;
		}
	}
}


##################################################################
#
#	Single molecule properties
#
##################################################################

sub molecule_mw {
	
	#<
	#? Calculate molecule molecular weight 
	#; Requires: Molecule, optional list of atoms, flag to suppress check for united atoms (not implemented)
	#  flag to supress 'no hydrogens' warning message
	#; Returns: Molecular weight as a scalar
	#. Calcluates where hydrogens are missing (for united atom forcefields).
	#>

	my $mol = $_[0];
	my $list = $_[1];
	my $nocheck= $_[2]; # Not implemented
	my $nowarn = $_[3]; # Not implemented
	
	my $mw = 0;

	atomic_masses();
	
	molecule_check_and_fix_connectivity($mol) if !$mol->{CONNECTIVITY_CHECKED};
	
	mol_calc_hcount($mol, $list);
	
	@$list = @{$mol->{ATOMS}} if !defined $list;
	
	foreach my $atom (@$list) {
		
		my $mass = $Silico::Atomic_Masses[$atom->{ELEMENT_NUM}];
		if (!defined $mass || $mass !~ /\d/) {
			silico_msg('w', "Undefined mass for atom \"$atom->{NAME}\" (number $atom->{NUM})!\n",
					"Please check the \@Silico::Atomic_Masses array in silico_definitions.pm.\n",
					"This molecule will be skipped (setting MW to -1).\n");
					
			return -1;
		}
		$mw += $mass;
		
		if ($atom->{HMISSING}) {
			$mw += $Silico::Atomic_Masses[1] * $atom->{HMISSING};
		}
	}
	
	if ($Silico::debug) {
		print "Molecular weight: $mw\n";
		print "Molecular formula: ";	
		my $formula = molecule_formula($mol, $list);
		foreach my $el (sort keys %$formula) {
			print $el.$formula->{$el}." ";
		}
		print "\n";
	}
	
	return $mw;
}


sub mol_count_heavy_atoms {
	
	#<
	#? Calculate the number of heavy (non-hydrogen) atoms in a molecule.
	#; Requires: molecule
	#; Returns: whole number
	#>
	
	my $mol = $_[0];
	
	my $count = 0;
	

	# The SDF noparse option counts these already
	return $mol->{NUMHEAVY} if $mol->{SOURCE_FILE_TYPE} eq 'sdf_noparse';
	
	my $i = 0;
	foreach my $atom (atoms($mol)) {
		
		++$i;
		
		if (!defined $atom->{ELEMENT}) {
			silico_msg('w', "Atom $i in molecule \"$mol->{NAME}\" has no element!\n",
					"Skipping.\n");
			next;
		}
		
		next if $atom->{ELEMENT} eq 'Du';
		next if $atom->{ELEMENT} eq 'Lp';
		next if $atom->{ELEMENT} eq 'H';
		++$count;
	}

	$mol->{NUMHEAVY} = $count;
		
	return $count;
}



sub molecule_formula {

	#<
	#? Calculate the molecular formula of a molecule
	#. Missing hydrogens are calculated and they are added to the formula as 'H*'
	#. Sets $mol->{FORMULA};
	#; Requires: Molecule, flag to suppress no-hydrogens warning (optional), atomlist (optional), flag to supress count of missing hydrogens.
	#; Returns: Molecular formula as a hash
	#>
	
	my $mol = $_[0] || croak();
	my $nowarn = $_[1];
	my $list = $_[2]; # Atom list
	my $nomissing = $_[3]; # Do not calculate number of missing hydrogens on each atom
	
	my $formula;

	# The SDF noparse option calculates the formula already
	return $mol->{FORMULA} if  $mol->{SOURCE_FILE_TYPE} && $mol->{SOURCE_FILE_TYPE} eq 'sdf_noparse';
	
	@$list = atoms($mol) if !$list;
	
	mol_calc_hcount($mol, $list);
	
	foreach my $atom (@$list) {
	
		++$formula->{$atom->{ELEMENT}};
		
		if (!$nomissing && $atom->{HMISSING}) {
			$formula->{'H*'} += $atom->{HMISSING};
		}
	}
	
	if ($Silico::debug) {
		silico_msg('g', "Molecular formula of $mol->{NAME}: ");
		
		foreach (sort keys(%$formula)) {
			silico_msg('g', "$_ ($formula->{$_}) ");
		}
		
		silico_msg('g', "\n");
	}
	
	my $name;
	if (!$list) {
		$name = $mol->{NAME} || "Mol";
	} else {
		$name = $list->[0]{SUBNAME} || "Mol";
	}
	
	$name =~ s/^\s*//;
	$name =~ s/\s*$//;
	
	if ($formula->{'H*'} && !$nowarn) {
		silico_msg('w', "Molecule '$name' is missing $formula->{'H*'} H atoms. These are included in the formula as 'H*'\n");
	}
	
	$mol->{FORMULA} = $formula;
	
	return $formula;
}

sub formula_string {

	#<
	#? Convert molecular formula to a string
	#; Requires: formula hash
	#; Returns: string
	#>

        my $f = $_[0];

        my $s = '';

        foreach (sort keys %$f) {
                $s .= $_.$f->{$_} if $f->{$_};
        }
        return sprintf "%-12s", $s;
}

sub molecule_formula_string {

	#<
	#? Return molecule formula as a string
	#; Requires: molecule - see other options for molecule_formuula
	#; Returns: string
	#>
	
	return formula_string(molecule_formula(@_));
}


sub mol_count_rot_bonds {

	#< 
	#? Count rotatable bonds 
	#. Rotatable bonds are defined as being single, connected to more than
	#  one atom (excluding hydrogen), not
	#  amide bonds and not in a ring with less than 10 atoms.
	#. This is close to the definition used in Veber et al. J. Med. Chem. 2002, 45, 2615
	#; Requires: molecule
	#; Sets: $mol->{SDF_DATA}{ROT_BONDS};
	#; Returns: num rotatable bonds
	#>
	
	my $mol = $_[0];

	my $termlist;

	# Find rings.
	molecule_find_rings($mol, 10);
	
	# Create a hash of terminal atoms
	foreach my $atom (atoms($mol)) {
	
		# Count number of nohydrogen atoms
		my $nonh = 0;
		foreach my $connum (@{$atom->{CONNECT}}) {
			my $con = $mol->{ATOMS}[$connum];
			++$nonh if $con->{ELEMENT} ne 'H';
		}
		
		next if $nonh > 1;
		
		++$termlist->{$atom->{NUM}};
	}
	
	my $i = -1;
	my $numRot = 0;
	foreach my $atom (atoms($mol)) {
	
		++$i;
		
		# Don't count bonds to atoms only connected to one atom
		next if $#{$atom->{CONNECT}} <= 0;
		
		# Skip if a terminal atom (ie only attatched to H and one heavy)
		next if $termlist->{$atom->{NUM}};
		
		my $j = -1;
		CON: foreach my $connum (@{$atom->{CONNECT}}) {
		
			++$j;
			my $con = $mol->{ATOMS}[$connum];
			
			# Count bonds only once
			next CON if $connum > $i;
			
			# Skip if a terminal atom (ie only attatched to H and one heavy)
			next CON if $termlist->{$con->{NUM}};
			
			# Don't count bonds to monovalent atoms
			next if $con->{ELEMENT} eq 'H';
			next if $con->{ELEMENT} eq 'F';
			next if $con->{ELEMENT} eq 'Cl';
			next if $con->{ELEMENT} eq 'Br';
			next if $con->{ELEMENT} eq 'I';

			# Only count single bonds
			next CON if $atom->{BORDERS}[$j] != 1;

			# Don't count amides though
			next if ($atom->{FG}{AMIDE_C} && $con->{ELEMENT} eq 'N');
			next if ($con->{FG}{AMIDE_C} && $atom->{ELEMENT} eq 'N');
			
			# Find if con is in the same ring as atom 
			foreach my $ring (@{$atom->{RINGS}}) {
				
				foreach my $ringatom (@$ring) {
					next CON if $ringatom == $connum;
				}
			}
						
			# Count the bond
			++$numRot; # This is a rotatable bond
		}
	}

	$mol->{SDF_DATA}{ROT_BONDS} = $numRot;
	return $numRot;
}

sub simple_atomtype {

	#<
	#? Simple polar/nonpolar atom type routine
	#>

	my $mol = $_[0];
	my $atoms = $_[1];

	foreach my $atom (@$atoms) {

		if (!($atom->{ELEMENT} eq 'C' || $atom->{ELEMENT} eq 'H')) {

			$atom->{POLAR} = 1;
			next;
		}

		my $polar = 0;
		foreach my $connum (atoms($mol)) {
			
			my $con = $mol->{ATOMS}[$connum];

			next if $con->{ELEMENT} eq 'C' || $con->{ELEMENT} eq 'H';

			$polar = 1;
			last;
		}

		if ($polar) {
			$atom->{POLAR} = 1;
			next;
		}

		$atom->{POLAR} = 0;
	}
}

sub count_polar_nonpolar {

        #<
        #? Count polar (N, O) and nonpolar (all other) heavy atoms
        #; Requires: molecule
        #; Sets: $mol->{SDF_DATA}{POLAR} and $mol->{SDF_DATA}{NONPOLAR}
        #; Returns: 1 or 0 if non-CHONSP+Halogen atom is encountered.
        #>

        my $mol = $_[0];

        my $nonpolar = 0;
        my $polar = 0;

        foreach my $atom (atoms($mol)) {

                my $el = $atom->{ELEMENT};

                next if $el eq 'H';

                if (!($el eq 'C' || $el eq 'O' || $el eq 'N' || $el eq 'S' || $el eq 'P' || $el eq 'F' || $el eq 'Cl' || $el eq 'Br' || $el eq 'I')) {

                        silico_msg('w', "Found undesirable atomic element '$el' in molecule: $mol->{NAME}\n");
                        return 0;
                }

                if ($el eq 'O' || $el eq 'N') {
                        ++$polar;
                } else {
                        ++$nonpolar;
                }
        }

        $mol->{SDF_DATA}{POLAR} = $polar;
        $mol->{SDF_DATA}{NONPOLAR} = $nonpolar;
        return 1;
}


##################################################################
#
#	Functional group identification
#
##################################################################

sub mol_label_functional_group {

	#<
	#? General functional group assigner
	#. Label atoms of important functional groups and generate approximate charge
	#  groups for carboxylates, amino groups, sulphate groups, guanidines and guanidiniums
	#. Identifies the following functional group atoms with the following
	#  $atom->{FG}{<group>}
	#
	#,	AA_CA	AA_CTER	AA_N	AA_NTER	
	#
	#	ACETAL_O
	#
	#	ACID_HALIDE_C	ACID_HALIDE_O
	#
	#	ALCOHOL_C	
	#	ALDEHYDE_C	ALDEHYDE_O	
	#	ALKENE_CH0	ALKENE_CH1	ALKENE_CH2	
	#	ALKOXIDE_C	ALKOXIDE_O	ALKYL_PHOSPHATE_P	
	#	ALKYNE	ALKYNE_ALPHA_CH0	ALKYNE_ALPHA_CH1	ALKYNE_ALPHA_CH2	ALKYNE_ALPHA_CH3
	#	ALKYNE_CH1	
	#	ALPHA_GUANIDINE_C	ALPHA_GUANIDINIUM_C	
	#	AMIDE_C AMIDE_O	AMIDO_C	
	#	AMIDO_N	AMINO_C	AMMONIUM_C	
	#	AROM_NITRILE_C	AROM_NITRILE_N
	#,	AROMATIC
	#		
	#,	CARBAMATE_C	CARBAMATE_N	CARBAMATE_O	
	#	CARBANION_ALPHA_C	CARBANION_C	
	#	CARBONYL_ALPHA_C	CARBONYL_C	
	#	CARBOXYLATE_C	CARBOXYLATE_O	
	#	CARBOXYLIC_C	CARBOXYLIC_O	
	#	CHLOROALKENE_C	CHLOROALKENE_Cl	
	#	DIALKYL_PHOSPHATE_ESTER_O	DIALKYL_PHOSPHATE_P	
	#	DICHLOROALKENE_C	
	#	DIOL_C	
	#	DISULFIDE_S	
	#	DPE_O
	#,	ESTER_C	ESTER_O	
	#	ETHER_C	ETHER_O	
	#	FORMYL_C	FORMYL_O		
	#,	GUANIDINE_C	GUANIDINE_N	
	#	GUANIDINIUM_C	GUANIDINIUM_N		
	#	
	#,	IMIDE_C	IMIDE_N	IMIDE_O	
	#	
	#,	KETONE_C	KETONE_O	
	#	
	#,	METHYL_C	
	#
	#	MONOALKYL_PHOSPHATE_ESTER_O	MPE_O
	#	
	#,	NITRILE_ANION_ALPHA_C	NITRILE_ANION_C	NITRILE_ANION_N	NITRILE_C	
	#	NITRILE_N	NITRO_C	NITRO_N	NITRO_O	
	#	
	#,	PHOSPHATE_ESTER_O	PHOSPHATE_O	PHOSPHATE_P	PHOSPHONATE_C	PHOSPHONATE_ESTER_C	PHOSPHONATE_O	
	#	PHOSPHONATE_P	PHOSPHONIUM_C	PHOSPHONIUM_P	P_AMIDE_N	P_AMINO_C	
	#	P_AMINO_N	P_AMMONIUM_C	P_AMMONIUM_N	
	#	
	#,	Q_AMMONIUM_N	
	#	
	#,	SN_ALPHA_C	SULFIDE_C	SULFIDE_S	SULFONAMIDE_C	SULFONAMIDE_N	SULFONAMIDE_N_ALPHA_C	
	#	SULFONAMIDE_O	SULFONAMIDE_S	SULFONATE_C	SULFONATE_O	SULFONATE_S	SULFONE_O	
	#	SULFONE_S	SULFOXIDE_C	SULFOXIDE_O	SULFOXIDE_S	S_AMIDE_N	S_AMINO_C	
	#	S_AMINO_N	S_AMMONIUM_C	S_AMMONIUM_N	S_CARBAMATE_N	
	#	
	#,	TFM_C	TFM_F	THIOLATE_C	THIOLATE_S	THIOL_C	THIOL_S	TRIOL_C	T_AMIDE_C	
	#	T_AMIDE_N	T_AMINO_C	T_AMINO_N	T_AMMONIUM_C	T_AMMONIUM_N	T_CARBAMATE_N	
	#	
	#,	UREA_C	UREA_N	UREA_O
	#
	#. List of all functional groups are stored in mol->{FG}
	#. Charge groups are stored as an array of arrays in mol->{CHARGE_GROUP}
	#.
	#. Note 1: The routine molecule_find_planar_atoms should be run first
	#  although it will be run within this routine if it has not already
	#. Note 2: This routine requires that bond orders are correct
	# 
	#>
	
	my $mol = $_[0];

	my $starttime = (times)[0];

	#if ($mol->{HAS_NO_COORDS}) {
		#silico_msg('n', "Molecule has no coordinates.\n",
				#"Can not determine functional groups.\n");
		#return;
	#}

	# Delete any existing values
	# And ensure that atom->{NUM} is sane
	my $i = 1;
	foreach my $atom (atoms($mol)) {
		delete $atom->{FG};
		delete $atom->{FG_ALPHA};
		delete $atom->{FG_BETA};
		$atom->{NUM} = $i;
		++$i;
	}
	
	delete $mol->{FG};
	delete $mol->{CHARGE_GROUPS};
	
	# Check to see if molecule has hydrogens
	require silico_check;
	mol_has_h($mol) if !defined $mol->{HAS_HYDROGENS};
	
	# Find planar atoms
	molecule_find_planar_atoms($mol) if (!defined $mol->{PLANAR_ATOMS});
	
	# Find rings up to 6 members
	molecule_find_rings($mol, 6);

	# Label carbonyls so that it can be added to neighbouring atoms
	label_fg_carbonyl($mol);

	# For each atom, identify its neighbours
	mol_find_neighbouring_atoms($mol);

	# Label functional groups
	label_fg_planar_carbon($mol);
	label_fg_aromatic($mol);
	label_fg_tfm($mol);
	label_fg_carbon($mol);
	label_fg_nitrogen($mol);
	label_fg_benzhydrylamine($mol);
	label_fg_phosphorous($mol);
	label_fg_sulphur($mol);
	label_fg_oxygen($mol);
	label_fg_alcohol($mol);
	label_fg_carbanion($mol);
	label_fg_halogen($mol);
	
	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
}

sub fgroup {

	#< 
	#? Make a functional group
	#>

	my $mol = shift;
	my $name = shift;
	my @atoms = @_;

	carp("Error: Undefined functional group name\n") if !defined $name;

	foreach my $atom (@atoms) {

		my $el = $atom->{ELEMENT} || 'X';
		my $label = ($name|| "Unk")."_".$el;
		$atom->{FG}{$label} = 1;
	}
	
	my $group->{NAME} = $name;
	@{$group->{ATOMS}} = @atoms;
	
	push @{$mol->{FG}}, $group;
}

sub label_fg_aromatic {

	#<
	#? Create AROMATIC functional groups with rings up to 6 atoms
	#; Requires: molecule
	#>
	
	my $mol = $_[0];
	
	molecule_find_rings($mol, 6) if $mol->{RINGSIZE_MAX} || $mol->{RINGSIZE_MAX} < 6;
	
	foreach my $ring (@{$mol->{ARINGS}}) {
		my @ringatoms;
		next if $#ringatoms > 5;
		foreach my $a (@$ring) {
			push @ringatoms, $mol->{ATOMS}[$a];
		}
		fgroup($mol, "AROMATIC", @ringatoms);
	}
}

sub label_fg_carbonyl {

	#<
	#? Label carbonyl, thioketone and formyl groups
	#; Requires: molecule
	#>
	
	my $mol = $_[0];
	
	my $numdouble;
	
	foreach my $atom (atoms($mol)) {
	
		next if !($atom->{ELEMENT} eq 'O' || $atom->{ELEMENT} eq 'S');
		
		$numdouble = 0;
		foreach (@{$atom->{CONNECT}}) {
			++$numdouble if $_ == 2;
		}
		next if $numdouble > 1;

		# Carbonyl and formyl
		if ($atom->{ELEMENT} eq 'O' && $#{$atom->{CONNECT}} == 0 && $atom->{BORDERS}[0] == 2 ) {

			foreach (@{$atom->{CONNECT}}) {

				my $con = $mol->{ATOMS}[$_];
		
				next if $con->{ELEMENT} ne 'C';
				fgroup($mol, 'CARBONYL', $atom, $con);
		
				# Formyl groups
				if (ch($con->{CONNECTED},'H') >= 1 && ch($con->{CONNECTED},'C') == 0) {
					fgroup($mol, 'FORMYL', $atom, $con);
				}
			}
		}
		
		# Thioketone
		if ($atom->{ELEMENT} eq 'S' && $#{$atom->{CONNECT}} == 0 && $atom->{BORDERS}[0] == 2) {
			
			foreach (@{$atom->{CONNECT}}) {

				my $con = $mol->{ATOMS}[$_];
				next if $con->{ELEMENT} ne 'C';
				fgroup($mol, 'THIOKETONE', $atom, $con);
			}
		}
	}
}

sub label_fg_planar_carbon {

	#<
	#? Label functional groups containing planar carbon atoms such as		
	# carboxylate / carboxylic acids and guanidine / guanidinium groups
	#. Note. Amides are labelled by label_fg_nitrogen
	#
	#.  Carboxylate carbon if bonded to two oxygen atoms that 
	#  have only one bonded atom each. Carbonyl oxygen
	#  Ester if connected to double bonded O and single bonded
	#  O with two carbon attachments if double bond to oxygen
	#>
	
	my $mol = $_[0];
	
	my @alist;
	
	PLANAR_CARBON: foreach my $atom (atoms($mol)) {
		
		next PLANAR_CARBON if !$atom->{PLANAR_ATOM};
		next PLANAR_CARBON if ($atom->{ELEMENT} ne 'C');

		# Charge groups (here) are obsolete and should be removed
		my $charge_group;

		my $con_o2 = 0; # Number of connected oxygens that are not bonded to any other atom
		my $i = 0;
		my $possible_acid1 = 0; # Oxygen with 1 hydrogen attached
		my $possible_acid2 = 0; # Oxygen with double bond
		
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			my $bo = $atom->{BORDERS}[$i];
			
			# Oxygen with only one connected atom
			++$con_o2 if ($con->{ELEMENT} eq 'O' && ($#{$con->{CONNECT}} == 0)); 
			
			# Acid 
			if ($con->{ELEMENT} eq 'O' && $bo == 1) {
				
				my $ccount = 0;
				my $hcount = 0;
				foreach my $connum1 (@{$con->{CONNECT}}) {
					my $con1 = $mol->{ATOMS}[$connum1];
					++$ccount if $con1->{ELEMENT} eq 'C';
					++$hcount if $con1->{ELEMENT} eq 'H';
				}
				
				$possible_acid1 = 1 if $hcount == 1;	
			}
			
			if ($con->{ELEMENT} eq 'O' && $bo == 2) {
				$possible_acid2 = 1;
			}

			++$i;
		}

		# Carboxylic acids/caboxylates
		# Atom is connected to two oxygen atoms and 1 carbon or two oxygens and no other heavy atoms (i.e. formic acid)
		if (ch($atom->{CONNECTED},'O') == 2 && ch($atom->{CONNECTED},'C') == 1 || ch($atom->{CONNECTED},'HEAVY') == 2) {	
			
			#print "-----\n";
			#print "pa1$ possible_acid1 pa2 $possible_acid2 con02 $con_o2\n";
			
			# Carboxylic acids
			if ($possible_acid1 == 1 && $possible_acid2 == 1) {
			
				@alist = ();
				foreach my $connum (@{$atom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$connum];
					push @alist, $con if $con->{ELEMENT} eq 'O';
				}
				fgroup ($mol, 'CARBOXYLIC', $atom, @alist);
				
				#print "here1\n";
				#atom_printout($atom);
				next PLANAR_CARBON;
			}

			# Carboxylate
			if ($con_o2 == 2 && $possible_acid2) {
			
				push @$charge_group, -1;
				push @$charge_group, $atom;
				
				@alist = ();
				foreach my $connum (@{$atom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$connum];
					next if ($con->{ELEMENT} ne 'O');
					push @alist, $con;
					push @$charge_group, $con;
				}
				push @{$mol->{CHARGE_GROUPS}}, $charge_group;
				fgroup ($mol, 'CARBOXYLATE', $atom, @alist);
				
				#atom_printout($atom);
				next PLANAR_CARBON;
			}
		}
		
		# Thioarboxylic acids/thiocaboxylates
		# Atom is connected to one oxygen atoms, 1 carbon and 1 sulfur
		if (ch($atom->{CONNECTED},'O') == 1 && ch($atom->{CONNECTED},'S') == 1 && ch ($atom->{CONNECTED},'C') == 1) {
		
			# Thiocarboxylic acids
			if ($possible_acid1 == 1 && $atom->{FG}{THIOKETONE_C}) {
		
				@alist = ();
				foreach my $connum (@{$atom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$connum];
					push @alist, $con if $con->{ELEMENT} eq 'O' || $con->{ELEMENT} eq 'S';
				}
				fgroup ($mol, 'THIOCARBOXYLIC', $atom, @alist);
				next PLANAR_CARBON;
			}

			# Thiocarboxylate groups
			# Two tautomers (only C(O)SH is observed according to Wikipedia)
			if ($con_o2 == 1 && $atom->{FG}{THIOKETONE_C}) {
			
				push @$charge_group, -1;
				push @$charge_group, $atom;
				
				@alist = ();
				foreach my $connum (@{$atom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$connum];
					next if ($con->{ELEMENT} ne 'O' && $con->{ELEMENT} ne 'S');
					push @alist, $con;
					push @$charge_group, $con;
				}
				push @{$mol->{CHARGE_GROUPS}}, $charge_group;
				fgroup ($mol, 'THIOCARBOXYLATE', $atom, @alist);
				next PLANAR_CARBON;
			}
		}
		
		# Ureas 
		if (ch($atom->{CONNECTED},'O') == 1 && $con_o2 == 1 && ch($atom->{CONNECTED},'N') == 2) {
			
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
				push @alist, $con if ($con->{ELEMENT} eq 'N');
			}
			fgroup ($mol, 'UREA', $atom, @alist);
			next PLANAR_CARBON;
		}
		
		# Carbamates
		if (ch($atom->{CONNECTED},'O') == 2 && $con_o2 >= 1 && ch($atom->{CONNECTED},'N') == 1) {
			
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
				push @alist, $con if ($con->{ELEMENT} eq 'N');
			}
			fgroup ($mol, 'CARBAMATE', $atom, @alist);
			next PLANAR_CARBON;
		}
		
		 # Carbonates
                if (ch($atom->{CONNECTED},'O') == 3 && $con_o2 >= 1) {

                        @alist = ();
                        foreach my $connum (@{$atom->{CONNECT}}) {
                                my $con = $mol->{ATOMS}[$connum];
                                push @alist, $con if ($con->{ELEMENT} eq 'O');
                        }
                        fgroup ($mol, 'CARBONATE', $atom, @alist);
                        next PLANAR_CARBON;
                }

		# Thioureas
		if ($atom->{FG}{THIOKETONE_C} && ch($atom->{CONNECTED},'N') == 2) {
	
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'S');
				push @alist, $con if ($con->{ELEMENT} eq 'N');
				
				if ($con->{ELEMENT} eq 'S') {
					
				}
			}
			
			fgroup ($mol, 'THIOUREA', $atom, @alist);
			next PLANAR_CARBON;
		}
		
		# Guanidines and Guanidiniums 
		# Provided we are not in an aromatic ring
		if (ch($atom->{CONNECTED},'N') == 3 && !$atom->{AROMATIC_RING}) {
		
			my $count = 0;
			foreach my $connum (@{$atom->{CONNECT}}) {
			
				my $con = $mol->{ATOMS}[$connum];
				$count += $#{$con->{CONNECT}}+1;
			}
			
			# If count = 8 we have hydrogens present and the unprotonated form
			# If count = 9 we have a guanidinium group
			# If count = anything else: we are not sure and we will assume a guanidinium group
			
			@alist = ();
			if ($count == 8) {
			
				# Neutral
				foreach my $connum (@{$atom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$connum];
					push @alist, $con;
					
				}
				fgroup($mol, 'GUANIDINE', $atom, @alist);
				
			} else {
			
				# Charged
				push @$charge_group, 1;
				push @$charge_group, $atom;
			
				foreach my $connum (@{$atom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$connum];
					push @alist, $con;
					push @$charge_group, $con;
					
				}
				fgroup($mol, 'GUANIDINIUM', $atom, @alist);
				push @{$mol->{CHARGE_GROUPS}}, $charge_group;	
			}
			
			next PLANAR_CARBON;
		}
		
		# Amidines (not in an aromatic ring)
		if (ch($atom->{CONNECTED},'N') == 2 && ch($atom->{CONNECTED},'C') == 1 && !$atom->{AROMATIC_RING}) {	
				
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				next if $con->{ELEMENT} eq 'C';
				push @alist, $con;
					
			}
			fgroup($mol, 'AMIDINE', $atom, @alist);
				
			next PLANAR_CARBON;
		}
		
		# Acid halides
		if (ch($atom->{CONNECTED},'O') == 1 && $con_o2 == 1 && 
		(ch($atom->{CONNECTED},'F') + ch($atom->{CONNECTED},'Cl') + ch($atom->{CONNECTED},'Br') + ch($atom->{CONNECTED},'I')) == 1 
		&& (ch($atom->{CONNECTED},'H') + ch($atom->{CONNECTED},'C')) == 1) {
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
			}
			fgroup($mol, 'ACID_HALIDE', $atom, @alist);
			next PLANAR_CARBON;
		}
		
		# Ketones
		if (ch($atom->{CONNECTED},'O') == 1 && $con_o2 == 1 && ch($atom->{CONNECTED},'C') == 2 && $atom->{FG}{CARBONYL_C}) {
		
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
				if ($con->{ELEMENT} eq 'C') {
					++$con->{FG}{ALPHA_KETONE_C};
				}
			}
			fgroup($mol, 'KETONE', $atom, @alist);
			next PLANAR_CARBON;
		}
		
		# Aldehydes
		if (ch($atom->{CONNECTED},'O') == 1 && $con_o2 == 1 && ch($atom->{CONNECTED},'H') >= 1 &&
		(ch($atom->{CONNECTED},'C') + ch($atom->{CONNECTED},'H')) == 2 && $atom->{FG}{CARBONYL_C}) {

			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
			}
			fgroup($mol, 'ALDEHYDE', $atom, @alist);
			next PLANAR_CARBON;
		}
	}
	
	return;
}



sub label_fg_carbanion {

	#<
	#? Label carbanions
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	# Carbanions
	CARBANION: foreach my $atom (atoms($mol)) {
	
		next CARBANION if ($atom->{ELEMENT} ne 'C');
		next if $atom->{FG}{GUANIDINIUM_C} ||$atom->{FG}{GUANIDINE_C};
		
		if (($atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE} == -1) ) {
		
			fgroup($mol, 'CARBANION', $atom);
		}
	}
}

sub label_fg_tfm {

	#<
	#? Label trifluoromethyl groups
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	my @alist;
	my $con;
		
	# Trifluoromethyl group
	TFM: foreach my $atom (atoms($mol)) {
		
		next TFM if ($atom->{ELEMENT} ne 'C');
		next TFM if (valence($atom) != 4);
		
		if (ch($atom->{CONNECTED},'F') == 3 && ch($atom->{CONNECTED},'C') == 1) {
			@alist = ();
			# Label attached fluorine atoms as TFM_F
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'F');
			}
			fgroup($mol, 'TFM', $atom, @alist);
		}
		next TFM;
	}
}
	
sub label_fg_carbon {

	#<
	#? Label carbon 'functional groups'
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	my @nlist;
	my $i;
	my $numbonds;
	my $numsingle;
	
	CARBON: foreach my $atom (atoms($mol)) {
	
		next if ($atom->{ELEMENT} ne 'C');
		
		$numsingle = 0;
		foreach my $bo (@{$atom->{BORDERS}}) {
			++$numsingle if $bo == 1;
		}
		
		$numbonds = $#{$atom->{CONNECT}} + 1;
		
		# All bonds are single
		if ($numsingle - $numbonds == 0) {
		
			if (ch($atom->{CONNECTED},'H') == 1 || ($mol->{HAS_HYDROGENS} == 0 && $numsingle == 3)) {
				fgroup($mol, 'METHYNE', $atom);
			} elsif (ch($atom->{CONNECTED},'H') == 2 || ($mol->{HAS_HYDROGENS} == 0 && $numsingle == 2)) {
				fgroup($mol, 'METHYLENE', $atom);
			} elsif (ch($atom->{CONNECTED},'H') == 3 || ($mol->{HAS_HYDROGENS} == 0 && $numsingle == 1)) {
				fgroup($mol, 'METHYL', $atom);
			} elsif (ch($atom->{CONNECTED},'H') == 0 || ($mol->{HAS_HYDROGENS} == 0 && $numsingle == 4)) {
				fgroup($mol, 'QUATERNARY', $atom);
			}
			
			next CARBON;
		}
		
		$i = -1;
		foreach my $con (connected($atom, $mol)) {
		
			++$i;
			next if $con->{ELEMENT} ne 'C';
			
			my $bo = $atom->{BORDERS}[$i];
			
			# Alkene
			if ( $bo == 2 && !$atom->{AROMATIC_RING}) {
		 
				next if $atom->{FG}{ALKENE_C};
				
				fgroup($mol, 'ALKENE', $con, $atom);
				
				# Mark alkenes with # of attached hydrogens
				$atom->{FG}{'ALKENE_CH'.ch($atom->{CONNECTED},'H')} = 1;
				$con->{FG}{'ALKENE_CH'.ch($con->{CONNECTED},'H')} = 1;
				
				# Chloroalkenes
				# Mark alkenes with # of attached chlorines
				if (ch($atom->{CONNECTED},'Cl') == 1) {
						$atom->{FG}{CHLOROALKENE_C} = 1;
				}
				if (ch($con->{CONNECTED},'Cl') == 1) {
					$con->{FG}{CHLOROALKENE_C} = 1;
				}
				if (ch($atom->{CONNECTED},'Cl') == 2) {
					$atom->{FG}{DICHLOROALKENE_C} = 1;
				}
			}
			
			# Alkyne
			if ($bo == 3) {
			
				next if $atom->{FG}{ALKYNE_C};
			 	if ($atom->{NUM} < $con->{NUM}) {
					fgroup($mol, 'ALKYNE', $con, $atom);
				
					$atom->{FG}{'ALKYNE_CH'.ch($atom->{CONNECTED},'H')} = 1;
					$con->{FG}{'ALKYNE_CH'.ch($con->{CONNECTED},'H')} = 1;
				}
			}
		}
		
		# Cyanates, Isocyanates, Thiocyanates and Isothiocyanates
		if ($numbonds == 2 && !$atom->{RINGS}[0] && ch($atom->{CONNECTED}, 'O') == 1 && ch($atom->{CONNECTED}, 'N') == 1) {
			
			@nlist = ();
			foreach my $con (connected($atom, $mol)) {
				push @nlist, $con if $con->{ELEMENT} eq 'N';
			}
			
			# Cyanate -OCN
			fgroup($mol, "CYANATE", $atom, connected($atom, $mol)) if ($#{$nlist[0]->{CONNECT}} == 0);
			
			# Isocyanate -NCO
			fgroup($mol, "ISOCYANATE", $atom, connected($atom, $mol)) if ($#{$nlist[0]->{CONNECT}} == 1);
		}
		
		if ($numbonds == 2 && !$atom->{RINGS}[0] && ch($atom->{CONNECTED}, 'S') == 1 && ch($atom->{CONNECTED}, 'N') == 1) {
		
			@nlist = ();
			foreach my $con (connected($atom, $mol)) {
				push @nlist, $con if $con->{ELEMENT} eq 'N';
			}
			
			# Thiocyanate -SCN
			fgroup($mol, "THIOCYANATE", $atom, connected($atom, $mol)) if ($#{$nlist[0]->{CONNECT}} == 0);
			
			# Isothiocyanate -NCS
			fgroup($mol, "ISOTHIOCYANATE", $atom, connected($atom, $mol)) if ($#{$nlist[0]->{CONNECT}} == 1);
		}
	}
}

sub label_fg_nitrogen {

	#<
	#? Label nitrogen containing functional groups
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	my @alist;
	my $flag;
	my $type;
	
	# Nitrogen containing functional groups
	NITROGEN: foreach my $atom (atoms($mol)) {
	
		next NITROGEN if ($atom->{ELEMENT} ne 'N');
		
		my $numbonds = $#{$atom->{CONNECT}} + 1;

		my $charge_group;
		my $con_o = 0;
		my $con_h = 0;
		my $con_c = 0;
		my $con_n = 0;
		my $con_n2 = 0;
		my $con_n3 = 0;
		my $numsingle = 0;
		my $numdouble = 0;
		my $numtriple = 0;
		my $numaromatic = 0;
		my $non_ch = 0;

		my $i = -1;
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			++$i;
			my $con = $mol->{ATOMS}[$connum];
			++$con_o if ($con->{ELEMENT} eq 'O');
			++$con_h if ($con->{ELEMENT} eq 'H');
			++$con_c if ($con->{ELEMENT} eq 'C');
			++$non_ch if ($con->{ELEMENT} ne 'C' && $con->{ELEMENT} ne 'H');
			
			my $bo = $atom->{BORDERS}[$i];

			if ($con->{ELEMENT} eq 'N') {
				++$con_n;
				++$con_n2 if $bo == 2;
				++$con_n3 if $bo == 3;
			}
			
			next if $con->{ELEMENT} eq 'H';
			# Count non-H bonds
			++$numsingle if $bo == 1;
			++$numdouble if $bo == 2;
			++$numtriple if $bo == 3;
			++$numaromatic if $bo == 4;
		}
		
		my $total_bo = $numsingle + 2*$numdouble + 3*$numtriple+1.5*$numaromatic; 
		
		# This section does not assume that hydrogens are present
		
		# Defined primary, secondary, tertiary  and quaternary nitrogens
		# as being those with 1, 2, 3, or non-hydrogen single bonds
		if ($numdouble == 0 && $numtriple==0 && $numaromatic == 0) {
			
			$type = "PRIMARY" 	if $numsingle == 1;
			$type = "SECONDARY" 	if $numsingle == 2;
			$type = "TERTIARY"	if $numsingle == 3;
			$type = "QUATERNARY"	if $numsingle == 4;
		}
		
		# If a nitrogen could not be identified as P, S, T or Q, then 
		# set its Type as the empty string.
		# This is so the substr() function in Amides & Carbamates (see
		# below) has something on which to operate.
		$type ||= '';

		# Amides, hydrazides, hydroxamic acids and esters
		# (see also hydrazine, azo, etc)
		if (ch($atom->{CONNECTED},'CO') == 1 && ch($atom->{CONNECTED},'CS') == 0) {
		
			my $name;
			my $acid = 0;
			my $carb;
		
			my @alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{FG}{CARBONYL_C}) {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {	
						push @alist, $con2 if $con2->{ELEMENT} eq 'O';
					}
					$carb = $con;
				}
				push @alist, $con if $con->{ELEMENT} eq 'N'  || $con->{ELEMENT} eq 'O';
				$acid = 1 if $con->{ELEMENT} eq 'O' && ch($con->{CONNECTED}, 'HEAVY') == 1;
			}
			
			# Exclude carbamates and thiocarbamates
			if (ch($carb->{CONNECTED}, "N") == 1 && ch($carb->{CONNECTED}, "O") == 1) {
			
				if ($con_o == 0 && $con_n == 0 && $numdouble == 0) {
					$name = 'AMIDE';
				} elsif ($con_o == 0 && $con_n == 1) {
					$name = 'HYDRAZIDE';
				} elsif ($con_o == 1 && $con_n == 0 && $acid) {
					$name = 'HYDROXAMIC_ACID';
				} elsif ($con_o == 1 && $con_n == 0) {
					$name = 'HYDROXAMIC_ESTER';
				} 

				if ($name) {
					fgroup ($mol, substr($type, 0, 1).'_'.$name, $atom, @alist);
					next NITROGEN;
				}
			}
		}
		
		# Imides (R-(C=O)-N-(C=O)-R)
		if (ch($atom->{CONNECTED},'CO') == 2) {
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{FG}{CARBONYL_C}) {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						push @alist, $con2 if $con2->{ELEMENT} eq 'O';
					}
				}
			}	
		
			fgroup ($mol, 'IMIDE', $atom, @alist);
			next NITROGEN;
		}
		
		# Monothioimides (R-(C=O)-N-(C=S)-R)
		if (ch($atom->{CONNECTED},'CO') == 1 && ch($atom->{CONNECTED},'CS') == 1) {
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if (ch($atom->{CONNECTED},'CO')|| ch($atom->{CONNECTED},'CS')) {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						push @alist, $con2 if $con2->{ELEMENT} eq 'O';
						push @alist, $con2 if $con2->{ELEMENT} eq 'S';
					}
				}
			}	
		
			fgroup ($mol, 'MONOTHIOIMIDE', $atom, @alist);
			
			next NITROGEN;
		}
		
		# Dithioimides (R-(C=S)-N-(C=S)-R)
		if (ch($atom->{CONNECTED},'CS') == 2) {
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if (ch($atom->{CONNECTED},'CS')) {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						push @alist, $con2 if $con2->{ELEMENT} eq 'S';
					}
				}
			}	
		
			fgroup ($mol, 'DITHIOIMIDE', $atom, @alist);
			
			next NITROGEN;
		}
		
		# Thioamides, thiocarbamates and thiohydroxamic acids
		if (ch($atom->{CONNECTED},'CS') == 1) {

			my $name;
			my @olist = ();
			my @alist = ();
			
			foreach my $con (connected($atom, $mol)) {
				
				if ($con->{FG}{THIOKETONE_C}) {
					foreach my $con2 (connected($con, $mol)) {
						push @alist, $con if $con2->{ELEMENT} eq 'S';
						push @olist, $con if $con2->{ELEMENT} eq 'O';
					}
				}
			}	

			if ($#alist == 0 && $con_o == 0 && $numdouble == 0) {
				$name = 'THIOAMIDE';	
			} elsif ($#alist == 0 && $con_o == 1){
				$name = 'THIOHYDROXAMIC_ACID';
			} elsif ($#olist == 0) {
				$name = 'THIOCARBAMATE';
			}

			fgroup ($mol, substr($type, 0, 1).'_'.$name, $atom, @alist, @olist) if $name;
			
			next NITROGEN;
		}
		
		# N-oxides (non protonated oxygen)
		if ($con_o == 1 && $total_bo == 4) {
		
			my @alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{ELEMENT} eq 'O' && $#{$con->{CONNECT}} == 0) {
					push @alist, $con;
				}
			}
			
			if ($alist[0]) {
				fgroup($mol, "N_OXIDE", $atom, @alist);
				next NITROGEN;
			}
		}
		
		# Guanidine
		if ($atom->{FG}{GUANIDINE_N}) {
			next NITROGEN;
		}
		
		# Guanidinium
		if ($atom->{FG}{GUANIDINIUM_N}) {
			next NITROGEN;
		}
		
		# Find connected nitrogens
		my $n2; # Connected nitrogen
		foreach my $con (connected($atom, $mol)) {
			$n2 = $con if ($con->{ELEMENT} eq 'N');
		}
		
		# Hydrazines, hydrazones, diazonium and azo groups
		$flag = 0;
		if ($n2 && ch($atom->{CONNECTED},'N') == 1  &&  ch($atom->{CONNECTED},'CO') == 0 && ch($n2->{CONNECTED},'CO') == 0 && ch($n2->{CONNECTED},'N') == 1 && !defined $atom->{PRINGS}[0] ) {
			
			next NITROGEN if $atom->{FG}{DIAZONIUM_N};
			next NITROGEN if $atom->{FG}{HYDRAZINE_N};
			next NITROGEN if $atom->{FG}{HYDRAZONE_N};
			next NITROGEN if $atom->{FG}{AZO_N};
			
			my $n2_numdouble = 0;
			
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{ELEMENT} eq 'N') {
				
					push @alist, $con;
					foreach (@{$con->{BORDERS}}) {
						++$n2_numdouble if $_ ==2;
					}
					
					foreach my $con2 (connected($con, $mol)) {
						next if $con2 == $atom;
						push @alist, $con2 if $con2->{ELEMENT} eq 'C';
					}
				}
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			
			my $bo = find_bondorder2($atom, $n2, $mol);
			
			if ($bo == 1 && $numdouble == 0 && $n2_numdouble == 0) {
				fgroup($mol, "HYDRAZINE", $atom, $atom, @alist);
			} elsif ($bo == 1 && ($numdouble == 1 || $n2_numdouble == 1)) {
				fgroup($mol, "HYDRAZONE", $atom, $atom, @alist);
			} elsif ($bo == 2 && ($numdouble == 1 && $n2_numdouble == 1)) {
				fgroup($mol, "AZO", $atom, $atom, @alist);
			} elsif ($bo == 3) {
				fgroup($mol, "DIAZONIUM", $atom, $atom, @alist);
			}
			
			next NITROGEN;
		}
		
		# Azide
		if (!$atom->{PLANAR_RING} && ch($atom->{CONNECTED},'N') == 2 && valence($atom) == 4) {
			fgroup($mol, "AZIDE", $atom, connected($atom, $mol));
			next NITROGEN;
		}	
		
		# Nitrile
		if (!$atom->{PLANAR_RING} && $numbonds == 1 && $numtriple == 1 && ch($atom->{CONNECTED},'C') == 1) {
			fgroup($mol, "NITRILE", $atom, $mol->{ATOMS}[$atom->{CONNECT}[0]]);
			next NITROGEN;
		}
		
		# Nitroamine
		if (!$atom->{PLANAR_RING} && ch($atom->{CONNECTED},'O') == 2 && ch($atom->{CONNECTED},'N') == 1 && $total_bo >= 4) {
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
				push @alist, $con if ($con->{ELEMENT} eq 'N');
			}
			fgroup($mol, "NITROAMINE", $atom, @alist);
			next NITROGEN;
		}
		
		# Nitro
		if (!$atom->{PLANAR_RING} && ch($atom->{CONNECTED},'O') == 2 && $total_bo >= 4) {
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if ($con->{ELEMENT} eq 'O');
			}
			fgroup($mol, "NITRO", $atom, @alist);
			next NITROGEN;
		}
		
		# Oxime & Nitroso
		if (!$atom->{PLANAR_RING} && ch($atom->{CONNECTED},'O') == 1 && $con_c == 1 && $total_bo == 3 && $numbonds == 2) {
		
			my $name;
			my $i = 0;
			foreach my $bo (@{$atom->{BORDERS}}) {
			
				my $con = $mol->{ATOMS}[$atom->{CONNECT}[$i]];
	
				if ($bo == 2) {
					$name = 'OXIME' if $con->{ELEMENT}  eq 'C';
					$name = 'NITROSO' if $con->{ELEMENT}  eq 'O';
				}
				++$i;
			}
			
			write_sdf($mol, "untyped_nitrogen") if !defined $name;
				
			fgroup($mol, $name, $atom, connected($atom, $mol));
			next NITROGEN;
		}
		
		# Nitrite
		if (!$atom->{PLANAR_RING} && ch($atom->{CONNECTED},'O') == 2 && $numbonds == 2 && $total_bo == 3) {
			fgroup($mol, "NITRITE", $atom, connected($atom, $mol));
			next NITROGEN;
		}
		
		# Heterocyclic amine
		if ($atom->{AROMATIC_RING} && (($mol->{HAS_HYDROGENS} > 0 && valence($atom) == 3 && ch($atom->{CONNECTED},'C') + ch($atom->{CONNECTED},'H') == 3) || (valence($atom) <= 3))) {
		
			$atom->{FG}{HET_AMINO_N} = 1;
		
			if (ch($atom->{CONNECTED},'HEAVY') == 2) {
				fgroup($mol, "S_HET_AMINE", $atom);
			} elsif (ch($atom->{CONNECTED},'HEAVY') == 3) {
				fgroup($mol, "T_HET_AMINE", $atom);
			} 
			next NITROGEN;
		}

		# All groups below can not be adjacent to a carbonyl or thiocarbonyl or non-carbon atom 
		next NITROGEN if ch($atom->{CONNECTED}, 'CO')  || ch($atom->{CONNECTED}, 'CS') || $non_ch;
		
		my $arom = 0;
		my $nname = '';
		foreach my $connum (@{$atom->{CONNECT}}) {
			my $con = $mol->{ATOMS}[$connum];
			++$arom if $con->{AROMATIC_RING};
		}
		if ($arom == 1)  {
			$nname = 'MONOARYL';
		} elsif ($arom == 2) {
			$nname = 'DIARYL';
		} elsif ($arom == 3) {
			$nname = 'TRIARYL';
		} elsif ($arom == 4) {
			$nname = 'TETRAARYL';
		}
		
		# Ammonium groups
		if (valence($atom) == 4) {
			
			push @$charge_group, 1;
			push @$charge_group, $atom;
			push @{$mol->{CHARGE_GROUPS}}, $charge_group;
			
			@alist = ();
			# label carbon attached to ammonium group as AMMONIUM_C
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				if ($con->{ELEMENT} eq 'C') {
					push @alist, $con;
				}
			}
			
			if (ch($atom->{CONNECTED},'C') == 1 && ch($atom->{CONNECTED},'H') == 3 && $numdouble == 0 && $numtriple == 0) {
				$nname = 'P_'.$nname;
			} elsif (ch($atom->{CONNECTED},'C') == 2 && ch($atom->{CONNECTED},'H') == 2) {
				$nname = 'S_'.$nname;
			} elsif (ch($atom->{CONNECTED},'C') == 3 && ch($atom->{CONNECTED},'H') == 1) {
				$nname = 'T_'.$nname;
			} elsif (ch($atom->{CONNECTED},'C')) {
				$nname = 'Q_'.$nname;
			}
			
			if ($nname) {
				$nname = $nname."AMMONIUM";
				fgroup($mol,$nname, $atom, @alist);
			}
			next NITROGEN;
		}

		# Neutral amine
		if (	!$atom->{FORMAL_CHARGE} && !$atom->{AROMATIC_RING} && (
			($mol->{HAS_HYDROGENS} > 0 && valence($atom) == 3 && ch($atom->{CONNECTED},'C') + ch($atom->{CONNECTED},'H') == 3) ||
			(valence($atom) <= 3))) {
		
			$atom->{FG}{AMINO_N} = 1;
			
			@alist = ();
			# label carbon attached to amino group as AMINO_C
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				if ($con->{ELEMENT} eq 'C') {
					push @alist, $con;
				}
			}
			
			if (ch($atom->{CONNECTED},'C') == 1 && ch($atom->{CONNECTED},'HEAVY') == 1  && $numdouble == 0 && $numtriple == 0) {
				$nname = 'P_'.$nname;
			} elsif (ch($atom->{CONNECTED},'C') == 2 && ch($atom->{CONNECTED},'HEAVY') == 2) {
				$nname = 'S_'.$nname;
			} elsif (ch($atom->{CONNECTED},'C') == 3 && ch($atom->{CONNECTED},'HEAVY') == 3) {
				$nname = 'T_'.$nname;
			} 
			
			if ($nname) {
				$nname = $nname."AMINE";
				fgroup($mol, $nname, $atom, @alist);
			}
			next NITROGEN;
		}
		
		# Amide ion (R2N-)
		if (($numbonds == 2 && $numsingle == 2 && $mol->{HAS_HYDROGENS}) || $atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE} == -1) {
		
			# Called "AMIDE_ION_N" to distinguish from "AMIDE_N"
			@alist = ();
			# label carbon attached to amino group as AMINO_C
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				if ($con->{ELEMENT} eq 'C') {
					push @alist, $con;
				}
			}
			fgroup($mol, "AMIDE_ION", $atom, @alist);
			
			next NITROGEN;
		}
	}
}

sub label_fg_benzhydrylamine {

	#<
	#? Label benzhydrylamines
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	my $con_aromatic;
	my $con_s_amide_n;
	
	BENZHYDRYLAMIDE: foreach my $atom (atoms($mol)) {
	
		next if $atom->{ELEMENT} ne 'C';
		next if defined $atom->{PRINGS}[0];
		
		$con_aromatic = 0;
		$con_s_amide_n = 0;
		
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			
			++$con_aromatic if defined $con->{PRINGS}[0];
			++$con_s_amide_n if defined $con->{FG}{S_AMIDE_N};
		}
		
		if ($con_aromatic == 2 && $con_s_amide_n == 1 && ch($atom->{CONNECTED},'H') == 1) {
			fgroup($mol, 'BENZHYDRYLAMIDE', $atom);
		}
	}
}	

sub label_fg_oxygen {

	#<
	#? Label oxygen containing functional groups excluding CARBONYL_O, CARBOXYLIC_O and CARBOXYLATE_O
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	my @alist;
	my $con_c;
	my $con_h;
	my @olist;
	my @clist;
	my $numbonds;
	my $numsingle;
		
	# Oxygen based functional groups
	OXYGEN: foreach my $atom (atoms($mol)) {
		
		next if ($atom->{ELEMENT} ne 'O');
				
		# Assigned elsewhere
		next if ($atom->{FG}{CARBONYL_O});
		next if ($atom->{FG}{CARBOXYLIC_O});
		next if ($atom->{FG}{CARBOXYLATE_O});
		next if ch($atom->{CONNECTED}, 'P');
		next if ch($atom->{CONNECTED}, 'S');
		next if ch($atom->{CONNECTED}, 'N');
				
		$numbonds = $#{$atom->{CONNECT}} + 1;
		$numsingle = 0;
		
		$con_c = 0;
		$con_h = 0;
		
		foreach my $con (connected($atom, $mol)) {
			++$con_c if ($con->{ELEMENT} eq 'C');
			++$con_h if ($con->{ELEMENT} eq 'H');
		}

		foreach (@{$atom->{BORDERS}}) {
			++$numsingle if $_ == 1;
		}
		
		# Esters
		my $flag = 0;
		if (ch($atom->{CONNECTED},'CO') == 1 && $con_c == 2) {
		
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{FG}{CARBONYL_C}) {
					
					# Set flag if carbonyl is not connected to C or H
					$flag = 1 if ch($con->{CONNECTED}, 'C') == 0 && ch($con->{CONNECTED}, 'H') == 0;
					
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						next if $con2 == $atom;
						push @alist, $con2 if $con2->{ELEMENT} eq 'O';
					}
				}
			}	
			
			# Add carbon attached to oxygen
			foreach my $con (connected($atom, $mol)) {
				next if $con->{FG}{CARBONYL_C};
				push @alist, $con;
				$atom->{FG}{ESTER_ALKOXY_O} = 1; # Label this for use by some programs
			}
			
			if (!$flag) {
				fgroup ($mol, 'ESTER', $atom, @alist);
				next OXYGEN;
			}	
		}
		
		# Anhydrides
		if (ch($atom->{CONNECTED},'CO') == 2) {
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{FG}{CARBONYL_C}) {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {

						next if $con2 == $atom;
						push @alist, $con2 if $con2->{ELEMENT} eq 'O';
					}
				}
			}	
		
			fgroup ($mol, 'ANHYDRIDE', $atom,@alist);
			#last;
			next OXYGEN;
		}
		
		# All groups below here can not be connected to a CO or CS
		next if ch($atom->{CONNECTED},'CO');
		next if ch($atom->{CONNECTED},'CS');
			
		# Alcohols
		# Modified so that hydrogens do not need to be present DKC 13/9/2011
		if (($con_h == 1 && $con_c == 1) || ($numbonds == 1 && $con_c == 1 && $mol->{HAS_HYDROGENS} == 0)) {
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			
			if ($alist[0]->{AROMATIC_RING}) {
				fgroup($mol, "ARYL_ALCOHOL", $atom, @alist);
			} else {
				fgroup($mol, "ALCOHOL", $atom, @alist);
			}
			next OXYGEN;
		}
		
		
		# Alkoxides (only if formal charge is present)
		if (!ch($atom->{CONNECTED},'CO') && $numbonds == 1 && $con_c == 1 && defined $atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE} == -1) {
			fgroup($mol, 'ALKOXIDE', $atom, connected($atom, $mol));
			next OXYGEN;
		}
		

		# Ethers, Acetals, Hemiacetals and Acetal esters
		if ($con_c == 2 && ch($atom->{CONNECTED}, 'CO') == 0) {

			next if $atom->{FG}{ACETAL_O};
			next if $atom->{FG}{HEMIACETAL_O};
			next if $atom->{FG}{ACETAL_ESTER_O};

			@alist = ();
			@olist = ();
			@clist = ();
			foreach my $con (connected($atom, $mol)) {
				push @alist, $con;
				foreach my $con2 (connected($con, $mol)) {
					next if $con2 == $atom;
					if ($con2->{ELEMENT} eq 'O') {
						push @olist, $con2;
						foreach my $con3 (connected($con2, $mol)) {
							next if $con3 == $con;
							push @clist, $con3 if $con3->{ELEMENT} eq 'C';
						}
					}
				}
			}

			# Ether
			if ($#olist == -1 && ch($atom->{CONNECTED}, 'C') == 2 && $numsingle == 2) {
				fgroup($mol, 'ETHER', $atom, @alist);
				next OXYGEN;
			}

			# Acetal
			if ($#olist == 0 && ch($olist[0]->{CONNECTED}, 'HEAVY') == 2 && ch($olist[0]->{CONNECTED}, 'CO') == 0) {
			
				fgroup($mol, 'ACETAL', $atom, $olist[0], @alist, @clist);
				next OXYGEN;
			}

			# Hemicetal
			if ($#olist == 0 && ch($olist[0]->{CONNECTED}, 'HEAVY') == 1 && ch($olist[0]->{CONNECTED}, 'CO') == 0) {
			
				fgroup($mol, 'HEMIACETAL', $atom, $olist[0],  @alist, @clist);
				next OXYGEN;
			}

			# Acetal ester
			if ($#olist == 0 && ch($olist[0]->{CONNECTED}, 'CO') == 1) {
				fgroup($mol, 'ACETAL_ESTER', $atom, $olist[0],  @alist, @clist);
				next OXYGEN;
			}	
		}

		
		# Peroxides
		if (ch($atom->{CONNECTED},'O') == 1 && $numbonds == 2 && ch($atom->{CONNECTED}, 'CO') == 0) {
		
			next OXYGEN if $atom->{FG}{PEROXIDE_O};
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{ELEMENT} eq 'O') {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						next if $con2 == $atom;
						push @alist, $con2 if $con2->{ELEMENT} eq 'C';
					}
				}
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			fgroup($mol, "PEROXIDE", $atom, @alist);
			next OXYGEN; 
		}
	}
}

sub label_fg_alcohol {

	#<
	#? Label alcohols, diols and triols
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	# Determination of diol and triol carbons
	# Diol parameters should be used for sugars (Damm et al. J. Com. Chem. p1955, 1997)
	ALCOHOL: foreach my $atom (atoms($mol)) {
	
		my $con_alcohol_c = 0;

		next if (!$atom->{FG}{ALCOHOL_C});

		foreach my $connum (@{$atom->{CONNECT}}) {

			my $con = $mol->{ATOMS}[$connum];

			if ($con->{FG}{ALCOHOL_C}) {
				++$con_alcohol_c;
			}
		}

		# Make assumption that we use diol parameters if we are in a non-aromatic ring.
		# Note that this is flawed: it should only pick up carbohydrates, but will
		# currently pick up rings such as cyclohexanol as well.
		if ($con_alcohol_c == 1 || ($atom->{RINGS}[0] && !defined $atom->{PRINGS}[0])) {
			$atom->{FG}{DIOL_C} = 1;
			next ALCOHOL;
		}

		if ($con_alcohol_c == 2) {
			$atom->{FG}{TRIOL_C} = 1;
			next ALCOHOL;
		}
	}
		
}

sub label_fg_sulphur {

	#<
	#? Label sulphur containing functional groups
	#; Requires: molecule
	#>

	my $mol = $_[0];

	my @alist;
	my @clist;
	my @hlist;
	my @olist;
	
	# Sulphur compounds
	# Note THIOKETONE groups have been identified previously
	
	SULPHUR: foreach my $atom (atoms($mol)) {

		next SULPHUR if ($atom->{ELEMENT} ne 'S');	
		
		my $con_carbonyl_c = 0;
		my $con_c = 0;
		my $con_h = 0;
		my $con_n = 0;
		my $con_o = 0;
		my $count_o_connect = 0;
		my $numsingle = 0;
		
		my $i = 0;
		foreach my $connum (@{$atom->{CONNECT}}) {
			my $con = $mol->{ATOMS}[$connum];
			my $bo = $atom->{BORDERS}[$i];
			++$con_c if ($con->{ELEMENT} eq 'C');
			++$con_h if ($con->{ELEMENT} eq 'H');
			++$con_n if ($con->{ELEMENT} eq 'N');
			++$con_o if ($con->{ELEMENT} eq 'O');
			++$con_carbonyl_c if ($con->{FG}{CARBONYL_C});
			++$numsingle if $bo == 1;
			++$i;
		}


		# Thiolates
		my $numbonds = $#{$atom->{CONNECT}} + 1;
		if ($numbonds == 1 && $numsingle == 1 && ch($atom->{CONNECTED},'C') == 1 && $mol->{HAS_HYDROGENS} > 0 && $con_carbonyl_c == 0) {
			
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			fgroup($mol, "THIOLATE", $atom, @alist);
			next SULPHUR;
		}
		
		# Sulfides
		if (ch($atom->{CONNECTED},'C') == 2 && $numbonds == 2 && !$atom->{PRINGS}[0] && $con_carbonyl_c == 0 && $numsingle ==2) {
		
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			fgroup($mol, "SULFIDE", $atom, @alist);
			next SULPHUR;
		}
		
		# Thiols
		if (((ch($atom->{CONNECTED},'H') == 1 && ch($atom->{CONNECTED},'C') == 1 && $numbonds == 2) || 
			($mol->{HAS_HYDROGENS} <= 0 && $numsingle == 1 && ch($atom->{CONNECTED},'C') == 1)) && $con_carbonyl_c == 0) { 
			
			@alist = ();
			foreach my $connum (@{$atom->{CONNECT}}) {
				my $con = $mol->{ATOMS}[$connum];
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			fgroup($mol, "THIOL", $atom, @alist);
			next SULPHUR;
		}
		
		# Disulfides
		if (ch($atom->{CONNECTED},'C') == 1 && ch($atom->{CONNECTED},'S') == 1 && $numbonds == 2 && $con_carbonyl_c == 0) {
		
			next SULPHUR if $atom->{FG}{DISULFIDE_S};
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{ELEMENT} eq 'S') {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						next if $con2 == $atom;
						push @alist, $con2 if $con2->{ELEMENT} eq 'C';
					}
				}
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			fgroup($mol, "DISULFIDE", $atom, @alist);
			next SULPHUR;
		}

		# Sulfenamides CSNR2
		if (ch($atom->{CONNECTED},'N') == 1 && ch($atom->{CONNECTED},'C') == 1 && $numbonds == 2 && $con_carbonyl_c == 0) {
			
			@alist = ();
			foreach my $con (connected($atom, $mol)) {
				if ($con->{ELEMENT} eq 'N') {
					push @alist, $con;
					foreach my $con2 (connected($con, $mol)) {
						next if $con2 == $atom;
						push @alist, $con2 if $con2->{ELEMENT} eq 'C';
					}
				}
				push @alist, $con if $con->{ELEMENT} eq 'C';
			}
			fgroup($mol, "SULFENAMIDE", $atom, @alist);
			next SULPHUR;
		}
	
		# All compounds below are some sort of sulfur oxide 
		next if ch($atom->{CONNECTED},'O') == 0;
				
		# Find carbons attached to oxygen
		@olist = ();
		@clist = ();
		@hlist = ();
		foreach my $con (connected($atom, $mol)) {
			next if $con->{ELEMENT} ne 'O';
			push @olist, $con;
			foreach my $con2 (connected($con, $mol)) {
				next if $con2 == $atom;
				++$count_o_connect;
				push @clist, $con2 if $con2->{ELEMENT} eq 'C';
				push @hlist, $con2 if $con2->{ELEMENT} eq 'H';
			}
		}
	
		# Sulfoxides CS(O)C
		if (ch($atom->{CONNECTED},'O') == 1 && ch($atom->{CONNECTED},'C') == 2 && $numbonds == 3) {
			fgroup($mol, 'SULFOXIDE', $atom, connected($atom, $mol));
			next SULPHUR;
		}
		
		# Sulfones CS(O)(O)C
		if (ch($atom->{CONNECTED},'O') == 2 && ch($atom->{CONNECTED},'C') == 2 && $numbonds == 4) {
			fgroup($mol, 'SULFONE', $atom, connected($atom, $mol));
			next SULPHUR;
		}
	}
}


sub label_fg_phosphorous {

	#<
	#? Label phosphorous containing functional groups
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	my @clist;
	my @hlist;
	my @olist;
	
	# Phosphorus compounds
	PHOSPHORUS: foreach my $atom (atoms($mol)) {
	
		next PHOSPHORUS if ($atom->{ELEMENT} ne 'P');
		
		my $numsingle = 0;
		my $numdouble = 0;
		foreach my $bo (@{$atom->{BORDERS}}) {
			++$numdouble if $bo == 2;
			++$numsingle if $bo == 1;
		}
		
		my $numbonds = $#{$atom->{CONNECT}} + 1;
		
		# Phosphines
		if ($numbonds <= 3 && 2*$numdouble+$numsingle <=3 && ch($atom->{CONNECTED},'C') <= 3 && ch($atom->{CONNECTED},'O') == 0 && ch($atom->{CONNECTED},'N') == 0 ){
			fgroup($mol, 'PHOSPHINE', $atom, connected($atom, $mol));
			next PHOSPHORUS;
		}
		
		# Phosphorane
		if (2*$numdouble+$numsingle ==5 && ch($atom->{CONNECTED},'HEAVY') <= 5 && ch($atom->{CONNECTED},'O') == 0) {
			fgroup($mol, 'PHOSPHORANE', $atom, connected($atom, $mol));
			next PHOSPHORUS;
		}
		
		# Phosphine oxide
		if ($numbonds == 4 && 2*$numdouble+$numsingle == 5 && ch($atom->{CONNECTED},'C') <= 3 && ch($atom->{CONNECTED},'O') == 1) {
			fgroup($mol, 'PHOSPHINE_OXIDE', $atom, connected($atom, $mol));
			next PHOSPHORUS;
		}
		
		
		# Phosphinic acids and esters
		if ($numbonds == 4 && ch($atom->{CONNECTED},'O') == 2 && ch($atom->{CONNECTED},'C') == 2) {
			fgroup($mol, 'PHOSPHINIC', $atom, connected($atom, $mol));
			next PHOSPHORUS;
		}
		
		# Phosphonic acids and esters
		if ($numbonds == 4 && ch($atom->{CONNECTED},'O') == 3 && ch($atom->{CONNECTED},'C') == 1) {
			fgroup($mol, 'PHOSPHONIC', $atom, connected($atom, $mol));
			next PHOSPHORUS;
		}
		
		# Find carbons attached to oxygen
		@olist = ();
		@clist = ();
		@hlist = ();
		my $count_o_connect;
		foreach my $con (connected($atom, $mol)) {
			next if $con->{ELEMENT} ne 'O';
			push @olist, $con;
			foreach my $con2 (connected($con, $mol)) {
				next if $con2 == $atom;
				++$count_o_connect;
				push @clist, $con2 if $con2->{ELEMENT} eq 'C';
				push @hlist, $con2 if $con2->{ELEMENT} eq 'H';
			}
		}
		
		# Phosphates and esters
		# Note these used to be MONOALKYL_PHOSPHATE AND DIALKYL_PHOSPHATEs
		if ($numbonds == 4) {
			
			if ($#clist == 1 && ch($atom->{CONNECTED},'O') == 4) {
				fgroup($mol, 'DIALKYL_PHOSPHATE', $atom, connected($atom, $mol), @clist);
			} elsif ($#clist == 0 && ch($atom->{CONNECTED},'O') == 4) {
				fgroup($mol, 'MONOALKYL_PHOSPHATE', $atom, connected($atom, $mol), @clist);
			} elsif (ch($atom->{CONNECTED},'O') >= 3) {
				fgroup($mol, 'PHOSPHATE', $atom, connected($atom, $mol));
			} 
			next PHOSPHORUS;
		}
	}
}

sub label_fg_halogen {

	#<
	#? Label halogen atoms
	#; Requires: molecule
	#>

	my $mol = $_[0];
	
	# Halogens
	HALOGEN: foreach my $atom (atoms($mol)) {
	
		my $el = $atom->{ELEMENT};
		
		#molecule_printout($mol, 'all') && atom_printout($atom) && die if !defined $el;
	
		next HALOGEN if !($el eq 'F' || $el eq 'Cl'  || $el eq 'Br' || $el eq 'I');
		
		# Ions
		if ($#{$atom->{CONNECT}} == -1) {
			fgroup($mol, 'HALIDE', $atom);
			next HALOGEN;
		}
		# Halogens in higher oxidation states
		if ($#{$atom->{CONNECT}} > 0) {
			fgroup($mol, 'HALOGEN_OXIDISED', $atom);
			next HALOGEN;
		}
		
		my $con_carbonyl_c = 0;
		my $con;
		my $con_c = 0;
		my $con_h = 0;
		my $con_n = 0;
		my $con_o = 0;
		my $con_s = 0;
		my $numsingle = 0;
		my $numdouble = 0;
		my $numtriple = 0;
		my $aromatic = 0;
		
		my $i = 0;
		foreach my $connum (@{$atom->{CONNECT}}) {
			my $con = $mol->{ATOMS}[$connum];
			my $bo = $atom->{BORDERS}[$i];
			
			if ($con->{ELEMENT} eq 'C') {
				++$con_c 
			} elsif ($con->{ELEMENT} eq 'H'){
				++$con_h
			} elsif ($con->{ELEMENT} eq 'N'){
				++$con_n 
			} elsif ($con->{ELEMENT} eq 'O'){
				++$con_o 
			} elsif ($con->{ELEMENT} eq 'S'){
				++$con_s 
			} else {
				#Non CHONS halogen
				fgroup($mol, "$con->{ELEMENT}_HALIDE", $atom);
				next HALOGEN;
			}
			
			++$con_carbonyl_c if ($con->{FG}{CARBONYL_C});
			++$numsingle if $bo == 1;
			++$numdouble if $bo == 2;
			++$numtriple if $bo == 3;
			++$aromatic if $con->{PLANAR_RING};
			++$i;
		}
		
		if ($con_h == 1) {
			fgroup($mol, 'HALOGEN_ACID', $atom);
			next HALOGEN;
		}
		if ($con_carbonyl_c == 1) {
			# Should already be labelled as an acid halide
			next HALOGEN;
		}
		if ($aromatic) {
			fgroup($mol, 'HALOGEN_AROMATIC', $atom);
			next HALOGEN;
		}
		if ($numtriple == 1) {
			fgroup($mol, 'ALKYNYL_HALOGEN', $atom);
			next HALOGEN;
		}
		if ($numdouble == 1) {
			fgroup($mol, 'ALKENYL_HALOGEN', $atom);
			next HALOGEN;
		}
		if ($con_c && $numtriple == 0 && $numdouble == 0) {
			
			my $type = "PRIMARY" 	if $numsingle == 1;
			$type = "SECONDARY" 	if $numsingle == 2;
			$type = "TERTIARY"	if $numsingle == 3;
			
			fgroup($mol, "$type\_ALKYL_HALOGEN", $atom);
			next HALOGEN;
		}
		if ($con_n == 1) {
			fgroup($mol, 'NITROGEN_HALIDE', $atom);
			next HALOGEN;
		}
		if ($con_o == 1) {
			fgroup($mol, 'OXYGEN_HALIDE', $atom);
			next HALOGEN;
		}
		if ($con_s) {
			fgroup($mol, 'SULFUR_HALIDE', $atom);
			next HALOGEN;
		}
	}
}

sub mol_find_neighbouring_atoms {

	#<
	#? Make a list of connected atom properties in the 'CONNECTED' hash for each atom
	#. Labels counts heteroatoms as 'HET'.  Counts nonhydrogens as 'HEAVY'. Carbonyls as 'CO'.
	#  Thioketo/nes as 'CS'
	#; Requires: Molecule, optional list of atoms
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $list = $_[1] || $mol->{ATOMS};
	
	foreach my $atom (@$list) {
	
		my $hash;
		%$hash = ();
	
		foreach my $connum (@{$atom->{CONNECT}}) {

			my $con = $mol->{ATOMS}[$connum];
			my $el = $con->{ELEMENT};
			
			++$hash->{$el};

			if ($el ne 'H') {
				++$hash->{HEAVY};
				++$hash->{HET} if $el ne 'C';
				++$hash->{CO}  if $con->{FG}{CARBONYL_C};
				++$hash->{CS}  if $con->{FG}{THIOKETONE_C};
			}
		}
		
		foreach my $bo (@{$atom->{BORDERS}}) {
			++$hash->{SP2} if $bo == 2 || $bo == 4;
			++$hash->{SP}  if $bo == 3;
		}

		# Replace any exitsting values
		$atom->{CONNECTED} = $hash;
	}
}

###################################################################
#
#   Heterocycles
#
##################################################################


sub mol_label_heterocycles {

	#<
	#? Label atoms in heterocyclic rings
	#. Requires that find_rings has been run and that bond orders are present
	#  and that aromatic rings are represented using aromatic bonds (ie
	#  make_aromatic_bonds has been run).  List of heterocycles is
	#  also put into $mol->{HETEROCYCLES}
	#; Requires: molecule, maximum ring size (default 8), optional prefix to add to heterocycle name 'eg HET'
	#; Returns: List of heterocycles in this molecule
	#>

	my $mol = $_[0];
	my $max_ringsize = $_[1] || 8;
	my $prefix = $_[2];

	my @alist;
	my $list;
	
	make_aromatic_bonds($mol);
	molecule_find_rings($mol, $max_ringsize) if $max_ringsize != 6;
	
	foreach my $ring (@{$mol->{RINGS}}) {
	
		@alist = ();
		
		my $label = identify_heterocycle($ring, $mol);
		$label = $prefix.$label if $prefix;

		next if !defined $label;
		
		push (@{$list->{$label}}, $ring);
			
		foreach my $atomnum (@$ring) {
			push @alist, $mol->{ATOMS}[$atomnum];
		}
		
		fgroup($mol, $label, @alist);
	}
	
	$mol->{HETEROCYCLES} = $list;
	
	return $list;		
}

sub identify_heterocycle {

	#<
	#? Identify a heterocycle
	#. This routine can name a number of common heterocycles.
	#. Requires that find_rings has been run and that bond orders are present
	#  and that aromatic rings are represented using aromatic bonds (ie
	#  make_aromatic_bonds has been run)

	#; Requires; ring, molecule
	#; Returns; heterocycle label or HETEROCYCLE <standardised_string> if unknown;
	#>

	my $ring = $_[0];
	my $mol = $_[1];
	
	my %heterocycles;
	my $hetname;
	my $std;
	
	(undef, $std) = standardise_ring_by_heteroatoms($ring, $mol);

	%heterocycles = qw/	
	
		C-C-C-			CYCLOPROPANE
		C=C-C-			CYCLOPROPENE
		N-C-C-			AZIRIDINE
		O-C-C-			EPOXIDE

		C-C-C-C			CYCLOBUTANE
		C=C-C-C			CYCLOBUTENE
		N-C-C-C-		AZETIDINE

		C-C-C-C-C-		CYCLOPENTANE
		C=C-C-C-C-		CYCLOPENTENE
		C=C-C=C-C-		CYCLOPENTADIENE
		
		N-C=C-C-C-		PYROLLINE-2
		N-C-C=C-C-		PYROLLINE-3
		N-C-C-C-C-		PYRROLIDINE
		
		N=C-C-C-C-		PYRROLINE-1
		N-C=C-C-C-		PYRROLINE-2
		N-C-C=C-C-		PYRROLINE-3
		
		O-O-C-C-C-		DIOXOLANE-12
		O-C-O-C-C-		DIOXOLANE-13
		
		O-C=N-C=C-		OXAZOLE
		O-N=C-C=C-		ISOXAZOLE
		O-N=N-C=C-		OXADIAZOLE-123
		
		S-C=N-C=C-		THIAZOLE
		S-N=C-C=C-		ISOTHIAZOLE
		S-C=N-N=C-		THIADIAZOLE-134

		N-C=C-C=C-		PYRROLE
		N-C-C=C-C=		PYRROLE-2H
		N-N=C-C=C-		PYRAZOLE
		N-N=C-C-C-		PYRAZOLINE-2
		N-N-C-C-C-		PYRAZOLIDINE
		N=C-N-C=C-		IMIDAZOLE
		N-C=N-C-C-		IMIDAZOLINE-2
		N-C-N-C-C-		IMIDAZOLIDINE
		N-C-N-C-C-		TETRAHYDROIMIDAZOLE
		
		N=N-N-C=C-		TRIAZOLE-123
		N-N=C-N=C-		TRIAZOLE-124
		N=N-N-N=C-		TETRAZOLE
		N-N=N-N=C-		TETRAZOLE
		
		S-C=C-C=C-		THIOPHENE
		
		O-C=C-C=C-		FURAN
		O-C=C-C-C-		DIHYDROFURAN-23
		O-C-C=C-C-		DIHYDROFURAN-25	
		O-C-C-C-C-		TETRAHYDROFURAN
		
		C-C-C-C-C-C-		CYCLOHEXANE
		C=C-C-C-C-C-		CYCLOHEXENE
		C=C-C=C-C-C-		CYCLOHEXADIENE-13
		C=C-C-C=C-C-		CYCLOHEXADIENE-14
		O-C=C-C=C-C-		PYRAN-2H
		O-C=C-C-C=C-		PYRAN-4H
		O-C=C-C-C-C-		DIHYDRO-34-PYRAN-2H
		O-C-C-C-C-C-		TETRAHYDROPYRAN
		O-O-C-C-C-C-		DIOXANE-12
		O-C-O-C-C-C-		DIOXANE-13
		O-C-C-O-C-C-		DIOXANE-14	

		S-C=C-C=C-		THIOPHENE	

		N-C-C-C-C-C-		PIPERIDINE
		N-C-C-N-C-C-		PIPERAZINE
		
		O-C-C-N-C-C-		MORPHOLINE
		S-C-C-N-C-C-		THIOMORPHOLINE
		S-C-C-S-C-C-		DITHIANE-14
		S-C-S-C-S-C-		TRITHIANE-135
		
		C=C=C=C=C=C=		BENZENE

		N=C=C=C=C=C=		PYRIDINE
		N=N=C=C=C=C=		PYRIDAZINE
		N=C=N=C=C=C=		PYRIMIDINE
		N=C=C=N=C=C=		PYRAZINE	
		N=N=N=C=C=C=		TRIAZINE-123
		N=N=C=N=C=C=		TRIAZINE-124
		N=C=N=C=N=C=		TRIAZINE-135
		

		C-C-C-C-C-C-C-		CYCLOHEPTANE
		C=C-C-C-C-C-C-		CYCLOHEPTENE
		N-C=C-C=C-C=C-		AZEPINE

		C-C-C-C-C-C-C-		CYCLOOCTANE
		C=C-C-C-C-C-C-		CYCLOOCTENE

		N.-C=C-C=C.-C=C-C=C-	INDOLIZINE
		N-C.=C=C=C=C=C.-C=C-	INDOLE
		N-C=C.-C=C-C=C-C.=C-	ISOINDOLE
		N-C.=C=C=C=C=C.-C-C=	INDOLE-3H
		N-C.=C=C=C=C=C.-C-C-	INDOLINE
		
		O-C.=C=C=C=C=C.-C=C-	BENZObFURAN
		S-C.=C=C=C=C=C.-C=C-	BENZObTHIOPHENE
	
		N-N-C.=C=C=C=C=C.-C=	INDAZOLE-1H
		N=C-N-C.=C=C=C=C=C.-	BENZIMIDAZOLE
		S-C=N-C.=C=C=C=C=C.-	BENZTHIAZOLE
		
		N=C=N=C.-N=C-N-C.=C=	PURINE
		
		N.-C=C-C=C-C.=C-C=C-C-	QUINOLIZINE-4H
		N=C.-C=C-C=C-C.=C=C=C=	QUINOLINE
		N=C=C.-C=C-C=C-C.=C=C=	ISOQUINOLINE
		N=N=C.=C=C=C=C=C.=C=C=	CINNOLINE
		N=N=C=C.-C=C-C=C-C.=C=	PHTHALAZINE
		N=C=N=C.-C=C-C=C-C.=C=	QUINAZOLINE
		N=C=C=N=C.-C=C-C=C-C.=	QUINOXALINE
		N=C.-N=C-C=C-C.=C=C=C=	NAPHTHYRIDINE-18
		N=C.-N=C-N=C-C.=N=C=C=	PTERIDINE
		
		N-C.=C=C=C=C=C.-C.=C=C=C=C=C.-	CARBAZOLE
		N=C.-C=C-C=C-C.=C-C.=C=C=C=C=C.-	ACRIDINE
		N=C.-C=C-C=C-C.=N-C.=C=C=C=C=C.-	PHENAZINE
		S-C.=C=C=C=C=C.-N-C.=C=C=C=C=C.-	PHENOTHIAZINE
		O-C.=C=C=C=C=C.-N-C.=C=C=C=C=C.-	PHENOXAZINE
		
		C.-C=C-C-C.=C=C=C=C=	INDENE
		
		C.-C-C-C-C-C.-C-C-C-C-	DECALIN
		C.=C=C=C=C=C.=C=C=C=C=  NAPHTHALENE
		
		N=C=C=C=C.=C=C=C=C=C.=	QUINOLINE
		N=C=C=N=C.=C=C=C=C=C.=	QUINOXALINE
		C.=C-C=C-C.=C-C=C-C=C-	AZULENE
		C.-C.=C=C=C=C=C.-C-C.=C=C=C=C=		FLUORENE
		C.=C-C.=C=C=C=C=C.-C=C.-C=C-C=C-	ANTHRACENE
			
	/;
	
	$hetname =  $heterocycles{$std} || "$std";

	return $hetname;
}

###################################################################
#
#   Routines for labelling individual atom properties/FG membership
#
##################################################################


sub mol_label_alpha_fg {

	#<
	#? Label atoms as being alpha to a given functional group atom
	#. Requires that mol_label_functional_group has been run
	#; Requires: Molecule
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	
	foreach my $atom (atoms($mol)) {
	
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			
			# Label atoms alpha to planar rings
			++$atom->{FG_ALPHA}{PLANAR_RING} if defined $con->{PRINGS}[0];
			++$atom->{FG_ALPHA}{PLANAR_ATOM} if defined $con->{PLANAR_ATOM};
			# Label atoms alpha to specific functional groups
			foreach my $fg (keys %{$con->{FG}}) {
				++$atom->{FG_ALPHA}{$fg};
			}
		}
	}
}

sub mol_label_beta_fg {

	#<
	#? Label atoms as being beta to a given functional group atom
	#. Requires that mol_label_functional_group and mol_label_alpha_fg have been run
	#; Requires: Molecule
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	
	foreach my $atom (atoms($mol)) {
	
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			
			# We are only really interested if we are connected through a carbon
			next if $con->{ELEMENT} ne 'C';
			
			# Label atoms beta to planar rings
			$atom->{FG_BETA}{PLANAR_RING} += $con->{FG_ALPHA}{PLANAR_RING} if defined $con->{FG_ALPHA}{PLANAR_RING};
			# Label atoms beta to specific functional groups
			foreach my $fg (keys %{$con->{FG_ALPHA}}) {
				$atom->{FG_BETA}{$fg} += $con->{FG_ALPHA}{$fg};
			}
		}
	}
}

sub label_waters {

	#<
	#? Label water molecule atoms $atom->{FG}{W};
	#; Requires: molecule, flag to force reading of all atoms in mol_has_h
	#; Sets: $mol->{NUMWATERS}. 
	#; Returns: number of water molecules
	#>

	my $mol = $_[0];
	my $read_all = $_[1];
	
	$read_all = 1 if !defined $read_all;

	my $numwaters = 0;
	my $h = mol_has_h($mol, $read_all);
	my $residues;

	# Shortcut rerunning 
	return $mol->{NUMWATERS} if defined $mol->{NUMWATERS}; 
	
	# Hydrogens are present but no bonds
	# Create atoms in water
	if ($h == -1) {

		# Generate bonds
		connect_atoms_in_waters($mol);
		$h = 1;
	}

	# No hydrogens at all. Assume all isolated oxygens are water
	if ($h == 0) {
	
		foreach my $atom (atoms($mol)) {
		
			next if $atom->{ELEMENT} ne 'O';
			next if $#{$atom->{CONNECT}} != -1;
			
			$atom->{FG}{W} = 1;
			++$numwaters;
			next;
		}
		
		$mol->{NUMWATERS} = $numwaters;
		return $numwaters;
	}

	# All bonds and hydrogens
	foreach my $atom (atoms($mol)) {
	
		next if $atom->{ELEMENT} ne 'O';
		
		my $anum0 = $atom->{CONNECT}[0];
		next if !defined $anum0;
		next if $mol->{ATOMS}[$anum0]{ELEMENT} ne 'H';
		
		my $anum1 = $atom->{CONNECT}[1];
		next if !defined $anum1;
		next if $mol->{ATOMS}[$anum1]{ELEMENT} ne 'H';
		
		$atom->{FG}{W} = 1;
		$mol->{ATOMS}[$anum0]{FG}{W} = 1;
		$mol->{ATOMS}[$anum1]{FG}{W} = 1;
		
		++$numwaters;
	}
	
	$mol->{NUMWATERS} = $numwaters;
	return $numwaters;
}

sub label_solvent {

	#<
	#? Label solvent molecule atoms $atom->{FG}{SOL};
	#; Requires: molecule, maximum number of atoms in solvent
	#; Returns: number of water molecules
	#>

	my $mol = $_[0];
	my $solvent_maxatoms = $_[1] || 10;
	
	my $numsol = 0;
	
	require silico_split;
	
	molecule_check_and_fix_connectivity($mol) if !$mol->{CONNECTIVITY_CHECKED};
	
	my $submols = molecule_label_submols($mol);

	foreach my $submol (@$submols) {
		
		next if $#{$submol} > $solvent_maxatoms;
		++$numsol;
		foreach my $atom (@$submol) {
		
			$atom->{FG}{SOL} = 1;
		}
	}
	
	return $numsol;
}

sub label_aa_backbone {

	#<
	#? Label alpha amino acid backbone atoms
	#. Add the labels AA_N, AA_CA, AA_C to the N, CA anc C atoms of amino acids.  
	#  Terminal amino acids are also labelled AA_CTER or AA_NTER.
	#; Requires: molecule
	#; Returns: numbers main chain-AAs, C-terminal AAs, N-terminal AAs
	#>
	
	my $mol = $_[0];
	my $force = $_[1];

	my $c;
	
	mol_label_functional_group($mol) if $force || !defined $mol->{FG};
	
	# remove any existing values
	foreach my $atom (atoms($mol)) {
		delete $atom->{FG}{AA_N};
		delete $atom->{FG}{AA_C};
		delete $atom->{FG}{AA_CA};
		delete $atom->{FG}{AA_O};
		delete $atom->{FG}{AA_CTER};
		delete $atom->{FG}{AA_NTER};
		delete $atom->{FG}{AA_CB};
	}
	
	# Find and label atom sequence (C') N CA C (N' (CA')) 
	# (Bracketed atoms are optional)
	
	my ($c1, $c2, $c3);
	$c1 = $c2 = $c3 = 0;
	foreach my $atom (atoms($mol)) {
	
		# Skip if we are not a carbon 
		next if $atom->{ELEMENT} ne 'C';

		my ($an, $ac, $atom_n, $atom_c, $atom_nter, $atom_cter, $atom_cb, $atom_np, $atom_cp, $atom_cap);

		foreach my $con (connected($atom, $mol)) {
		
			my $fg = $con->{FG};

			# Two of the following reqired if atom is a CA
			$atom_n =    $con if ($fg->{S_AMIDE_N} || $fg->{T_AMIDE_N});
			$atom_c =    $con if ($fg->{S_AMIDE_C} || $fg->{T_AMIDE_C});
			$atom_nter = $con if ($fg->{P_AMMONIUM_N} || $fg->{P_AMINE_N}) ;
			$atom_cter = $con if ($fg->{CARBOXYLATE_C} || $fg->{CARBOXYLIC_C});
			
			#atom_printout($con) if $con->{NAME} =~ 'OXT';
			#atom_printout($con) if $con->{ELEMENT} eq 'N';
		}
		
		# Normal residue
		my $flag = 0;
		if ($atom_n && $atom_c) {
			
			$atom_n->{FG}{AA_N} = 1;
			$atom->{FG}{AA_CA} = 1;		
			$atom_c->{FG}{AA_C} = 1;
			$flag = 1;
			++$c1;
			#print "Normal residue\n";
		}

		# N terminal residue
		elsif ($atom_nter && $atom_c) {
			
			$atom_nter->{FG}{AA_N} = 1;
			$atom->{FG}{AA_CA} = 1;
			$atom_c->{FG}{AA_C} = 1;

			$atom_nter->{FG}{AA_NTER} = 1;
			$atom->{FG}{AA_NTER} = 1;
			$atom_c->{FG}{AA_NTER} = 1;
			$flag = 2;
			++$c2;
			#print "N-terrminal residue\n";
		}

		# C terminal residue
		elsif ($atom_n && $atom_cter) {
			
			$atom_n->{FG}{AA_N} = 1;
			$atom->{FG}{AA_CA} = 1;
			$atom_cter->{FG}{AA_C} = 1;

			$atom_n->{FG}{AA_CTER} = 1;
			$atom->{FG}{AA_CTER} = 1;
			$atom_cter->{FG}{AA_CTER} = 1;
			$flag = 3;
			++$c3;
			#print "C-terminal residue\n";
		}

		next if !$flag;
		
		# Find N' and CA' 
		foreach my $con (connected(($atom_c || $atom_cter), $mol)) {
				
			next if $con->{ELEMENT} ne 'N';
			my $fg = $con->{FG};
			if ($fg->{S_AMIDE_N} || $fg->{T_AMIDE_N} || $fg->{AMMONIUM_N} || $fg->{P_AMINE_N}) {
					
				$atom_np = $con;
					
				foreach my $con2 (connected($con, $mol)) {
						
					next if $con2->{ELEMENT} ne 'C';
					my $fg = $con2->{FG};
					next if !($fg->{QUATERNARY_C} || $fg->{METHYNE_C} || $fg->{METHYLENE_C});
					$atom_cap = $con2;
					#print "@{$con2->{FG}}\n";
					#atom_printout($con2);
					last;
					
				}
				last;
			}
		}
		
		# Find C'	
		foreach my $con2 (connected(($atom_n || $atom_nter), $mol)) {
				
			next if $con2->{ELEMENT} ne 'C';
			my $fg = $con2->{FG};
			$atom_cp = $con2 if ($fg->{S_AMIDE_C} || $fg->{T_AMIDE_C} || $fg->{CARBOXYLATE_C});
			last;
		}	
		
		# Label attached atoms (carbonyl oxygen, CB, HA, HB) if we have found CA
		my $cb_count = 0;
		if (($atom_n ||$atom_nter) && ($atom_c || $atom_cter)) {
			$c = $atom_c || $atom_cter;
			foreach my $con (connected($c, $mol)) {
				next if $con->{ELEMENT} ne 'O';
				$con->{FG}{AA_O} = 1;
			}
			foreach my $con (connected($atom, $mol)) {
				
				if ($con->{ELEMENT} eq 'C' && !$con->{FG}{AA_C}) {
					$con->{FG}{AA_CB} = 1;
					$atom_cb = $con;
					++$cb_count;
					foreach my $con2 (connected($con, $mol)) {
						next if $con2->{ELEMENT} ne 'H';
						$con2->{FG}{AA_HB} = 1;
					}
				}

				# Label CA hydrogens
				if ($con->{ELEMENT} eq 'H') {
					$con->{FG}{AA_HA} = 1;
				}
			}
		}
		
		# Label terminal hydrogens
		if ($atom_n && ($atom_c || $atom_cter)) {
			foreach my $con (connected($atom_n, $mol)) {
				next if $con->{ELEMENT} ne 'H';
				$con->{FG}{AA_HN} = 1;
			}
		}
		
		# Make dihedral list
		#print join "\t", $atom_cp->{NAME}, ($atom_n || $atom_nter)->{NAME}, $atom->{NAME}, ($atom_c || $atom_cter)->{NAME}, $atom_np->{NAME}, $atom_cap->{NAME}, "\n";
		#print join "\t", $atom_cp->{NUM}, ($atom_n || $atom_nter)->{NUM}, $atom->{NUM}, ($atom_c || $atom_cter)->{NUM}, $atom_np->{NUM}, $atom_cap->{NUM}, "\n";
		
		# Chirality 
		if (($atom_n || $atom_nter) && ($atom_c ||$atom_cter) && ($atom_cb && $cb_count == 1)) {
			my $c = aa_chirality(($atom_n || $atom_nter), $atom, ($atom_c ||$atom_cter), $atom_cb);
			$atom->{CHIRALITY} = $c;
			#print "c $c\n";
		}
	}

	silico_msg('g', "Found $c1 main chain AAs.  $c2 N-terminal and $c3 C-terminal AAs.\n");
	return ($c1, $c2, $c3);
}


sub is_aa_backbone {

	#<
	#? Return true if atom is in backbone residue (C O CA N)
	#; Requires: atom
	#; Returns: true or false
	#>
	
	foreach my $t (qw "AA_N AA_C AA_O AA_CA") {
		return 1 if $_[0]->{FG}{$t};
	}
	return 0;
}


sub label_amides {
	
	#<
	#? Label amides AMIDE_C, AMIDE_O, AMIDE_N
	#. Bondorders must be correct
	#. Note: This algorithm handles carbamates OK, but ureas somewhat imprecisely.
	#; Requires: molecule, list of atoms to check (uses entire molecule if this is not defined)
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $atom_list = $_[1];
	
	my $count = 0;
	
	
	if (!defined $atom_list) {
		@$atom_list = @{$mol->{ATOMS}};
	}
	
	foreach my $atom (@$atom_list) {
	
		next if $atom->{ELEMENT} ne 'C';
	
		my $i = -1;
		my $flag1 = 0;
		my $flag2 = 0;
		
		# Find if we have =O and -N 
		foreach my $connum (@{$atom->{CONNECT}}) {
			
			my $con = $mol->{ATOMS}[$connum];
			++$i;
			
			if ($con->{ELEMENT} eq 'O' && $atom->{BORDERS}[$i] == 2) {
				++$flag1;
			}
			
			if ($con->{ELEMENT} eq 'N' && $atom->{BORDERS}[$i] == 1) {
				++$flag2;
			}
		}
		
		# Skip if we didn't have =O and -N
		next if !($flag1 == 1 && $flag2 == 1);
		
		++$count;
		# Label atoms
		$atom->{FG}{AMIDE_C} = 1;
		$i = -1;
		
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			++$i;
			
			if ($con->{ELEMENT} eq 'O' && $atom->{BORDERS}[$i] == 2) {
				$con->{FG}{AMIDE_O} = 1;
			}
			
			if ($con->{ELEMENT} eq 'N' && $atom->{BORDERS}[$i] == 1) {
				$con->{FG}{AMIDE_N} = 1;
			}
		}
	}  
	
	silico_msg('g', "Found $count amides\n"); 
}


sub mol_label_united_atom {

	#<
	#? Calculate missing hydrogens for carbon atoms
	#  and add UNITED_ATOM and MISSING_HYDROGEN labels
	#. Used for united atom forcefields 
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	
	my $con;
	my $flag;
	my $valence;
	
	foreach my $atom (atoms($mol)) {
		delete $atom->{MISSING_HYDROGENS};
		delete $atom->{UNITED_ATOM};
	}
	
	# First find united atom carbons
	foreach my $atom (atoms($mol)) {
	
		next if $atom->{ELEMENT} ne 'C';
	
		$valence = valence($atom);

		if ($valence < 4) {
			$atom->{MISSING_HYDROGENS} = 4-$valence;
			$atom->{UNITED_ATOM} = 1;
		}
	}
	
	# Find all atoms that are attached to UA carbons and mark them
	# with a UNITED_ATOM label
	$flag = 1;
	while ($flag == 1) {
	
		$flag = 0;
		foreach my $atom (atoms($mol)) {
			
			#Skip atoms already labelled
			next if !$atom->{UNITED_ATOM};
			
			#Skip terminal carbon atoms
			next if !$atom->{ELEMENT} eq 'C';
						
			foreach my $con (connected($atom, $mol)) {
			
				next if $con->{UNITED_ATOM};
				$con->{UNITED_ATOM} = 1;
				$flag = 1;
			}
		}
	}	
}




sub mol_get_functional_group {

	#<
	#? Get array of functional groups given a functional group name
	#; Requires: molecule, group name (string)
	#; Returns: pointer to array of functional groups (each fg is a hash)
	#>

	my $mol = $_[0];
	my $fgname = $_[1];
	
	my $fglist;
	
	foreach my $fg (@{$mol->{FG}}) {
	
		push @$fglist, $fg if $fg->{NAME} eq $fgname;
	}
	
	return $fglist;
}
	
		


##################################################################
#
#       Breadth-first Tree Searching
#
##################################################################


sub add_bondorders_to_atom_tree {
	
	#<
	#? Add bondorders to the branches in an atom tree.
	#; Requires: molecule, atom tree
	#; Returns: modified atom tree
	#>
	
	my $mol = $_[0];
	my $tree = $_[1];
	
	my $branch;
	my $modified_tree;
	
	foreach my $branch (@$tree) {
		
		my $modified_branch;
		
		for (my $i = 0; $i < $#{$branch}; ++$i) {
			
			my $atomoffset = $branch->[$i];
			my $conoffset = $branch->[$i+1];
			my $atom = $mol->{ATOMS}[$atomoffset];
			
			#
			# Following code looks suspect - should not use $j outside loop...
			#
			my $j;
			for ($j = 0; $j <= $#{$atom->{CONNECT}}; ++$j) {
				last if $conoffset == $atom->{CONNECT}[$j];				
			}
			
			my $bo = $atom->{BORDERS}[$j];
			
			push @$modified_branch, $atomoffset;
			push @$modified_branch, $bo;
		}
		
		push @$modified_branch, $branch->[$#{$branch}];
		
		push @$modified_tree, $modified_branch;	
	}
	
	return $modified_tree;
}

sub find_atom_paths {

	#<
	#? Find the path connecting each atom to a given root atom
	#; Requires: molecule, starting atom, optional maximum path length, optional exclusion of H
	#; Returns: Nothing (information added to $atom->{PATH})
	#>
	
	my $mol = $_[0];
	my $start = $_[1];
	my $depth = $_[2];
	my $noh = $_[3];
	
	my @receive;
	my @send;
	my $sendnum;
	
	my $i = 0;
	
	foreach my $atom (atoms($mol)) {
	
		if ($atom == $start) {
			@send = ($i);
			last;
		}
		
		++$i;
	}
	
	my $j = 0;
	while (defined $send[0]) {
	
		@receive = ();
		
		foreach my $sendnum (@send) {
					
			my $sendatom = $mol->{ATOMS}[$sendnum];
			
			CON: foreach my $connum (@{$sendatom->{CONNECT}}) {
			
				my $con;
				my $path_atom;
				
				$con = $mol->{ATOMS}[$connum];
				
				next if ($con->{ELEMENT} eq 'H' && $noh);
				
				# Stop going round and round a ring
				foreach my $path_atom (@{$sendatom->{PATH}}) {
					next CON if $con == $mol->{ATOMS}[$path_atom];
				}
				
				@{$con->{PATH}} = @{$sendatom->{PATH}} if (defined $sendatom->{PATH}[0]);
				
				push @{$con->{PATH}}, $sendnum;
				
				push @receive, $connum;
			}
		}
		
		++$j;
		last if (defined $depth && $j == $depth);
		
		@send = @receive;
	}
	
	# Debugging
		
	if (get_sflag('debug')) {
	
		my $i;
		foreach my $atom (atoms($mol)) {
	
			my $patomnum;
		
			++$i;
			print "Atom $i of ".($#{$mol->{ATOMS}}+1).": ";
		
			foreach my $patomnum (@{$atom->{PATH}}) {
			
				my $patom = $mol->{ATOMS}[$patomnum];
			
				print $patom->{SUBNAME};
				print $patom->{SUBID};
				print "-";
			}
		
			print $atom->{SUBNAME};
			print $atom->{SUBID};
			print "\n";
		}
	}
	
	# End debugging
}

sub construct_atom_relationships {

	#<
	#? Construct the parent-child relationships between pairs of atoms,
	#  relative to an arbitrary common ancestor, based on the contents
	#  of $atom->{PATH}.
	#; Requires: molecule, optional flag to terminate at substructure boundaries
	#;
	#>
	
	my $mol = $_[0];
	my $term = $_[1];
	
	my $relationships;
	
	foreach my $atom (atoms($mol)) {
		
		my $predecessor;
		my $relate;
		
		# We needn't worry about any atom with no path
		next if $#{$atom->{PATH}} < 0;
		
		$predecessor = $mol->{ATOMS}[$atom->{PATH}[$#{$atom->{PATH}}]];
		
		if (		(!$term)
		
			||	(	(  !defined $atom->{SUBID}
					|| !defined $predecessor->{SUBID}
					|| $atom->{SUBID} eq $predecessor->{SUBID})
				&&	(  !defined $atom->{SUBNAME}
					|| !defined $predecessor->{SUBNAME}
					|| $atom->{SUBNAME} eq $predecessor->{SUBNAME})
				&&	(  !defined $atom->{CHAIN}
					|| !defined $predecessor->{CHAIN}
					|| $atom->{CHAIN} eq $predecessor->{CHAIN})
				&&	(  !defined $atom->{SEGID}
					|| !defined $predecessor->{SEGID}
					|| $atom->{SEGID} eq $predecessor->{SEGID})
				)
		    ) {
		
			my $i = -1;
			foreach my $atom2 (atoms($mol)) {
				++$i;
				$relate->{PARENT} = $i if $atom2 == $predecessor;
				$relate->{CHILD} = $i if $atom2 == $atom;
			}
			
			push @$relationships, $relate;
		}
	}
	
	# Debugging
	if (get_sflag('debug')) {
		
		my $i = 0;
		foreach my $r (@$relationships) {
			
			++$i;
			printf "Relationship $i:\tParent: $r->{PARENT}\tChild: $r->{CHILD}\n";
		}
	}
	
	# End debugging
	
	return $relationships;
}

sub construct_residue_relationships {

	#<
	#?
	#;
	#;
	#>
	
	my $mol = $_[0];
	
	my $relationships;
	
	foreach my $atom (atoms($mol)) {
		
		my $predecessor;
		my $relate;
		
		next if $atom->{ELEMENT} eq 'H';
		# We needn't worry about any atom with no path
		next if $#{$atom->{PATH}} < 0;
		
		$predecessor = $mol->{ATOMS}[$atom->{PATH}[$#{$atom->{PATH}}]];
		
		if ($atom->{SUBID} != $predecessor->{SUBID}) {
			
			$relate->{PARENT} = $predecessor->{SUBNAME}.$predecessor->{SUBID};
			$relate->{CHILD} = $atom->{SUBNAME}.$atom->{SUBID};
			
			push @$relationships, $relate;
		}
	}
	
	# Debugging
	
	if (get_sflag('debug')) {
		
		my $i = 0;
		foreach my $r (@$relationships) {
			
			++$i;
			printf "Relationship $i:\tParent: $r->{PARENT}\tChild: $r->{CHILD}\n";
		}
	}
	
	# End debugging
	
	return $relationships;
}

sub construct_family_tree {

	#<
	#? Construct a family tree from parent-child relationships with an
	#  arbitrary ultimate ancestor. Can take any type of correctly formatted
	#  input. So far, atoms and residues use this.
	#; Requires: array of relationships (each one a hash containing PARENT and CHILD entries)
	#; Returns: branches
	#>
	
	my $relationships = $_[0];
	
	my $branches;
	my $common_ancestor;
	my $extant_branch;
	my $extant_branches;
	my $initbranch;
	my $new_branches;
	my @receive;
	my $relate1;
	my @send;
	my $sender;
	
	foreach my $relate1 (@$relationships) {
		
		my $flag = 0;
		my $relate2;
		
		foreach my $relate2 (@$relationships) {
		
			$flag = 1 if $relate1->{PARENT} eq $relate2->{CHILD};
		
		}
		
		$common_ancestor = $relate1->{PARENT} if $flag == 0;
	}
	
	@$initbranch = ($common_ancestor);
	@$extant_branches = ($initbranch);
	@send = ($common_ancestor);
	
	while (defined $send[0]) {
		
		@receive = ();
		
		foreach my $sender (@send) {
			
			foreach my $relate1 (@$relationships) {
				
				my $holding_array;
				
				if ($relate1->{PARENT} eq $sender) {
					
					foreach my $extant_branch (@$extant_branches) {
						
						next if $extant_branch->[$#{$extant_branch}] ne $relate1->{PARENT};
						
						foreach (@$extant_branch) {push @$holding_array, $_}
					}
					
					push @$holding_array, $relate1->{CHILD};
					
					$new_branches->[$#{$new_branches}+1] = $holding_array;
					
					push @receive, $relate1->{CHILD};
				}
			}
		}
		
		@send = @receive;
		$extant_branches = $new_branches;
	}
	
	$branches = prune_family_tree($new_branches);
	return $branches;
}

sub prune_family_tree {

	#<
	#? Taking a list of branches in a family tree, throw away all those
	#  which are incomplete.
	#; Requires: list of branches
	#; Returns: list of branches
	#>
	
	my $branches_in = $_[0];
	
	my $branch_in1;
	my $branch_in2;
	my $branches_out;
	
	BRANCH: for (my $i = 0; $i <= $#{$branches_in}; ++$i) {
	
		$branch_in1 = $branches_in->[$i];
	
		foreach my $branch_in2 (@$branches_in) {
		
			for (my $j = 0; $j < $#{$branch_in2}; ++$j) {
			
				my $thing = $branch_in2->[$j];
				
				next BRANCH if $branch_in1->[$#{$branch_in1}] eq $thing;
			}
		}
		
		push @$branches_out, $branch_in1;
	}
	
	# Debugging
	
	if (get_sflag('debug')) {
			
		foreach my $branch (@$branches_out) {
			
			print join "-", @$branch;
			print "\n";
		}
	}
	
	# End debugging
	
	return $branches_out;
}




##################################################################
#
#      Filtering molecules by properties
#
##################################################################

sub carbanion_filter_neutraliser {

        #<
        #? Some incorrectly formatted SDF files can contain incorrect carbanions on aromatic rings
        #; Requires: molecule
        #; Returns: 1 if carbanion is detected (except for isocyanides -N+%C-) and all metal containing compounds. The carbanion is also corrected.
        #>

        my $mol = $_[0];

	my $count = 0;
	my @metals = qw(Ag Al Bi Ce Co Cr Cu Er Eu Hf In Ir Mg Mn Mo Ni Pb Pd Pr Pt Rh Ru Ru Sn Ta V W Yb Zr);

        ATOM: foreach my $atom (atoms($mol)) {
		next if !$atom->{FORMAL_CHARGE} || $atom->{FORMAL_CHARGE} != -1;
		next if $atom->{ELEMENT} ne 'C';
		
		my $charge_is_bad = 1;
		
		# Isonitrile
		my $i = 0;
		foreach my $bo (@{$atom->{BORDERS}}) {
			if ($bo == 3 && $mol->{ATOMS}[$atom->{CONNECT}[$i]]{ELEMENT} eq 'N') {
				#print "*** Isonitrile ***\n";
			 	$charge_is_bad = 0;
			}
			++$i;
		}
		
		# Metal compound
		molecule_formula($mol, 1) if !$mol->{FORMULA};
		
		foreach (@metals) {
			if ($mol->{FORMULA}{$_}) {
				print "*** $_ compound ***\n";
				$charge_is_bad = 0;
				last ATOM; # 
			}
		}
		
		if ($charge_is_bad) {
			$atom->{CHARGE} = 0;
			++$count;
		}
        }
	
	#if ($count) {
	#
		#++$Silico::counter;
		#my $ens;
		#push @$ens, $mol;
		#write_mol_any($ens, "out_$Silico::counter");
	#}

        return $count;
}

sub chonspx_filter {

        #<
        #? Find if molecule contains other than CHONSP+Halogen elements
	#. H* is included to take implicit hydrogens into account
        #; Requires: molecule
        #; Returns: Element name (true) if non CHONSPX elements are present or 0 otherwise
        #>

        my $mol = $_[0];

	my $formula = molecule_formula($mol, 1);
	my %good = qw(C 1 H 1 H* 1 O 1 N 1 S 1 P 1 F 1 Cl 1 Br 1 I 1 );

	foreach my $key (keys %$formula) {
		next if $good{$key};
		#print "Element $key is bad\n";
		return $key;
	}
	return 0;
}

sub reos_filter {

        #<
        #? Implements a simple filter for the Rapid Elimination of Swill
        #; Requires: molecule
        #; Returns: true or false
        #>

        my $mol = $_[0];

        my $max_allowed;
        my $mol_fg_hash;

        %$max_allowed = qw(ALKYNE 2 NITRILE 2 ACID_HALIDE 0 ALDEHYDE 1 PRIMARY_ALKYL_HALOGEN 0 THIOCYANATE 0 ISOTHIOCYANATE 0 CYANATE 0 ISOCYANATE 0);

        # Label functional groups (perhaps again)
        mol_label_functional_group($mol);

        # Count functional groups
        foreach (@{$mol->{FG}}) {
                ++$mol_fg_hash->{$_->{NAME}};
        }

        foreach my $fg (keys %$mol_fg_hash) {;

                next if !defined $max_allowed->{$fg};
                next if $mol_fg_hash->{$fg} <= $max_allowed->{$fg};

                print "Rejected by REOS filter: $fg\n";
                return 0;
        }

        # Skip compounds with more than 4 halogens
        return 0 if count_halogens($mol) > 4;

        return 1;
}


sub count_halogens {

        my $mol = $_[0];

        my $count = 0;

        foreach my $atom (atoms ($mol)) {

                my $el = $atom->{ELEMENT};
                ++$count if $el eq 'F' || $el eq 'Cl' || $el eq 'Br' || $el eq 'I';
        }

        return $count;
}


##################################################################
#
#	Atom type errors
#
##################################################################

sub warn_atom_type_error {

	#<
	#? Add atom type warning to mol->{WARN};
	#. Could be made simpler
	#; Requires: Molecule, atom, numsingle, numdouble, $numtriple, numatomatic, numbonds
	#; Returns: Nothing
	#>

        my ($mol, $atom, $numsingle, $numdouble, $numtriple,$numaromatic, $numbonds) = @_;
        
        if ($Silico::warnings) {
                print_atom_type_error( $atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
        }
        
        my $string =  "Unknown $atom->{ELEMENT} atom type. Use -warnings flag for more information";
        
        ++$mol->{WARN}{$string};
}

sub print_atom_type_error {

	my ($atom, $numsingle, $numdouble, $numtriple,$numaromatic, $numbonds) = @_;

	my $dash = 0;

	#print "$Silico::debug\n";
	#print "$Silico::warnings\n";

	if (!$Silico::debug) {
		my $string =  "Unknown $atom->{ELEMENT} atom type. ";
		$string .= "Num: $atom->{NUM}. " if $atom->{NUM};
		$string .= "Name: $atom->{NAME}. " if $atom->{NAME};
		$string .= "Subname: $atom->{SUBNAME}. " if $atom->{SUBNAME};
		$string .= "Subid: $atom->{SUBID}. " if $atom->{SUBID};	
		$string .= " Use -debug flag for extra information\n";
		
		silico_msg ('c', $string);
	
	} else {

		print "\n";
		print heading("Bonding\n");
		print "Single\t\t$numsingle\n";
		print "Double\t\t$numdouble\n";
		print "Triple\t\t$numtriple\n";
		print "Aromatic\t$numaromatic\n";
		print "\t\t-\n";
		print "TOTAL\t\t$numbonds\n";
		print "\n";
		print heading("Functional groups\n");
		foreach my $fg (keys %{$atom->{FG}}) {
			printf "%-15s %-s\n", $fg, ($atom->{FG}{$fg} || "UNDEFINED");
		}
		unless ($numbonds == 0) {
			print "\n";
			print heading("Connected atoms\n");
		
			foreach my $el (keys %{$atom->{CONNECTED}}) {
				next if $el eq 'C';
				next if $el eq 'HET';
				next if $el eq 'HEAVY';
				next if $el eq 'H';
				print "$el\t\t$atom->{CONNECTED}{$el}\n";
				$dash = 1;
			}
			if ($dash) {
				print "\t\t-\n";
				$dash = 0;
			}
			if (defined $atom->{CONNECTED}{HET} || defined $atom->{CONNECTED}{C}) {
				print "HET\t\t$atom->{CONNECTED}{HET}\n" if defined $atom->{CONNECTED}{HET};
				print "C\t\t$atom->{CONNECTED}{C}\n" if defined $atom->{CONNECTED}{C};
				$dash = 1;
			}
			if ($dash) {
				print "\t\t-\n";
				$dash = 0;
			}
			if (defined $atom->{CONNECTED}{HEAVY} || defined $atom->{CONNECTED}{H}) {
				print "HEAVY\t\t$atom->{CONNECTED}{HEAVY}\n" if defined $atom->{CONNECTED}{HEAVY};
				print "H\t\t$atom->{CONNECTED}{H}\n" if defined $atom->{CONNECTED}{H};
				$dash = 1;
			}
			if ($dash) {
				print "\t\t-\n";
				$dash = 0;
			}
			print "TOTAL\t\t$numbonds\n";
			print heading("Connected functional groups\n");
		
			foreach my $alpha (keys %{$atom->{FG_ALPHA}}) {
				printf "%-15s %-s\n", $alpha, ($atom->{FG_ALPHA}{$alpha} || "UNDEFINED");
			}
		}
		print "\n";
	}
}


return 1;
