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
#! silico_rings.pm
#? Silico routines to identify rings in molecules.
#. $Revision: 1.15-42-g45ba5cf $
#>

use strict;
package Silico;

##############################
# Routines to identify rings #
##############################

sub molecule_find_rings {

	#<
	#? Find rings and planar rings in a molecule.
	#
	#. Rings are arrays of atom numbers.  This routine returns all rings found up to
	#  max_ring_size.  It does not find the minimum set of smallest rings.  So naphthalene
	#  will return three rings, two of size six and one of size 10.
	#. Sets the following arrays for each molecule:
	#; @{$molecule->{RINGS}} All rings in the molecule
	#; @{$molecule->{PRINGS}} All planar rings
	#; @{$molecule->{ARINGS}} All aromatic rings (where a planar ring is not aromatic)
	#. The maximum ring size is stored in $mol->{RINGSIZE_MAX}.
	#. The following arrays are set for each atom.  Where the atom is not a member of any
	#  ring, the array is empty.
	#; @{$atom->{RINGS}} All rings that that atom belongs to
	#; @{$atom->{PRINGS}} All planar rings that that atom belongs to
	#. In addition, the following flags are set
	#; $atom->{AROMATIC_RING} flag for all 5 and 6 membered aromatic rings.
	#.
	#; $atom->{NOTRING} flag for all non-ring atoms
	#; Requires: ensemble, maximum ring size (optional, default 8), use
	#  non-recursive subroutine to find rings (optional), force (optional)
	#; Returns: rings, planar rings, aromatic rings
	#>

	my $mol = $_[0];
	my $max_ring_size = $_[1] || 8;	#Maximum search distance from each atom
	my $notrec = $_[2];
	my $force = $_[3];

	my %ring_count;
	my $all_rings;
	my $aromatic_rings;
	my $hash;
	my $planar_rings;
	my $ring_count;
	my $starttime = (times)[0];
		
	# We only want to run molecule_find_rings if at least one of the following criteria is met:
	# - the programmer has requested a forced regeneration of rings
	# - the molecule is known to have bad rings (e.g., through deletion of bonds)
	# - the ring size to search for exceeds whatever ring size has been searched for
	#   previously

	my $do = 1 if ($force || $mol->{BAD_RINGS} || !defined $mol->{RINGSIZE_MAX} || $max_ring_size > $mol->{RINGSIZE_MAX});
	
	# If none of these criteria are met, return whatever was previously known.
	if (!$do) {
		return $mol->{RINGS}, $mol->{PRINGS}, $mol->{ARINGS};
	}
	
	silico_msg('w', "Calling find_rings with a large number of atoms: $mol->{NUMATOMS}. This may be slow\n") if ($mol->{NUMATOMS} || 0) > 100000;
	#print "find_rings\n\n";
	#carp();

	if ($Silico::TIMING || $Silico::debug_rings) {
		silico_msg('v', "Now entering ".subname()." with maximum ring size of $max_ring_size\n");
	}

	my $rings;
	@$rings = ();

	# Get rid of old rings if any
	delete_rings($mol) if $mol->{RINGS}[0];
	delete $mol->{BAD_RINGS};

	# Find the rings for each atom
	my $i = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {

		++$i;

		if ($Silico::TIMING && $i %10000 == 0) {
			print "Atom: $i Time: ".calc_human_readable_time((times)[0]-$starttime)."\n";
		}

		# Skip water if we know about it
		next if $atom->{FG}{W};

		my $nonhcount = 0;
		# Count number of nonhydrogen connections
		foreach	my $connum (@{$atom->{CONNECT}}) {
			++$nonhcount if $mol->{ATOMS}[$connum]{ELEMENT} ne 'H';
		}

		# We can't be in a ring if we have 0 or 1 nonh connections
		if ($nonhcount <= 1) {
			$atom->{NOTRING} = 1;
			next;
		}

		# We have already found all rings if we have one ring and two nonh connections
		if ($nonhcount == 2  && defined $atom->{RINGS}[0]) {
			next;
		}

		if (($nonhcount > 4 && ($atom->{ELEMENT} eq 'C' || $atom->{ELEMENT} eq 'N')) || ($nonhcount > 3 && $atom->{ELEMENT} eq 'O')) {
			
			my $nm;
			my $i = 0;
			my $nb = $#{$atom->{CONNECT}}+1;
			my $s = atominfo($atom);
			#silico_msg('w', "Too many bonds to atom: $s num bonds: $nb\n");
			++$mol->{WARN}{FIND_RINGS_TOO_MANY_BONDS};
		}

		if ($notrec) {
	
			# Warning. This step will be time consuming
			foreach my $atom2 (atoms($mol)) {
				delete $atom2->{RING_PATH};
			}
			
			$rings = find_atom_rings_not_recursive($mol, $i, $max_ring_size);
		} else {
			$rings = find_atom_rings($mol, $i, $max_ring_size);
		}
		
		# Mark atom if it is not in a ring
		$atom->{NOTRING} = 1 if !defined $rings;

		# Save a list of all rings in $all_rings and remove the duplicates
		foreach my $ring (@$rings) {

			$ring = standardise_ring($ring);
			
			$hash = join ('',@$ring);
			next if (++$ring_count{$hash} > 1);
			push @$all_rings, $ring;

			# Add ring $atom->{RINGS} for all atoms in ring
			foreach my $atomnum (@$ring) {
				push @{$mol->{ATOMS}[$atomnum]{RINGS}}, $ring;
			}
		}
	}

	silico_msg('t', "Found ".($#{$all_rings}+1)." rings in time ".calc_human_readable_time((times)[0]-$starttime)."\n");

	# Identify planar rings (all flat rings) and aromatic rings
	# and add list to each atom in ring
	$i = 0;
	foreach my $ring (@$all_rings) {

		++$i;
		if (($Silico::TIMING || $Silico::debug_rings) && $i %10 == 0) {
			silico_msg('v', "Finding planar rings. Ring: $i Time:".calc_human_readable_time((times)[0]-$starttime)."\n");
		}
	
		my ($planar, $aromatic) = is_ring_planar($mol, $ring);
	
		if ($planar) {
			foreach my $atomnum (@$ring) {
				$mol->{ATOMS}[$atomnum]{PLANAR_RING} = 1;
				push @{$mol->{ATOMS}[$atomnum]{PRINGS}}, $ring;
			}
		
			push @$planar_rings, $ring;
		}
		if ($aromatic) {
			foreach my $atomnum (@$ring) {
				$mol->{ATOMS}[$atomnum]{AROMATIC_RING} = 1;
				push @{$mol->{ATOMS}[$atomnum]{ARINGS}}, $ring;
			}
		
			push @$aromatic_rings, $ring;
		}
	}

	$mol->{RINGS} = $all_rings;
	$mol->{PRINGS} = $planar_rings;
	$mol->{ARINGS} = $aromatic_rings;
	$mol->{RINGSIZE_MAX} = $max_ring_size;
	
	silico_msg('t', "Finished in find_rings in time ".calc_human_readable_time((times)[0]-$starttime)."\n");
	
	if ($mol->{WARN}{FIND_RINGS_TOO_MANY_BONDS}) {
	
		silico_msg('w', "Found $mol->{WARN}{FIND_RINGS_TOO_MANY_BONDS} atoms with more than four nonhydrogen atoms connected\n"); 
	}

	return $all_rings, $planar_rings, $aromatic_rings;
}


sub find_atom_rings {

	#<
	#? Find all rings that belong to a particular atom
	#. Calls traverse.
	#; Requires: molecule, starting atom, maximum search depth.
	#; Returns: a pointer to array of arrays containing ring atoms
	#  or undef if no rings are found.
	#>

	my $mol = $_[0];
	
	use vars qw(	$far_atoms
			$traverse_start
			$traverse_depth
			$traverse_ringcount
			$traverse_rings
		);
	
	@$traverse_rings = ();
	
	$far_atoms = $mol->{ATOMS};
	$traverse_start = $_[1];	# Starting atom for ring search
	$traverse_depth = $_[2]+1;	# Maximum search depth
	$traverse_ringcount = 0;	# Number of rings found
		
	traverse($traverse_start);

	return $traverse_rings;
}

sub find_atom_rings_not_recursive {

	#<
	#? Subroutine to find rings, based on a non-recursive algorithm
	#. Algorithm by Balducci and Pearlman, JCICS 1994, 34, 822-831
	#. Slightly modified to suit Silico conventions
	#; Requires: Molecule, starting atom, maximum ring size
	#; Returns: List of rings containing this atom
	#>
	
	my $mol = $_[0];
	my $i = $_[1];
	my $max_ring_size = $_[2];
	
	my @receive;
	my $plength;
	my @ring;
	my $rings;
	my @send;
	@send = ($i);
	
	while (defined $send[0]) {
	
		# Clear the list of receivers
		@receive = ();
	
		# For each sender number...
		foreach my $sendnum (@send) {
			
			# Get the sender
			my $sendatom = $mol->{ATOMS}[$sendnum];
			
			# Set the length of the path to this sender so far
			$plength = $#{$sendatom->{RING_PATH}} if (defined $sendatom->{RING_PATH}[0]);
			
			# This test here should find if the ring is too big
			# The path is too long if the length of the path is
			# at least half the ring size
			# Why: This means that if the ring were this size,
			# we would be more than halfway round. We should have
			# already met the path coming up the other side by now
			# This is not a perfect test. If $max_ring_size is odd,
			# this will not exclude rings with a size of $max_ring_size + 1.
			# These rings are tested for later.
			next if (defined $plength && $plength+1 >= $max_ring_size/2);
			
			# For each connected atom...
			foreach my $connum (@{$sendatom->{CONNECT}}) {
				
				# ...skip if this atom is the one the signal just came from
				next if (defined $sendatom->{RING_PATH}[0] && $connum == $sendatom->{RING_PATH}[$plength]);
	
				# get the atom table
				my $con = $mol->{ATOMS}[$connum];
				
				# Skip this atom if it is definitely not in a ring
				# (this is not formally part of the original algorithm)
				next if $con->{NOTRING};
				
				# Skip this atom if it has less than two connected heavy atoms
				# (this is not formally part of the original algorithm)
				next if ($#{$con->{CONNECT}} < 1);
				
				my $connected_heavy_atoms = 0;
				foreach my $connum2 (@{$con->{CONNECT}}) {
				
					my $con2 = $mol->{ATOMS}[$connum2];
					
					++$connected_heavy_atoms if ($con2->{ELEMENT} ne 'H');
				}
				
				next if ($connected_heavy_atoms < 2);
				
				# The ring path of the connected atom may already be defined
				# and it not be the one we came from. If this is the case,
				# we have a ring
				
				# However, we also need to exclude the case where it forks
				# further down the path - these should be picked up on another sweep
				
				if (defined $con->{RING_PATH}[0]) {
				
					if (defined $con->{RING_PATH}[1] && defined $sendatom->{RING_PATH}[1]
						&& $con->{RING_PATH}[1] == $sendatom->{RING_PATH}[1]) {
					
						next;
						
					} elsif (defined $con->{RING_PATH}[1] or defined $sendatom->{RING_PATH}[1]) {
								
						@ring = (@{$sendatom->{RING_PATH}}, $sendnum, $connum, reverse @{$con->{RING_PATH}});
						@ring = @ring[0..$#ring-1];
						
						# This test is here to discount those rings of size
						# $max_ring_size + 1.
						# This is only a problem if $max_ring_size is odd.
						@{$rings->[$#{$rings}+1]} = @ring if ($#ring < $max_ring_size);
						next;
					}
				}
				
				# Set the ring path of this atom to the one of its predecessor...
				@{$con->{RING_PATH}} = @{$sendatom->{RING_PATH}} if (defined $sendatom->{RING_PATH}[0]);
				
				# ...and add the predecessor to finish it off
				push @{$con->{RING_PATH}}, $sendnum;
				
				# otherwise, add to the list of receivers
				push @receive, $connum;
			}
		}
		
		@send = @receive;
	}
	
	return $rings;
}

sub traverse {

	#<
	#? Recursive subroutine to find rings
	#. Called from find_atom_rings
	#; Requires: path
	#; Returns: nothing
	#>
	
	my @path = @_;

	my $atom = $path[0];
		
	# Reached an atom which we already know is not in a ring
	if ($far_atoms->[$atom]{NOTRING}) {
		return;
	}
	
	my $connected_atoms = $far_atoms->[$atom]{CONNECT};

	# No or one connected atom
	return if ($#{$connected_atoms} <= 0);
	
	# Path has looped back on itself
	foreach (@path[1..$#path-1]) {
		return if ($_ == $atom);
	}
	
	# We have found a ring.  Save it.
	if ($atom == $traverse_start && $#path > 0) {
	
		# Save if the ring is big enough
		save_ring(@path) if $#path > 2;
		return;
	}
	
	# Reached maximum depth. Go no further
	return if ($#path+1 == $traverse_depth);

	foreach (@$connected_atoms) {
	
		next if $#{$far_atoms->[$_]{CONNECT}} <= 0; # Connected atom has one or fewer bonded atoms
		next if ($_ == ($path[1] || -1)); # Skip if connected atom is the parent (making sure it is defined)
		traverse($_, @path);
	}
}

sub save_ring {
	
	my @path = @_;

	# Delete the last atom from the ring (@path contains first atom twice)
	my @ring_prime = @path[reverse (1..$#path)];

	# Have we already found this ring?
	foreach my $ring (@$traverse_rings) {
		return if @$ring == @ring_prime;
	}

	# Save the ring atoms in an array (@path contains first atom twice)
	@{$traverse_rings->[$traverse_ringcount]} = @path[0..$#path-1];
	++$traverse_ringcount;
}


sub delete_rings {

	#<
	#? Remove all rings from a molecule
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];

	#print "\n\ndelete_rings\n";
	#carp();

	my $atom;
	
	delete $mol->{ARINGS};
	delete $mol->{RINGS};
	delete $mol->{PRINGS};
	delete $mol->{KEKULE_BONDS};
	delete $mol->{AROMATIC_BONDS};
	delete $mol->{RINGSIZE_MAX};
	delete $mol->{BAD_RINGS};
	
	foreach my $atom (atoms($mol)) {
		
		delete $atom->{ARINGS};
		delete $atom->{PRINGS};
		delete $atom->{RINGS};
		delete $atom->{PLANAR_RING};
		delete $atom->{AROMATIC_RING};
		delete $atom->{NOTRING};
	}
}

#################################
# Routines to standardise rings #
#################################

sub standardise_ring {

	#<
	#? Order a ring (list of atom numbers) that is in a standard 'orientation'.
	#. Rearrange atoms in a ring list so that atom[0] has the lowest number
	#  and atom[1] has a lower number than atom[-1].
	#; Requires: ptr to array of atom numbers.
	#; Returns: ptr to normalised array of atom numbers.
	#>

	my @array = @{$_[0]};
	my $minatom = 99999;

	if ($#array < 0) {
		return ();
	}

	# Find lowest atom number
	foreach (@array) {
		$minatom = $_ if ($_ < $minatom);
	}
	
	# Move the lowest numbered atom to the start of the array
	while ($array[0] != $minatom) {
		unshift @array, (pop @array);
	}
	
	# Reverse array if required
	if ($array[1] > $array[-1]) {
		@array = reverse (@array);
		unshift @array, (pop @array);
	}

	return \@array;
}

sub standardise_ring_by_heteroatoms {

	#<
	#? Order a list of ring atoms so that is in a standard 'orientation'
	#  organised by hetero ring atoms and bondorders.
	#. Rearrange atoms in a ring list so that the highest atomic number has the lowest position
	#  in the ring with a subsequent requirement that the highest bondorder is at the lowest position
	#  in the ring.
	#; Requires: ptr to array of atom numbers.
	#; Returns: ptr to normalised array of atom numbers, string made of atom elements and bondorders -
	#  for example pyridine will be returned as 'N4C4C4C4C4C4'.
	#>

	my $array = $_[0];
	my $mol = $_[1];

	my @best_array;
	my $best_score = 0; # best score

	if ($#{$array} < 0) {
		return ();
	}

	# Find the best score in one ring direction
	for (0..$#{$array}) {
		
		# Cycle array;
		unshift @$array, (pop @$array);

		my $score = score_ring_het($array, $mol);

		# Keep best array
		if ($score gt $best_score) {
			@best_array = @$array;
			$best_score = $score;
		}
	}
	
	# Reverse the array
	@$array = reverse (@$array);

	# Repeat the process
	for (0..$#{$array}) {
		
		# Cycle array;
		unshift @$array, (pop @$array);

		# Score array
		my $score = score_ring_het($array, $mol);
		
		# Keep best array
		if ($score gt $best_score) {
			@best_array = @$array;
			$best_score = $score;
		}
	}
	
	# Produce as string that represents the atoms and bondorders in ring
	my $string = '';
	foreach my $i (0..$#best_array) {

		my $atom = $mol->{ATOMS}[$best_array[$i]];

		# Add element
		$string.= $atom->{ELEMENT};

		# Find  out if the atom is a bridgehead atom (ie one that connects two atoms
		# in the same large ring
		my $count = 0;
		CON: foreach my $connum (@{$atom->{CONNECT}} ) {

			foreach my $ratom (@best_array) {

				next if $ratom != $connum;
				++$count;
				next CON;
			}
		}

		$string .= "." if $count >= 3;
		
		my $j = $i+1;
		$j = 0 if $j > $#best_array;
		if (find_bondorder($best_array[$i], $best_array[$j], $mol) == 1) {
		
			# Single bond
			$string .= '-';
		} else {
			
			# Double or aromatic bond
			$string .= '=';
		}

	}

	return \@best_array, $string;
}

sub score_ring_het {
	
	#<
	#? Score a heterocyclic ring so that it can be sorted into a canonical order
	#. Score is a string made of atom element_numbers printed as C = 006, O = 008 etc.
	#  followed by a string of bondorders.  Using this system, benzene is represented
	#  as '006006006006006006444444'.
	#. With added bridgeheads!
	#; Requires: array of ring atom numbers, molecule
	#; Returns: score string
	#>

	my @ring = @{$_[0]};
	my $mol = $_[1];

	my $score = '';
	my $ringsize = $#ring;

	# Elements
	foreach my $i (0..$ringsize) {
			
		my $atom = $mol->{ATOMS}[$ring[$i]];
		$score .= sprintf "%03d", $atom->{ELEMENT_NUM};
	}

	# Bridgeheads (These are found by finding that three ring atoms are connected to the current atom)
	foreach my $i (0..$ringsize) {

		my $atom = $mol->{ATOMS}[$ring[$i]];

		my $count = 0;
		CON: foreach my $connum (@{$atom->{CONNECT}} ) {

			foreach my $ratom (@ring) {

				next if $ratom != $connum;
				++$count;
				next CON;
			}
		}
		if ($count >= 3) {
			$score .= '1';
		} else {
			$score .= '0';
		}
	}

	# Add bondorders
	foreach my $i (0..$ringsize) {
		
		my $j = $i + 1;
		$j = 0 if $j > $ringsize;
		$score.= find_bondorder($ring[$i], $ring[$j],$mol);
	}

	return $score;
}

##############################################
# Routines to find planar and aromatic rings #
##############################################

sub is_ring_planar {

	#<
	#? Robust ring planarity and aromaticity 
	#. This routine must be able to cope with: 
	#  molecules with no hydrogens and/or no bondorders,
	#  molecules with very poor geometry (eg 2D molecules) and molecules with distorted
	#  geometries (such as from high temperature MD simulations).
	#. This routine is expected to fail if the molecule has
	#  distorted geometry, bad bondorders and missing atoms
	#; Requires: molecule, ring (array of atoms)
	#; Returns: flag for planarity, flag for aromaticity
	#>

	my $mol = $_[0];
	my $ring = $_[1];
	
	my $aromatic = 0;
	my $planar = 0;
	my $starttime = (times)[0];

	mol_has_h($mol) if !defined $mol->{HAS_HYDROGENS};

	# Molecules that are flagged as having bad geometry
	if ($mol->{HAS_BAD_GEOM} || $mol->{HAS_NO_COORDS}) {
	
		# Molecules with 'good' bondorders
		if (is_ring_aromatic_by_bondorder($mol, $ring)) {
			
			return (1, 1);
		}
		
	} else {
		
		if (is_ring_aromatic_by_bondorder($mol, $ring)) {
			
			return (1, 1);
		}
		
		# Molecules with good geometry but poor bondorders
		if (is_ring_planar_by_dihedral($mol, $ring)) {

			$planar = 1;
			
			# Molecules with good geometry but possibly poor connectivity (incl. no H's)
			if (is_ring_aromatic_by_connectivity($mol, $ring, 0)) {

				$aromatic = 1;
			}
			
			return (1, $aromatic);
		}
		
		# Molecules with poor geometry and bondorders,
		# but with hydrogens (on carbon) present
		if ($mol->{HAS_HYDROGENS} == 2 && is_ring_aromatic_by_connectivity($mol, $ring, 1)) {
			return (1, 1);
		}
	}
	
	# If you've reached here, this ring is probably neither planar nor aromatic.
	return (0, 0);
}

sub is_ring_planar_by_dihedral {

	#<
	#? Finds planar rings by checking dihedral angles
	#. Tests dihedrals in rings
	#. Attempts to determine planarity (as distinct from
	#  aromaticity)
	#. Returns 0 if any single in-ring dihedral is larger
	#  than the hardcoded tolerance:
	#           18 for 5-membered rings
	#           26 for all other rings
	#; Requires: molecule (pointer to hash), $atoms
	#  (pointer to array of atom numbers)
	#; Returns: 2 if atoms are few enough in number that they must lie in a plane (<= 3)
	#, 1 if there are enough atoms that they could be non-planar, but are planar (within tolerance)
	#, 0 otherwise (i.e., not planar within tolerance)
	#>
	
	my $mol = $_[0];
	my $atoms = $_[1];

	my $radian = 57.29577951308232088;
	my $tolerance;
	my $torsion;
	
	if ($#{$atoms} == 4) {
					
		$tolerance = 18;	# Tolerance for 5 membered rings
	} else {
	
		$tolerance = 26;	# Default tolerance for planarity in degrees
					# This angle is large to accommodate molecules
					# like porphyrins (rings >= 6 members)
	}
	
	# Three or less atoms necessarily lie in a single plane.
	return 2 if ($#{$atoms} < 3);

	# Make double length array
	my @a = (@$atoms, @$atoms);

	for (my $i=0; $i <= $#{$atoms}; ++$i) {
		
		my $string;
		
		my $torsion = $radian * dihedral(	$mol->{ATOMS}[$a[$i]],
							$mol->{ATOMS}[$a[$i+1]],
							$mol->{ATOMS}[$a[$i+2]],
							$mol->{ATOMS}[$a[$i+3]]	);

		# Account for torsions that are near 180 (eg porphyrin)
		$torsion = $torsion - 180 if ($torsion > 90);
		$torsion = $torsion + 180 if ($torsion < -90);

		$string = sprintf "(".subname().") torsion %3d %3d %3d %3d %6.2f tolerance %6.2f\n",$a[$i], $a[$i+1], $a[$i+2], $a[$i+3],abs($torsion), $tolerance;
		
		return 0 if (abs($torsion) > $tolerance);
	}
	
	return 1;
}


sub is_ring_aromatic_by_bondorder {

	#<
	#? Find aromatic rings ie those with all aromatic bonds or 4n+2 pi electrons
	#. Subroutine assumes that it is operating on a planar ring
	#; Requires: molecule , ring
	#; Returns: 1 if aromatic, 0 otherwise
	#>

	my $mol = $_[0];
	my $ring = $_[1];
	
	my @borders;
	my $rcount;
	
	#print "is_ring_aromatic_by_bondorder\n";
	
	$rcount = $#{$ring}+1; # Number of atoms in ring
	
	# Return 0 if we dont have 4n+2 atoms or 5 atoms
	return 0 if !($rcount == 5 || ($rcount-2)%4 == 0);
	
	my $flag;
	my $count_O_sp3 = 0;
	my $count_N_sp3 = 0;
	
	foreach my $atomnum (@$ring) {
	
		my $numsingle = 0;
		my $numdouble = 0;
	
		my $atom = $mol->{ATOMS}[$atomnum];
		my $el = $atom->{ELEMENT};
		
		foreach my $bo (@{$atom->{BORDERS}}) {
			++$numsingle if $bo == 1;
			++$numdouble if $bo == 2 || $bo == 4;
		}
		
		$flag = 1;
		if ($el eq 'C') {
			
			next if $numdouble == 1 or $numdouble == 2;
			$flag = 0;
			#print "here1\n";
			last;
		}
		
		$flag = 1;
		if ($el eq 'O' || $el eq 'S' || $el eq 'Se') {
		
			++$count_O_sp3;
			$flag == 0 && last if $numsingle != 2;
			next if $rcount == 5 && $count_O_sp3 == 1;
			next if $rcount == 5 && $count_O_sp3 <= 2;
			$flag = 0;
			last; 
		}
		
		if ($el eq 'N' || $el eq 'P'  || $el eq 'As') {
		
			++$count_N_sp3 if $numsingle+$numdouble == 3;
			$flag == 0 && last if $numsingle+$numdouble == 4; # Quaternary!
			
			next if $rcount == 5 && $count_N_sp3 == 1;
			next if $rcount == 5 && $count_N_sp3 <= 2;
			next if $rcount == 6 && ($numdouble == 1 or $numdouble == 2);
			$flag = 0;
			#print "here3\n";
			last; 
		}
	}
	
	#print "Flag $flag\n";

	return $flag;
}


sub is_ring_aromatic_by_connectivity {

	#<
	#? Find 4n+2 member aromatic rings (ie those with all aromatic bonds) 
	#  by connectivity where bondorders and hydrogens are present or missing
	#. Subroutine assumes that it is operating on a planar ring
	#; Requires: Molecule, ring, flag to indicate hydrogens are present
	#; Returns: 1 if aromatic, 0 otherwise
	#>

	my $mol = $_[0];
	my $ring = $_[1];
	my $h_present = $_[2];
	
	# Calculate the ring size
	my $size = $#{$ring} + 1;
	
	# If the ring is not 4n + 2 we are not aromatic
	return 0 if ($size-2)/4 - int(($size-2)/4) != 0;

	foreach my $anum (@$ring) {
	
		my $atom = $mol->{ATOMS}[$anum];
		
		if ($atom->{ELEMENT} eq 'C') {
		
			# Not aromatic if hydrogens present and 
			# any carbon has other than three connected atoms
			return 0 if $h_present && $#{$atom->{CONNECT}} != 2;

			# Not aromatic if hydrogens not present and 
			# any carbon has more than three connected atoms
			return 0 if $h_present && $#{$atom->{CONNECT}} > 2;
			next;
		} 

		if ($atom->{ELEMENT} eq 'N') {
	
			# OK if three connected atoms and charged
			next if $#{$atom->{CONNECT}} == 3 && $atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE} == 1;
			# Not aromatic if nitrogen has more than two connected atoms 
			return 0 if $#{$atom->{CONNECT}} > 1;
			next;
		}

		if ($atom->{ELEMENT} eq 'O' || $atom->{ELEMENT} eq 'S') {

			# Not aromatic atom explicitly has 
			# a positive charge (unlikely!)
			return 0 if !($atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE} == 1);
			next;
		}

		# Other elements are probably not aromatic...
		return  0;
	}
	
	return 1;
}


###########################
# Planar atom subroutines #
###########################


sub molecule_find_planar_atoms {

	#<
	#? Call geometry or bondorder based routines to find planar atoms
	#. Uses geometry routine by default
	#. Note that the geometry and bondorder routines do not return the same
	#  sets of atoms.  Also the results will vary if hydrogens are missing
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];

	if ($mol->{HAS_NO_COORDS} && $mol->{HAS_BAD_BONDORDERS}) {
		return molecule_find_planar_atoms_connections($mol);
	} elsif (!$mol->{HAS_BAD_GEOM}) {
		return molecule_find_planar_atoms_geom($mol);
	} elsif (!$mol->{HAS_BAD_BONDORDERS}) {
		return molecule_find_planar_atoms_bondorders($mol);
	} else {
		silico_msg('w', "Molecule has neither good geometry nor good bondorders! No planar atoms found.\n");
		return;
	}
}

sub molecule_find_planar_atoms_geom {

	#<
	#? Find planar atoms using molecule geometry or bondorder
	#. A planar atom is defined as one with having a double bond,
	#  aromatic bond or three substituents and
	#  an improper dihedral less than 15 degrees from planar.
	#. Does not find planar tetracoordinate atoms.
	#; Returns: an array of planar atoms and sets $mol->{PLANAR_ATOMS}
	#>

	my $mol = $_[0];

	my $radian = 57.29577951308232088;
	my $tol = 30; # Deviation from planarity (was originally 25)

	ATOM: for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {
	
		my $atom = $mol->{ATOMS}[$i];

		next if (!defined $atom->{CONNECT});

		my @con = @{$atom->{CONNECT}};

		# Not 3 connected atoms
		next if ($#con != 2);
		
		foreach my $bo (@{$atom->{BORDERS}}) {

			# Check for double or aromatic bond
			next if ($bo != 2 && $bo != 4);
		
			$atom->{PLANAR_ATOM} = 1;
			push @{$mol->{PLANAR_ATOMS}}, $atom;
			next ATOM;
		}

		my $atom1 = $mol->{ATOMS}[$con[0]];
		my $atom2 = $atom;
		my $atom3 = $mol->{ATOMS}[$con[1]];
		my $atom4 = $mol->{ATOMS}[$con[2]];
		my $dihedral = 180 - $radian * abs(dihedral($atom1, $atom2, $atom3, $atom4));

		if ($dihedral <= $tol) {
			$atom->{PLANAR_ATOM} = 1;
			push @{$mol->{PLANAR_ATOMS}}, $atom;
			next ATOM;
		}
	}
	
	return $mol->{PLANAR_ATOMS};
}

sub molecule_find_planar_atoms_bondorders {

	#<
	#? Find planar atoms using bondorders only
	#. All atoms which have double bonds are planar.  This
	#  includes S=O and P=O which may cause some problems
	#; Returns: an array of planar atoms and sets $mol->{PLANAR_ATOMS}
	#>

	my $mol = $_[0];

	ATOM: foreach my $atom (atoms($mol)) {

		next if (!defined $atom->{CONNECT});
		
		next if ($atom->{CONNECT} <= 0);

		foreach my $bo (@{$atom->{BORDERS}}) {

			if ($bo == 2 || $bo == 4) {
				$atom->{PLANAR_ATOM} = 1;
				push @{$mol->{PLANAR_ATOMS}}, $atom;
				next ATOM;
			}
		}
	}
	
	return $mol->{PLANAR_ATOMS};
}

sub molecule_find_planar_atoms_connections {

	#<
	#? Find planar atoms using connected atoms only (assumes no coordinates)
	#; Returns: an array of planar atoms and sets $mol->{PLANAR_ATOMS}
	#>

	my $mol = $_[0];
	
	silico_msg('d', "Unfinished subroutine\n");

	ATOM: foreach my $atom (atoms($mol)) {
	
		my $ncon = $#{$atom->{CONNECT}}+1;

		if ($atom->{ELEMENT} eq 'C') {
		
			if ($ncon == 3) {
				atom_printout($atom);
				die "Incomplete subroutine";
				push @{$mol->{PLANAR_ATOMS}}, $atom;
				next ATOM;
			}
		}
		
		if ($atom->{ELEMENT} eq 'O') {
			if ($ncon == 2) {
				atom_printout($atom);
				die "Incomplete subroutine";
				push @{$mol->{PLANAR_ATOMS}}, $atom;
				next ATOM;
			}
			
		}
	}
	
	return $mol->{PLANAR_ATOMS};
}

#######################################
# Routines to convert ring bondorders #
#######################################


sub make_aromatic_bonds {

	#<
	#? Make aromatic bonds (bondorder 4) in all 6-membered aromatic rings
	#. Rings must contain only C and N
	#; Requires: molecule, kekule flag
	#; Returns: nothing
	#>

	my $mol = $_[0];
	
	(undef, undef, my $aromatic_rings) = molecule_find_rings($mol, 6);

	RING:foreach my $ring (@$aromatic_rings) {

		next if $#{$ring} != 5;

		foreach my $i (0..$#{$ring}) {
		
			my $atom1 = $mol->{ATOMS}[$ring->[$i-1]];
			my $atom2 = $mol->{ATOMS}[$ring->[$i]];
	
			bond_modify_order($mol, $atom1, $atom2, 4);
		}
	}

	delete $mol->{KEKULE_BONDS};
	$mol->{AROMATIC_BONDS} = 1;
}

sub convert_aromatic_bonds_kekule {

	#<
	#? Convert rings with aromatic bonds to kekule representation
	#. Converts systems up to 18 atoms
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	
	my ($rings, $prings, $aromatic_rings) = molecule_find_rings($mol, 18);

	return if !defined $aromatic_rings;
	
	# Sort biggest rings to smallest
	@$aromatic_rings = sort {$#{$b} <=> $#{$a}} @$aromatic_rings;

	RING: foreach my $ring (@$aromatic_rings) {

		my $n = $#{$ring} +1;

		# Ring must have 2n+2 atoms;
		my $m = ($n-2)/2;
		next if $m != int($m);

		# Skip rings which are not made of aromatic bonds
		my @bos = find_bondorders_ring($ring, $mol);
		
		foreach (@bos) {
		
			silico_msg('d', "Found a RING record with nonbonded atoms. Has the molecule been sorted?\n") if $_ == -1;
			next RING if $_ != 4;
		}
		
		my $bo = 1;
		foreach my $i (0..$#{$ring}) {
		
			my $atom1 = $mol->{ATOMS}[$ring->[$i-1]];
			my $atom2 = $mol->{ATOMS}[$ring->[$i]];
	
			bond_modify_order($mol, $atom1, $atom2, $bo);

			# Flip bond
			$bo = !($bo-1)+1;
		}
	}

	$mol->{KEKULE_BONDS}= 1;
	delete $mol->{AROMATIC_BONDS};
}




###########################################
# Routines for ring bondorders & elements #
###########################################

sub find_bondorders_ring {

	#<
	#? Return array containing bondorders of bonds in ring
	#; Requires: ring (array of atom numbers), molecule
	#; Returns: array of bondorders in ring
	#>

	my $ring = $_[0];
	my $mol = $_[1];

	my @bolist;

	my $i = -1;
	foreach my $atomnum (@$ring) {

		++$i;
		my $atom = $mol->{ATOMS}[$atomnum];
		my $prevatomnum = $ring->[$i-1]; # Relies on index -1 being the last element of array
		my $prevatom = $mol->{ATOMS}[$prevatomnum];
		
		# Find the bond order of connection to prevatom
		my $bo = find_bondorder($atomnum, $prevatomnum, $mol);
		
		push @bolist, $bo;
	}
	
	return @bolist;
}

sub find_elements_ring {

	#<
	#? Return array containing elements of atoms in ring
	#  and a hash counting the number of times each element occurs
	#.
	#; Requires: ring, molecule
	#; Returns: pointer to array of elements, hash containing a count of each element
	#>

	my $ring = $_[0];
	my $mol = $_[1];

	my $elist;
	my $hash;

	foreach my $atomnum (@$ring) {

		my $atom = $mol->{ATOMS}[$atomnum];

		push @$elist, $atom->{ELEMENT};
		++$hash->{$atom->{ELEMENT}};
	}

	return $elist, $hash;
}

sub shared_rings {
	
	#<
	#? Counts the number of common rings a pair of atoms share.
	#; Requires: Molecule, atom 1, atom 2
	#; Returns: Number of shared rings
	#>
	
	my $mol = $_[0];
	my $atom1 = $_[1];
	my $atom2 = $_[2];
	
	my $count = 0;
	foreach my $ring (@{$atom1->{RINGS}}) {
		
		my $count2 = 0;
		foreach my $ringatomnum (@$ring) {
			next if $mol->{ATOMS}[$ringatomnum] != $atom2;
			++$count2;
		}
		
		++$count if ($count2 > 0);
	}
	
	return $count;
}

sub ring_printout {

	my $ring = $_[0];
	my $mol = $_[1];
	
	print "ring: $ring ";
	
	foreach (@$ring) {
		
		my $atom = $mol->{ATOMS}[$_];
		print "$atom->{ELEMENT} ";
	}
	
	print "\n";
}

return 1;
