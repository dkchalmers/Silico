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
#! silico_hb.pm
#? Silico hydrogen bond routines
#. $Revision: 1.1.2.2 $
#>

use strict;
package Silico;

##################################################################
#
#       Hydrogen bonds
#
##################################################################


sub mol_find_donors_acceptors {

	#<
	#? Count molecule H-bond donors and acceptors and label them
	#. Requires that hydrogens are present on molecule
	#. Donors are defined as N, O or S with attatched H.
	#. Acceptors: are N, O, but not in a planar 5-memebered ring. 
	#  N is counted provided valence < 4 except for amides; 
	#  S is counted provided valence < 3;
	#. Atoms are labelled as HBACCEPTOR and/or HBDONOR
	#; Requires: molecule
	#; Returns: number of donors (scalar), number of acceptors (scalar)
	#: Sets:  $mol->{HBDONORS} and $mol->{HBACCEPTORS} (atomlists)
	#>

	my $mol = $_[0];

	my $numDonor = 0;
	my $numAccept = 0;

	# Find rings.
	molecule_find_rings($mol, 10);
	mol_label_functional_group($mol) if !$mol->{CHARGE_GROUP};
	
	@{$mol->{HBDONORS}} = ();
	@{$mol->{HBACCEPTORS}} = ();
	
	foreach my $atom (atoms($mol)) {
	
		my $el = $atom->{ELEMENT};
		
		next if $el eq 'C' || $el eq 'H';
	
		# Count number of connected hydrogens
		# (was "Count number of nonhydrogen connected atoms"?)
		my $conH = 0;
		foreach my $con (connected($atom, $mol)) {

			++$conH if $con->{ELEMENT} eq 'H';
		}
		
		if ($el eq 'N') {
		
			if ($conH) {
				++$numDonor;
				$atom->{HBDONOR} = 1;
				push @{$mol->{HBDONORS}}, $atom;
			}
			
			next if valence($atom) >= 4;
			next if $atom->{FG}{AMIDE_N};
			
			if ($atom->{PRINGS} && $atom->{PRINGS} == 4) {
				next;
			}

			++$numAccept;
			$atom->{HBACCEPTOR} = 1;
			push @{$mol->{HBACCEPTORS}}, $atom;
		}

		if ($el eq 'O') {
			next if $atom->{PRINGS} && $atom->{PRINGS} == 4;
			
			if ($conH) {
				++$numDonor;
				$atom->{HBDONOR} = 1;
				push @{$mol->{HBDONORS}}, $atom;
			}
			if ($atom->{PRINGS} && $atom->{PRINGS} == 4) {
				next;
			}

			++$numAccept;
			$atom->{HBACCEPTOR} = 1;
			push @{$mol->{HBACCEPTORS}}, $atom;
		}

		if ($el eq 'S' ) {
			
			if ($conH) {
				++$numDonor;
				$atom->{HBDONOR} = 1;
				push @{$mol->{HBDONORS}}, $atom;
			}
			
			next if valence($atom) > 2;
			
			if ($atom->{PRINGS} && $atom->{PRINGS} == 4) {
				next;
			}

			++$numAccept;
			$atom->{HBACCEPTOR} = 1;
			push @{$mol->{HBACCEPTORS}}, $atom;
		}
	}
	
	silico_msg('g', "Found $numDonor H-bond donors and $numAccept H-bond acceptors\n");
	return $numDonor, $numAccept;
}

sub is_hbond_mcdonald_thornton {

	#<
	#? Identify a hydrogen bond (D-H:::A) according to
	#  the definition used by McDonald and Thornton
	#, J. Mol. Biol. 1994, 238, 777-793
	#; Requires: donor molecule, donor atom, acceptor molecule, acceptor atom
	#  (donor cutoff distance, H cutoff distance, dha angle max, HAR angle max - 
	#  defaults 3.9, 2.5, 90, 90)
	#; Returns: array of the hbonds. Each record contains (dmol, datom, hydrogen, amol, aatom)
	#>
	
	my $dmol = $_[0];
	my $donor = $_[1];
	my $amol = $_[2];
	my $acceptor = $_[3];
	my $dcutoff = $_[4];
	my $hcutoff = $_[5];
	my $dha_cutoff = $_[6];
	my $har_cutoff = $_[7];
	
	my $hnum;
	my $list;
	
	$dcutoff ||= 3.9;
	$hcutoff ||= 2.5;
	$dha_cutoff ||= 90;
	$har_cutoff ||= 90;
	
	@$list = ();
	
	# Check that the donor and acceptor are close enough together
	my $ddist = distance($donor, $acceptor);
	
	return $list if $ddist > $dcutoff;
	
	H: foreach my $h (connected($donor, $dmol)) {
	
		next if $h->{ELEMENT} ne 'H';
		my $hbond;
		
		my $hdist = distance($h, $acceptor);
		next if $hdist > $hcutoff;
	
		my $dha_angle = bndangle($donor, $h, $acceptor);
		next if $dha_angle < $dha_cutoff;
	
		foreach my $r (connected($acceptor, $amol)) {
			
			my $har_angle = bndangle($h, $acceptor, $r);
			next H if $har_angle < $har_cutoff;
		}
		
		@$hbond = ($dmol, $donor, $h, $amol, $acceptor); 
		
		push @$list, $hbond;
	}
	
	return $list;
}

sub find_intramolecular_hbonds_h {

	#<
	#? Find hydrogen bonds between atoms within the one molecule
	#  using hydrogen positional information
	#. mol_find_donors_acceptors must be run prior to
	#  commencing this routine
	#; Requires: molecule
	#; Optional: D-A cutoff, H-A cutoff
	#, D-H-A angle cutoff, H-A-R angle cutoff
	#, flag to include waters
	#; Returns: list of hydrogen bonds
	#; Sets: $atom->{HBONDS} - a list of hydrogen bonds belonging to each atom
	#>
	
	my $mol = $_[0];
	my $dcutoff = $_[1];
	my $hcutoff = $_[2];
	my $dha_cutoff = $_[3];
	my $har_cutoff = $_[4];
	my $water = $_[5];
	
	
	my $dha_angle ;
	my $hblist;
	my $hash;
	
	mol_find_donors_acceptors($mol) if !$mol->{HBDONORS};
	
	foreach my $acceptor (@{$mol->{HBACCEPTORS}}) {
	
		next if ($acceptor->{FG}{W} && !$water);
	
		foreach my $donor (@{$mol->{HBDONORS}}) {
	
			next if $acceptor == $donor;
			next if ($donor->{FG}{W} && !$water);

			# Only count pairs once
			next if $hash->{"$donor->{NUM}_$acceptor->{NUM}"};
		
			#Returns: array of the several D-H-A triplets that are found each contains (dmol, datom, hydrogen, amol, aatom)
			my $list = is_hbond_mcdonald_thornton( $mol, $donor, $mol, $acceptor, $dcutoff, $hcutoff, $dha_cutoff, $har_cutoff);
			
			next if $#{$list} == -1;

			foreach my $hbond (@$list) {
			
				push @$hblist, $hbond;
				push @{$donor->{HBONDS}}, $hbond;
				push @{$acceptor->{HBONDS}}, $hbond;
				#atom_printout($donor);
			}

			++$hash->{"$donor->{NUM}_$acceptor->{NUM}"};
		}
	}
	
	return $hblist;
}

sub find_intermolecular_hbonds_h {

	#<
	#? Find hydrogen bonds between two molecules
	#  using hydrogen atom information
	#. mol_label_donors_acceptors must be run
	#  before commencing this routine
	#; Requires: Two molecules
	#; Optional: D-A distance, H-A distance
	#, D-H-A angle, H-A-R angle
	#; Returns: List of hydrogen bonds (D-H-A)
	#>
	
	my $mol1 = $_[0];
	my $mol2 = $_[1];
	my $dcutoff = $_[2];
	my $hcutoff = $_[3];
	my $dha_cutoff = $_[4];
	my $har_cutoff = $_[5];
	
	my $hblist;
	my $mola;
	my $molb;	
	
	mol_find_donors_acceptors($mol1) if !$mol1->{HBDONORS};
	mol_find_donors_acceptors($mol2) if !$mol2->{HBDONORS};

	foreach my $i (0..1) {
	
		if ($i) {
			$mola = $mol1;
			$molb = $mol2;
		} else {
			$mola = $mol2;
			$molb = $mol1;
		}
	
		foreach my $donor (@{$mola->{HBDONORS}}) {
	
			foreach my $acceptor (@{$molb->{HBACCEPTORS}}) {
		
				#Returns: array of the several D-H-A triplets that are found each contains (dmol, datom, hydrogen, amol, aatom)
				my $list = is_hbond_mcdonald_thornton(
					$mola, $donor, $molb, $acceptor,
					$dcutoff, $hcutoff, $dha_cutoff, $har_cutoff);
			
				foreach my $hbond (@$list) {
					push @$hblist, $hbond;
				}
			}
		}
	}
	
	return $hblist;
}

sub find_intermolecular_hbonds_h_simple {

	#<
	#? Find hydrogen bonds between two molecules using 
	#  information about the position of the hydrogen atom
	#. Uses simple cutoff between donor hydrogen and acceptor
	#. You must call mol_find_donors_acceptors on each molecule first
	#>

	my $mol1 = $_[0];
	my $mol2 = $_[1];
	my $cutoff = $_[2] || 4.5;

	my $list1_accept;
	my $list1_donorh;
	my $list2_accept;
	my $list2_donorh;
	my $cutoff_sq;
	my $hblist;
	my $con;
	my $conatom;
	
	$cutoff_sq = $cutoff*$cutoff;

	# Make atom lists to speed things up
	foreach my $atom1 (atoms($mol1)) {

		if ($atom1->{HBACCEPTOR}) {
		
			push @$list1_accept, $atom1;
		}
		
		if ($atom1->{HBDONOR}) {
		
			foreach $con (@{$atom1->{CONNECT}}) {
				$conatom = $mol1->{ATOMS}[$con];
				
				next if ($conatom->{ELEMENT} ne 'H');
				
				push @$list1_donorh, $conatom;
			}
		}
	}
	
	foreach my $atom2 (atoms($mol2)) {

		if ($atom2->{HBACCEPTOR}) {
		
			push @$list2_accept, $atom2;
		}
		
		if ($atom2->{HBDONOR}) {
		
			foreach $con (@{$atom2->{CONNECT}}) {
				$conatom = $mol2->{ATOMS}[$con];
				
				next if ($conatom->{ELEMENT} ne 'H');
				
				push @$list2_donorh, $conatom;
			}
		}
	}
	
	foreach my $atom1 (@$list1_accept) {

		foreach my $atom2 (@$list2_donorh) {
		
			next if distance_sq($atom1, $atom2) > $cutoff_sq;
			next if distance($atom1, $atom2) > $cutoff;
			
			my $pair;
			@$pair = ($atom1, $atom2);
			push @$hblist, $pair;
		}
	}
	
	foreach my $atom1 (@$list1_donorh) {

		foreach my $atom2 (@$list2_accept) {
		
			next if distance_sq($atom1, $atom2) > $cutoff_sq;
			next if distance($atom1, $atom2) > $cutoff;
		
			my $pair;
			@$pair = ($atom1, $atom2);
			push @$hblist, $pair;
		}
	}
	return $hblist;
}

sub find_intermolecular_hbonds_noh {

	#<
	#? Find hydrogen bonds between two molecules _without_ using 
	#  information about the position of the hydrogen atom
	#. Uses simple cutoff between donor and acceptor heteroatoms
	#. You must call mol_find_donors_acceptors on each molecule first
	#>

	my $mol1 = $_[0];
	my $mol2 = $_[1];
	my $cutoff = $_[2] || 4.5;

	my $list1;
	my $list2;
	my $cutoff_sq;
	my $hblist;
	
	$cutoff_sq = $cutoff*$cutoff;

	# Make atom lists to speed things up
	foreach my $atom1 (atoms($mol1)) {

		next if !($atom1->{HBDONOR} || $atom1->{HBACCEPTOR});
		
		push @$list1, $atom1;
	}
	
	foreach my $atom2 (atoms($mol2)) {

		next if !($atom2->{HBDONOR} || $atom2->{HBACCEPTOR});
		
		push @$list2, $atom2;
	}
	
	foreach my $atom1 (@$list1) {

		foreach my $atom2 (@$list2) {
		
			my $pair;
		
			next if distance_sq($atom1, $atom2) > $cutoff_sq;
			next if distance($atom1, $atom2) > $cutoff;
		
			if ($atom2->{HBACCEPTOR} && $atom1->{HBDONOR}) {
				@$pair = ($atom1, $atom2);
				push @$hblist, $pair;
				next;
			}
			if ($atom1->{HBACCEPTOR} && $atom2->{HBDONOR}) {
				@$pair = ($atom1, $atom2);
				push @$hblist, $pair;
				next;
			}
		}
	}

	return $hblist;
}

return 1;
