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
#! silico_molecules.pm
#? Silico routines which operate on any atom, molecule, or ensemble.
#. $Revision: 1.15-67-gf7b58c9 $
#>

use strict;
package Silico;

####################################################
####################################################
###                                              ###
###     SUBROUTINES DEALING WITH AN ENSEMBLE     ###
###                                              ###
####################################################
####################################################

#############
# ENSEMBLES #
#############

sub ensemble_cluster {

	#<
	#? Clusters molecules based on RMS distance from the lowest energy
	#  structure.
	#
	#. Assumes that molecules are sorted by energy.  Sets the value
	#  $mol->{CLUSTER_COUNT} to the number of molecules in that cluster.
	#; Requires: ensemble, threshold for clustering, flag to keep all molecules, flag to rename molecules with cluster membership
	#; Returns: ensemble containing the first structure from each cluster, or all memebers of each cluster
	#>

	my $ensemble = $_[0];
	my $threshold = $_[1];
	my $keep_all = $_[2];
	my $rename = $_[3];

	my $all;
	my $clusters; 
	my $parents;

	if (!defined $threshold) {
		silico_msg('d', "No threshold was specified!\n");
	}

	my $i = 1;
	foreach my $mol1 (@$ensemble) {
		$mol1->{NUM} = $i;
		delete $mol1->{CLUSTER_COUNT};
		++$i;
	}

	my $cluster_count = -1;
	$i = -1;
	foreach my $mol1 (@$ensemble) {

		++$i;
		# Skip if this molecule belongs to a cluster already
		next if (defined $mol1->{CLUSTER_COUNT});

		++$cluster_count;

		$mol1->{CLUSTER_COUNT} = $cluster_count+1;
		$mol1->{RMS} = 0;

		# RMS from cluster 1
		my $rms = rms($mol1, $ensemble->[0]);
		$mol1->{RMS2} = $rms; 	

		push @{$clusters->[$cluster_count]}, $mol1;

		foreach my $mol2 (@$ensemble) {

			# Skip if this molecule belongs to a cluster already
			next if (defined $mol2->{CLUSTER_COUNT});
			
			my $rms = rms($mol1, $mol2);

			if ($rms < $threshold && $rms >= 0) {
				$mol2->{CLUSTER_COUNT} = $cluster_count+1;
				$mol2->{RMS} = sprintf "%.2f", $rms;
				push @{$clusters->[$cluster_count]}, $mol2;
			}
		}
	}

	# Print out the RMS matrix for debugging
	if ($Silico::debug) {
		silico_msg('g', heading("RMS matrix\n"));
	
		my $i = 0;
		foreach my $mol1 (@$ensemble) {
			silico_msg('g', "$i\t");
			foreach my $mol2 (@$ensemble) {
				silico_msg('g', sprintf("%4.2f  ", rms($mol1, $mol2)));
			}
			silico_msg('g', "\n");
			++$i;
		}
		silico_msg('g', "\n");
	}

	# Print out the membership of each cluster
	silico_msg('c', 
		heading("Cluster membership and energies\n"),
		"Using threshold $threshold Angstroms RMS\n",
		"\n");

	printf "          %10s %15s\n", "Energy", "RMS to clust 1";

	$i = 0;
	foreach my $cluster (@$clusters) {

		++$i;

		if (!$keep_all) {
			push @$parents, $cluster->[0];
		} else {
			foreach (@$cluster) {
				push @$all, $_;
			}
		}

		my $en;
		if (defined $cluster->[0]{ENERGY}) {
			$en = sprintf  "%6.2f", $cluster->[0]{ENERGY};
		} else {
			$en = '-';
		}

		printf "cluster (%s) %6s      %6.2f       ", $cluster->[0]{CLUSTER_COUNT}, $en, $cluster->[0]{RMS2};

		my $j = 0;
		foreach my $mol1 (@$cluster) {

			++$j;
			silico_msg ('c', "  ".($mol1->{NUM}));
			my $name = $mol1->{NAME};
			$name ||= '';
			if (length $name > 25) {
				$name = substr ($name, 0, 22) if length $name > 10;
				$name .= "...";
			}
			$name = sprintf $name."_%d_%d", $i, $j if $rename;
			silico_msg ('c', " ($name, $mol1->{RMS})");
			$mol1->{NAME} = $name if $rename;
		}
		silico_msg('c', "\n");
	}

	# Return all members of each closter
	return $all if $keep_all;

	# or return the parent of each cluster
	return $parents;
}

sub ensemble_connect_atoms {

	#<
	#? Connect all atoms in an ensemble by distance
	#  (only within each molecule though)
	#; Calls: Connect_atoms_by_distance
	#; Requires: Ensemble
	#; Returns: Nothing
	#>

	my $ensemble = $_[0];
	my $mol;

	foreach my $molecule (@$ensemble) {
		connect_atoms_by_distance($mol);
	}

	return undef;
}

sub ensemble_consolidate{

        #<
	#? Combine all molecules in an ensemble into a single molecule
        #. Rewritten version. DKC Septemer 2025
        #. Works with gromacs .itp files 
        #; Requires: Ensemble, flag to prevent use of deep_copy
        #; Returns: New ensemble with all atoms in the molecule $newens->[0]
        #>
        
        my $ens = $_[0];

        croak("Thing to subroutine is not an ensemble\n") if ref $ens ne 'ARRAY';
	return () if !defined $ens->[0]; # Empty ensemble
        
        my @fields = qw(CONNECT GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS GMX_EXCLUSIONS);

	my $mol0 = $ens->[0];
        
	my $newmol;
	$newmol->{NUMATOMS} = 0;
	@{$newmol->{ATOMS}} = ();
        
	my $i = 0;
	foreach my $mol (@$ens) {
	
		my $offset = $i;
		
		foreach my $atom (atoms($mol)) {
		
			++$i;
		
			foreach my $field (@fields) {
				
				next if !defined $atom->{$field};
				
				my $j = 0;
				foreach (@{$atom->{$field}}) {
					
					$atom->{$field}[$j] += $offset;
					++$j;
				}
			}
			
			push @{$newmol->{ATOMS}}, $atom;
			$atom->{NUM} = $i;
		}
        }
	
	$newmol->{NUMATOMS} = $i;
	$newmol->{NUMBONDS} = molecule_count_numbonds($newmol);
	
        my $newens;
        $newens->[0] = $newmol;

        return $newens;
}


sub ensemble_consolidate_old {

	#<
	#? Combines all molecules in an ensemble into a single molecule
	#. Works with gromacs .itp files
	#; Requires: Ensemble, flag to prevent use of deep_copy
	#; Returns: New ensemble with all atoms in the molecule $newens->[0]
	#>
	
	my $ens;
	my $flag = $_[1];

	if ($flag) {
		$ens = $_[0];
	} else {
		$ens = deep_copy($_[0]);
	}
	
	(ref ($ens) eq 'ARRAY' && ref ($ens->[0]) eq 'HASH') || silico_msg('d', "Error in data passed to subroutine\n");
	
	my @fields =  qw(ATOMS GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS GMX_EXCLUSIONS);

	my $newmol = molecule_copy_noatoms($ens->[0]);

	my $count = 0;
	foreach my $mol (@$ens) {
	
		++$count;
		# Renumber atoms from subsequent molecules
		moleculeFmo2($mol, $newmol->{NUMATOMS}+1) if $count > 1;

		# Combine ATOMS and GROMACS itp ANGLES etc
		foreach( @fields ) {
			next if !defined $mol->{$_};
			push @{$newmol->{$_}}, @{$mol->{$_}};
		}

		$newmol->{NUMATOMS} = $#{$newmol->{ATOMS}}+1;
		$newmol->{NUMBONDS} += $mol->{NUMBONDS};
	}
	
	# Give the ensemble the name of the first molecule in original ensemble
	$newmol->{NAME} = $ens->[0]{NAME};
	$newmol->{SOURCE_FILE_NAME} = $ens->[0]{SOURCE_FILE_NAME};
	
	my $newens;
	$newens->[0] = $newmol;
		
	return $newens;
}

sub ensemble_copy {

	#<
	#? Creates a copy of an ensemble.
	#; Calls: Molecule_copy
	#; Requires: Ensemble
	#; Returns: Duplicate ensemble as an arrayref
	#>

	my $molecules = $_[0];
	my $newmolecules;
	
	for (my $i = 0; $i <= $#{$molecules}; ++$i) {
		$newmolecules->[$i] = deep_copy( $molecules->[$i]);
	}
	return $newmolecules;
}

sub ensemble_renumber {

	#<
	#? Renumber all atoms in an ensemble.
	#; Calls: Molecule_renumber2
	#; Requires: Ensemble
	#; Returns: Nothing
	#>

	my $molecules = $_[0];

	my $mol;
	
	foreach my $mol (@$molecules) {
		molecule_renumber2($mol);
	}

	return undef;
}


###################################################
###################################################
###                                             ###
###     SUBROUTINES DEALING WITH A MOLECULE     ###
###                                             ###
###################################################
###################################################


#####################
# MOLECULE CREATION #
#####################

sub create_molecule {
	
	#<
	#? Create a new molecule.
	#; Requires: Optional molecule name
	#; Returns: molecule
	#>
	
	my $name = $_[0] || 'Unnamed';
	my $mol;
	
	$mol->{NAME} = $name;
	$mol->{NUMATOMS} = 0;
	$mol->{NUMBONDS} = 0;
	@{$mol->{ATOMS}} = ();
	
	return $mol;
}

sub deep_copy_orig {

	#<
	#? Recursive deep copy routine for data structures
	#. Taken from a tutorial by Randall Schwartz
	#; Requires: Data structure
	#; Returns: copy of data structure
	#>
	
	my $this = shift;
	
	if (not ref $this) {
		return $this;
	} elsif (ref $this eq "ARRAY") {
		[map deep_copy($_), @$this];
	} elsif (ref $this eq "HASH") {
		+{map { $_ => deep_copy($this->{$_}) } keys %$this};
	} else {
		silico_msg('d', "\"$_\" is an unknown data type!\n");
	}
}


sub deep_copy {

        #<
        #? Alternative deep copy routine for data structures
        #; Requires: Data structure
        #; Returns: copy of data structure
        #>

	use Storable qw(dclone);

	if (not ref $_[0]) {
		return $_[0];
	} else {
		return dclone($_[0]);
	}
}


sub molecule_copy_header {
	
	#<
	#? Copy a molecule's header data to another molecule.
	#; Requires: two molecules
	#; Returns: nothing
	#>
	
	my $sourcemol = $_[0];
	my $destmol = $_[1];
	
	my @l = qw(HAS_NO_COORDS SOURCE_FILE_NAME SOURCE_FILE_TYPE OUTPUT_FILE_BASE 
		CELL CELL_ALPHA CELL_BETA CELL_GAMMA);
	
	foreach (@l) {
		$destmol->{$_} = deep_copy($sourcemol->{$_}) if defined $sourcemol->{$_};
	}
}

sub molecule_copy_coordinates {
	
	#<
	#? Copy coordinates from one molecule to another (similar) molecule
	#; Requires: two molecules
	#; Returns: nothing
	#>
	
	my $sourcemol = $_[0];
	my $destmol = $_[1];
	
	my $satom;
	my $datom;
	my $i = 0;
	
	my @c = qw(X Y Z);
	
	silico_msg('d', "Molecules have different numbers of atoms") if $sourcemol->{NUMATOMS} != $destmol->{NUMATOMS};
	
	foreach my $satom (atoms($sourcemol)) {
		
		$datom = $destmol->{ATOMS}[$i];
		
		foreach (@c) {
			$datom->{$_} = $satom->{$_};
		}
		++$i;
	}
}

sub molecule_copy_noatoms {

	#<
	#? Create a copy of a molecule, but without any atoms
	#. Note: this routine fails for cyclic references.  These could be implemented
	#  by keeping a list of references in a hash.
	#; Requires: Molecule
	#; Returns: The same molecule, but containing no $mol->{ATOMS} records
	#  or anything in $mol->{NUMATOMS}
	#>

	my $mol = $_[0];
	my $newmolecule;
	
	my @skip = qw(ARINGS ATOMS BINS BINSIZE BONDS NUMATOMS NUMBINS NUMBONDS RINGS PRINGS 
			PLANAR_ATOMS RESIDUES GMX_BONDS GMX_ANGLES GMX_DIHEDRALS GMX_PAIRS GMX_EXCLUSIONS WARN);
	KEY: foreach my $key (keys %$mol) {
				
		# Skip values that will no longer be vaild
		
		foreach (@skip) {
			next KEY if $key eq $_;
		}
					
		$newmolecule->{$key} = deep_copy($mol->{$key});
	}
	
	$newmolecule->{NUMATOMS} = 0;
	$newmolecule->{NUMBONDS} = 0; # This should be true
	
	return $newmolecule;
}


##########################
# SORTING MOLECULAR DATA #
##########################


sub molecule_pack {

	#<
	#? Remove any undefined atoms from a molecule.
	#. Bonds are updated. Atom->{NUM} is renumbered
	#; Requires: molecule, nowarnings flag.
	#; Returns: packed molecule
	#>

	my $mol = $_[0];
	my $nowarn = $_[1];
	
	my $ctable;
	my $newatoms;
		
	# Generate translation table
	my $i = -1;
	my $j = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {

		++$i;
		
		# Skip undefined atoms that arise from deletions
		next if (!defined $atom);

		$ctable->[$i] = $j;
		$newatoms->[$j] = $atom;

		++$j;
	}
	
	$mol->{NUMATOMS} = $j;

	# Translate using connection table
	$i = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {

		++$i;

		next if (!defined $atom);

		foreach ($j=0; $j<=$#{$atom->{CONNECT}}; ++$j) {

			if (!defined $ctable->[$atom->{CONNECT}[$j]]) {
				silico_msg('w', "Connection to undefined atom $atom->{NUM}\n") if !$nowarn;
				next;
			}

			$atom->{CONNECT}[$j] = $ctable->[$atom->{CONNECT}[$j]];
		}
	}

	# Renumber all gromacs GMX_BONDS, GMX_ANGLES, etc.
	pack_gromacs_angles($mol, $ctable);
	
	$mol->{ATOMS} = $newatoms;
	$mol->{NUMBONDS} = molecule_count_numbonds($mol);

	# Rings are probably bad.  Mark them bad
	$mol->{BAD_RINGS} = 1;

	# Renumber atom->{NUM}
	molecule_renumber($mol);
	
	return $mol;
}


sub mol_rename_atoms {

        #<
        #? Name atoms using one of several methods:
	#; Consecutive; Atom numbering increases by atom
	#; Consecutive-H; Atom numbering increases by atom, but hydrogens are named after the parent heavy atom
	#; Element; Atom numbering is independent for each element
	#; Element-H; As for element, but hydrogens are named after the parent heavy atom
	#.
	#; Requires: molecule, method (consecutive, element, consecutive-h), flag to only number if existing names are not unique
	#; Returns: Nothing
        #>

        my $mol = $_[0];
        my $method = lc($_[1] || "element");
	my $hparent = 1 if $method =~ 'h';

	my $count = 0;
	my $el_hash;

        foreach my $atom (atoms($mol)) {

                my $el = $atom->{ELEMENT};
		++$el_hash->{$el};

                next if ($el eq 'H' && $hparent);

                ++$count;
		
                my $n = $count if $method =~ 'consecutive';
		$n = $el_hash->{$el} if $method =~ 'element';
		
		my $num;
		if ($n <=  99) {
			$num = $n;
		} else {
			$num = enc_string($n-100);
		}
		
		$atom->{NAME} = "$el$num";
               
                if ($hparent) {
		
			# Count attached hydrogens
			my $htotal = 0;
                        foreach my $con (connected($atom, $mol)) {
                                next if $con->{ELEMENT} ne 'H';
                                ++$htotal;
                        }
			
			# Rename hydrogens
                        foreach my $con (connected($atom, $mol)) {

                                next if $con->{ELEMENT} ne 'H';

				# Try to keep atom names under 4 chars
				my $x = length($num) + length($atom->{ELEMENT});
				++$x if $htotal;

				my $base;
				if ($el ne 'C' && $x <= 4) {
					$base = 'H'.$atom->{ELEMENT};
				} else {
					$base = 'H';
				}

				if ($htotal == 1) {
					$con->{NAME} = $base.$num;
				} else {
					my $hlabel = 'A';
					foreach my $con (connected($atom, $mol)) {
						next if $con->{ELEMENT} ne 'H';
						$con->{NAME} = $base.$num.$hlabel;
						 ++$hlabel;
					}
                                }
                        }
                }
        }

	my %namehash = ();
	foreach my $atom (atoms($mol)) {
		silico_msg('w', "Algorithm has unfortunately generated a duplicate atom name '$atom->{NAME}\n") if $namehash{$atom->{NAME}};
		++$namehash{$atom->{NAME}};
	}

	return 1;
}

sub enc_string {

	# Encode numbers as characters
	# Taken from perlmonks

  	my $num = shift;
  	my $end = "";
	
	my @alphanum = ("A" .. "Z");
	my $divisor = scalar @alphanum;

  	my $a = int($num / $divisor),
	my $b = $num % $divisor;
 	$end = $alphanum[$b] . $end;

 	until ($a < $divisor) {
		$a = int($a / $divisor),
		$b = $a % $divisor;
   	 	$end = $alphanum[$b] . $end;
 	}

  	$end = $alphanum[$a] . $end;

 	return $end;
}

sub dec {

	my $num = shift;

	my @alphanum = ("A" .. "Z");
	my $divisor = scalar @alphanum;

	# strip leading 0's
	$num =~ s/$0+//g;

	my ($y, $result) = (0, 0);

	foreach (split(//, reverse($num))) {

		my $found = 0;

		foreach my $item (@alphanum) {
			if ($item eq $_) {
				last;
			}
			$found++;
		}

		my $temp = $found * ($divisor ** $y);
		$result += $temp;
		$y++;
	}

	return $result;
}

sub mol_fix_atom_names {
	
	#<
	#? Removes initial and final spaces from atom names. Other spaces
	#  are replaced by underscores.
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	
	foreach my $atom (@{$mol->{ATOMS}}) {
		
		$atom->{NAME} =~ s/^\s+//g;
		$atom->{NAME} =~ s/\s+$//g;
		$atom->{NAME} =~ s/\s/_/g;
	}
}

sub molecule_renumber {

	#<
	#? Renumber atoms in a molecule (Note!. Only renumbers $atom->{NUM})
	#. Renumbers from 1 or from an optional second argument.
	#; Requires: molecule, (optional starting number).
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $start = $_[1] || 1;
	
	my $i = 0;
	
	foreach my $atom (@{$mol->{ATOMS}}) {

		next if !defined $atom;
		$atom->{NUM} = $i+$start;
		++$i;
	}
}


sub molecule_renumber2 {

	#<
	#? Renumber atoms in a molecule using an offset.  
	#. Useful to renumber one molecule before combining with another
	#. Renumbers $atom->NUM, CONNECT, GMX_BONDS, GMX_ANGLES, GMX_DIHEDRALS, GMX_PAIRS, and GMX_EXCLUSIONS
	#; Requires: molecule, Initial value for atom->{NUM} (default 1).
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $start = $_[1] || 1;

	my @glist = qw(GMX_BONDS GMX_ANGLES GMX_DIHEDRALS GMX_PAIRS GMX_EXCLUSIONS);

	my $diff = $start -1; # Starting point for array numbering
	
	my $i=0;
	foreach my $atom (@{$mol->{ATOMS}}) {

		next if !defined $atom;
		$atom->{NUM} = $i+$start;
		
		foreach my $connum (@{$atom->{CONNECT}}) {
			$connum += $diff;  # This seems odd and needs to be checked!
		}
		++$i;
	}
	
	if (defined($mol->{GMX_BONDS})) {
	
		foreach my $g (@glist) {
			
			foreach my $a (@{$mol->{$g}}) {
			
				foreach (@{$a->{ATOMS}}) {
					$_ += $diff;  # CHECK THIS!
				}
			}
		}
	}

	return undef;
}

sub molecule_renumber2_new {


	#
	#  Warning - does not work correctly
	#
	

	#<
	#? Renumber atoms in a molecule using an offset.  
	#. Useful to renumber one molecule before combining with another
	#. Renumbers $atom->NUM, CONNECT, GMX_BONDS, GMX_ANGLES, GMX_DIHEDRALS, GMX_PAIRS, and GMX_EXCLUSIONS
	#; Requires: molecule, Initial value for atom->{NUM} (default 1).
	#; Returns: nothing
	#>

	my $mol = $_[0];
	#my $start = $_[1] || 1;
	my $start = 1;

	my $i = 0;
	foreach my $atom (atoms($mol)) {

		$atom->{NEW_NUM} = $i + $start;
		++$i;
	}
	
	molecule_reorder($mol);

	return undef;
}


sub molecule_reorder {

	#<
	#? Sort atom order within a molecule using a specified order
	#. Sort molecules according to specified field or $atom->{NEW_NUM} by default;
	#; Requires: molecule.
	#>

	my $mol = $_[0];
	my $field = $_[1] || "NEW_NUM";
	
	my $ctable;
	my $warn;

	# Save the old ordering
	my $i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {

		if (!defined $atom->{$field}) {
			++$warn;
			atom_printout($atom);
		}
		
		$atom->{OLDINDEX} = $i;
		++$i;
	}
	
	silico_msg('d', "Field '$field' was not defined for $warn atoms (of total $mol->{NUMATOMS})\n") if $warn;

	# Sort the atoms
	@{$mol->{ATOMS}} = sort {$a->{$field} <=> $b->{$field}} @{$mol->{ATOMS}};

	# Make translation table
	$i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
		$ctable->[$mol->{ATOMS}[$i]{OLDINDEX}] = $i;
		++$i;
	}
	
	# Translate connection table
	$i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {

		# Translate connection table
		for (my $j=0; $j<=$#{$atom->{CONNECT}}; ++$j) {

			#. Replace the old atom number with
			# the translated one
			$atom->{CONNECT}[$j] = $ctable->[$atom->{CONNECT}[$j]];
		}
		++$i;
	}
	
	#
	# Needs testing
	#
	if (defined($mol->{GMX_BONDS})) {
	
		silico_msg('w', "Please test this code before using");
	
		my @glist = qw(GMX_BONDS GMX_ANGLES GMX_DIHEDRALS GMX_PAIRS GMX_EXCLUSIONS);
		
		foreach my $g (@glist) {
			
			foreach my $a (@{$mol->{$g}}) {
			
				foreach (@{$a->{ATOMS}}) {
					
					$_ =  $ctable->[$_];
				}
			}
		}
	}
	
	# Mark rings as bad
	$mol->{BAD_RINGS} = 1;

	# Molecule now needs to be renumbered
	molecule_renumber($mol);

	return $ctable;
}

sub mol_renumber_substructures {

        #<
        #? Renumbers substructures so all residues are numbered
        #  sequentially from 1
        #. Increments count on SUBID, SUBNAME, SEGID, or CHAIN
        #; Requires: Molecule, starting number (optional, default 1)
        #; Returns: Highest new substructure number
        #>

        my $mol = $_[0];
        my $subid = $_[1] || 1;

        $mol->{ATOMS}[0]{NEWSUBID} = $subid;

        my $prevatom;
        my $i = -1;
        foreach  (@{$mol->{ATOMS}}) {

                ++$i;
                next if $i == 0;

                my $atom = $mol->{ATOMS}[$i];
                $prevatom = $mol->{ATOMS}[$i-1];

                if (defined $atom->{SUBID} && ($atom->{SUBID}||0) ne ($prevatom->{SUBID}||0) ||
                        ($atom->{SUBNAME} ne $prevatom->{SUBNAME} ) ||
                        ($atom->{SEGID}||'') ne ($prevatom->{SEGID}||'') ||
                        defined $atom->{CHAIN} && ($atom->{CHAIN}||'') ne ($prevatom->{CHAIN}||'')) {

                        ++$subid;
                }

                $atom->{NEWSUBID} = $subid;

        }

        foreach my $atom (@{$mol->{ATOMS}}) {

                $atom->{SUBID} = $atom->{NEWSUBID};
                delete $atom->{NEWSUBID};
        }

        return $subid;
}


#################################
# QUERYING MOLECULE PROPERTIES  #
#################################

sub atoms {
	
	#<
	#? Return atoms in molecule as an array
	#. Note that this routine is actually quite expensive to call if 
	#  there are a large number of atoms in the molecule.  If it is being
	#  called many times it is much better to use @{$mol->{ATOMS}}
	#; Requires: molecule
	#; Returns: array of atoms
	
	# This often seems to trip us up.
	if (ref($_[0]) ne 'HASH') {

		silico_msg('d', "Argument provided to subroutine \"atoms\" is undefined!\n") if (!defined $_[0]);
		silico_msg('d', "Molecule \"$_[0]\" is not a hash!\n");
	}

	# The 'ATOMS' array can be undefined if there are no atoms
	if (defined $_[0]->{NUMATOMS} && $_[0]->{NUMATOMS} == 0) {
		return ();
	}

	# This is a problem though
	if (!defined $_[0]->{ATOMS} ) {
		silico_msg('d', "No defined {ATOMS} key in molecule \"$_[0]\"!\n");
	}

	return @{$_[0]->{ATOMS}};
}

sub ens {

	#<
        #? Convert a list of molecules or ensembles to an ensemble
        #; Requires: molecule or ensemble
        #; Returns: an ensemble
        #>

	my $ens;
	
	foreach my $t (@_) {

		if (ref $t eq 'HASH') {
			push @$ens, $t;
		} elsif (ref $t eq 'ARRAY') {
			foreach (@{$t}) {
				push @$ens, $_;
			}
		} else {
			silico_msg('d', "Illegal thing '".(ref $t)."' passed to subroutine 'ens'");
		}
	}

	return $ens;
}

sub mols {

	#<
        #? Return a list of molecules from a set of molecules or ensembles
        #; Requires: molecule or ensemble
        #; Returns: array of molecules
        #>

	my $ens;
	foreach (@_) {

		if (ref $_ eq 'HASH') {
			push @$ens, $_;
		} elsif (ref $_[0] eq 'ARRAY') {
			push @$ens, @$_;
		} else {
			silico_msg('d', "Illegal thing '".(ref $_)."' passed to subroutine 'mols'");
		}
	}

	return @$ens;
}

sub connected {
	
	#<
	#? Return atoms connected to an atom as an array
	#; Requires: atom, molecule
	#; array of atoms
	#>

	my $atom = $_[0];
	my $mol= $_[1];
	
	silico_msg('d', "Atom required for connected argument 1\n") if (!defined $atom->{NUM});
	silico_msg('d', "Molecule required for connected argument 2\n") if (!defined $mol->{ATOMS});
	
	my @array;
	
	foreach my $connum (@{$atom->{CONNECT}}) {
	
		next if !defined  $mol->{ATOMS}[$connum];
		push @array, $mol->{ATOMS}[$connum];
	}
	
	return @array;
}


sub mol_calc_charge {
	
	#<
	#? Calculates the total of the partial atomic charges for a molecule and checks that it is integral
	#. Also can calculate the charge on each residue and on each segment
	#. Will return a warning if there are missing charges or if the total charge is not integral
	#; Requires: Molecule, flag to calculate residue totals, flag to calculate segment totals
	#; Returns:	ERR if there are atoms with no set charge
	#		Total charge otherwise
	#>
	
	my $mol = $_[0];
	my $r = $_[1];
	my $s = $_[2];
	
	my $charge;
	my @missing;
	my $res;
	my %res_charge;
	my $rkey;
	my $seg;
	my %seg_charge;
	my $skey;
	my $total_charge = 0;
	
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		if (!defined $atom->{CHARGE}) {
			push @missing, "Number $atom->{NUM}, Name '$atom->{NAME}', Residue '$atom->{SUBNAME} $atom->{SUBID}'";
			next;
		}
		
		$charge = $atom->{CHARGE};
		$total_charge += $charge;
		
		if ($r) {
		
			$atom->{SUBID} ||= '1';
			$atom->{SUB} ||= 'SUB';
			$res = '';
			
			if ($atom->{CHAIN}) {
				$res .= $atom->{CHAIN}."\t";
			}
			
			if ($atom->{SEGID}) {
				$res .= $atom->{SEGID}."\t";
			}
			
			$res .= $atom->{SUBNAME}." ".$atom->{SUBID};
			$res_charge{$res} += $charge;
		}
		
		if ($s) {
			
			$seg = $atom->{SEGID} || 'no_segid';
			$seg_charge{$seg} += $charge;
		}
	}
	
	# Printout charges sorted by residue
	if ($r) {
	
		silico_msg('c', heading("Residue charges\n"));
		
		foreach my $rkey (sort mol_calc_charge_sort_res keys(%res_charge)) {
			
			my $string;
			
			$res_charge{$rkey} = sprintf "%.5f", $res_charge{$rkey};
			$res_charge{$rkey} = 0.00000 if $res_charge{$rkey} eq "-0.00000";
			
			$string = $rkey;
			if ($res_charge{$rkey} < 0) {
				$string .= " " x (15 - length $rkey);
			} else {
				$string .= " " x (16 - length $rkey);
			}
			$string .= "$res_charge{$rkey}";
			
			silico_msg('c', "$string\n");
		}
		
		silico_msg('c', "\n");
	}
	
	# Printout charges sorted by SEGID
	if ($s) {
	
		silico_msg('c', heading("Segment charges\n"));
		
		foreach my $skey (sort keys(%seg_charge)) {
			
			my $string;
			
			$seg_charge{$skey} = sprintf "%.5f", $seg_charge{$skey};
			$seg_charge{$skey} = 0.00000 if $seg_charge{$skey} eq "-0.00000";
			
			$string = $skey;
			if ($seg_charge{$skey} < 0) {
				$string .= " " x (15 - length $skey);
			} else {
				$string .= " " x (16 - length $skey);
			}
			$string .= "$seg_charge{$skey}";
			
			silico_msg('c', "$string\n");
		}
		
		silico_msg('c', "\n");
	}
	
	if (defined $missing[0]) {
	
		my $c = $#missing+1;
		
		print "\n";
		silico_msg('w', "$c atoms have no defined charge!\n",
				"The first 10 of these are the following:\n");
			
		for (my $i = 0; $i < 10; ++$i) {
			silico_msg('q', "$missing[$i]\n");
		}
		
		silico_msg('q', "\n\n");
		return 'ERR';
	}
	
	$total_charge = sprintf ("%.5f", $total_charge);
	$total_charge = 0.00000 if $total_charge eq "-0.00000";
	
	return $total_charge;
}

sub mol_calc_charge_sort_res {
	
	# Get the residue numbers: the last set of digits in each string
	my $a ||= '';
	$a =~ m/(\d+$)/;
	my $a2 = $1 || 0;
	
	my $b ||= '';
	$b =~ m/(\d+$)/;
	
	my $b2 = $1 || 0;
	
	return $a2 <=> $b2;
}

sub mol_check_unit_cell {
	
	#<
	#? Make sure a molecule's unit cell is present
	#; Requires: Molecule
	#; Returns: 1 if OK, 0 otherwise
	#>
	
	my $mol = $_[0];
	
	if (defined $mol->{CELL}[0]) {
		if (check_data_type($mol->{CELL}[0], "DECIMAL >0")) {
			silico_msg('g', "a axis: $mol->{CELL}[0]\n");
		} else {
			silico_msg('w', "Bad unit cell a axis!\n",
					"a axis: \"$mol->{CELL}[0]\"\n");
			return 0;
		}
	} else {
		silico_msg('w', "Unit cell a axis is not defined!\n");
		return 0;
	}

	if (defined $mol->{CELL}[1]) {
		if (check_data_type($mol->{CELL}[1], "DECIMAL >0")) {
			silico_msg('g', "b axis: $mol->{CELL}[1]\n");
		} else {
			silico_msg('w', "Bad unit cell b axis!\n",
					"b axis: \"$mol->{CELL}[1]\"\n");
			return 0;
		}
	} else {
		silico_msg('w', "Unit cell b axis is not defined!\n");
		return 0;
	}

	if (defined $mol->{CELL}[2]) {
		if (check_data_type($mol->{CELL}[2], "DECIMAL >0")) {
			silico_msg('g', "c axis: $mol->{CELL}[2]\n");
		} else {
			silico_msg('w', "Bad unit cell c axis!",
					"c axis: \"$mol->{CELL}[2]\"\n");
			return 0;
		}
	} else {
		silico_msg('w', "Unit cell c axis is not defined!\n");
		return 0;
	}
	
	if (defined $mol->{CELL_ALPHA}) {
		if (check_data_type($mol->{CELL_ALPHA}, "DECIMAL >0 <180")) {
			silico_msg('g', "alpha angle: $mol->{CELL}[0]\n");
		} else {
			silico_msg('w', "Bad unit cell alpha angle!\n",
					"alpha angle: \"$mol->{CELL_ALPHA}\"\n");
			return 0;
		}
	} else {
		silico_msg('w', "Unit cell alpha angle is not defined!\n");
		return 0;
	}
	
	if (defined $mol->{CELL_BETA}) {
		if (check_data_type($mol->{CELL_BETA}, "DECIMAL >0 <180")) {
			silico_msg('g', "beta angle: $mol->{CELL_BETA}\n");
		} else {
			silico_msg('w', "Bad unit cell beta angle!\n",
					"beta angle: \"$mol->{CELL_BETA}\"\n");
			return 0;
		}
	} else {
		silico_msg('w', "Unit cell beta angle is not defined!\n");
		return 0;
	}
	
	if (defined $mol->{CELL_GAMMA}) {
		if (check_data_type($mol->{CELL_GAMMA}, "DECIMAL >0 <180")) {
			silico_msg('g', "gamma angle: $mol->{CELL_GAMMA}\n");
		} else {
			silico_msg('w', "Bad unit cell gamma angle!\n",
					"gamma angle: \"$mol->{CELL_GAMMA}\"\n");
			return 0;
		}
	} else {
		silico_msg('w', "Unit cell gamma angle is not defined!\n");
		return 0;
	}
	
	return 1;
}

sub mol_copy_unit_cell {
	
	#<
	#? Copies as much as is defined of a unit cell record from one
	#  Silico molecule to another.
	#; Requires: two molecules
	#; Returns: nothing
	#>
	
	my $mol1 = $_[0];
	my $mol2 = $_[1];
	
	if (defined $mol1->{CELL}[0]) {
		$mol2->{CELL} = deep_copy($mol1->{CELL});
	}
	
	if (defined $mol1->{CELL_ALPHA}) {
		$mol2->{CELL_ALPHA} = deep_copy($mol1->{CELL_ALPHA});
	}
	if (defined $mol1->{CELL_BETA}) {
		$mol2->{CELL_BETA} = deep_copy($mol1->{CELL_BETA});
	}
	if (defined $mol1->{CELL_GAMMA}) {
		$mol2->{CELL_GAMMA} = deep_copy($mol1->{CELL_GAMMA});
	}
		
	if (defined $mol1->{SPACE_GROUP}) {
		$mol2->{SPACE_GROUP} = deep_copy($mol1->{SPACE_GROUP});
	}
}

sub molecule_get_chains {

	#<
	#? Return an array of chains
	#; Requires: molecule
	#; Returns: array of atom lists
	#>

	my $mol = $_[0];
	
	my $oldchain = '';
	my $oldsegid = '';
	my $chaincount = -1;
	my $chainlist;
	
	foreach my $atom (@{$mol->{ATOMS}}) {
			
		if ((defined $atom->{CHAIN} && $atom->{CHAIN} ne $oldchain)
			|| (defined $atom->{SEGID} && $atom->{SEGID} ne $oldsegid)) {
			
			$oldchain = $atom->{CHAIN} || '';
			$oldsegid = $atom->{SEGID} || '';
			++$chaincount;
		}

		push @{$chainlist->[$chaincount]}, $atom;
	}
	
	return $chainlist;
}




#########################
# QUERYING CONNECTIVITY #
#########################

sub mol_bonded_interactions {

	#<
	#? Creates tables containing bonds, angles, dihedrals and improper dihedrals
	#  in a molecule
	#; Calls:   Mol_find_all_bonds
	#           Mol_find_all_angles
	#           Mol_find_all_dihedrals
	#           Mol_find_planar_impropers
	#           Mol_find_double_bond_impropers
	#; Requires: Molecule, flag to ignore waters, flag to skip hydrogen dihedrals, flag to find only unique dihedrals
	#; Returns: Nothing
	#>

	my $mol = $_[0];
	my $skip_water = $_[1];
	my $skip_h = $_[2];
	my $unique = $_[3];
	
	mol_find_all_bonds($mol, $skip_water);
	mol_find_all_angles($mol, $skip_water);
	mol_find_all_dihedrals($mol, $skip_water, $skip_h, $unique);
	
	if ($mol->{HAS_NO_COORDS}) {
		silico_msg('w', "Molecule has no coordinates!\n",
				"Re-identification of improper dihedral angles has been skipped.\n",
				"Existing improper dihedral information may be preserved.\n");
	} else {
		mol_find_planar_impropers($mol, $skip_water, $unique);
	}
}

sub molecule_count_numbonds {

	#<
	#? Return the number of bonds.
	#. Note: Each bond is specified twice in CONNECT records,
	#  once for each of the two atoms in the bond. So the total
	#  number of CONNECT records is divided by two, giving the
	#  number of bonds.
	#; Requires: Molecule
	#; Returns: Number of bonds
	#>
	
	my $mol = $_[0];

	my $count = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		next if !defined $atom;
			
		if (defined $atom->{CONNECT} || defined $atom->{CONNECT}[0]) {
			$count += (1 + $#{$atom->{CONNECT}});
		}
	}
	
	return $count/2;
}

sub create_angle {
	
	#<
	#? Create an ANGLE record between three atoms (whether they are actually
	#  connected or not). This is used, for example, in reading in PSF files.
	#; Requires: molecule, three atom offsets
	#; Returns: 1 if successful
	#>
	
	my $mol = $_[0];
	my $a1 = $_[1];
	my $a2 = $_[2];
	my $a3 = $_[3];
	
	if (!defined $a1 || !defined $a2 || !defined $a3) {
		
		silico_msg('e', "One or more atoms have not been supplied!\n",
				"Three atom offsets are required.\n",
				"No angle has been created.\n");
		return undef;
	}
	
	my $angle;
	
	@$angle = ($a1, $a2, $a3);
	
	push @{$mol->{ANGLES}}, $angle;
	return 1;
}

sub mol_find_all_angles {

	#<
	#? Find all atom-atom-atom angles in a molecule
	#. Angle array is stored in $mol->{ANGLES}
	#. Label_waters must be run first if water is to be ignored
	#; Requires: Molecule, (optional) flag to ignore waters
	#; Returns: Arrayref of angles, each of which is an arrayref of 3 atom numbers
	#>

	my $mol = $_[0];
	my $flag = $_[1];
	
	# Clear mol->ANGLES
	@{$mol->{ANGLES}}  = ();
	
	my $i = -1;
	my $j = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		++$i;
	
		next if $flag && $atom->{FG}{W};
	
		foreach my $connum1 (@{$atom->{CONNECT}}) {
		
			foreach my $connum2 (@{$atom->{CONNECT}}) {
			
				next if $connum2 <= $connum1;
				
				my $angle;
				
				++$j;
				
				$angle->[0] = $connum1;
				$angle->[1] = $i;
				$angle->[2] = $connum2;
				
				$mol->{ANGLES}[$j] = $angle;
			}
		}
	}
	
	# Sort the array of angles
	@{$mol->{ANGLES}} = sort {
		return 1 if $a->[1] > $b->[1];
		return -1 if $a->[1] < $b->[1];
		return 1 if $a->[0] > $b->[0];
		return -1 if $a->[0] < $b->[0];
		return 1 if $a->[2] > $b->[2];
		return -1 if $a->[2] < $b->[2];
		return 0;
					
		} @{$mol->{ANGLES}};
	
	return $mol->{ANGLES};;
}

sub mol_find_all_bonds {

	#<
	#? Find all bonds in a molecule
	#. Bond array is stored in $mol->{BONDS}
	#. Label_waters must be run first if water is to be ignored
	#; Requires: Molecule, (optional) flag to ignore waters
	#; Returns: Arrayref of bonds, each of which is an arrayref of 2 atom numbers
	#>
	
	my $mol = $_[0];
	my $flag = $_[1];
	
	my $i = -1;
	my $j = -1;

	# Clear mol->BONDS
	@{$mol->{BONDS}}  = ();
	
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		++$i;
	
		next if $flag && $atom->{FG}{W};
	
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $bond;
			
			next if $connum <= $i;
			++$j;
			
			$bond->[0] = $i;
			$bond->[1] = $connum;
			$mol->{BONDS}[$j] = $bond;
		}
	}
	
	@{$mol->{BONDS}} = sort {
		return 1 if $a->[0] > $b->[0];
		return -1 if $a->[0] < $b->[0];
		return 1 if $a->[1] > $b->[1];
		return -1 if $a->[1] < $b->[1];
		return 0;
	} @{$mol->{BONDS}};

	return $mol->{BONDS};
}

sub mol_find_all_dihedrals {

	#<
	#? Find all atom-atom-atom-atom dihedral (torsion) angles in a molecule
	#. Dihedral array is stored in $mol->{DIHEDRALS}
	#. Label_waters must be run first if water is to be ignored
	#; Requires: Molecule, (optional) flag to ignore waters, (optional) flag to ignore dihedrals
        #  containing hydrogen, flag to only find unique dihedrals (ie only one dihedral per bond)
	#; Returns: Arrayref of sorted dihedrals, each of which is an arrayref of 4 atom numbers
	#>
	
	my $mol = $_[0];
	my $flag = $_[1];
	my $skiph = $_[2];
	my $unique = $_[3];
	
	my $atoms = $mol->{ATOMS};
	my $hash;
	my $j = -1;

	# Clear mol->DIHEDALS
	@{$mol->{DIHEDRALS}}  = ();

	# AtomB
	my $atomnumb = -1;
	foreach my $atomb (@{$mol->{ATOMS}}) {

		++$atomnumb;
		
		# Skip water molecules if flag is set
		next if $flag && $atomb->{FG}{W};

		# Skip if atomb is terminal
		next if $#{$atomb->{CONNECT}} <= 0;

		#AtomC
		foreach my $atomnumc (@{$atomb->{CONNECT}}) {
		
			next if $atomnumc < $atomnumb; # Select each dihedral only once
			my $atomc = $atoms->[$atomnumc];

			# (AtomA) All atoms connected to AtomB
			foreach my $atomnuma (@{$atomb->{CONNECT}}) {
			
				next if $atomnumc == $atomnuma;
				
				my $atoma = $atoms->[$atomnuma];
				next if $skiph && $atoma->{ELEMENT} eq 'H';
			
				# AtomD (All atoms connected to atomC)
				foreach my $atomnumd (@{$atomc->{CONNECT}}) {
				
					next if $atomnumd == $atomnumb;
					next if $skiph && $atoms->[$atomnumd]{ELEMENT} eq 'H';

					if ($unique) {
						my $h = $atomnumb."_".$atomnumc;	
						last if $hash->{$h};
						++$hash->{$h};
					}
										
					my $tors;
					
					++$j;
					
					@$tors = ($atomnuma, $atomnumb, $atomnumc, $atomnumd);
					$mol->{DIHEDRALS}[$j] = $tors;
				}
			}
		}
	}
	
	@{$mol->{DIHEDRALS}} = sort {
					return 1 if $a->[1] > $b->[1];
					return -1 if $a->[1] < $b->[1];
					return 1 if $a->[2] > $b->[2];
					return -1 if $a->[2] < $b->[2];
					return 1 if $a->[0] > $b->[0];
					return -1 if $a->[0] < $b->[0];
					return 1 if $a->[3] > $b->[3];
					return -1 if $a->[3] < $b->[3];
					return 0;
				} @{$mol->{DIHEDRALS}};
	
	my $numlist = deep_copy($mol->{DIHEDRALS});
	return $numlist;
}

sub mol_find_chiral_impropers {

	#<
	#? Find improper torsion angles for chiral centres (some force fields require this)
	#. Impropers are stored in $mol->{C_IMPROPERS}
	#. Returns improper torsions where an atom has more than two attached atoms
	#  Only one improper is produced in cases where an atom has four attached atoms
	#  using the first three connected atoms.
	#; Requires: Molecule
	#; Returns: Arrayref of improper dihedrals, each of which is an arrayref of 4 atom numbers
	#>
	
	my $mol = $_[0];
	
	my $i = -1;
	my $numlist;
	my $tornum = -1;
	
	# Chiral impropers
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		++$i;
		next if ($#{$atom->{CONNECT}} < 2); # Skip if we have 0, 1 or 2 attatched atoms
		
		my $tor;
		++$tornum;
		
		@$tor = ($i, $atom->{CONNECT}[0], $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
		
		#my $angle = bndangle($mol->{ATOMS}[$tor->[0]], $mol->{ATOMS}[$tor->[1]],$mol->{ATOMS}[$tor->[2]], $mol->{ATOMS}[$tor->[3]]);
				
		$mol->{C_IMPROPERS}[$tornum] = $tor;
		push @$numlist, $tor;
	}
	
	return $numlist;
}



sub create_planar_improper {
	
	#<
	#? Create a planar improper dihedral, e.g. when reading in a PSF file
	#; Requires: molecule, list of four atom offsets (NB. not atom numbers)
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $ao1 = $_[1];
	my $ao2 = $_[2];
	my $ao3 = $_[3];
	my $ao4 = $_[4];
	
	my $imptor;
	@$imptor = ($ao1, $ao2, $ao3, $ao4);
	
	push @{$mol->{P_IMPROPERS}}, $imptor;
	
	return undef;
}

sub mol_find_planar_impropers {
	
	#<
	#? Finds improper dihedrals centred on planar atoms
	#. Impropers are stored in $mol->{P_IMPROPERS}
	#.
	#.      J
	#.      |
	#.      I
	#.     / \
	#.    K   L
	#.
	#. Orders improper dihedrals suitable for the OPLS force field
	#. Atom I is always position 2, while the unique atom is usually
	#  position 1
	#
	#. Note this subroutine needs a good clean up to remove duplication
	#
	#; Requires: Molecule, flag to ignore waters
	#; Returns: Arrayref of improper dihedrals, each of which is an arrayref of 4 atom numbers
	#>

	#
	# Not likely to work because AMIDE_N etc needs to be converted to S_AMIDE_N etc
	#
	
	my $mol = $_[0];
	my $wflag = $_[1];
	
	my $torlist;
	my $tornum = -1;
	
	molecule_find_rings ($mol);
	
	@{$mol->{P_IMPROPERS}} = ();
	
	my $i = -1;
	ATOM: foreach my $atom (@{$mol->{ATOMS}}) {
		
		++$i;
		
		# Ignore waters
		# This goes after the "++$i" so the count is correct
		next ATOM if $wflag && $atom->{FG}{W};
		
		# Ignore atoms that aren't connected to exactly 3 others
		next ATOM if ($#{$atom->{CONNECT}} != 2);
		
		# Get the atoms this atom is connected to
		my $con0 = $mol->{ATOMS}[$atom->{CONNECT}[0]];
		my $con1 = $mol->{ATOMS}[$atom->{CONNECT}[1]];
		my $con2 = $mol->{ATOMS}[$atom->{CONNECT}[2]];
		
		# Create improper entries for all amide carbons
		# O, C, R, N and N, C, R, O
		# O, C, H, N and N, C, H, O
				
		my $tor1;
		my $tor2;
		if ($atom->{FG}{AMIDE_C}) {
			
			if (($con0->{FG}{AMIDE_O} && $con1->{FG}{AMIDE_N}) || ($con1->{FG}{AMIDE_O} && $con0->{FG}{AMIDE_N})) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif (($con0->{FG}{AMIDE_O} && $con2->{FG}{AMIDE_N}) || ($con2->{FG}{AMIDE_O} && $con0->{FG}{AMIDE_N})) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif (($con1->{FG}{AMIDE_O} && $con2->{FG}{AMIDE_N}) || ($con2->{FG}{AMIDE_O} && $con1->{FG}{AMIDE_N})) {
				
				@$tor1 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
		}
		
		# Create improper entries for all amide nitrogens
		# H, N, C(O), R
		# H, N, C(O), H
		# R1, N, C(O), R2
		
		if ($atom->{FG}{AMIDE_N}) {
			
			if ($con0->{FG}{AMIDE_C}) {
						
				@$tor1 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			
			} elsif ($con1->{FG}{AMIDE_C}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif ($con2->{FG}{AMIDE_C}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
		}
		
		# Create improper entries for all carboxylate carbons
		# O, C, R, O
		
		if ($atom->{FG}{CARBOXYLATE_C}) {
			
			if ($con0->{FG}{CARBOXYLATE_O} && $con1->{FG}{CARBOXYLATE_O}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
			
			if ($con0->{FG}{CARBOXYLATE_O} && $con2->{FG}{CARBOXYLATE_O}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
			
			if ($con1->{FG}{CARBOXYLATE_O} && $con2->{FG}{CARBOXYLATE_O}) {
				
				@$tor1 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
		}
		
		# Create improper entries for all other carbonyls
		# X, C, O, Y
		
		if ($atom->{FG}{CARBONYL_C}) {
			
			if ($con0->{FG}{CARBONYL_O}) {
				
				@$tor1 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif ($con1->{FG}{CARBONYL_O}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif ($con2->{FG}{CARBONYL_O}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
		}
		
		# Create improper entries for all atoms in a single aromatic ring
		# X, Ar, Y, Z and Y, Ar, X, Z where Z is outside the ring
		
		# Ensure that only atoms in precisely one aromatic ring are included
		if ($#{$atom->{PRINGS}} == 0) {
			
			# Clear the temporary holding places for in-ring and out-of-ring atoms
			
			my @inring;
			my $outring;
			
			# Compare each atom connected to this atom
			# with the atoms in the ring to which this atom belongs
			CON: foreach my $connum (@{$atom->{CONNECT}}) {
				
				my $flag = 0;
				
				foreach my $ratom (@{$atom->{PRINGS}[0]}) {
					
					# If the atom is in this ring, add it to the "inring" array
					if ($connum == $ratom) {
						push @inring, $connum;
						$flag = 1;
					}
				}
				
				$outring = $connum if ($flag == 0);
			}
			
			# Get the full atom table for each of the in-ring atoms
			my $inring0 = $mol->{ATOMS}[$inring[0]];
			my $inring1 = $mol->{ATOMS}[$inring[1]];
			my $outringa = $mol->{ATOMS}[$outring];
			
			# Sort by element
			@$tor1 = ($inring[0], $i, $inring[1], $outring);
			++$tornum;
			$mol->{P_IMPROPERS}[$tornum] = $tor1;
			push @$torlist, $tor1;
			
			@$tor2 = ($inring[1], $i, $inring[0], $outring);
			++$tornum;
			$mol->{P_IMPROPERS}[$tornum] = $tor2;
			push @$torlist, $tor2;
				
			next ATOM;
		}
		
		# Create improper entries for all atoms in more than one aromatic ring
		# Reverse sorted by element
		#
		# COMMENTED OUT BPR 19/12/2006
		# (not done in Tinker, so is questionable)
		# Instead, these atoms are skipped, according to the following instruction.
		
		next ATOM if (defined $atom->{PRINGS}[1]);
			
		# Skip guanidinium C (per OPLS force field)
		# BPR 18/7/09
		next if ($atom->{FG}{GUANIDINIUM_C} && !defined $atom->{PRINGS}[0]);
		next if ($atom->{FG}{GUANIDINIUM_N} && !defined $atom->{PRINGS}[0]);
		
		# Create improper entries for all alkene carbons
		if ($atom->{FG}{ALKENE_CH0} || $atom->{FG}{ALKENE_CH1} || $atom->{FG}{ALKENE_CH2}) {
			
			if ($con0->{FG}{ALKENE_CH0} || $con0->{FG}{ALKENE_CH1} || $con0->{FG}{ALKENE_CH2}) {
				
				@$tor1 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[0], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif ($con1->{FG}{ALKENE_CH0} || $con1->{FG}{ALKENE_CH1} || $con1->{FG}{ALKENE_CH2}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[2], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
				
			} elsif ($con2->{FG}{ALKENE_CH0} || $con2->{FG}{ALKENE_CH1} || $con2->{FG}{ALKENE_CH2}) {
				
				@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[1]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor1;
				push @$torlist, $tor1;
				
				@$tor2 = ($atom->{CONNECT}[1], $i, $atom->{CONNECT}[2], $atom->{CONNECT}[0]);
				++$tornum;
				$mol->{P_IMPROPERS}[$tornum] = $tor2;
				push @$torlist, $tor2;
				
				next ATOM;
			}
		}
		
		# Check for other double or aromatic systems
		my $flag = 0;
		my $j = -1;
		CON: foreach my $connum (@{$atom->{CONNECT}}) {
			
			++$j;
			
			my $bo = $atom->{BORDERS}[$j];
			
			if ($bo == 2 || $bo == 4) {
				$flag = 1;
				last;
			}
		}
		
		
		# We do not have a double or aromatic bond
		next ATOM if !$flag;
		
		# Now, add impropers for all planar or aromatic systems which are left
		@$tor1 = ($atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]);
		++$tornum;
		$mol->{P_IMPROPERS}[$tornum] = $tor1;
		push @$torlist, $tor1;
		
		silico_msg('n', "Improper torsion of unknown type created for atoms $atom->{CONNECT}[0], $i, $atom->{CONNECT}[1], $atom->{CONNECT}[2]\n");
	}
	
	return $torlist;
}

sub mol_transfer_coordinates {

	#<
	#? Transfer coordinates from molecule 2 to to molecule 1
	#; Requires: molecule 1, molecule 2
	#; Returns: new molecule with properties of molecule 1 and coordinates of molecule2
	#>
	
	my $mol1 = $_[0];
	my $mol2 = $_[1];
	
	my $newmol;
	
	# Make a complete copy of psfmol
	$newmol = deep_copy($mol1);
	
	# Transfer coordinates to molecule1
	my $i = 0;
	foreach my $atom (atoms($newmol)) {
		$atom->{X} = $mol2->{ATOMS}[$i]{X};
		$atom->{Y} = $mol2->{ATOMS}[$i]{Y};
		$atom->{Z} = $mol2->{ATOMS}[$i]{Z};
		++$i;
	}
	
	delete $newmol->{HAS_NO_COORDS} if !$mol2->{HAS_NO_COORDS};
	return $newmol;
}


################################################
################################################
###                                          ###
###     SUBROUTINES DEALING WITH BONDING     ###
###                                          ###
################################################
################################################


sub molecule_check_and_fix_connectivity {

	#<
	#? Check that a molecule has bonds and that correct bond orders.  Generate them
	#  if they are missing
	#. A reasonable number of bonds is deemed to be (numatoms-1)/2;
	#. Valid options BOND, NOBOND
	#; Requires: molecule, $options (optional)
	#; Returns: 1 or undef on error
	#>

	my $mol = $_[0];
	my $options = $_[1] || '';

	my $bondcount;
	my $quiet = 1 if $options =~/\bQUIET\b/i;
	return 1 if $options =~ /\bNOBOND\b/;
	return 1 if get_lflag('nobond');
	
	# Flag to mark that this has been run so that we can avoid rerunning this subroutine
	$mol->{CONNECTIVITY_CHECKED} = 1; 
	
	if ($mol->{HAS_NO_COORDS} ) {
	
		if (!$mol->{NUMBONDS} && $mol->{NUMATOMS} > 1) {
			silico_msg('w', "Molecule '".($mol->{NAME} || '')."' does not have coordinates or bonds. Can not generate bonds.\n");
			#carp;
			return undef;
		}

		# We have bonds already - we have no choice but to use these
		return 1;
	} 
	
	if ($mol->{HAS_BAD_GEOM}) {
		#silico_msg('w', "Molecule has poor geometry. Skipping bond creation\n");
		return 0;
	}
	
	if (!defined $mol->{NUMBONDS}) {
	
		# No bonds
		silico_msg('c', "Number of bonds is not defined. Creating bonds.\n");
		$mol->{NUMBONDS} = 0;
		$bondcount = molecule_generate_bonds($mol);
		silico_msg('c', "Generated $bondcount bonds\n") if $bondcount && !$quiet;
		return 1;
	} 
	
	if ($mol->{HAS_BAD_BONDS}) {
	
		# Crummy bonds
		silico_msg('n',  "Molecule \"$mol->{NAME}\" is marked as having bad bonds. Creating bonds.\n") if !$quiet;
		$bondcount = molecule_generate_bonds($mol);
		silico_msg('c', "Generated $bondcount bonds\n") if $bondcount && !$quiet;
		return 1;
	} 
	
	if ($options =~ /\bBOND\b/) {

		# Option to force generation of bonds
		silico_msg('c', "Force is set. Creating bonds\n");
		$bondcount = molecule_generate_bonds($mol);
		silico_msg('c', "Generated $bondcount bonds\n") if $bondcount && !$quiet;
		return 1;
	} 
	
	# No bonds to carbon or no bonds to oxygen (e.g.water)
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		next if ($atom->{ELEMENT} ne 'C' && $atom->{ELEMENT} ne 'O');
		next if $#{$atom->{CONNECT}} >= 0;

		if (!$quiet) {
			silico_msg('c', "Found a carbon with 0 bonds\n") if $atom->{ELEMENT} eq 'C';
			silico_msg('c', "Found an oxygen with 0 bonds\n") if $atom->{ELEMENT} eq 'O';
		}
		$bondcount = molecule_generate_bonds($mol);
		silico_msg('c', "Generated $bondcount bonds\n") if $bondcount && !$quiet;
		return 1;
	}
	
	if ($mol->{NUMBONDS} < ($mol->{NUMATOMS}-1)/2) {

		# Too few bonds
		silico_msg('c', "Too few bonds found: $mol->{NUMBONDS} bonds for $mol->{NUMATOMS} atoms.\n") if !$quiet;
		$bondcount = molecule_generate_bonds($mol);
		silico_msg('c', "Generated $bondcount bonds\n") if $bondcount && !$quiet;
		return 1;
	} 
	
	if ($mol->{HAS_BAD_BONDORDERS}) {
	
		# Needs bondorders
		silico_msg('n', "Molecule \"$mol->{NAME}\" is marked as having bad bondorders.\n") if !$quiet;
		molecule_generate_bondorders($mol);
	}
	
	return 1;
}

sub molecule_generate_bonds {

	#<
	#? Create new bonds for structure with no connection table
	#  for example a pdb file
	#. Uses a fast residue-based approach to generate bonds in waters
	#  if the molecule contains more than 500 atoms
	#; Requires: molecule
	#; Returns: bondcount
	#>

	my $mol = $_[0];
	my $starttime = (times)[0];
	
	if ($mol->{HAS_BAD_GEOM}) {
		silico_msg('w', "Molecule has poor geometry. Skipping bond creation\n");
		return 0
	}
	my $bondcount = 0;
	
	silico_msg('c', subname()."\n") if $Silico::debug || $Silico::TIMING;

	if ($mol->{NUMATOMS} > 500) {
		connect_atoms_in_waters($mol);
		$bondcount = connect_atoms_by_distance($mol, 1);
	} else {
		$bondcount = connect_atoms_by_distance($mol);
	}
	
	molecule_generate_bondorders($mol);
	
	# Clear the HAS_BAD_BONDS flag
	delete $mol->{HAS_BAD_BONDS};
	
	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");

	return $bondcount;
}

sub molecule_generate_bondorders {

	#<
	#? Calculate bondorders for bonds
	#; Requires: molecule, flag to reset all bondorders to single
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $makesingle = $_[1];

	my $bondorder;
	my $starttime = (times)[0];
	
	# Convert all existing bonds to single if flag set
	if ($makesingle) {
		foreach my $atom (@{$mol->{ATOMS}}) {
			foreach my $bondorder (@{$atom->{BORDERS}}) {
				$bondorder = 1 if $bondorder != 0;
			}
		}
	}
	
	make_aromatic_bonds($mol); # Make aromatic bond types in aromatic rings
	fix_dative_bonds($mol, 2);	# Fix O=S and O=P bonds
	fix_multiple_bonds($mol);
	
	# Reset bondorders flag
	delete $mol->{HAS_BAD_BONDORDERS};

	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
}

sub fix_dative_bonds {

	# <
	#? Generate the correct bond orders for dative bonds
	#. Different file formats require different bonds orders for O-S and O-P bonds.
	#. For example Merck format requires that these be single bonds.  Most others
	#  seem to require double bonds.
	#. Fixes the valence of dative bonds to S and P (to specified order 1 or 2)
	#; Requires: molecule, order of dative bonds
	#; Returns: count of modified bonds
	#>
	
	my $mol = $_[0];
	my $order = $_[1];
	
	my $conatom;
	my $connum;
	my $count = 0;
	
	my $maxvalence = maximum_valences();

	foreach my $atom (@{$mol->{ATOMS}}) {
	
		next if !($atom->{ELEMENT} eq 'S' ||$atom->{ELEMENT} eq 'P');
			
		foreach my $connum (@{$atom->{CONNECT}}) {
			
			$conatom = $mol->{ATOMS}[$connum];
			
			# Do not keep increasing order past maxvalence
			last if $order >= 2 && (valence($atom) >= $maxvalence->{$atom->{ELEMENT}});
			
			if ($conatom->{ELEMENT} eq 'O' && $#{$conatom->{CONNECT}} == 0) {
			
				bond_modify_order($mol, $atom, $conatom, $order);
				++$count;
			}
		}
	}

	return $count;
}

sub fix_multiple_bonds {
	
	my $mol = $_[0];
	
	# Mol_renumber is required by bond_modify_order which is called by 
	# subroutines below
	molecule_renumber($mol);
	mol_score_hybridisation($mol);
	# Now that the probable hybridisation states of each atom have been 
	# determined, the next stage can be done. This involves a determination
	# of the likelihood of any one bond being double, aromatic, triple...
	# taking into account the respective hybridisation states of the two
	# atoms involved in the bond (to the extent these are known), and the bond length.
	
	fix_double_bonds_carbonyl($mol);
	fix_double_bonds_carboxyl($mol);
	mol_fix_multiple_bonds_general($mol);
}

sub fix_double_bonds_carbonyl {

	#<
	#? Assign carbonyl bond orders
	#. Uses molecule geometry to generate bond orders
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		next if $atom->{ELEMENT} ne 'C';
		next if !$atom->{PLANAR_ATOM};
		
		my $con_o1 = 0;
		
		foreach my $connum (@{$atom->{CONNECT}}) {
			my $con = $mol->{ATOMS}[$connum];
			
			++$con_o1 if ($con->{ELEMENT} eq 'O' && valence($con) == 1);
		}
		
		if ($con_o1 ==1) {
			change_bondorder($atom, 'O', $mol, 2);
		}
	}
	
	return undef;
}

sub fix_double_bonds_carboxyl {

	#<
	#? Assign ester, carboxylic acid, carboxylate, amide and guanidinium bond orders
	#. Uses molecule geometry to assign bond orders
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	
	if (!defined $mol->{PLANAR_ATOMS}) {
		molecule_find_planar_atoms($mol);
	}
	
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		# Amide, caboxylate and guanidinium carbons are all planar
		next if ($atom->{ELEMENT} ne 'C' || !$atom->{PLANAR_ATOM});
		
		my $con_o = 0;
		my $con_n = 0;
		my $con_o1 = 0;
		
		# Amide carbon if doubly bonded to O and also to N
		# Carboxylate carbon if bonded to two oxygen atoms
		#  which are have only one bond each
		
		foreach my $connum (@{$atom->{CONNECT}}) {
			
			my $con = $mol->{ATOMS}[$connum];
			
			++$con_o if ($con->{ELEMENT} eq 'O');
			++$con_o1 if ($con->{ELEMENT} eq 'O' && valence($con) == 1);
			++$con_n if ($con->{ELEMENT} eq 'N');
		}
				
		# We can't increase the bond order of a carbon if it already
		# greater than 3
		next if (valence($atom) > 3);
				
		# Carboxylate
		if ($con_o1 == 2) {
			change_bondorder($atom, 'O', $mol, 10);
		}
		# Ester/Carboxylic acid
		if ($con_o == 2 && $con_o1 == 1) {
			change_bondorder($atom, 'O', $mol, 10);
		}
		# Amide
		if ($con_o && $con_n) {
			change_bondorder($atom, 'O', $mol, 2);
		}
		
		# Guanidiniums
		# Make double bond to first doubly substituted nitrogen
		if ($con_n == 3) {
			my $j = 0;
			foreach my $connum (@{$atom->{CONNECT}}) {
				
				my $con = $mol->{ATOMS}[$connum];
					
				if ($#{$con->{CONNECT}} > 0 || $j ==2) {
					bond_modify_order_check($mol, $atom, $con, 2);
					last;
				}
				++$j;
			}
		}
	}
	
	return undef;
}

sub fix_double_bonds_general {

	#<
	#? Fix general planar double bond orders
	#. Uses molecule geometry to generate assignment
	#. Contains a distance check
	#; Requires: molecule, ccflag (only generate double C=C double bonds)
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $ccflag = $_[1]; # Set this flag to only generate C=C double bonds
	
	my $tol = 10; # Tolerance for planarity
	
	# A          D
	#  \        /
	#   \      /
	#    B----C
	#   /      \
	#  /        \
	
	silico_msg('g', subname()."\n");
	my $maxvalence = maximum_valences();
	
	# modify maximum valences
	$maxvalence->{'Na'} = 0;
	$maxvalence->{'N'} = 3;

	my $atoms = $mol->{ATOMS};

	# Loop over all atoms in molecule
	# Atom B
	ATOMB: foreach my $atomb (@$atoms) {
		
		# Carbon-carbon double bonds only flag
		next if ($ccflag && $atomb->{ELEMENT} ne "C");

		# Check that atomb has more than one connected atom
		next if ($#{$atomb->{CONNECT}} < 1);

		my $numb = $atomb->{NUM}-1;
		
		# Do not continue if we are already at max valence
		my $v1 = valence($atomb);
		next if $v1 >= ($maxvalence->{$atomb->{ELEMENT}}  || 4);
		
		# Atom C
		ATOMC: foreach my $numc (@{$atomb->{CONNECT}}) {
		
			# Only calculate planarity of bonds in one direction
			next if $numb > $numc;
			
			my $atomc = $atoms->[$numc];
			
			# Carbon-carbon double bonds only flag
			next if ($ccflag && $atomc->{ELEMENT} ne "C");
			
			# Do not continue if we are already at max valence
			my $v2 = valence($atomc);
			next if $v2 >= ($maxvalence->{$atomc->{ELEMENT}} || 4);
						
			# Check that atomc has more than one connected atom
			next if ($#{$atomc->{CONNECT}} < 1);
			
			my $bo = find_bondorder2($atomb, $atomc, $mol);
			
			# Skip if this bond is already aromatic
			next if $bo == 4;
			
			# Skip if it is already double
			next if $bo == 2;
			
			# Distance check;
			# Make bond only if bondlength is less than 1.49 angstroms
			if (distance($atomb, $atomc) > 1.49) {
				next;
			}
			
			# Do not add double bonds to an atom in a ring that already has one
			# This is necessary to stop making ketenes in 5-membered rings when
			# hydrogens are missing
			if ($atomb->{RINGS}[0]) {
				foreach my $bo (@{$atomb->{BORDERS}}){
					next ATOMB if $bo > 1;
				}
			}

			if ($atomc->{RINGS}[0]) {
				foreach my $bo (@{$atomc->{BORDERS}}){
					next ATOMC if $bo > 1;
				}
			}

			#Atom A
			my $i = 0;
			my $sum = 0;
			ATOMA: foreach my $numa (@{$atomb->{CONNECT}}) {

				next if $numa == $numc;
				my $atoma = $atoms->[$numa];
								
				# Atom D
				ATOMD: foreach my $numd (@{$atomc->{CONNECT}}) {

					next if $numd == $numb;
					#next if $numd >= $numa; # Check each pair once
					
					my $atomd = $atoms->[$numd];
					my $torsion = (180/$Silico::Pi) * abs(dihedral($atoma, $atomb, $atomc, $atomd));
					my $diff = $torsion if ($torsion <= 90);
					$diff = 180-$torsion if ($torsion > 90);
					
					++$i;
					$sum += $diff;
				}
			}
			
			my $ave = $sum/$i;
			
			if ($ave > $tol) {
				next ATOMC;
			}
			
			bond_modify_order_check($mol, $atomb, $atomc, 2);
		}
	}
	
	return undef;
}

sub mol_score_hybridisation {
	
	#<
	#? Determine the hybridisation state of each atom
	#  in a molecule.
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	
	foreach my $atom (@{$mol->{ATOMS}}) {
		atom_score_hybridisation($atom, $mol);
	}
	
	if ($mol->{WARN}{ATOM_TOO_MANY_BONDS}) {
	
		silico_msg('w', "Found $mol->{WARN}{ATOM_TOO_MANY_BONDS} atoms with > 4 atoms connected.\n".
			"Bond types will for those atoms will be determined using length only. Use -warnings flag to print them out\n");
		
		if ($Silico::warnings) {
		
			foreach my $atom (@{$mol->{ATOM_TOO_MANY_BONDS}}) {
			
				print (heading("Atom"));
			
				print "Subid $atom->{SUBID}  Subname $atom->{SUBNAME}  Num $atom->{NUM}  Name $atom->{NAME}\n";
			
				print (heading("Connected atoms"));
		
				foreach (@{$atom->{CONNECT}}) {
					my $a = $mol->{ATOMS}[$_];
					print "Subid $a->{SUBID}  Subname $a->{SUBNAME}  Num $a->{NUM}  Name $a->{NAME}\n";
				}
			
				print "\n";
			}
		}
	}
	
	return undef;
}

sub atom_score_hybridisation {
	
	#<
	#? Attempt to determine the hybridisation scheme of an atom.
	#. Assumption is that all hydrogens are present.
	#; Requires: atom, molecule
	#; Returns: nothing
	#>
	
	my $atom = $_[0];
	my $mol = $_[1];
	
	my @angles;
	
	use vars qw($Tet);
	
	# Can't determine anything if no atoms are connected or we have no coordinates
	return if !defined $atom->{CONNECT}[0];

	my $numcon = $#{$atom->{CONNECT}} + 1;
	
	# Can't determine anything if only one atom is connected
	return if $numcon == 1;
	
	# If two or more atoms are connected, the geometry can be
	# computed. Only the most common geometries are described,
	# however.
	
	my $atoms = $mol->{ATOMS};
	
	if ($atom->{ANGLES}) {
	
		# Angles may be defined from read_gromacs_itp
		@angles = $atom->{ANGLES};
	
	} else {
	
		return if $mol->{HAS_NO_COORDINATES}; # We can't calculate anything
			
		for (my $i = 0; $i <= $#{$atom->{CONNECT}}; ++$i) {
		
			my $con = $atoms->[$atom->{CONNECT}[$i]];
		
			for (my $j = $i+1; $j <= $#{$atom->{CONNECT}}; ++$j) {
			
				my $con2 = $atoms->[$atom->{CONNECT}[$j]];
				my $angle = bndangle($con, $atom, $con2);
				push @angles, $angle;
			}
		}
	}
	
	my $mean_angle = mean(@angles);
	
	my $score120tet = ($mean_angle - $Tet)/(120 - $Tet);
	my $score180120 = ($mean_angle - 120)/(180-120);
	
	if ($numcon >= 2 && $numcon <= 4) {
		if ($score180120 > 0) {
			$atom->{GEOMETRY_SCORE} = 2 + $score180120;
		} elsif ($score120tet > 0) {
			$atom->{GEOMETRY_SCORE} = 1 + $score120tet;
		} else {
			$atom->{GEOMETRY_SCORE} = 1;
		}
		
	} elsif ($atom->{ELEMENT} eq "C" || $atom->{ELEMENT} eq "N" || $atom->{ELEMENT} eq "O") {
				
		++$mol->{WARN}{ATOM_TOO_MANY_BONDS};
		push @{$mol->{ATOM_TOO_MANY_BONDS}}, $atom;
	}
	
	return undef;
}

sub mol_fix_multiple_bonds_general {
	
	#<
	#? Determine the bondorder score of any bond.
	#
	#  Note: This subroutine is probably deficient in setting bondorders
	#  of 5-membered planar rings with poor bond lengths and no hydrogens.
	#
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	
	my @bonds2;
	my $d;
	my @forbidden = qw( H Du He Ne Ar Kr Xe Rn );
	my $lscore;
	my $maxvalence = maximum_valences();
	my $score;

	# Determine whether the molecule has hydrogens. If so, we can
	# be much more relaxed about relying on bond length, since a lot more info
	# will come from an atom's chemical environment.
	# These values may well need additional tweaking.
	my $n;
	if (mol_has_h($mol) > 1) {
		$n = 3;
	} else {
		$n = 12;
	}
	
	my $bonds = mol_find_all_bonds($mol);
	
	foreach my $bond (@$bonds) {
		
		my $g1 = 0;
		my $g2 = 0;
		my $gc = 0;
		
		my $anum1 = $bond->[0];
		my $anum2 = $bond->[1];
		$bond->[3] = 0;

		my $atom1 = $mol->{ATOMS}[$anum1];
		my $atom2 = $mol->{ATOMS}[$anum2];
		my $e1 = $atom1->{ELEMENT};
		my $e2 = $atom2->{ELEMENT};
		
		if ($e1 eq 'H' || $e2 eq 'H') {
			$bond->[2] = 1;
			next;
		}
		
		my $bl1 = bl("$e1-$e2") || 0;
		my $bl2 = bl("$e1=$e2") || 0;
		my $bl3 = bl("$e1\%$e2") || 0;
		my $d = distance($atom1, $atom2) || 0;
		
		# Be somewhat forgiving on the length.
		if ($d < $bl3) {
			$lscore = 2.95;
		} elsif ($d < $bl2) {
			$lscore = 1.95 + (($bl2 - $d)/($bl2 - $bl3));
		} elsif ($d < $bl1) {
			$lscore = 0.95 + (($bl1 - $d)/($bl1 - $bl2));
		} else {
			$lscore = 0.95;
		}
		
		# Work out the average score. $n is a variable used to 
		# alter the weighting given to the length. Defaults to 2.
		if (defined $atom1->{GEOMETRY_SCORE}) {
			$g1 = $atom1->{GEOMETRY_SCORE};
			++$gc;
		}
		
		if (defined $atom2->{GEOMETRY_SCORE}) {
			$g2 = $atom2->{GEOMETRY_SCORE};
			++$gc;
		}
		
		$score = ($g1 + $g2 + ($lscore * $n))/($gc + $n);
				
		$bond->[2] = $score;

		# Flag to indicate that both atoms are part of a planar ring
		
		$bond->[3] = 1 if defined $atom1->{PRINGS}[0] && defined $atom2->{PRINGS}[0];
		
		#print "$atom1->{NAME} $atom2->{NAME} g1 $g1 g2 $g2 lscore $lscore gc $gc score $score\n";
	}
	
	# Sort bonds by (1) belonging to a planar ring (2) their score. Highest first.
	@bonds2 = sort bsort @$bonds;
	
	BOND: foreach my $bond (@bonds2) {
		
		my $anum1 = $bond->[0];
		my $anum2 = $bond->[1];
		my $atom1 = $mol->{ATOMS}[$anum1];
		my $atom2 = $mol->{ATOMS}[$anum2];

		#print "$atom1->{ELEMENT} $atom1->{NUM} \t $atom2->{ELEMENT} $atom2->{NUM}\n";
		
		# Skip atoms of known unhelpful elements (i.e., hydrogen,
		# noble gases, dummy atoms)
		foreach (@forbidden) {
			next BOND if $atom1->{ELEMENT} eq $_;
			next BOND if $atom2->{ELEMENT} eq $_;
		}
		
		# Skip any bond of which either atom has no known maximum valence.
		# This should, in principle, remove certain "uninitialised value" errors.
		if (!defined $maxvalence->{$atom1->{ELEMENT}}) {
			silico_msg('w', "No known maximum valence for atom $atom1->{NUM} (\"$atom1->{NAME}\")!\n",
					"This bond will be skipped.\n");
			next BOND;
		} elsif (!defined $maxvalence->{$atom2->{ELEMENT}}) {
			silico_msg('w', "No known maximum valence for atom $atom2->{NUM} (\"$atom2->{NAME}\")!\n",
					"This bond will be skipped.\n");
			next BOND;
		}
		
		my $cbo = find_bondorder2($atom1, $atom2, $mol);
		
		# Do not modify aromatic rings that are already aromatic.
		next if $cbo == 4 && $atom1->{AROMATIC_RING} && $atom2->{AROMATIC_RING};
		
		$cbo = 1.5 if $cbo == 4;
		
		my $av1 = $maxvalence->{$atom1->{ELEMENT}};
		if (!defined $av1) {
			silico_msg('w',"Could not determine maximum valence of $atom1->{ELEMENT}\n");
			$av1 = 0;
		}
		$av1 = $av1- valence($atom1);

		my $av2 = $maxvalence->{$atom2->{ELEMENT}};
		if (!defined $av2) {
			silico_msg('w',"Could not determine maximum valence of $atom2->{ELEMENT}\n");
			$av2 = 0;
		}
		$av2 = $av2- valence($atom2);
		
		# Select the lowest of the two "available valences".
		my $av = ($av1 < $av2 ? $av1 : $av2);
				
		# Bond score. For some reason, this seems to work better
		# if the cutoffs are more forgiving towards lower-order
		# bonds.
		my $score;
		if ($bond->[2] > 2.5) {
			$score = 3;
		} elsif ($bond->[2] > 1.5) {
			$score = 2;
		} else {
			$score = 1;
		}
		
		# Target bond order
		my $tbo = ($score > ($cbo + $av) ? ($cbo + $av) : $score);
		
		#if ($tbo != $cbo) {
			#print "xx atom1 $atom1->{NAME} atom2 $atom2->{NAME} score $score av $av cbo $cbo tbo $tbo\n", "changing to bond order: $tbo\n";
		#}
		
		$tbo = 4 if $tbo == 1.5;
				
		# If the target bond order is greater than 1, change the bond order.
		# But, call bond_modify_order_check (instead of bond_modify_order) so
		# the resulting bond orders are sensible given the atom environments.
		if ($tbo > 1) {
			bond_modify_order_check($mol, $atom1, $atom2, $tbo); 
		}
	}
	
	return undef;
}

sub bsort {

	return 1 if $b->[3] > $a->[3];
	return 0 if $b->[3] < $a->[3];

	$b->[2] <=> $a->[2];
}


#################
# BOND CREATION #
#################

sub connect_atoms_by_distance {

        #<
        #? Bond atoms together that are closer than their combined vdW radii.
        #  A faster routine that uses atom bins to substantially speed things up.
        #  All existing bonds are deleted before new bonds are calculated.  All new bonds
        #  are created a single bonds.
        #. Note1: This routine utilises an overall maximum bond length of 2.3 Angstroms.
        #  I-I bonds fall outside this range
        #. Note2: Does not bond dummy atoms
        #. Note3: If a molecule contains 'Alternate atoms' eg from a pdb file - then
        #  atoms with specifier 'A' are bonded to atoms without an alternate specifier. Alternate atoms
        #  with specifier 'B' or greater are bonded only to themselves.
        #; Requires: molecule, flag to skip making bonds in water molecules
        #; Returns: number of bonds in molecule
        #>

        my $mol       = $_[0];
        my $waterflag = $_[1];
      
        my $bond = bondlengths();
        my $starttime = (times)[0];
	
	mol_check_cell($mol);
	
	my $cx = $mol->{CELL}[0];
	my $cy = $mol->{CELL}[1];
	my $cz = $mol->{CELL}[2];
	
	my $cx2 = $cx/2;
	my $cy2 = $cy/2;
	my $cz2 = $cz/2;

        my $maxbondlength = 2.3;
	
	# Defined here for speed. 'my' statements within loops
	# slow things down
	my ($atom1, $atom2);
	my ($radius1, $e1, $num1, $bin, $n);
	my ($x1, $y1, $z1);
	my ($xd, $yd, $zd, $xd2, $yd2, $zd2, $radius, $distsq);
	my ($bx, $by, $bz);

        # Ensure that atom numbers are correct
        # Delete any previous bonds but
	# don't delete bonds in waters if waterflag is set
        my $i = 0;
        foreach my $atom (@{$mol->{ATOMS}}) {
                ++$i;
                $atom->{NUM} = $i;
                next if $waterflag && $atom->{FG}{W};
                undef $atom->{CONNECT};
                undef $atom->{BORDERS};
        }

        # Just select nonwater atoms
	my $atoms;
	foreach my $atom (@{$mol->{ATOMS}}) {
		next if $waterflag && $atom->{FG}{W};
		next if $atom->{ELEMENT} eq 'Du';
		push @$atoms, $atom;
	}
	
	my $wlarge = 1 if $#{$atoms} > 100000;
	silico_msg('w', "Creating bonds in a large number of atoms: ".($#{$atoms}+1)."\n") if $wlarge;

     	calc_binsize($mol, $maxbondlength);
	my $bins = molecule_bin_atoms($mol, $atoms); 
	my $binsize = $mol->{BINSIZE} || croak();
	my $numbins =  $mol->{NUMBINS} || croak();

	my $count = 0;
	my $oldtime = (times)[0];
        my $numeval = 0;
	my $numbonds = 0;
	
        foreach $atom1 (@$atoms) {
	
		++$count;

		if ($count % 50000 == 0) {
			silico_msg('t', "Count: $count Time: ". calc_human_readable_time((times)[0] - $oldtime)."\n");
			print "Bonded $count atoms\n" if $wlarge;
		}

                $radius1 = $Vdw_Radii[ $atom1->{ELEMENT_NUM} ];
                $e1      = "$atom1->{ELEMENT}-"; # Note appended dash
		$num1    = $atom1->{NUM};

              	$bin = atom_calc_bin($atom1, $mol);

		$x1 = $atom1->{X};
                $y1 = $atom1->{Y};
                $z1 = $atom1->{Z};
		
		# Make lists of bins within +/- one bin of the current bin

		my $blist;
		foreach $i (0..2) {
			@{$blist->[$i]} = ($bin->[$i]-1..$bin->[$i]+1);
			foreach $n (@{$blist->[$i]}) {
				$n = $n % $numbins->[$i]; # Modulo
			}
		}
		
		foreach $bx (@{$blist->[0]}) {
		
                        foreach $by (@{$blist->[1]}) {
			
                               	 foreach $bz (@{$blist->[2]}) {
					
					foreach $atom2 (@{$bins->[$bx][$by][$bz]}) {
			
						# Atoms in bins are in numerical order
						last if $atom2->{NUM} >= $num1;
						
						#++$numeval;
						
               			             	if ( $radius = $bond->{$e1.$atom2->{ELEMENT}} ) {
							# Do nothing, radius is set
               			      		} else {
                			 		$radius = (($radius1 + $Vdw_Radii[ $atom2->{ELEMENT_NUM} ]) / 2);
                				}

                 			 	$xd = $x1 - $atom2->{X};
						# Note - could optimise this next statement
						# This calculation only needs to be done if the
						# distance is close to the cell length
						# or bins are a long way apart
						$xd2 = abs($xd - floor(($xd+$cx2)/$cx)*$cx);  
						next if $xd2 > $radius;
						
						$yd = $y1 - $atom2->{Y};
						$yd2 = abs($yd - floor(($yd+$cy2)/$cy)*$cy);
						next if $yd2 > $radius;
					
						$zd = $z1 - $atom2->{Z};
						$zd2 = abs($zd - floor(($zd+$cz2)/$cz)*$cz);
						next if $zd2 > $radius;
						
						$distsq = $xd2*$xd2+$yd2*$yd2+$zd2*$zd2;
						
						next if ($distsq > $radius*$radius);
						
						# Do not generate bond if atoms have different alternate designations
                           			next if ( 
							defined $atom1->{ALT}
                           				&& defined $atom2->{ALT}
                            				&& $atom1->{ALT} ne $atom2->{ALT}
                            				&& (uc($atom1->{ALT}) gt 'A' || uc($atom2->{ALT}) gt 'A'));
							
						# Bond is too short.  Don't make it it will only cause problems later
                      			  	if ($distsq < $radius*$radius*0.5) {
							++$mol->{WARN}{SHORT_BONDS};
							my $pair;
							@$pair = ($atom1, $atom2);
							push @{$mol->{SHORT_BOND_LIST}}, $pair;
							next;
						}
						
						# Dont make bonds between hydrogen atoms
						next if ($e1 eq 'H-' && $atom2->{ELEMENT} eq 'H');
						
	                			# Make bond
						bond_create($mol, $atom1->{NUM}-1, $atom2->{NUM}-1, 1);
						++$numbonds;
					}
				}
                        }
                }
        }
	
	
	$numbonds += rescue_unbonded_hydrogens($mol);
	if ($mol->{WARN}{SHORT_BONDS}) {
	
		silico_msg('w', "Found $mol->{WARN}{SHORT_BONDS} pairs of very close atoms (r < 0.1 Ang).\n". 
			"These atoms were not bonded. Use -warnings flag to print out atoms\n");
	
		if ($Silico::warnings) {
			foreach my $pair (@{$mol->{SHORT_BOND_LIST}}) {
		
				atom_printout($pair->[0]);
				print "---\n";
				atom_printout($pair->[1]);
				print "===========================\n";
			}
		}
	}
	
	if ($mol->{WARN}{OVERBONDED_ATOMS}) {
	
		silico_msg('w', "Found $mol->{WARN}{OVERBONDED_ATOMS} atoms with excess bonds. Use -warnings flag to print out atoms\n");
	
		if ($Silico::warnings) {
			foreach my $atom (@{$mol->{OVERBONDED_ATOMS_LIST}}) {
		
				atom_printout($atom);
				print "===========================\n";
			}
		}
	}
	
	if (0) {
		my $ratio = $numbonds/$numeval;
		silico_msg('t', "Time in " . subname() . ": " . calc_human_readable_time((times)[0] - $starttime) . "\n", 
			sprintf "Atoms bonded, not bonded, ratio:  %d,  %d,  %6.4f\n", $numbonds, $numeval, $ratio);
	}

        my $bondcount = molecule_count_numbonds($mol);
        $mol->{NUMBONDS} = $bondcount;
	
	silico_msg('c', "Finished bonding\n") if $wlarge;
        return $bondcount;
}


sub connect_atoms_by_distance_old {

        #<
        #? Bond atoms together that are closer than their combined vdW radii.
        #  A faster routine that uses atom bins to substantially speed things up.
        #  All existing bonds are deleted before new bonds are calculated.  All new bonds
        #  are created a single bonds.
        #. Note1: This routine utilises an overall maximum bond length of 2.3 Angstroms.
        #  I-I bonds fall outside this range
        #. Note2: Does not bond dummy atoms
        #. Note3: If a molecule contains 'Alternate atoms' eg from a pdb file - then
        #  atoms with specifier 'A' are bonded to atoms without an alternate specifier. Alternate atoms
        #  with specifier 'B' or greater are bonded only to themselves.
        #; Requires: molecule, flag to skip making bonds in water molecules
        #; Returns: number of bonds in molecule
        #>

        my $mol       = $_[0];
        my $waterflag = $_[1];
      
        my $bond = bondlengths();
        my $starttime = (times)[0];
	
	# Can not use cell information if angles are not 90 deg
	# Delete and warn if this is not true
	if (
		(defined $mol->{CELL_ALPHA} && defined $mol->{CELL_BETA} && defined $mol->{CELL_GAMMA}) &&
		!($mol->{CELL_ALPHA} == 90 && $mol->{CELL_BETA} == 90 && $mol->{CELL_GAMMA} == 90)) {
		
		silico_msg('w', "Can not use cell dimensions because cell angles are not all 90 deg\n",
			"\tDeleting cell information\n");
		undef $mol->{CELL};
		undef $mol->{CELL_ALPHA};
		undef $mol->{CELL_BETA};
		undef $mol->{CELL_GAMMA};
	}

	# Find molecule dimensions if no cell is defined or the cell is obviously rubbish
	# which can be true with a PDB file taken from electron diffraction data
	if (!defined $mol->{CELL}  || $mol->{CELL}[0] < 2) {
		@{$mol->{CELL}} = molecule_dimensions($mol, 10);
		silico_msg('w', "No good cell data\n","\tSetting system dimensions to @{$mol->{CELL}}\n");
	}
	
	my $cx = $mol->{CELL}[0];
	my $cy = $mol->{CELL}[1];
	my $cz = $mol->{CELL}[2];
	
	my $cx2 = $cx/2;
	my $cy2 = $cy/2;
	my $cz2 = $cz/2;

        my $maxbondlength = 2.3;
	
	# Defined here for speed. 'my' statements within loops
	# slow things down
	my ($atom1, $atom2);
	my ($radius1, $e1, $num1, $bin, $n);
	my ($x1, $y1, $z1);
	my ($xd, $yd, $zd, $xd2, $yd2, $zd2, $radius, $distsq);
	my ($bx, $by, $bz);

        # Ensure that atom numbers are correct
        # Delete any previous bonds but
	# don't delete bonds in waters if waterflag is set
        my $i = 0;
        foreach my $atom (@{$mol->{ATOMS}}) {
                ++$i;
                $atom->{NUM} = $i;
                next if $waterflag && $atom->{FG}{W};
                undef $atom->{CONNECT};
                undef $atom->{BORDERS};
        }

        # Just select nonwater atoms
	my $atoms;
	foreach my $atom (@{$mol->{ATOMS}}) {
		next if $waterflag && $atom->{FG}{W};
		next if $atom->{ELEMENT} eq 'Du';
		push @$atoms, $atom;
	}
	
	my $wlarge = 1 if $#{$atoms} > 100000;
	silico_msg('w', "Creating bonds in a large number of atoms: ".($#{$atoms}+1)."\n") if $wlarge;

     	calc_binsize($mol, $maxbondlength);
	my $bins = molecule_bin_atoms($mol, $atoms); 
	my $binsize = $mol->{BINSIZE} || croak();
	my $numbins =  $mol->{NUMBINS} || croak();

	my $count = 0;
	my $oldtime = (times)[0];
        my $numeval = 0;
	my $numbonds = 0;
	
        foreach $atom1 (@$atoms) {
	
		++$count;

		if ($count % 50000 == 0) {
			silico_msg('t', "Count: $count Time: ". calc_human_readable_time((times)[0] - $oldtime)."\n");
			print "Bonded $count atoms\n" if $wlarge;
		}

                $radius1 = $Vdw_Radii[ $atom1->{ELEMENT_NUM} ];
                $e1      = "$atom1->{ELEMENT}-"; # Note appended dash
		$num1    = $atom1->{NUM};

              	$bin = atom_calc_bin($atom1, $binsize, $numbins);

		$x1 = $atom1->{X};
                $y1 = $atom1->{Y};
                $z1 = $atom1->{Z};
		
		# Make lists of bins within +/- one bin of the current bin

		my $blist;
		foreach $i (0..2) {
			@{$blist->[$i]} = ($bin->[$i]-1..$bin->[$i]+1);
			foreach $n (@{$blist->[$i]}) {
				$n = $n % $numbins->[$i]; # Modulo
			}
		}
		
		foreach $bx (@{$blist->[0]}) {
		
                        foreach $by (@{$blist->[1]}) {
			
                               	 foreach $bz (@{$blist->[2]}) {
					
					foreach $atom2 (@{$bins->[$bx][$by][$bz]}) {
			
						# Atoms in bins are in numerical order
						last if $atom2->{NUM} >= $num1;
						
						#++$numeval;
						
               			             	if ( $radius = $bond->{$e1.$atom2->{ELEMENT}} ) {
							# Do nothing, radius is set
               			      		} else {
                			 		$radius = (($radius1 + $Vdw_Radii[ $atom2->{ELEMENT_NUM} ]) / 2);
                				}

                 			 	$xd = $x1 - $atom2->{X};
						# Note - could optimise this next statement
						# This calculation only needs to be done if the
						# distance is close to the cell length
						# or bins are a long way apart
						$xd2 = abs($xd - floor(($xd+$cx2)/$cx)*$cx);  
						next if $xd2 > $radius;
						
						$yd = $y1 - $atom2->{Y};
						$yd2 = abs($yd - floor(($yd+$cy2)/$cy)*$cy);
						next if $yd2 > $radius;
					
						$zd = $z1 - $atom2->{Z};
						$zd2 = abs($zd - floor(($zd+$cz2)/$cz)*$cz);
						next if $zd2 > $radius;
						
						$distsq = $xd2*$xd2+$yd2*$yd2+$zd2*$zd2;
						
						next if ($distsq > $radius*$radius);
						
						# Do not generate bond if atoms have different alternate designations
                           			next if ( 
							defined $atom1->{ALT}
                           				&& defined $atom2->{ALT}
                            				&& $atom1->{ALT} ne $atom2->{ALT}
                            				&& (uc($atom1->{ALT}) gt 'A' || uc($atom2->{ALT}) gt 'A'));
							
						# Bond is too short.  Don't make it it will only cause problems later
                      			  	if ($distsq < $radius*$radius*0.5) {
							++$mol->{WARN}{SHORT_BONDS};
							my $pair;
							@$pair = ($atom1, $atom2);
							push @{$mol->{SHORT_BOND_LIST}}, $pair;
							next;
						}
						
						# Dont make bonds between hydrogen atoms
						next if ($e1 eq 'H-' && $atom2->{ELEMENT} eq 'H');
						
	                			# Make bond
						bond_create($mol, $atom1->{NUM}-1, $atom2->{NUM}-1, 1);
						++$numbonds;
					}
				}
                        }
                }
        }
	
	$numbonds += rescue_unbonded_hydrogens($mol);

	if ($mol->{WARN}{SHORT_BONDS}) {
	
		silico_msg('w', "Found $mol->{WARN}{SHORT_BONDS} pairs of very close atoms (r < 0.1 Ang).\n". 
			"These atoms were not bonded. Use -warnings flag to print out atoms\n");
	
		if ($Silico::warnings) {
			foreach my $pair (@{$mol->{SHORT_BOND_LIST}}) {
		
				atom_printout($pair->[0]);
				print "---\n";
				atom_printout($pair->[1]);
				print "===========================\n";
			}
		}
	}
	
	if ($mol->{WARN}{OVERBONDED_ATOMS}) {
	
		silico_msg('w', "Found $mol->{WARN}{OVERBONDED_ATOMS} atoms with excess bonds. Use -warnings flag to print out atoms\n");
	
		if ($Silico::warnings) {
			foreach my $atom (@{$mol->{OVERBONDED_ATOMS_LIST}}) {
		
				atom_printout($atom);
				print "===========================\n";
			}
		}
	}
	
	if (0) {
		my $ratio = $numbonds/$numeval;
		silico_msg('t', "Time in " . subname() . ": " . calc_human_readable_time((times)[0] - $starttime) . "\n", 
			sprintf "Atoms bonded, not bonded, ratio:  %d,  %d,  %6.4f\n", $numbonds, $numeval, $ratio);
	}

        my $bondcount = molecule_count_numbonds($mol);
        $mol->{NUMBONDS} = $bondcount;
	
	silico_msg('c', "Finished bonding\n") if $wlarge;
        return $bondcount;
}

sub rescue_unbonded_hydrogens {

	#<
	#? Fix unbonded hydrogens
	#; Requires: molecule
	#; Returns: total number of bonds added
	#>

	my $mol = $_[0];
	my $print = get_sflag('warnings');
	
	my $bondcount = 0;
	
	# Fix hydrogens with zero connected atoms
        foreach my $atom1 (@{$mol->{ATOMS}}) {
		
		# Hydrogens
                if ($atom1->{ELEMENT} eq 'H' && !defined $atom1->{CONNECT}[0]) {

			my $at;
                        my $min = 999;

                        foreach my $atom2 (@{$mol->{ATOMS}}) {

                                next if $atom1 == $atom2;

				my $d2;
				$d2 = general_distance_sq($atom1, $atom2, $mol);
                                next if $d2 > $min;
                                $min = $d2;
                                $at  = $atom2;
                        }
				
                        if ($min < 4) {
                                
				bond_create($mol, $atom1->{NUM} - 1, $at->{NUM} - 1, 1);
					
				++$bondcount;
				
				if ($print && $bondcount < 10) {
					my $d = sprintf "%6.4f", sqrt($min);
					my $dn = sprintf "%6.4f", sqrt(general_distance_sq($atom1, $at, $mol));
					print "Creating bond to hydrogen with length $d ($dn). '$atom1->{SUBNAME}' '$atom1->{SUBID}' '$atom1->{NAME}' '$atom1->{ELEMENT}' - '$at->{SUBNAME}' '$at->{SUBID}' '$at->{NAME}' '$at->{ELEMENT}'\n";
					print "Warnings truncated\n" if $bondcount == 9;
				}
                        }
		}
	}

	return $bondcount;
}


sub connect_atoms_by_distance_simple {

	#<
	#? Bond atoms together that are closer than their combined vdW radii.
	#. Very simple connect atoms routine for testing
	#; Requires: molecule
	#; Returns: number of bonds created
	#>
	
	my $mol = $_[0];

	my $starttime = (times)[0];

	my $atoms = $mol->{ATOMS};
	my $bond2  = bondlengths_2();
	
	# Remove existing bonds
	bond_delete_all($mol);

	my $i = -1;
	my $bondcount = 0;
	foreach my $atom1 (@$atoms) {
	
		++$i;
	
		my $x1 = $atom1->{X};
		my $y1 = $atom1->{Y};
		my $z1 = $atom1->{Z};
		
		my $e1  = $atom1->{ELEMENT};
		my $radius1 = $Vdw_Radii[ $atom1->{ELEMENT_NUM} ];

		my $j = -1;
		foreach my $atom2 (@$atoms) {

			++$j;
		
			last if $j >= $i;	
			
			# Calculate bonding distance
                    	my $e2  = $atom2->{ELEMENT};
              		my $key = $e1 . "-" . $e2;
			my $radius;    
               		if (defined $bond2->{$key}) {
               			$radius = sqrt($bond2->{$key});
               		} else {
                		$radius = ($radius1 + $Vdw_Radii[$atom2->{ELEMENT_NUM}]) / 2;
                	}
			
			my $a = $x1-$atom2->{X};
			next if $a > $radius;
			
			my $b = $y1-$atom2->{Y};
			next if $b > $radius;
			
			my $c = $z1-$atom2->{Z};
			next if $c > $radius;
			
			next if sqrt($a**2 + $b**2 + $c**2) > $radius;
			
			$bondcount += bond_create($mol, $i, $j, 1);
		}
	}

	silico_msg('c', "Finished in ".subname().". Made $bondcount bonds\n");
	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
	
	return $bondcount;
}

sub connect_atoms_by_distance_simple_bv_periodic {

	#<
	#? Bond atoms together that are closer than their combined vdW radii.
	#. Very simple connect atoms routine for testing
	#; Requires: molecule
	#; Returns: number of bonds created
	#>
	
	my $mol = $_[0];

	my $starttime = (times)[0];
	my $atoms = $mol->{ATOMS};
	my $mat = $mol->{BV_MAT} || croak();
	my $matinv = $mol->{BV_MATINV} || croak();
	my $bond2  = bondlengths_2();
	
	# Remove existing bonds
	bond_delete_all($mol);
	
	mol_check_cell($mol);

	my $i = -1;
	my $bondcount = 0;
	foreach my $atom1 (@$atoms) {
	
		++$i;
	
		my $x1 = $atom1->{X};
		my $y1 = $atom1->{Y};
		my $z1 = $atom1->{Z};
		
		my $e1  = $atom1->{ELEMENT};
		my $radius1 = $Vdw_Radii[ $atom1->{ELEMENT_NUM} ];

		my $j = -1;
		
		my @r;
		my ($a, $b, $c);
		foreach my $atom2 (@$atoms) {

			++$j;
		
			last if $j >= $i;	
			
			# Calculate bonding distance
                    	my $e2  = $atom2->{ELEMENT};
              		my $key = $e1 . "-" . $e2;
			my $radius;    
               		if (defined $bond2->{$key}) {
               			$radius = sqrt($bond2->{$key});
               		} else {
                		$radius = ($radius1 + $Vdw_Radii[$atom2->{ELEMENT_NUM}]) / 2;
                	}
			
			($a, $b, $c) = distance_periodic_points_xyz_mat($x1, $y1, $z1, $atom2->{X}, $atom2->{Y}, $atom2->{Z}, $mat, $matinv);
			
			next if $a > $radius;
			
			next if $b > $radius;
		
			next if $c > $radius;
			
			next if sqrt($a**2 + $b**2 + $c**2) > $radius;
			
			$bondcount += bond_create($mol, $i, $j, 1);
		}
	}

	silico_msg('c', "Finished in ".subname().". Made $bondcount bonds\n");
	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
	
	return $bondcount;
}


sub connect_atoms_in_waters {

	#<
	#? Connect atoms in water molecules quickly by using the fact that water residues have 3 atoms
	#; Requires: molecule
	#; Returns: number of bonds made

	my $mol = $_[0];

	my $bondcount = 0;
	my $count = 0;
	my $residues;
	my $res;
	my $starttime = (times)[0];
	
	# Renumber atoms to ensure that atoms are numbered in order
	molecule_renumber($mol);

	$residues = molecule_get_residues($mol);
	
	foreach my $res (@$residues) {

		# Water molecules all have three atoms.
		next if $#{$res} != 2;

		# These three atoms must be an oxygen and two hydrogens.
		if ($res->[0]{ELEMENT} eq 'O' && $res->[1]{ELEMENT} eq 'H' && $res->[2]{ELEMENT} eq 'H') {

			# Create bonds
			$bondcount += bond_create($mol, $res->[0]{NUM}-1, $res->[1]{NUM}-1);
			$bondcount += bond_create($mol, $res->[0]{NUM}-1, $res->[2]{NUM}-1);
		
		} elsif ($res->[0]{ELEMENT} eq 'H' && $res->[1]{ELEMENT} eq 'O' && $res->[2]{ELEMENT} eq 'H') {
		
			# Create bonds
			$bondcount += bond_create($mol, $res->[0]{NUM}-1, $res->[1]{NUM}-1);
			$bondcount += bond_create($mol, $res->[1]{NUM}-1, $res->[2]{NUM}-1);
		
		} elsif ($res->[0]{ELEMENT} eq 'H' && $res->[1]{ELEMENT} eq 'H' && $res->[2]{ELEMENT} eq 'O') {
		
			# Create bonds
			$bondcount += bond_create($mol, $res->[0]{NUM}-1, $res->[2]{NUM}-1);
			$bondcount += bond_create($mol, $res->[1]{NUM}-1, $res->[2]{NUM}-1);
		
		} else {
			next;
		}
		
		# Label water atoms
		$res->[0]{FG}{W} = 1;
		$res->[1]{FG}{W} = 1;
		$res->[2]{FG}{W} = 1;
	
		++$count;
	}

	silico_msg('c', "Created $bondcount bonds in $count water molecules.\n") if $count;
	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
	
	return $bondcount;
}


sub bond_create {

	#<
	#? Create a bond between two atoms (takes atom offsets numbered from ZERO as arguments).
	#. Creates a bond in both directions. ie makes connections atom2->atom1 and atom1->atom2
	#. Silico internal bond types:
	#.	1 single, 2 double, 3 triple, 4 aromatic,
	#.	5 single or double, 6 single or aromatic, 7 double or aromatic,
	#.	8 any (undefined), 0 zero order.
	#; Requires: molecule, atomnum1, atomnum2, bondorder, flag to not increment mol2->NUMBONDS.
	#>

	my $mol = $_[0];
	my $atomnum1 = $_[1];
	my $atomnum2 = $_[2];
	my $order = $_[3] || 1;
	my $flag = $_[4];
	
	my $atom1 = $mol->{ATOMS}[$atomnum1];
	my $atom2 = $mol->{ATOMS}[$atomnum2];
	
	if ($atomnum1 == $atomnum2) {
		silico_msg('e', "Attempt to create a bond from an atom to itself!\n",
				"Atom number: ".($atomnum1+1)."\n");
		return 0;
	}
	
	if ($atomnum1 > $#{$mol->{ATOMS}} || $atomnum2 > $#{$mol->{ATOMS}}) {
		silico_msg('e', "Attempt to create a bond to a nonexistent atom!\n",
				"Atom 1: ".($atomnum1+1).", atom 2: ".($atomnum2+1).", number of atoms: ".($#{$mol->{ATOMS}}+1)."\n");
		return 0;
	}
	
	push @{$atom1->{CONNECT}}, $atomnum2;
	push @{$atom1->{BORDERS}}, $order;
	
	push @{$atom2->{CONNECT}}, $atomnum1;
	push @{$atom2->{BORDERS}}, $order;
	
	++$mol->{NUMBONDS} if !$flag;
	
	return 1;
}

sub bond_create_atom {

	#<
	#? Create a bond between two atoms (takes atoms as arguments).
	#. See bond_create
	#>

	my $mol = $_[0];
	my $atom1 = $_[1] || croak();
	my $atom2 = $_[2] || croak();
	my $order = $_[3];
	my $flag = $_[4];

	# Check that atoms are numbered continuously
	if ($mol->{ATOMS}[$atom1->{NUM}-1] != $atom1 || $mol->{ATOMS}[$atom2->{NUM}-1] != $atom2) {
		
		atom_printout($atom1);
		atom_printout($atom2);

		molecule_printout($mol);
		
		silico_msg('d', "Bond_create_atom requires that \$atom->{NUM} is numbered consecutively starting from 0\n");
	}

	return bond_create($mol, $atom1->{NUM}-1, $atom2->{NUM}-1, $order, $flag);
}

#################
# BOND DELETION #
#################

sub bond_delete_all {
	
	#<
	#? Delete all bonds
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];

	# Delete all existing bonds
	foreach my $atom (atoms($mol)) {

		delete $atom->{CONNECT};
		delete $atom->{BORDERS};
	}

	$mol->{NUMBONDS} = 0;
	delete $mol->{RINGS};
}

sub bond_delete {

	#<
	#? Delete a bond between two atoms.
	#; Requires: molecule, atomnum1, atomnum2
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $num1 = $_[1];
	my $num2 = $_[2];
	
	my $atom1 = $mol->{ATOMS}[$num1];
	my $atom2 = $mol->{ATOMS}[$num2];
	
	my $newcon1;
	my $newcon2;
	my $newborder1;
	my $newborder2;
	my $del = 0;
	
	my $i = 0;
	foreach my $connum (@{$atom1->{CONNECT}}) {
		
		if ($connum == $num2) {
			++$del;
		} else {
			
			push @$newcon1, $connum;
			push @$newborder1, $atom1->{BORDERS}[$i];
		}
		++$i;
	}
	$atom1->{CONNECT} = $newcon1;
	$atom1->{BORDERS} = $newborder1;
	
	$i =0;
	foreach my $connum (@{$atom2->{CONNECT}}) {
		
		if ($connum == $num1) {
			++$del;
		} else {
			push @$newcon2, $connum;
			push @$newborder2, $atom2->{BORDERS}[$i];
		}
		++$i;
	}
	$atom2->{CONNECT} = $newcon2;
	$atom2->{BORDERS} = $newborder2;
	
	$del == 0 && silico_msg('w', "Bond to be deleted ($num1, $num2) was not found!\n");
	$del == 1 && silico_msg('w', "Only half a bond ($num1, $num2) was deleted!\n");
		
	--$mol->{NUMBONDS};
	
	return;
}

sub bond_delete_atom {

	#<
	#? Delete a bond between two atoms.
	#. Same as bond_delete but using atoms instead of atom numbers.
	#  this version is to be preferred
	#; Requires: molecule, atom1, atom2
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $atom1 = $_[1];
	my $atom2 = $_[2];
	
	my $newcon1;
	my $newcon2;
	my $newborder1;
	my $newborder2;
	my $del = 0;
	
	my $i = 0;
	foreach my $connum (@{$atom1->{CONNECT}}) {
		
		my $con = $mol->{ATOMS}[$connum];
		
		if ($con == $atom2) {
			++$del;
		} else {
			
			push @$newcon1, $connum;
			push @$newborder1, $atom1->{BORDERS}[$i];
		}
		++$i;
	}
	$atom1->{CONNECT} = $newcon1;
	$atom1->{BORDERS} = $newborder1;
	
	$i =0;
	foreach my $connum (@{$atom2->{CONNECT}}) {
	
		my $con = $mol->{ATOMS}[$connum];
		
		if ($con == $atom1) {
			++$del;
		} else {
			push @$newcon2, $connum;
			push @$newborder2, $atom2->{BORDERS}[$i];
		}
		++$i;
	}
	$atom2->{CONNECT} = $newcon2;
	$atom2->{BORDERS} = $newborder2;
	
	$del == 0 && silico_msg('w', "Bond to be deleted ($atom1->{NUM}, $atom2->{NUM}) was not found!\n");
	$del == 1 && silico_msg('w', "Only half a bond ($atom1->{NUM}, $atom2->{NUM}) was deleted!\n");
		
	--$mol->{NUMBONDS};
	
	return;
}

#####################
# BOND MODIFICATION #
#####################

sub bond_modify_order {

	#<
	#? Modify the order of an existing bond
	#  Note: this requires that atom->{NUM} is numbered correctly
	#. Modifies bond in both directions.
	#. Silico internal bond types:
	#	1 single, 2 double, 3 triple, 4 aromatic,
	#	5 single or double, 6 single or aromatic, 7 double or aromatic,
	#	8 any (undefined), 0 zero order.
	#; Requires: Molecule, atom 1, atom 2, bondorder (default 1)
	#>

	my $mol = $_[0]; # Not used!
	my $atom1 = $_[1] || croak();
	my $atom2 = $_[2] || croak();
	my $order = $_[3] || 1;

	silico_msg('d', "Atoms have undefined atom numbers\n") if !defined $atom1->{NUM} || !defined $atom2->{NUM};

	my $atomnum1 = $atom1->{NUM}-1;
	my $atomnum2 = $atom2->{NUM}-1;

	#carp();

	#print "$atomnum1 $atom1->{NAME} $atomnum2 $atom2->{NAME} $order\n";

	my $i = 0;
	my $flag1;
	my $flag2;
	foreach my $atomnum (@{$atom1->{CONNECT}}) {
		if ($atomnum == $atomnum2) {
			$atom1->{BORDERS}[$i] = $order;
			$flag1 = 1;
			last;
		}
		++$i;
	}
	
	$i = 0;
	foreach my $atomnum (@{$atom2->{CONNECT}}) {
		if ($atomnum == $atomnum1) {
			$atom2->{BORDERS}[$i] = $order;
			$flag2 = 1;
			last;
		}
		++$i;
	}
	
	if (!$flag1 || !$flag2) {
		silico_msg('w', "Ill-defined or missing bond between atoms ".($atomnum1+1)." and ".($atomnum2+1)."!\n",
				"Bond orders have not been properly modified.\n");
	}
	
	return 1;
}

sub bond_modify_order_check {

	#<
	#? Modify the order of a bond but check that we are not breaking
	#  bonding rules before doing it
	#. Calls bond_modify_order
	#; Requires: Molecule, atom 1, atom 2, order
	#; Returns: 1 if bond was modified 0 if not
	#>

	my $mol = $_[0];
	my $atom1 = $_[1] || croak();
	my $atom2 = $_[2] || croak();
	my $order = $_[3];

	my $current_bond = find_bondorder2($atom1, $atom2, $mol);
	my $v1 = valence($atom1);
	my $v2 = valence($atom2);

	#print "a1 $atom1->{ELEMENT} a2 $atom2->{ELEMENT}  current $current_bond v1 $v1 v2 $v2 order $order\n";

	# If we have carbon or nitrogen and already have one double bond
	# then atom1 and atom2 MUST be linear to add an additional double bond.

	my $i = 0;
	foreach my $atom ($atom1, $atom2) {
		++$i;
		my $doublecount = 0;
		if ($order == 2 && ($atom->{ELEMENT} eq 'C' || $atom->{ELEMENT} eq 'N')) {

			my $j = 0;
			foreach (@{$atom->{BORDERS}}) {
				++$doublecount if ($_ == 2);
			}
			
			# check that we have more than one connected atom before measuring angle
			if ($doublecount) {
			
				return (0) if ($#{$atom->{CONNECT}} < 1);
			
				my $con1 = $mol->{ATOMS}[$atom->{CONNECT}[0]];
				my $con2 = $mol->{ATOMS}[$atom->{CONNECT}[1]];
						
				if (bndangle($con1, $atom, $con2) < 170) {
					return 0;
				}
			}
		}
	}

	my $c1 = checkvalence($atom1, $order);
	my $c2 = checkvalence($atom2, $order);

	#print "bmoc: $atom1->{NAME} c1 $c1 $atom2->{NAME} c2 $c2\n";
	
	if ($c1 && $c2) {
		my $val = bond_modify_order($mol, $atom1, $atom2, $order);
		return 1;
	}
		
	return 0;
}

sub checkvalence {

	#<
	#? Check if a proposed valence is less than the maximum valence for an atom
	#  taking into account atomic charge (where present) or whether an
	#  atom is an aromatic ring and might be charged
	#; Requires: atom, proposed valence
	#; Returns: 1 if OK 0 if not
	#>

	my $atom = $_[0];
	my $valence = $_[1];

	my $maxvalence = maximum_valences();

	my $el = $atom->{ELEMENT};
	my $mv = $maxvalence->{$atom->{ELEMENT}};

	return 1 if !defined $mv; # Unusual element or dummy atom

	# Increase maximum valence if 1) we have a positive charge or 2)
	# the atom is in a planar ring and is already at the maximum valence
	if (
		($atom->{CHARGE} && $atom->{CHARGE} == 1) ||
		($el ne 'C' && $atom->{PLANAR_RING} && ($#{$atom->{BONDS}}+1)==$mv)
	) {
		++$mv;
	}

	# Reduce the maximum valence if we have a negative charge
	if ($atom->{CHARGE} && $atom->{CHARGE} == -1) {
		--$mv;
	}

	return 1 if $valence <= $mv;
	return 0;
}


#################
# BOND QUERYING #
#################


sub find_bondorder {

	#<
	#? Returns the bondorder between two atoms
	#; Requires: atomnum1, atomnum2, molecule
	#; Returns: bondorder or -1 if not connected
	#>
	
	#$_[0] atomnum1
	#$_[1] atomnum2
	my $mol = $_[2] || die;
	
	my $atom1 = $mol->{ATOMS}[$_[0]];

	# Loop over all atoms connected to atom1 and find if
	# we are connected to atom2
	my $i = -1;
	foreach my $connum (@{$atom1->{CONNECT}}) {

		++$i;
		next if $connum != $_[1];
			
		return  $atom1->{BORDERS}[$i];
	}

	return  -1;
}

sub find_bondorder2 {

	#<
	#? Returns the bondorder between two atoms
	#. Same as find_bondorder but uses atoms as arguments
	#  rather than atom numbers.  
	#; Requires: atom1, atom2
	#; Returns: bondorder or -1 if not connected
	#>
	
	#my $atom1 = $_[0];
	#my $atom2 = $_[1];
	#my $mol = $_[2];
	
	die if !defined $_[2];
	
	# Loop over atoms connected to atom1 and find if
	# we are connected to atom2
	
	my $i = -1;
	foreach my $connum (@{$_[0]->{CONNECT}}) {
	
		++$i;
		next if $_[2]->{ATOMS}[$connum] != $_[1];
		return $_[0]->{BORDERS}[$i];
		last;
	}
	
	return -1;
}



##########################
# SUPPORTING INFORMATION #
##########################

sub bl {
	
	#<
	#? Return an expected bond length based on an input string.
	#  Single bonds A-B, double A=B, triple A%B
	#; Requires: String
	#; Returns: number, or the undefined value if no bond length
	#>
	
	standard_bond_lengths() if !defined $Silico::Bond_Lengths;
	
	my $string = $_[0];
	if (defined $Silico::Bond_Lengths{"$string"}) {
		return $Silico::Bond_Lengths{"$string"};
	} else {
		return undef;
	}
}


sub bondlengths {

	#<
	#? Routine that returns a hash of maximum reasonable bondlengths for common
	#  atom pairs
	#. This routine is used by the various 'connect_atom' routines.  For single bonds, the bond lengths
	#  are maximum allowable lengths and are generally 10% longer than
	#  standard bond lengths. For double and triple bonds they are 5% larger.
	#. Note. Some bondlengths are longer than the maximum bondlength in the connect routine.
	#;  Requires: nothing
	#;  Returns: Hash maximum bondlengths.  Keys have the form El_El where 'El' is
	#   the atom element
	#>
	
	my $bond;
	
	%$bond = qw (
		 H-H	0.814	 H-C 	1.2	 H-N 	1.15	 H-O 	1.15	 H-F 	1.012	 H-Si 	1.595	 
		 H-Na	1.9	 H-P 	1.518	 H-S 	1.452	 H-Cl 	1.397	 H-Br 	1.562	 H-I	1.771
		 C-H	1.2	 C-C 	1.694	 C-N 	1.617	 C-O 	1.573	 C-F 	1.551	 C-Si 	2.134	 
		 C-P 	2.057	 C-S 	1.991	 C-Cl 	1.936	 C-Br 	2.101	 C-I	2.31
		 N-H	1.15	 N-C 	1.617	 N-N 	1.54	 N-O 	1.496	 N-F 	1.474	 N-Si 	2.057	 
		 N-P 	1.98	 N-S 	1.914	 N-Cl 	1.859	 N-Br 	2.024	 N-I	2.233
		 O-H	1.15	 O-C 	1.573	 O-N 	1.496	 O-O 	1.452	 O-F 	1.43	 O-Si 	2.013	 
		 O-P 	1.936	 O-S 	1.87	 O-Cl 	1.815	 O-Br 	1.98	 O-I	2.189
		 F-H	1.012	 F-C 	1.551	 F-N 	1.474	 F-O 	1.43	 F-F 	1.408	 F-Si 	1.991	 
		 F-P 	1.914	 F-S 	1.848	 F-Cl 	1.793	 F-Br 	1.958	 F-I	2.167
		 Si-H	1.595	 Si-C 	2.134	 Si-N 	2.057	 Si-O 	2.013	 Si-F 	1.991	 Si-Si 	2.574	 
		 Si-P 	2.497	 Si-S 	2.431	 Si-Cl 	2.376	 Si-Br 	2.541	 Si-I	2.75
		 Na-H	1.9	 P-H	1.518	 P-C 	2.057	 P-N 	1.98	 P-O 	1.936	 P-F 	1.914	 P-Si 	2.497	 
		 P-P 	2.42	 P-S 	2.354	 P-Cl 	2.299	 P-Br 	2.464	 P-I	2.673
		 S-H	1.452	 S-C 	1.991	 S-N 	1.914	 S-O 	1.87	 S-F 	1.848	 S-Si 	2.431	 
		 S-P 	2.354	 S-S 	2.288	 S-Cl 	2.233	 S-Br 	2.398	 S-I	2.607
		 Cl-H	1.397	 Cl-C 	1.936	 Cl-N 	1.859	 Cl-O 	1.815	 Cl-F 	1.793	 Cl-Si 	2.376	 
		 Cl-P 	2.299	 Cl-S 	2.233	 Cl-Cl 	2.2	 Cl-Br 	2.343	 Cl-I	2.552
		 Br-H	1.562	 Br-C 	2.101	 Br-N 	2.024	 Br-O 	1.98	 Br-F 	1.958	 Br-Si 	2.541	 
		 Br-P 	2.464	 Br-S 	2.398	 Br-Cl 	2.343	 Br-Br 	2.508	 Br-I	2.717
		 I-H	1.771	 I-C 	2.31	 I-N 	2.233	 I-O 	2.189	 I-F 	2.167	 I-Si 	2.75	 
		 I-P 	2.673	 I-S 	2.607	 I-Cl 	2.552	 I-Br 	2.717	 I-I	2.926

		C=C	1.41	C=N	1.33	C=O	1.28	N=O	1.21	N=C	1.33	O=C	1.28	O=N	1.21

		C%C	1.27	C%N	1.21	C%O	1.18	N%O	1.13	N%C	1.21	O%C	1.13	O%N	1.13

		);
		
	return $bond;
}



sub bondlengths_2 {

	#<
	#? Routine that returns a hash of squared bondlengths.  See bondlengths.
	#;  Requires: nothing
	#;  Returns: Hash maximum bondlengths squared.  Keys have the form El_El where 'El' is
	#   the atom element
	#>
	
	my $bond = bondlengths();
	my $bond2;
	
	foreach (keys %$bond) {
		
		$bond2->{$_} = $bond->{$_}**2;
	}
	
	return $bond2;
}





#################
# ATOM ADDITION #
#################


sub copy_atom {

	#<
	#? Return a new copy of an atom.
	#. Warning.  This copy is only one level deep and does not include hashes.
	#; Requires: atom.
	#; Returns: new atom.
	#>

	my $atom = $_[0];
	my $newatom;
	
	foreach my $key (keys %$atom) {
						
		if ($#{$atom->{$key}} == -1) {

			# Copy scalars
			$newatom->{$key} = $atom->{$key};
		} else {
					
			# Copy arrays
			@{$newatom->{$key}} = @{$atom->{$key}};
		}
	}
	
	return $newatom;
}

sub mol_add_atom {

	#<
	#? Create a new atom and add it to the specified molecule.
	#; Requires: molecule, atom name, atom element, x, y, z,
	#  subname, subid.
	#; Returns: the new atom
	#>

	my $mol = $_[0];

	my $atom;

	++$mol->{NUMATOMS};
	my $num = $mol->{NUMATOMS} -1;

	# Set values or defaults
	$atom->{NAME} = ($_[1] || "X");
	$atom->{ELEMENT} = ($_[2] || "Du");
	$atom->{ELEMENT_NUM} = element_symbol2number($atom->{ELEMENT});
	$atom->{X} = ($_[3] || 0);
	$atom->{Y} = ($_[4] || 0);
	$atom->{Z} = ($_[5] || 0);
	$atom->{NUM} = $mol->{NUMATOMS};
	$atom->{SUBNAME} = ($_[6] || "UNK");
	$atom->{SUBID} = ($_[7] || 9999);
	$atom->{TYPE} = '';

	$mol->{ATOMS}[$num] = $atom;
	return $atom;
}


#################
# ATOM DELETION #
#################


sub molecule_delete_atom {

	#<
	#? Deletes a single atom by offset.
	#. Warning: You will need to run 'molecule_pack' after deleting
	#  a set of atoms or all hell will break loose.
	#; Requires: molecule, number of atom to be deleted.
	#; Returns: 1
	#>
	
	my $mol = $_[0];
	my $offset = $_[1];

	die "Offset undefined\n" if !defined $_[1];
	
	my $atom = $mol->{ATOMS}[$offset];
	
	# Delete any bonds to this atom
	foreach my $connum (@{$atom->{CONNECT}}) {
		
		my $con = $mol->{ATOMS}[$connum];
		
		# Skip if connected atom has already been deleted
		next if !defined $con;
		
		my $new_connection;
		my $new_borders;
		
		# Loop through atoms connected to $con
		my $i = -1;
		foreach my $connection (@{$con->{CONNECT}}) {
			++$i;

			next if $connection == $offset;

			push @$new_connection, $connection;
			push @$new_borders, ${$con->{BORDERS}}[$i];
		}
		
		$con->{CONNECT} = $new_connection;
		$con->{BORDERS} = $new_borders;
		
		# Adjust number of bonds
		--$mol->{NUMBONDS};
	}
		
	# Adjust number of atoms
	--$mol->{NUMATOMS};
	
	# Mark atom as undefined
	$mol->{ATOMS}[$offset] = undef;
	
	return 1;
}

sub molecule_delete_atoms {
	
	#<
	#? Delete a list of atoms, identified via atom offsets.
	#; Requires: molecule, list
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $list = $_[1];
	
	my @angles = ();
	my @bonds = ();
	my %deleted;
	my @dihedrals = ();
	my $improper;
	my @impropers = ();
	
	foreach my $offset (@$list) {
		
		molecule_delete_atom($mol, $offset);
		$deleted{$offset} = 1;
	}
	
	# Take out bonds, angles, dihedrals and impropers
	foreach my $bond (@{$mol->{BONDS}}) {
		next if defined $bond->[0] && $deleted{$bond->[0]};
		next if defined $bond->[1] && $deleted{$bond->[1]};
		push @bonds, $bond;
	}
	
	foreach my $angle (@{$mol->{ANGLES}}) {
		next if defined $angle->[0] && $deleted{$angle->[0]};
		next if defined $angle->[1] && $deleted{$angle->[1]};
		next if defined $angle->[2] && $deleted{$angle->[2]};
		push @angles, $angle;
	}
	
	foreach my $dihedral (@{$mol->{DIHEDRALS}}) {
		next if defined $dihedral->[0] && $deleted{$dihedral->[0]};
		next if defined $dihedral->[1] && $deleted{$dihedral->[1]};
		next if defined $dihedral->[2] && $deleted{$dihedral->[2]};
		next if defined $dihedral->[3] && $deleted{$dihedral->[3]};
		push @dihedrals, $dihedral;
	}
	
	foreach my $improper (@{$mol->{P_IMPROPERS}}) {
		next if defined $improper->[0] && $deleted{$improper->[0]};
		next if defined $improper->[1] && $deleted{$improper->[1]};
		next if defined $improper->[2] && $deleted{$improper->[2]};
		next if defined $improper->[3] && $deleted{$improper->[3]};
		push @impropers, $improper;
	}
	
	@{$mol->{BONDS}} = @bonds;
	@{$mol->{ANGLES}} = @angles;
	@{$mol->{DIHEDRALS}} = @dihedrals;
	@{$mol->{P_IMPROPERS}} = @impropers;
}

sub molecule_delete_dummy {
	
	#<
	#? Delete all dummy atoms (eg lone pairs) from a molecule
	#; Requires: Molecule, optional list identifying a subset
	#; Returns: number of deleted atoms
	#>
	
	my $mol = $_[0];
	my $list = $_[1];
	
	my $delcount = 0;
	my $dellist;
		
	my $i = 0;
	@$dellist = ();
		
	foreach my $atom (@{$mol->{ATOMS}}) {
		
		my $flag;
		
		if (defined $list) {
			foreach my $atom2 (@$list) {
				$flag = 1 if $atom == $atom2;
			}
		} else {
			$flag = 1;
		}
		
		# Delete dummy atoms
		if  ($atom->{ELEMENT_NUM} == 0 && $flag) {
			push @$dellist, $i;
			++$delcount;
		}
		
		++$i;
	}
		
	# Only do any actual deletion if dummy atoms have been found.
	if ($delcount > 0) {
		molecule_delete_atoms($mol, $dellist);
		silico_msg('c', "$delcount dummy atoms deleted\n");
		
		# Pack and renumber remaining atoms
		molecule_pack($mol);	
		$mol->{NUMBONDS} = molecule_count_numbonds($mol);
	}
		
	return $delcount;
}

sub mol_delete_duplicate_atoms {
	
	#<
	#? Delete an apparently duplicated atom.
	#  Looks for atoms that are close in space (within 0.001 A). Then
	#  checks to make sure they are the same element. No further checks
	#  are done.
	#  Hydrogen atoms are deleted where these are attached to deleted
	#  heavy atoms, unless a flag is sent to keep them. This flag has
	#  no effect where the hydrogen atoms are deleted by reason of being
	#  duplicates in their own right.
	#; Requires: Molecule, tolerance, optional flag to keep hydrogens
	#  when these are attached to deleted atoms
	#; Returns: Number of atoms deleted
	#>
	
	my $mol = $_[0];
	my $tol = $_[1] || 0.001;
	my $keeph = $_[2];
	
	my $delcount;
	my $dellist;
	
	@$dellist = ();
	
	for (my $i = 0; $i < $#{$mol->{ATOMS}}; ++$i) {
		
		my $atom1 = $mol->{ATOMS}[$i];
		my $j;
		
		for ($j = $i+1; $j <= $#{$mol->{ATOMS}}; ++$j) {
			
			my $atom2 = $mol->{ATOMS}[$j];
			my $connum;
			
			if (distance($atom1, $atom2) < $tol) {
				
				if ($atom2->{ELEMENT} eq $atom1->{ELEMENT}) {
				
					silico_msg('c', "Deleting atom ".($j+1)." (\"$atom2->{NAME}\")\n");
				
					push @$dellist, $j;
					++$delcount;
					
					if (!$keeph) {
						
						foreach my $connum (@{$atom2->{CONNECT}}) {
							
							my $con = $mol->{ATOMS}[$connum];
							
							if ($con->{ELEMENT} eq 'H') {
								
								silico_msg('c', "Deleting atom ".($connum+1)." (\"$con->{NAME}\")\n");
								
								push @$dellist, $connum;
								++$delcount;
							}
						}
					}
					
				} else {
					silico_msg('w', "Atoms $atom1->{NUM} (\"$atom1->{NAME}\") and $atom2->{NUM} (\"$atom2->{NAME}\") are in the same position!\n",
							"They are different elements and neither will be deleted.\n");
				}
			}
		}
	}
	
	# Only do any actual deletion if dummy atoms have been found.
	if ($delcount > 0) {
		molecule_delete_atoms($mol, $dellist);
		silico_msg('c', "$delcount dummy atoms deleted\n");
		
		# Pack and renumber remaining atoms
		molecule_pack($mol);	
		$mol->{NUMBONDS} = molecule_count_numbonds($mol);
	}
	
	return $delcount;
}


#################
# ATOM QUERYING #
#################


sub bonded {

	#<
	#? Count the number of bonds to an atom and the number of bonded C, H, O, N atoms
	#. Dummy atoms are not included
	#; Requires: atom, molecule
	#; Returns: hash containing the element (key) and element count (value) of all connected atoms.
	#  The hash also contains the key value NOT_H. The values are zeroed for C, H, N, O, P and NOT_H.
	#  Other elements may return 'undef'
	#>

	my $atom = $_[0];
	my $mol = $_[1];

	my $h;
	%$h = ();
	
	foreach (qw(C H N O P NOT_H)) {
		$h->{$_} = 0;
	}

	foreach my $con (connected($atom, $mol)) {

		++$h->{$con->{ELEMENT}};
		++$h->{NOT_H} if $con->{ELEMENT} ne 'H';
	}

	return $h;
}

sub find_nearest_atoms {

	#<
	#? Find the nearest atoms to a nominated atom and return list of atoms
	#; Requires: atom, molecule, cutoff (default 5 angstroms)
	#; Returns: pointer to an array of atoms
	#>
	
	my $atom1 = $_[0];
	my $mol = $_[1];
	my $cutoff = $_[2] || 5;

	my $cutoff2 = $cutoff*$cutoff;

	my $list;
	@$list = ();

	my $x1 = $atom1->{X};
	my $y1 = $atom1->{Y};
	my $z1 = $atom1->{Z};

	foreach my $atom2 (@{$mol->{ATOMS}}) {

		next if $atom1 == $atom2; # Ignore self

		my $a = $x1-$atom2->{X};
		next if $a > $cutoff2;

		my $b = $y1-$atom2->{Y};
		next if $b > $cutoff2;

		my $c = $z1-$atom2->{Z};
		next if $c > $cutoff2;

		my $d2 = $a**2 + $b**2 +$c**2;

		if ($d2 < $cutoff2) {

			push @$list, $atom2;
		}
	}
	return @$list;
}



sub valence {

	#<
	#? Return the valence of an atom
	#. Note that THIS ROUTINE CAN RETURN NON-INTEGRAL VALUES for bonds involving
	#  aromatic bond types (eg carboxylates or guanidines represented using aromatic bonds).
	#  The valence is corrected for aromatic bridgehead atoms which return a valence
	#  of 4 rather than 4.5 (which would result from three aromatic bonds).
	#. The minimum possible valence is returned for poorly defined bond types (>=5)
	#; Requires: atom, molecule (optional but necessary to detect lone pairs)
	#; Returns: valence
	#>

	my $atom = $_[0];
	my $mol = $_[1]; # Optional but necessary to detect lone pairs
	
	my $aromcount = 0;
	my $i = 0;
	my $val = 0;
	
	foreach my $bo (@{$atom->{BORDERS}}) {
	
		# Skip lone pairs and dummy atoms
		if ($mol) {
			next if $mol->{ATOMS}[$atom->{CONNECT}[$i]]{ELEMENT_NUM} == 0;
		}
	
		if ($bo <= 3) {
			$val += $bo;
			next;
		}
		if ($bo == 4 || $bo == 7) {
			$val += 1.5;
			++$aromcount;
			next;
		}
		$val += 1;
		++$i;
	}
	
	# Return valence of 4 if we have three aromatic bonds
	# eg bridgehead carbon
	
	# This causes incorrect functioning of valence
	# All 0.5 valences are converted to integral values!
	#if (((2*$val) & 1) == 1) {
	#	$val -= 0.5;
	#}
	
	$val = 4 if $atom->{ELEMENT} eq 'C' && $aromcount == 3;
	
	return $val;
}

sub atom_guess_element {

	#<
	#? Determine the  atom element using available information
	#; Requires: atom, 'quiet' flag, flag to 
	#  force preference for 'CHONSP' elements (ie
	#  H will be chosen over Hg)
	#; Returns:  element and element number. Also sets
	#  $atom->{ELEMENT} and $atom->{ELEMENT_NUM}
	#>

	my $atom = $_[0];
	# Assignments to local variables not used to speed up this
	# frequently called subroutine
	#my $mol = $_[1];
	#my $quiet = $_[2];
	#my $chonsp = $_[3];
	
	#
	# Use the element field if it was read
	#
	
	if ($atom->{ELEMENT}) {
		
		my $el = $atom->{ELEMENT};
		$el =~ s/[^A-Za-z]//g; # Remove any non-char
		$el = ucfirst lc $el;

		# Return values if we have a real element
		my $n = $Silico::Atomic_Elements{$el};
		if (defined $n) {

			$atom->{ELEMENT_NUM} = $n;
			$atom->{ELEMENT} = $el;
			return ($el, $n);
		}
		
		# Don't warn me about lone pairs
		if ($el eq 'Lp') {
			$atom->{ELEMENT_NUM} = 0;
			return ('LP',0);
		}
		
		# Otherwise warn and obtain element from atom name
		silico_msg('w', "Atom $atom->{NUM} (element field '$el') does not fit PDB format.\n",
				"Using atom name instead.\n") if ($el ne '' && !$_[2]);
	}
	
	#
	# Process atom and substructure names
	#
	
	# Assume that the element makes up the first one or two characters
	
	my $afield = $atom->{NAME} || '';
	$afield =~ s/[\W ]//g;		# Remove nonword characters and spaces
	$afield =~ s/^\d*//g;		# Remove leading digits
	$afield =~ s/\d.*//g;		# Remove trailing digits and everything following (reduces spurious recognition as heavier elements)
	$afield = "Du" if $afield eq '';
	$afield = substr ($afield,0,2) if length($afield) > 2;
	
	my $aname = ucfirst (lc $afield);
	
	# Substructure name (RESFIELD is used in read_mol2)
	my $sname = $atom->{SUBNAME} || $atom->{RESFIELD} || '';
	$sname =~ s/ //g;
	
	# Compare only first three characters of substructure name
	$sname = substr($sname."   ", 0, 3); # Compare only first three characters
	
	#
	# If residue name is a standard one (or the chonsp flag is set) then assume only CHONSP elements, LP or Q atom
	#
	
	if ($_[3] || $Silico::Amino_Nucleic_Acids{$sname}) {
	
		# Return element if the first character matches
		# element name
		
		if ($aname =~ /^H/) {
			$atom->{ELEMENT} = 'H';
			$atom->{ELEMENT_NUM} = 1;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /^C/) {
			$atom->{ELEMENT} = 'C';
			$atom->{ELEMENT_NUM} = 6;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /^N/) {
			$atom->{ELEMENT} = 'N';
			$atom->{ELEMENT_NUM} = 7;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /^O/) {
			$atom->{ELEMENT} = 'O';
			$atom->{ELEMENT_NUM} = 8;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /^S/) {
			$atom->{ELEMENT} = 'S';
			$atom->{ELEMENT_NUM} = 16;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /^P/) {
			$atom->{ELEMENT} = 'P';
			$atom->{ELEMENT_NUM} = 15;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		
		# Pseudo atoms
		if ($aname =~ /^Q/) {
			silico_msg('w', "Atom of name \"$atom->{NAME}\" found!\n",
					"Assuming pseudo atom (Du).\n",
					"Residue: \"$atom->{SUBNAME} $atom->{SUBID}\"\n") if !$_[2];
			$atom->{ELEMENT} = 'Du';
			$atom->{ELEMENT_NUM} = 0;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /^Lp/) {
			silico_msg('w', "Atom of name \"$atom->{NAME}\" found!\n",
					"Assuming lone pair (Du).\n",
					"Residue: \"$atom->{SUBNAME} $atom->{SUBID}\"\n") if !$_[2];
			$atom->{ELEMENT} = 'Du';
			$atom->{ELEMENT_NUM} = 0;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		
		# Last ditch try
		if ($aname =~ /C/) {
			$atom->{ELEMENT} = 'C';
			$atom->{ELEMENT_NUM} = 6;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /N/) {
			$atom->{ELEMENT} = 'N';
			$atom->{ELEMENT_NUM} = 7;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /O/) {
			$atom->{ELEMENT} = 'O';
			$atom->{ELEMENT_NUM} = 8;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /S/) {
			$atom->{ELEMENT} = 'S';
			$atom->{ELEMENT_NUM} = 16;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		if ($aname =~ /H/) {
			$atom->{ELEMENT} = 'H';
			$atom->{ELEMENT_NUM} = 1;
			return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
		}
		silico_msg('w', "Found standard amino acid \"$sname\" with nonstandard atom name \"$aname\"\n");
	}
	
	#
	# Determine atomic element using the first two characters of atom name
	#
	
	# Assume that CA is C-alpha and not calcium if residue name does not start with CA
	# And other elements
	if ($aname eq 'Ca') { 
		if ($sname =~ /^CA/i) {
			$atom->{ELEMENT} = 'Ca';
			$atom->{ELEMENT_NUM} = 20;
			return ('Ca',20);
		} else {
			$atom->{ELEMENT} = 'C';
			$atom->{ELEMENT_NUM} = 6;
			return ('C',6);
		}
	}
	if ($aname eq 'Cd') { 
		if ($sname =~ /^CD/i) {
			$atom->{ELEMENT} = 'Cd';
			$atom->{ELEMENT_NUM} = 48;
			return ('Ca',20);
		} else {
			$atom->{ELEMENT} = 'C';
			$atom->{ELEMENT_NUM} = 6;
			return ('C',6);
		}
	}
	if ($aname eq 'Ce') { 
		if ($sname =~ /^CE/i) {
			$atom->{ELEMENT} = 'Ce';
			$atom->{ELEMENT_NUM} = 58;
			return ('Ce',58);
		} else {
			$atom->{ELEMENT} = 'C';
			$atom->{ELEMENT_NUM} = 6;
			return ('C',6);
		}
	}
	
	if ($aname eq 'Hg') { 
		if ($sname =~ /^HG/i) {
			$atom->{ELEMENT} = 'Hg';
			$atom->{ELEMENT_NUM} = 80;
			return ('Hg',60);
		} else {
			$atom->{ELEMENT} = 'H';
			$atom->{ELEMENT_NUM} = 1;
			return ('H',1);
		}
	}
	if ($aname eq 'He') { 
		if ($sname =~ /^HE/i) {
			$atom->{ELEMENT} = 'He';
			$atom->{ELEMENT_NUM} = 2;
			return ('He',2);
		} else {
			$atom->{ELEMENT} = 'H';
			$atom->{ELEMENT_NUM} = 1;
			return ('H',1);
		}
	}
	if ($aname eq 'Ho') { 
		if ($sname =~ /^HO/i) {
			$atom->{ELEMENT} = 'Ho';
			$atom->{ELEMENT_NUM} = 2;
			return ('Ho',67);
		} else {
			$atom->{ELEMENT} = 'H';
			$atom->{ELEMENT_NUM} = 1;
			return ('H',1);
		}
	}
	# Assume that NA is Nitrogen unless residue name is NA or SOD
	# or if molecule has connectivity and atom has no bonds
	if ($aname eq 'Na') {
		
		if ($sname =~ /^NA/i  || $sname =~ /^SOD/i) {
			$atom->{ELEMENT} = 'Na';
			$atom->{ELEMENT_NUM} = 11;
			return ('Na',11);
		} elsif (defined $atom->{CONNECT} && $atom->{CONNECT} == 0) {
			$atom->{ELEMENT} = 'Na';
			$atom->{ELEMENT_NUM} = 11;
			return ('Na',11);
		} else {
			$atom->{ELEMENT} = 'N';
			$atom->{ELEMENT_NUM} = 7;
			return ('N',7);
		}
	}

	if ($aname eq 'Ne') { 
		if ($sname =~ /^ND/i) {
			$atom->{ELEMENT} = 'Ne';
			$atom->{ELEMENT_NUM} = 10;
			return ('Ne',60);
		} else {
			$atom->{ELEMENT} = 'N';
			$atom->{ELEMENT_NUM} = 7;
			return ('N',7);
		}
	}
	if ($aname eq 'Nd') { 
		if ($sname =~ /^ND/i) {
			$atom->{ELEMENT} = 'Nd';
			$atom->{ELEMENT_NUM} = 60;
			return ('Nd',60);
		} else {
			$atom->{ELEMENT} = 'N';
			$atom->{ELEMENT_NUM} = 7;
			return ('N',7);
		}
	}
	

	# Match two characters of aname
	if (defined $Silico::Atomic_Elements{$aname}) {
		$atom->{ELEMENT} = $aname;
		$atom->{ELEMENT_NUM} = $Silico::Atomic_Elements{$aname};
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	
	# Maybe the atomic element is only the 1st character
	$aname = substr ($afield,0,1);
	$aname = ucfirst lc $aname;
	if (defined $Silico::Atomic_Elements{$aname}) {
		$atom->{ELEMENT} = $aname;
		$atom->{ELEMENT_NUM} = $Silico::Atomic_Elements{$aname};
		return ($aname, $Silico::Atomic_Elements{$aname});
	}
	
	# Pseudo atoms
	if ($aname =~ /^Q/ || $aname =~ /Du/ || $aname =~ /Lp/) {
		silico_msg('n', "Found atom of name \"$atom->{NAME}\". Assuming pseudo atom (Du).\n") if $aname =~ /^Q/;
		silico_msg('n', "Found atom of name \"$atom->{NAME}\". Assuming dummy atom (Du).\n") if $aname =~ /Du/;
		silico_msg('n', "Found atom of name \"$atom->{NAME}\". Assuming lone pair (Du).\n") if $aname =~ /Lp/;
		$atom->{ELEMENT} = 'Du';
		$atom->{ELEMENT_NUM} = 0;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	
	# 
	# More desperate guesses
	#
	
	# Maybe the atomic element is only the 2nd character
	$aname = substr ($afield,1,1);
	$aname = ucfirst lc $aname;
	if (defined $Silico::Atomic_Elements{$aname}) {
		$atom->{ELEMENT} = $aname;
		$atom->{ELEMENT_NUM} = $Silico::Atomic_Elements{$aname};
		my $msg = "Guessed atom element $atom->{ELEMENT} from atom name $aname\n";
		++$_[1]->{WARN}{$msg};
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	
	$aname = $atom->{NAME};
	$aname =~ s/ //g;
	
	# Make last ditch attempts for CHONSP and some other common elements
	# First character of name
	if ($aname =~ /^H/) {
		$atom->{ELEMENT} = 'H';
		$atom->{ELEMENT_NUM} = 1;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ /^C/) {
		$atom->{ELEMENT} = 'C';
		$atom->{ELEMENT_NUM} = 6;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ /^N/) {
		$atom->{ELEMENT} = 'N';
		$atom->{ELEMENT_NUM} = 7;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ /^O/) {
		$atom->{ELEMENT} = 'O';
		$atom->{ELEMENT_NUM} = 8;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ /^S/) {
		$atom->{ELEMENT} = 'S';
		$atom->{ELEMENT_NUM} = 16;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ /^P/) {
		$atom->{ELEMENT} = 'P';
		$atom->{ELEMENT_NUM} = 15;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}

	# Make last ditch attempts for CHONSP. Any character in name.
	if ($aname =~ 'H') {
		$atom->{ELEMENT} = 'H';
		$atom->{ELEMENT_NUM} = 1;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ 'C') {
		$atom->{ELEMENT} = 'C';
		$atom->{ELEMENT_NUM} = 6;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ 'N') {
		$atom->{ELEMENT} = 'N';
		$atom->{ELEMENT_NUM} = 7;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ 'O') {
		$atom->{ELEMENT} = 'O';
		$atom->{ELEMENT_NUM} = 8;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ 'S') {
		$atom->{ELEMENT} = 'S';
		$atom->{ELEMENT_NUM} = 16;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ 'P') {
		$atom->{ELEMENT} = 'P';
		$atom->{ELEMENT_NUM} = 15;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	
	# Common ions
	if ($aname =~ 'NA') {
		$atom->{ELEMENT} = 'Na';
		$atom->{ELEMENT_NUM} = 11;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}
	if ($aname =~ 'MG') {
		$atom->{ELEMENT} = 'Mg';
		$atom->{ELEMENT_NUM} = 12;
		return ($atom->{ELEMENT}, $atom->{ELEMENT_NUM});
	}

	# Oh well.  Call it dummy
	silico_msg('w', "Could not derive element of atom $atom->{NUM}!\n",
			"Atom name: \"$aname\"\n",
			"Setting atom type to Dummy Atom (Du).\n");
	$atom->{ELEMENT} = 'Du';
	$atom->{ELEMENT_NUM} = 0;
	return ("Du", 0);
}


#####################
# ATOM MODIFICATION #
#####################


sub change_bondorder {

	#<
	#? Change bond order between 'atom' and any connected atom
	#  of the specified element
	#. If order = 10 then this routine alternates single,double
	#  (for carboxylate groups)
	#. Note: calls bond_modify_order and requires that atom->{NUM} is correct
	#; Requires: atom, desired element of connected atom, molecule, bondorder
	#; Returns: nothing
	#>

	my $atom = $_[0];
	my $element = $_[1];
	my $mol = $_[2];
	my $order = $_[3];
	
	my $flip = 0;
	
	if ($order == 10) {
		$order = 1;
		$flip = 1;
	}
	
	foreach my $connum (@{$atom->{CONNECT}}) {
				
		my $conatom = $mol->{ATOMS}[$connum];
				
		if ($conatom->{ELEMENT} eq $element) {
		
			if ($flip) {
				--$order;
				$order = !$order;
				++$order;
			}
			
			bond_modify_order_check($mol, $atom, $conatom, $order);
		}
	}

	return;
}

sub make_atom_names_unique {

	#<
	#? Make all atom names in a molecule unique
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	
	my $nhash;
	my $ehash;
	my $count = 0;
	
	foreach my $atom (atoms($mol)) {
	
		if ($nhash->{$atom->{NAME}||'X'}) {
			
			while (1) {
			
				my $el = $atom->{ELEMENT}||'Du';
				++$ehash->{$el};
				
				my $name = $el.$ehash->{$el};
				
				if (!$nhash->{$name}) {
				
					$atom->{NAME} = $name;
					++$count;
					#print "n $name\n";
					last;
				}
			}
		}
		
		++$nhash->{$atom->{NAME}};
	}
	
	#mol_warn($mol,"Changed $count duplicate atom names");
	silico_msg('c', "Changed $count duplicate atom names");
}


#################
# ATOM PRINTING #
#################

sub atominfo {

	#<
        #? Convert atom info into a string for printing & catch errors
        #; Requires: atom
        #; Returns: string
        #>

	my $atom = $_[0];
	
	my $s;
	my @w = qw(8 2 4 5 5 2);
	
	my $i = 0;
	foreach (qw (NUM ELEMENT SUBNAME SUBID NAME CHAIN)) {

		$s .= ucfirst(lc($_));
		$s .= ": ";
		my $v = ($atom->{$_} || '-');
		$v =~ s/ //g;
		$s .= sprintf "%-$w[$i]s ", $v;
		++$i;
	}
	
 	
	return $s;
}

sub print_con {

	my $atom = $_[0];
	my $con = $_[1];
	
	print "\t";
	print_current($atom);
	print "\tConnected to ";
	print_current($con);
}

sub print_current {

	my $atom = $_[0];
	
	print "Atom: $atom->{NUM} (name \"$atom->{NAME}\"";
	print " Res: '$atom->{SUBNAME}'" if defined $atom->{SUBNAME};
	print " Subid: '$atom->{SUBID}'" if defined $atom->{SUBID};
	print "\n";
}

sub print_two {

	my $atom1 = $_[0];
	my $atom2 = $_[1];
	
	print "\tAtom1 $atom1->{NAME} Res $atom1->{SUNAME} $atom1->{SUBID}\n";
	print "\tAtom2 $atom2->{NAME} Res $atom2->{SUNAME} $atom2->{SUBID}\n";

}

sub compare_atom_order {

	#<
	#? Compare two atoms and see if they are in the correct
	#  order for output to mol2, merck files etc.
	#. This routine should be the same as the one used in
	#  pdb_molecule_sort_atoms
	#>
	
	my $a = $_[0];
	my $b = $_[1];
		
	# Sort by SEGID
	if (defined $a->{SEGID} || defined $b->{SEGID}) {
	
		return -1 if !defined $a->{SEGID};
		return 1 if !defined $b->{SEGID};
		return -1 if $a->{SEGID} lt $b->{SEGID};
		return 1 if $a->{SEGID} gt $b->{SEGID};
	}
			
	# Sort by chain
	if (defined $a->{CHAIN} || defined $b->{CHAIN}) {
	
		return -1 if !defined $a->{CHAIN};
		return 1 if !defined $b->{CHAIN};
		return -1 if $a->{CHAIN} lt $b->{CHAIN};
		return 1 if $a->{CHAIN} gt $b->{CHAIN};
	}

	# Sort by substructure number
	if (defined $a->{SUBID} || defined $b->{SUBID}) {
		return -1 if $a->{SUBID} < $b->{SUBID};
		return 1 if $a->{SUBID} > $b->{SUBID};
	}

	# Sort by atom number
	$a->{NUM} <=> $b->{NUM};
}


##########################
# SUPPORTING INFORMATION #
##########################


sub element_number2symbol {

	#<
	#? Convert an element number into the corresponding atomic symbol.
	#; Requires: element number.
	#; Returns: element number or 0 Du.
	#>

	my $number = ($_[0]);
	my $quiet = $_[1];

	if ($number < $#Silico::Atomic_Elements && $number >= 0) {
		return $Silico::Atomic_Elements[$number];
	}
	#silico_msg('w', "Could not determine the atomic symbol for element number $number!\n") if !$quiet;
	return 0;
}

sub element_symbol2number {

	#<
	#? Convert an element symbol into the corresponding atomic number.
	#; Requires: element symbol.
	#; Returns: element number or 0 if it fails.
	#>
	
	my $symbol = ucfirst lc ($_[0]);
	my $val;
		
	$val = $Silico::Atomic_Elements{$symbol};

	return $val if defined $val;
	
	return 0 if $symbol eq 'Lp'; # Lone pair
	return 0 if $symbol eq 'A';  # Query atom
	return 0 if $symbol eq 'Q';  # Query atom
	
	silico_msg('w', "Could not determine the atomic number for element '$symbol'!\n");
	return 0;
}


sub ensemble_insight_rename {

	#<
	#? Generate insight-safe names for an ensemble.
	#  Calls: insight_rename.
	#; Requires: ensemble.
	#>

	my $ligands = $_[0];
	
	foreach my $mol (@$ligands) {

		$mol->{NAME} = insight_rename($mol->{NAME});
	}
}

sub insight_rename {

	#<
	#? Generate insight-safe names
	#
	#. Insight has a number of reqirements for molecule names. For example
	#  they can not contain spaces, punctuation, start with an underscore
	#  or a digit  and are of limited length.
	#  This routine generates unique, Insight-safe names by using
	#  the hash %insight_namecount to store previously used names.
	#; Requires: name, pointer to a hash to store names in.
	#; Returns: modified and unique name.
	#>

	my $name = $_[0];
	my $insight_namecount = $_[1];

	# Noname!
	$name = "Mol" if !$name;

	# Remove any non-word characters (alphanumeric + '_');
	my @f = split /\W/, $name;
	
	my $i = 0;
	foreach (@f) {
		$f[$i] = ucfirst(lc $_);
		++$i;
	}
	
	$name = join '', @f;

	# Mask leading digits or underscore
	$name =~ s/^(?=\d)/M/;
	$name =~ s/^_/M_/;
	
	# Make sure that we have _some_ name left
	$name = "Mol" if !$name;

	# Check that the name is not longer than 16 characters
	# to give a total of *19* chars with added number
	# (Despite the fact that Insight thinks it has a limit of 20 chars)
	# If it is truncate it to first 12 and last 3 chars
	
	if (length ($name) > 16 ) {
		
		$name = substr ($name, 0, 12)."_".substr ($name, -3, 3);
	}
	
	# Convert to upper case to avoid insight bugs
	#$name = uc $name;

	# Remember the names we have used
	++$insight_namecount->{$name};
	
	# Add a counter if necessary
	if ($insight_namecount->{$name} > 1) {
		
		my $count;
		if ($insight_namecount->{$name} < 100) {
			$count = "_".sprintf "%02d", $insight_namecount->{$name};
		} else {
			$count = sprintf "%04d", $insight_namecount->{$name};
		}

		$name = $name.$count;
	}

	return $name;
}

return 1;
