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
#! silico_rbox.pm
#? Subroutines for random_box and bilayer_builder
#. $Revision: 1.1.2.17 $
#>

use strict;
package Silico;

sub check_number_components {

	#<
	#? Check that number of components is set for each component. If it is not, then ask
	#>

	my $comp = $_[0];
	my $num_components = $_[1];
	my $nummols = $_[2] || 0;

	my $i = 0;
	foreach my $mol (@$comp) {

		if (!defined $num_components->[$i]) {

			my $ans = prompt("Molecule ".($i+1)." ($mol->{NAME}): Number of molecules to add (".($nummols||0)."): ");

			while (1) {
				next if !defined $ans;
				if ($ans !~ /^\d+$|\./) {
					chomp;
					print "Using default value (";
					if ($nummols) {
						print "$nummols";
					} else {
						print "0";
					}

					print ")\n";
					$ans = $nummols || 0;
				}
				last if $ans;
			}

			$num_components->[$i] = $ans;
		}
		++$i;
	}

	return $num_components;
}

sub make_molecule_copies {

        #<
        #? Generate multiple copies of molecules to be packed
        #; Requires: ensemble of molecules to be packed, number to pack, starting chain, starting segid, number of bonds to randomise
	#; Returns: molecules
        #>

        my $comp = $_[0];
        my $num_components = $_[1] || 0;
        my $chain = $_[2] || 'A';
        my $segid = $_[3] || 'A';
	my $random = $_[4];
	my $oformat = $_[5] || '';

        my $atomcount = 0;
        my $molcount = 0;
	my $mollist;
        my $subcount = 1;

	my $maxsub;
        $maxsub = 9999 if $oformat eq 'pdb' || get_lflag('pdb-compatible');
        $maxsub = 99999 if $oformat eq 'gro';

        print "Using output format: $oformat.  Maximum residue number: $maxsub\n" if $maxsub;

	my $i = 0;
        foreach my $mol (@$comp) {

		mol_rotate_random_bond($mol, $random) if $random;

                for (my $j = 0; $j < int($num_components->[$i]); ++$j) {

                        if ($maxsub && $subcount > $maxsub) {
                               
			        $subcount = 1;
                                pdb_increment_chain($chain);
                                silico_msg ('c', "$molcount: Reached $maxsub residues. Incrementing chain to $chain\n");
                        }
			
                        my $newmol = copy_minimal_molecule($mol);
                        push @$mollist, $newmol;

                        $newmol->{NUM} = $molcount;
			
			#print "subcount $subcount\n";

                        foreach my $atom (atoms($newmol)) {
                                ++$atomcount;
                                $atom->{NUM} = $atomcount;
                                $atom->{SUBID} = $subcount;
                                $atom->{SEGID} = $segid;
                                $atom->{CHAIN} = $chain;
                        }

                        ++$subcount;
                        ++$molcount;
                }

                ++$i;
        }

        return $mollist;
}


sub mol_rotate_random_bond {

        my $mol = $_[0];
        my $numbonds = $_[1] || 1;

	my $best_overlaps;
        my $num;

        use vars qw($Pi);

        for $num (1..$numbonds) {

                $best_overlaps = mol_check_atom_overlap($mol, 1);
                record_best_coords($mol);

                my @bond = select_rotate_bond($mol);
                my ($rlist1, $rlist2) = mol_partition_on_bond2($mol, @bond, $mol->{ATOMS}[0]);

                my $count = 0;

                while (1) {

                        ++$count;

                        # Choose a random angle up to 360 degrees
                        my $angle = rand(2 * $Silico::Pi);

                        # Perform a rotation by this angle
                        rotate_torsion_angle2($mol, @bond, $rlist2, $angle);

                        my $overlaps = mol_check_atom_overlap($mol, 1);

                        if ($overlaps <= $best_overlaps) {

                                $best_overlaps = $overlaps;
                                last;
                        }

                        if ($count == 10) {

                                revert_to_best_coords($mol);
                                last;
                        }
                }
        }
}

sub copy_minimal_molecule {

	#<
	#? Make a very minimal copy of a molecule
	#. Note to save space, records which are not modified by the program are duplicates and NOT copies
	#>

	my $mol = $_[0];

	my $newmol;

	my @list_copy = qw (X Y Z CONNECT);
	my @list_dup  = qw (BORDERS ELEMENT ELEMENT_NUM NAME SUBNAME);

	# Records which are NOT modified (duplicates)
	# Records which are modified (copies)
	my $i = 0;
	foreach my $atom (atoms($mol)) {
	
		foreach (@list_copy) {
			next if !defined $mol->{ATOMS}[$i]{$_};
			$newmol->{ATOMS}[$i]{$_} = deep_copy($mol->{ATOMS}[$i]{$_});
		}
		foreach (@list_dup) {
			next if !defined $mol->{ATOMS}[$i]{$_};
			$newmol->{ATOMS}[$i]{$_} = $mol->{ATOMS}[$i]{$_};
		}
		++$i;
	}

	$newmol->{NUMATOMS} = $mol->{NUMATOMS};
	$newmol->{NUMBONDS} = $mol->{NUMBONDS};
	return $newmol;
}

sub ensemble_consolidate_nocopy {

	#<
	#? Combines all molecules in an ensemble into a single molecule
	#  This is the same as the version in silico_molecules but without the
	#  deep copy - which is much faster.
	#; Requires: Ensemble
	#; Returns: New ensemble with all atoms in the molecule $newens->[0]
	#>
	
	my $ens = $_[0];
	
	my $count = 0;
	my $atomcount = 0;
	my $newmol;
	my $numbonds = 0;

	foreach my $mol (@$ens) {
	
		next if !defined $mol;
	
		++$count;
		$numbonds += $mol->{NUMBONDS};

		# Copy molecule data
		if ($count == 1) {
			$newmol->{ATOMS} = $mol->{ATOMS};
			$newmol->{CELL} = $mol->{CELL};
			$newmol->{CELL_ALPHA} = $mol->{CELL_ALPHA};
			$newmol->{CELL_BETA}  = $mol->{CELL_BETA};
			$newmol->{CELL_GAMMA} = $mol->{CELL_GAMMA};
			$atomcount = $mol->{NUMATOMS};
			next;
		}
		
		# Copy atoms renumber atoms from subsequent molecules
		molecule_renumber2($mol, $atomcount+1);
		foreach my $atom (@{$mol->{ATOMS}}) {
			push @{$newmol->{ATOMS}},  $atom;
			++$atomcount;
		}
	}
	
	$newmol->{NUMATOMS} = $atomcount;
	$newmol->{NUMBONDS} = $numbonds;
	
	# Give the ensemble the name of the first molecule in original ensemble
	$newmol->{NAME} = $ens->[0]{NAME};
	$newmol->{SOURCE_FILE_NAME} = $ens->[0]{SOURCE_FILE_NAME};
	$newmol->{CELL} = $ens->[0]{CELL};
	$newmol->{BINS} = $ens->[0]{BINS};
	$newmol->{NUMBINS} = $ens->[0]{NUMBINS};
	$newmol->{BINSIZE} = $ens->[0]{BINSIZE};
	
	return $newmol;
}

sub pack_molecules {

	#<
        #? Pack molecules in a random or semirandom (e.g. bilayer) way
	#; Requires: array of molecules to pack, molecule to embed, molecule cell
	#; Also requires subroutine 'randomise_and_test_molecule_position'
        #>

	my $pack = $_[0];		# New molecules to pack
	my $embed = $_[1];		# Existing molecule
	my $cell = $_[2];		# Unit cell
	
	my $scale_change = get_lflag('incr') || 10000; # Change scaling when reaching this number of attempts
	my $scale_init = get_flag('init-scale-factor','l') || 1; 	# Initial scale factor for atom-atom packing distance
	my $scale_min = get_flag('min-scale-factor','l') || 0.7;	# Minimum scale factor for atom-atom packing distance

	my $freq = 10;
	my $maxclash = 1;
	my $nummols = 0;
	my $prev_numatoms = 0;
	my $scale = $scale_init;
	my $scale_incr = 0.01;
	my $total_tries = 0;
	
	# Get molecule to be embedded or take the first molecule from
	# the list and use it
	my $exist;
	if ($embed) {
   		
		# Mark each atom in embedded molecule so that vDW scaling 
		# is not used.
		foreach my $atom (atoms($embed)) {
			$atom->{DONT_SCALE} = 1;
		}
		$exist = $embed;
	} else {
   		$exist = shift @$pack;
		++$nummols;
		++$total_tries;
		randomise_and_test_molecule_position($exist, $cell);
		mol_move_into_cell($exist, $cell);
	} 

	$exist->{CELL} = $cell; # mol->{CELL} is used by multiple routines
	my ($binsize, $numbins) = calc_binsize($exist);
	molecule_bin_atoms($exist, undef, $binsize, $numbins);
	
	my $packed;
	for (my $i=0; $i <= $#{$pack}; ++$i) {

		my $mol = $pack->[$i];
		my $count_insert = 0;

		if ($mol->{NUMATOMS} < $prev_numatoms && $scale != $scale_init) {
			silico_msg('c', "Packing smaller molecule. Resetting scale factor to $scale_init\n");
			$scale = $scale_init;
		}
		
		while (1) {
			
			++$count_insert;
			++$total_tries;
	
			randomise_and_test_molecule_position($mol, $cell);
			
			my ($touch, undef) = clash_periodic_binned($mol, $exist, $maxclash, $scale, 1.5, get_lflag("ignore-hydrogens")); 
		
			if ($total_tries % $freq == 0) {
				$freq *= 5 if $total_tries >= $freq * 10;
				printit($mol, $total_tries, $nummols, $count_insert, 0, 0);
				write_pdb($packed, "random_box_intermediate.pdb") if get_sflag('write');
			}

			if ($touch == 0) {
				++$nummols;
				# Print final line to screen
				last;
			}

			 if (($count_insert == $scale_change) && ($scale > $scale_min)) {
				$scale = sprintf "%.2f", $scale - $scale_incr;
				printit($mol, $total_tries, $nummols, $count_insert, 1, $scale);
				$count_insert = 0;
			} 
		}
		
		mol_move_into_cell($mol, $cell);

		# Add molecule to existing molecules
		push @$packed, $mol;
		my $ens;
		@$ens = ($exist, $mol);

		# Consolidate ensemble without using deep_copy (faster and not required)
		$exist = ensemble_consolidate_nocopy($ens);
	
		# Bin the atoms that we have just addded
		molecule_bin_atoms($exist, $mol->{ATOMS}, $binsize, $numbins); 
		
		$prev_numatoms = $mol->{NUMATOMS};
	}

	printf ("total_tries: %-10d mols_added: %-10d\n", $total_tries, $nummols);
	
	return $exist;
}



sub printit {

	#<
        #? Print routine for pack_molecules
	#>

	my ($mol, $total_tries, $nummols, $count_insert, $reduceflag, $scale) = @_;

	my $name = $mol->{ATOMS}[0]{SUBNAME} || 'Mol';
	printf ("total_tries: %-10d mols_added: %-10d name: %-6s mol_tries: %-6d", $total_tries, $nummols, $name, $count_insert);
	print " Reducing scale factor to $scale" if $reduceflag;			
	print "\n";
}

sub orient_lipid {

	my $mol = $_[0];

	my @vec = lipophilicity_dipole($mol);

	foreach my $atom (atoms($mol)) {

		# Translate to the origin of the lipophilicity_dipole
		my $i = 0;
		foreach (qw(X Y Z)) {
			$atom->{$_} -= $vec[$i];
			++$i;
		}
	}

	calc_rotation_vector_onto_vector()

}


sub lipophilicity_dipole {
	
	#<
        #? Calculate the polar vector of a molecule 
	#  along the Y axis
	#: Requires: molecule
	#>

	my $mol = $_[0];
	
	my $polar;
	my $nonpolar;
	my $pcount = 0;
	my $ncount = 0;

	my $scale_ion_and_hbond_donors = 2;
	
	molecule_formal_charge($mol);
	mol_calc_hcount ($mol);
	
	foreach my $atom (atoms($mol)) {
	
		my $e = $atom->{ELEMENT};
		
		if ($atom->{FORMAL_CHARGE}) {
		
			# Charged atoms are worth more
			$polar->{X} += $atom->{X} * $scale_ion_and_hbond_donors;
			$polar->{Y} += $atom->{Y} * $scale_ion_and_hbond_donors;
			$polar->{Z} += $atom->{Z} * $scale_ion_and_hbond_donors;
			$pcount += 5;
			
		} elsif ( $e eq 'O' || $e eq 'N') {
		
			if ($atom->{HCOUNT}) {
			
				# OH and NH atoms are worth more
				$polar->[0] += $atom->{X} * $scale_ion_and_hbond_donors;
				$polar->[1] += $atom->{Y} * $scale_ion_and_hbond_donors;
				$polar->[2] += $atom->{Z} * $scale_ion_and_hbond_donors;
				$pcount += 5;

			} else {
			
				# Atoms without attached H are worth 1
				$polar->[0] += $atom->{X};
				$polar->[1] += $atom->{Y};
				$polar->[2] += $atom->{Z};
				$pcount += 1;
			}
			
		} else {
		
			# Carbons etc are worth 1
			$nonpolar->[0] += $atom->{X};
			$nonpolar->[1] += $atom->{Y};
			$nonpolar->[2] += $atom->{Z};
			$ncount += 1;
		}
	}
	
	if ($pcount == 0 || $ncount == 0) {
		silico_msg('w', "Can not calculate polar vector for molecule\n");
		return 0;
	}
	
	my $fakemol = create_molecule("vector");
	my $fa = $fakemol->{ATOMS};
	$fa->[0]{ELEMENT} = "C";
	$fa->[1]{ELEMENT} = "P";
	$fa->{NUMATOMS} =2;

	my $i = 0;
	foreach ('X', 'Y', 'Z') {
		$polar->{$_} /= $pcount;
		$nonpolar->{$_} /= $ncount;
		$fa->[0]{$_} = $nonpolar->{$i};
		$fa->[1]{$_} = $polar->{$i};
		++$i;
	}
	
	#write_mol_any($fakemol);
	
	return @$nonpolar, @$polar;
}



#
# Packing routines
#

sub test_all_shapes { 

	#<
        #? Test whether a molecule complies with a set of geometric restrictions provided by flags
	#; Requires: molecule, cell
	#; Returns: true or false
        #>

	my $mol = $_[0];
        my $cell = $_[1];

        my $inner = ($Silico::sphere||0)-($Silico::thick||0);

	my $p = $Silico::plane && test_plane($mol, $cell, $Silico::invert);
	my $s = $Silico::sphere && test_shell($mol, $cell, $inner, $Silico::sphere, $Silico::invert);
	my $c = $Silico::cylinder && test_cylinder($mol, $cell, $iSilico::invert);
	my $h = $Silico::holes && test_cylinder_arb($mol, $cell, $Silico::invert);
		
	my $shape_val =  $p || $s || $c;
	
	++$Silico::count;
	if (0 && $Silico::count < 100) {
	
		$p ||= 0;
		$s ||= 0;
		$c ||= 0;
		$h ||= 0;
		$shape_val ||= 0;
		
		#print "p $p s $s c $c h $h shape_val $shape_val\n";
	}

	if ($Silico::invert) {
		return 1 if $shape_val || $h;
	} else {
		return 1 if ($shape_val || (!$Silico::plane && !$Silico::sphere && !$Silico::cylinder)) && ($h || !$Silico::holes);
	}

	return 0;
}


sub test_plane {

	#<
	#? Test if molecule lies completely within (or completely outside) a (set of) plane(s)
	#. Plane centres are set by by Silico:origin_points
	#. Plane direction default xz or xy by flag
	#. If invert NOT set
	#. Returns 1 if molecule lies completely inside any one plane.
	#. If invert IS set
	#. Returns 1 if molecule lies completely outside ALL planes
	#; Requires: molecule, cell, invert flag
	#; Returns: true or false
	#>

	my $mol = $_[0];
	my $cell = $_[1];
	my $invert = $_[2];

	my $o = $Silico::origin_points;
	my $plane = get_flag('plane-size', 'l'); # Width of plane
	my $xy = get_lflag('xy');
	
	my $cy = $cell->[1];
	my $cz = $cell->[2];
	my $cy2 = $cy/2;
	my $cz2 = $cz/2;
	
	# Loop over planes
	ORIGIN: for (my $i = 0; $i < $#{$o}; $i = $i + 3) {
	
		my $dum; # 
		#$dum->{X} = $cell->[0]*$o->[$i];
		$dum->{Y} = $cell->[1]*$o->[$i+1];
		$dum->{Z} = $cell->[2]*$o->[$i+2];

		foreach my $atom (@{$mol->{ATOMS}}) {

			my $val;
			if ($xy) {
				$val = abs($atom->{Z} - $dum->{Z});
				# Make sure distance is not greater than half cell dimension
				$val = $val - floor(($val+$cz2)/$cz)*$cz;
				#die if $val > $cz2;
				
			} else {
				$val = abs($atom->{Y} - $dum->{Y});
				# Make sure distance is not greater than half cell dimension
				$val = $val - floor(($val+$cy2)/$cy)*$cy;
				#die if $val > $cy2;
			}
			
			my $inside = $val < $plane/2;
			
			if ($invert) {
				if ($inside) {
					return 0;
				} else {
					next;
				}
			} else {
				if ($inside) {
					next;
				} else {
					next ORIGIN;
				}
			}
		}

		return 1 if !$invert;
	}

	return 1 if $invert;
	return 0;
}

sub test_shell {

	#<
	#? Test if all atoms of a molecule molecule lies within a shell with radii inner & outer
	#  from the centre of the cell
	#; Requires: molecule, cell, inner radius, outer radius, invert flag
	#; Returns: true or false
	#<

	my $mol = $_[0];
	my $cell = $_[1];
	my $inner = $_[2] || 0;
	my $sphere = $_[3] || 9999999;
	my $invert = $_[4];

	my $o = $Silico::origin_points;

	ORIGIN: for (my $i=0; $i < $#{$o}; $i = $i + 3) {
	
		# Convert fractional coordinates to Cartesian
		my $dum; 
		$dum->{X} = $cell->[0]*$o->[$i];
		$dum->{Y} = $cell->[1]*$o->[$i+1];
		$dum->{Z} = $cell->[2]*$o->[$i+2];
		
		foreach my $atom (@{$mol->{ATOMS}}) {

			my $d = sqrt(distance_periodic_sq($atom, $dum, $cell));
			my $d2 = sqrt(distance_periodic_sq($atom, $dum, $cell));
			
			die if $d != $d2;

			my $inside = $d > $inner && $d < $sphere;

			if ($invert) {
				if ($inside) {
					return 0;
				} else {
					next;
				}
			} else {
				if ($inside) {
					next;
				} else {
					next ORIGIN;
				}
			}
		}

		return 1 if !$invert;
	}

	return 1 if $invert;
	return 0;
}

sub test_cylinder {

	#< 
	#? Test if all atoms of a molecule molecule lie inside a cylinder or cylindrical shell
	#; Requires: molecule, cell
	#; Returns: true or false
	#<

	my $mol = $_[0];
	my $cell = $_[1];
	my $invert = $_[2];
	
	my $cylinder = get_lflag('cylinder');
	my $cylinder2 = $cylinder**2;
	my $dum;
	my $inner2;
	my $length = get_lflag('length') || 0;
	my $o = $Silico::origin_points;
	my $thick = get_lflag('thick')||0;
	my ($x, $y, $z);
	my $xy = get_lflag('xy');

	$inner2 = (($cylinder||0)-($thick||$cylinder))**2;
	
	ORIGIN: for (my $i=0; $i < $#{$o}; $i = $i + 3) {
	
		$dum->{X} = $cell->[0]*$o->[$i];
		$dum->{Y} = $cell->[1]*$o->[$i+1];
		$dum->{Z} = $cell->[2]*$o->[$i+2];

		foreach my $atom (@{$mol->{ATOMS}}) {
	
			($x, $y, $z) = distance_periodic_xyz($atom, $dum, $cell);
			
			my $d2;
			
			if ($xy) {
				$d2 = $x**2 + $y**2;
			} else {
				$d2 = $x**2 + $z**2;
			}
			
			# Are we inside the cylinder?
			my $inside = ($d2 > $inner2 && $d2 < $cylinder2) && ($length == 0 || $z < $length/2);


			#
			# Test this code it differs from test_cylinder_arb!
			#
			if ($invert) {
			
				# If invert is set, all atoms must be outside all cylinders
				if ($inside) {
					return 0;
				} else {
					next;
				}
			} else {
			
				# If invert not set all atoms must be inside one cylinder
				if ($inside) {
					next;
				} else {
					next ORIGIN;
				}
			}
		}
		
		return 1 if !$invert; # Atoms were all inside one cylinder
	}
	
	return 1 if $invert; # All atoms were outside all cyinders if invert set
	return 0; # In no case were all atoms inside one cylinder if invert NOT set
}

sub test_cylinder_arb {

	#<
	#? Test if all atoms of a molecule are INSIDE or OUTSIDE a set of cylinders lying along arbitrary axes
	#. If invert is FALSE then function returns TRUE if molecule is OUTSIDE ALL cylinders
	#. If invert is TRUE then function returns TRUE if molecule is INSIDE ANY cylinder
	#. Note the sense of this test is the opposite of the others. Probably should be modified.
	#; Requires: Molecule, cell, invert flag
	#; Returns: True or false
	#>
	
	#
	# Note: cylinder length is not implemented
	#
	
	use vars qw($hole_vectors);
	
	my $mol = $_[0];
	my $cell = $_[1];
	my $invert = $_[2];

	my $dum; # point on axis closest to atom (dummy atom)
	my $o = $Silico::origin_points;
	my $radius2 = get_lflag('hr')**2;
	my ($vx, $vy, $vz); # Cylinder axis vector
	
	silico_msg('d', "Hole vectors not defined\n") if !defined $hole_vectors->[0];
	silico_msg('d', "Routine will not work with multiple origins\n") if defined $o->[3];
		
	my $i;
	ORIGIN: for ($i=0; $i < $#{$o}; $i = $i + 3) {
		
		CYL: foreach my $vec (@$hole_vectors) {
				
		($vx, $vy, $vz) = unit_vector($vec->[0], $vec->[1], $vec->[2]);
			
			foreach my $atom (@{$mol->{ATOMS}}) {
			
				my $a;
				$a->{X} = $atom->{X}-$cell->[0] * $o->[$i];
				$a->{Y} = $atom->{Y}-$cell->[1] * $o->[$i+1];
				$a->{Z} = $atom->{Z}-$cell->[2] * $o->[$i+2];
		
				# Point on axis closest to atom
				my $len = ($a->{X} * $vx + $a->{Y} * $vy + $a->{Z} * $vz);
				
				# Compare only to positive part of vector
				if ($len < 0) {
					if ($invert) {
						# Atom is outside. This cylinder fails. Test next one
						next CYL;
					} else {
						# Atom is inside. Good. Test next atom
						next;
					}
				}
				
				$dum->{X} = $len * $vx;
				$dum->{Y} = $len * $vy;
				$dum->{Z} = $len * $vz;
				
				if (distance_periodic_sq($a, $dum, $cell) < $radius2) {
					# Inside
					if ($invert) {
						# Atom is inside. Good. Test next atom
						next;
					
					} else {
						# Fail. Can not be inside any cylinder
						return 0;
					}
				} else {
					# Outside
					if ($invert) {
						# Atom is outside. This cylinder fails. Test next one
						next CYL;
					
					} else {
						# Atom is outside. Good. Test next atom
						next;	
					}
				}
			}
			
			# Molecule lies entirely INSIDE (invert) or OUTSIDE (!invert) a cylinder
			return 1 if $invert;
		}
	}
	
	return 0 if  $invert; 
	return 1;
}

sub find_maxchain_maxsegid {

	my $mol = $_[0];
	
	my $maxsegid = '';
	my $maxchain= 'A';
	
	# Find maximum existing segid and chain
	foreach my $atom (atoms($mol)) {
	
		my $s = $atom->{SEGID} || '';
		$s =~ s/ //g;
		
		$maxsegid = $s if ((!defined($maxsegid)) || $s gt $maxsegid);
		
		if (defined $atom->{CHAIN}) {
			$maxchain = $atom->{CHAIN} if !$maxchain || $atom->{CHAIN} gt $maxchain;
		}
	}
	
	print "Maximum chain found: $maxchain\n" if $maxchain;
	print "Maximum SEGID found: $maxsegid\n" if $maxsegid;
	
	return $maxchain, $maxsegid;
	
}

sub define_pdb_chain_segid {

	#<
	#? Make sure that pdb CHAIN and SEGID are defined and sane
	#
	#. Needs a tidy up
	#
	#; Requires: Molecule
	#; Returns: Maximum pdb chain
	#>
	
	my $mol = $_[0];
	
	my $chain;
	my $maxchain;
	my $maxsegid;
	my $segid;
	
	($maxchain, $maxsegid) = find_maxchain_maxsegid($mol);
	
	if ($maxchain) {
		$chain = $maxchain;
		pdb_increment_chain($chain);
		
	} else {
		$chain = 'A';
	}
	
	if ($maxsegid) {
		$segid = ++$maxsegid;
	} else {
		$segid = "SA";
	}
	print "Initial chain: $chain  SEGID: $segid\n";
	
	my $oldsubid = $mol->{ATOMS}[0]{SUBID};
	my $oldchain = 'XXXXXXXXXX';
	# Make sure that segid and chain are defined
	foreach my $atom (atoms($mol)) {
	
		if ($atom->{SUBID} < $oldsubid && $atom->{CHAIN} eq $oldchain) {
		
			pdb_increment_chain($chain);
			
			silico_msg ('c', "Residue number restart.  Increasing chain to $chain.\n");
			
			if ($chain eq 'A') {
				++$segid;
				silico_msg ('c', "Chain restart.  Increasing SEGID to $segid.\n");
				
			}
		}
		
		$atom->{CHAIN} ||= $chain;
		$atom->{SEGID} ||= $segid;
		$oldsubid = $atom->{SUBID};
		$oldchain = $atom->{CHAIN}
	}
	
	return $chain, $segid;
}

sub origins {

	#<
	#? Origin points for holes 
	#  o    Single point at 0, 0, 0
	#  a	Single point in centre of cell
	#  b 	Single point in centre of cell (X, Y) and 0 in Z
	#  c    4 hexagonally arranged holes in xy plane
	#>
	
	my $type = $_[0] || '';
	
	my $all_origins;
	
	# Type, number of xyz triples, positions (fractional cell)
	
	my @data = qw (

		o       1       0   0   0 
		a	1 	0.5 0.5 0.5
		b	1 	0.5 0.5 0.0
		c       4       0.25 0.75 0.5 0 0.25 0.5 0.75 0.75 0.5 0.5 0.25 0.5
	);
	
	while ($data[0]) {
	
		my $type = shift @data;
		
		my $num = shift @data || die;
		
		for (1..($num*3)) {
			push @{$all_origins->{$type}}, shift @data;
		}
	};
	
	$Silico::origin_points =  $all_origins->{$type};
	
	silico_msg('d', "Invalid origin point type '$type'\n") if !defined $Silico::origin_points;

	print "Origin points: ";
	my $i = 0;
	foreach (@$Silico::origin_points) {
        	print "(" if $i%3==0;
        	print $_;
        	print ")" if ($i+1)%3==0;
        	print " ";
        	++$i;
	}
	print "\n";

	return $Silico::origin_points;
}


sub holes {

	# Points on a sphere taken from R. H. Hardin, N. J. A. Sloane and W. D. Smith
	# http://www2.research.att.com/~njas/electrons/
	
	my $num = $_[0]; # Number of holes to make
	
	die "Number of holes $num is too large\n" if $num > 25;
	die "Number of holes $num no good\n" if $num > 2 && $num < 4;
	
	
	# Vectors evenly distributed on a sphere: x,y,z triples
	my $string = "

0 0 1 x
0 0 1
0 0 -1
x
x	
0.720811059899 -0.605940036634 -0.336553246800
-0.843553324940 -0.170505595776 -0.509259884334
-0.156171632120 -0.175543065384 0.972005685949
0.278913892764 0.951988698026 -0.126192548332
x
0.855565360836 -0.156802100197 0.493377152604
-0.658633419068 -0.744679382224 -0.107956643971
-0.196931925292 0.901481477843 -0.385420500112
0.443792542701 -0.268572623162 -0.854936795986
-0.443792549436 0.268572625140 0.854936791869
x
-0.477973841915 0.243304922042 0.844004574250
-0.169605180501 0.917224510901 -0.360462590776
0.477973836299 -0.243304924988 -0.844004576582
0.861844008545 0.315439234995 0.397143543271
-0.861844004834 -0.315439241394 -0.397143546240
0.169605180866 -0.917224508446 0.360462596851
x
-0.570584671282 -0.374416093280 -0.730921146219
-0.500743916118 -0.546814264721 0.671006475654
0.442735410739 -0.827918915986 -0.344290029816
0.079015190293 0.743083893065 0.664517063416
0.844210175839 -0.137265964975 0.518138238185
-0.795376086465 0.596517091729 -0.107444126667
0.500743902593 0.546814263304 -0.671006486901
x
-0.358740453187 -0.111796186592 -0.926718349830
0.230506662846 -0.878717440372 0.417998012398
-0.687412312891 0.726243033921 -0.005947080580
0.704082308404 -0.438460053911 -0.558588295722
-0.832316101440 -0.552053522700 0.049867979287
-0.270751592548 0.183947828588 0.944910986014
0.398985382895 0.806565799153 -0.436190641660
0.815646098192 0.264270542800 0.514667390372
x
0.914109572223 -0.182781178690 -0.361931942064
0.293329304506 0.734642489361 -0.611766566546
-0.480176899428 -0.046026929940 0.875963279468
-0.705684904851 0.704780196051 -0.072757750931
0.370605109670 0.769162968265 0.520615194684
-0.904030464226 -0.412626217894 -0.111662545460
-0.162180419233 -0.247163999394 -0.955304908927
0.063327560246 -0.997971078243 -0.006583851785
0.610701141906 -0.322016246902 0.723429092590
x
0.978696890330 0.074682616274 0.191245663177
0.537258145625 0.448413180814 -0.714338368164
-0.227939324473 -0.303819959434 -0.925060590777
0.274577116268 0.833436432027 0.479573895237
-0.599426405232 0.240685139624 0.763386303437
-0.424664555168 0.830194107787 -0.361161679833
-0.402701180119 -0.893328907767 0.199487398294
0.552788606831 -0.770301636525 -0.317899583084
0.290107593166 -0.385278374104 0.876012647646
-0.978696887344 -0.074682599351 -0.191245685067
x
0.153486836562 -0.831354332797 0.534127105044
0.092812115769 0.691598091278 -0.716294626049
0.686120068086 0.724987503180 0.060269166267
0.101393837471 0.257848797505 0.960850293931
-0.143059218646 -0.243142754178 -0.959382958495
-0.909929380017 0.200934944687 -0.362841110384
-0.405338453688 0.872713317547 0.272162090194
0.896918545883 -0.184616420020 0.401813264476
0.731466092268 -0.415052523977 -0.541007170195
-0.439821168531 -0.864743799130 -0.242436592901
-0.773718984882 -0.203685975092 0.599892453681
x
-0.519657468039 0.854102387855 -0.021569120792
0.519657456743 -0.854102395012 0.021569109539
-0.118190661273 0.429081803333 -0.895499734024
-0.768920201670 0.071819589916 0.635298095360
0.521506258421 0.836678341438 -0.167333724625
0.119333282457 0.615878166459 0.778751341429
0.915718081273 0.043626880273 0.399445980012
0.768920197428 -0.071819593191 -0.635298100124
0.118190681568 -0.429081793598 0.895499736009
-0.521506259512 -0.836678336747 0.167333744677
-0.915718079511 -0.043626889027 -0.399445983095
-0.119333287884 -0.615878161696 -0.778751344364
x
-0.754267365796 0.173499508397 -0.633228759202
0.330321246216 -0.372227441281 0.867372242037
-0.533035489131 -0.458931921778 0.710812674690
-0.686807855664 -0.708384153523 -0.162747843106
-0.321357401597 0.941504078299 -0.101486407881
0.617869468935 0.786068318474 -0.018273424674
0.019352413289 0.576222568464 0.817063666854
0.311368579326 -0.948900371701 0.051358469549
-0.882505894394 0.313694805173 0.350398224264
-0.038448968078 -0.536051111410 -0.843309482225
0.184148321180 0.468947287620 -0.863815858410
0.936881686583 0.003852303897 0.349625321023
0.813567705110 -0.236010693892 -0.531419365069
x
0.183504500440 -0.774327848934 -0.605592668948
0.153461674965 0.930458382726 0.332711154506
-0.708967620726 -0.280759148895 -0.646946066588
-0.976939230344 -0.107848809998 0.184277981311
0.372001336525 0.124628732488 0.919827529846
-0.280635908687 -0.474838246429 0.834129562169
-0.374884405704 -0.924316474887 0.071419441413
0.609952952660 -0.713378465035 0.345034144926
-0.545724577900 0.497759969260 0.674106592519
-0.009781260146 0.151939133540 -0.988341452459
0.836141855120 -0.174860987699 -0.519894636535
-0.609952964510 0.713378454079 -0.345034146631
0.930390223148 0.274617378657 0.242815419632
0.421433390542 0.757547917835 -0.498512837869
x
0.216001983700 -0.898615365314 0.381881615504
-0.444561363553 -0.418726663881 0.791854263732
-0.635186145240 -0.770065881729 -0.059473512527
0.500548727571 -0.192538217024 0.844026069688
-0.221229771469 0.531169628316 -0.817872981685
-0.891941312897 0.079142510044 -0.445170930601
-0.169315707509 0.513429353508 0.841262438331
-0.263965919566 -0.440254059158 -0.858194824444
0.598524035730 0.006919955250 -0.801074960833
-0.899631198936 0.223429384229 0.375157321885
0.985305483557 -0.166404497880 0.038505157549
0.627905121857 0.627079786181 0.460983838881
-0.350119335865 0.936470382954 0.020968369105
0.527048413266 0.773124394001 -0.352843650184
0.420616988873 -0.804160710609 -0.420008214424
x
0.733571208539 -0.509367037098 -0.449909439244
-0.527936513082 0.197513112343 0.825997341768
0.210490445990 0.934016014999 -0.288631002964
-0.394388830898 0.876999734693 0.274461136430
0.901500079093 0.427890594987 -0.064863287901
-0.126116414468 -0.264945505828 -0.955980401966
-0.807197676003 -0.422051151370 -0.412679945579
0.264549795600 -0.885039709623 0.383038011218
-0.616921323093 -0.639391147814 0.458897636963
-0.969209372877 0.238904864368 0.059646100528
-0.080651424152 -0.925976042468 -0.368868156063
0.393234069013 0.625539858854 0.673844827799
0.874906967098 -0.244744309460 0.417897142739
0.175293851258 -0.219349658857 0.959769656152
-0.509883948542 0.523218324223 -0.682833028065
0.478759086918 0.286782058057 -0.829786591763
x
0.912696377884 -0.226318895321 0.340242677244
-0.937516993909 0.174499515532 0.301017948322
0.217875182394 0.779593349993 -0.587166597776
-0.097521799429 -0.689116202969 -0.718058742333
-0.529244112066 -0.429928221895 0.731479592238
-0.868219130327 -0.253956581808 -0.426264702141
-0.488636597163 -0.872424666927 0.010463099477
-0.090780894434 0.980271295409 0.175576241570
0.241453202690 -0.865087838370 0.439685549927
-0.334318092896 0.477784149911 0.812375355889
0.278050056752 -0.144343667896 0.949659450265
-0.278049987827 0.144343774391 -0.949659454259
0.635522963862 -0.713687640763 -0.294551377228
0.573721357879 -0.050347268156 -0.817501655107
0.881411428344 0.431341360377 -0.192505908512
0.556850158511 0.603621752363 0.570577497838
-0.673293126781 0.653755790140 -0.345368979915
x
-0.115864955267 0.878504081289 0.463471564715
-0.828166364901 -0.445191611216 0.340506830107
0.253495768810 0.923311572233 -0.288505867837
-0.394060866631 -0.892001634136 -0.221470354875
-0.301636361322 -0.023759230509 0.953126961372
0.641266153580 -0.370493887364 0.671946426211
0.389101298212 0.376221841256 0.840866996552
0.423305773042 -0.905396870145 -0.032691436801
-0.095201322600 -0.770256663988 0.630588122118
-0.389101207805 -0.376221819797 -0.840867047989
0.404032202999 -0.499183509195 -0.766533628152
-0.131340369752 0.383991458853 -0.913947628041
-0.929433392676 -0.008826753533 -0.368884340954
0.964071703007 -0.248886818223 -0.092849895935
0.827037608987 0.531769491246 0.182290980308
-0.805777212123 0.401941679568 0.434932144880
-0.581600320724 0.752341956713 -0.309390767645
0.669871864060 0.292136717514 -0.682589059406
x
-0.493849369713 0.150082991549 0.856497458071
0.160995035918 -0.251402381535 0.954398994639
0.216767340524 0.650955278007 0.727508863257
0.514906117305 -0.792513457502 -0.326793681150
-0.313466895221 -0.942439288426 -0.116390262620
0.027820519392 0.999490279267 0.015658874500
0.785474997111 0.078620022781 0.613879402596
-0.501855939280 -0.654402699938 0.565595016353
-0.028731240581 -0.561058899325 -0.827277116390
-0.790845401663 -0.362223441679 -0.493313013174
0.631801474150 0.016389327300 -0.774956958297
0.402342628259 -0.798230082272 0.448273516106
0.785284752008 0.617575749685 -0.043909585115
-0.626007935033 0.705139068047 0.333005945877
-0.979233960994 -0.074186490464 0.188672240324
0.217208195399 0.709925667115 -0.669944883570
-0.272097559178 0.164027352376 -0.948186661983
0.974863117801 -0.210944356125 -0.071724334571
-0.711269733310 0.555278833773 -0.430999748540
x
-0.434938946482 0.685132899717 -0.584312435739
-0.474060013098 0.860976920054 0.184352507755
-0.829256490167 -0.550289985613 0.097542838019
0.464421656080 0.422384590848 -0.778398216068
0.240386556712 0.941384293978 -0.236664138397
0.910939837519 -0.098176548471 -0.400686882429
-0.185116778405 0.035545557996 -0.982073465510
0.266800060066 -0.684012097813 0.678929435209
0.339677264891 0.071129244772 0.937848594526
0.347081475440 -0.527235773025 -0.775600985720
-0.403753197462 0.355548772404 0.842952208598
-0.410800740586 -0.675911546523 -0.611871173379
0.268568798654 0.770374026033 0.578268674928
-0.865131507158 0.028438831904 -0.500738163278
0.884860772610 -0.185065380644 0.427518675596
0.843316814297 0.529181929554 0.093718921021
-0.464421687775 -0.422384607415 0.778398188168
0.596562711719 -0.798812839696 -0.077530498020
-0.938578491452 0.227080195093 0.259817244193
-0.156558083055 -0.985288478562 0.068528684817
x
-0.549492281923 -0.499531375822 -0.669721312694
-0.229274217564 0.798245119516 0.556990181538
0.051584885924 -0.126471373803 -0.990628079126
-0.880573755176 -0.472011665596 0.042365661057
0.874734231412 0.472842236169 -0.106114297292
0.464585318035 0.484632949097 -0.741141947885
0.517653358822 0.668752788677 0.533670973294
-0.214762501849 0.160053900575 0.963462410637
0.178941144896 -0.784300408657 -0.594014255422
-0.517841089726 -0.483670070722 0.705623035677
-0.594929339937 0.796921583984 -0.104761965723
0.586167742727 -0.809073263040 0.042518612639
0.211591983303 -0.690878714671 0.691314280351
0.180352667723 0.967957510562 -0.174731717191
0.511467833181 -0.015518531629 0.859162284319
-0.318668087641 0.525871531163 -0.788612568142
-0.893219122289 0.152006707555 -0.423147208943
-0.254923991529 -0.965694963273 0.049467124975
-0.864507633341 0.245460494780 0.438606540533
0.946614484046 -0.182258255989 0.265900257086
0.793779200124 -0.243314387822 -0.557416083486
x
0.268370018872 -0.725061606257 0.634242225103
-0.164287143257 0.957224109600 -0.238184253388
0.957434682736 0.195909546256 -0.211986504242
-0.591872052500 0.738850387217 0.322160796464
-0.467436167418 -0.502376750030 0.727407059644
0.704035966340 -0.377129383126 -0.601753094285
-0.922930714224 0.021572888337 0.384361167693
-0.846507760234 -0.520886911130 -0.110006534704
0.467436178435 0.502376737773 -0.727407061029
0.118493974706 -0.139937653967 -0.983044572215
-0.291271467276 -0.949098924590 0.119883959279
0.461155213711 -0.879338228320 -0.118744048614
0.123384179729 0.777704965935 0.616401922574
0.536066715387 0.843426372897 -0.035558826703
-0.633932658146 -0.209698378663 -0.744416533214
-0.123384160799 -0.777704956507 -0.616401938259
0.764649497998 0.350186218064 0.540999776238
0.875374065824 -0.377879086249 0.301542104952
-0.875374082613 0.377879061851 -0.301542086787
0.306758933840 -0.059207804241 0.949943889094
-0.302873999755 0.507579158749 -0.806616847008
-0.363289222149 0.245610228245 0.898719398284
x
0.652755932306 -0.746414814390 -0.129516862607
-0.571782213729 -0.244225346813 -0.783210750716
-0.709688239042 -0.669143489929 -0.220430472601
0.074260622701 0.149277887820 -0.986002774907
-0.016875357639 -0.997959076091 -0.061586563074
-0.009413501650 0.687898371154 0.725745972741
0.772596564939 0.523362459372 0.359424935090
0.263157839075 0.959418508514 0.101311782404
0.400751973491 0.108379743006 0.909753640855
0.839823682437 -0.276754944100 0.467014864147
0.009413450197 -0.687898401663 -0.725745944490
0.677548194000 0.538161884565 -0.501308518588
0.609302219506 -0.266442386686 -0.746832819232
-0.975664874445 0.152681067858 -0.157373899655
0.030264912093 0.815154628228 -0.578452216846
0.984938230299 -0.008800677927 -0.172682455866
0.298116671784 -0.690842928425 0.658682395582
-0.280144459258 -0.117059244474 0.952793899661
-0.509780373863 0.859977207332 0.023731272475
-0.684554448284 0.440828515136 0.580564748824
-0.393141091044 -0.758368362464 0.519920676013
-0.874060593742 -0.248142537028 0.417664171058
-0.587825140081 0.476911936064 -0.653465079349
x
0.134985506170 -0.985952791173 0.098366695086
0.650088694587 -0.462793906466 -0.602666150790
-0.206347275695 0.662332382064 0.720233724205
-0.435854930194 0.350113596741 -0.829126618318
0.985232180539 0.063955280317 -0.158830955895
-0.614137826481 -0.392074584561 -0.684917695951
0.353418655052 -0.587198098366 0.728212638957
-0.944785803171 0.119466198906 -0.305135401821
0.021302190861 -0.808690754556 -0.587848177815
0.090322707483 0.786736539321 -0.610645090220
0.019144141922 0.992112273707 0.123881952644
0.685499579269 0.715694962504 -0.133682637122
-0.492664112126 -0.044270667252 0.869092734203
-0.582890874715 -0.806420953907 -0.099616631511
0.048922860209 -0.155267715919 -0.986660270885
0.542444409187 0.622999428166 0.563582980089
0.763772109691 -0.640056077552 0.083548680702
-0.913538674752 -0.294880095995 0.280165698682
-0.347061443006 -0.702169846181 0.621695956149
0.842177009528 -0.028740243511 0.538434659941
0.207815962624 0.070701024559 0.975609497086
-0.578727030552 0.786414267629 -0.215934304310
-0.804218513225 0.456634532848 0.380417515892
0.575100475666 0.281355246534 -0.768178799588
x
0.532212965784 0.104012540712 -0.840196852187
0.250139840140 -0.721834058886 0.645279514479
-0.633506688611 -0.753135632571 -0.177358378534
-0.218922463231 0.850392517067 0.478440719431
-0.413149806364 -0.780294228218 0.469519067675
0.734236034340 -0.661803691272 0.151371463942
0.936469940841 0.054189066241 0.346536570943
-0.333112408327 -0.242859001882 0.911073887576
0.103315488089 -0.991677012629 -0.076828461807
-0.926539655526 0.372983070696 -0.049070313963
0.801488967040 0.543143312217 -0.250221458126
-0.916836769339 -0.261947685782 0.301320009791
-0.537755795155 0.438027927278 -0.720382009562
-0.213559057919 -0.664684481487 -0.715951862105
-0.124138610474 -0.051815651202 -0.990911067493
0.273979937314 0.960273112918 -0.053014550417
-0.443315490443 0.878048334487 -0.180284492503
0.024170841904 0.368771526312 0.929205753202
0.568365159741 0.637098308775 0.520640749604
0.930930389545 -0.168157726945 -0.324178328533
-0.825458835938 -0.182178469442 -0.534255290516
0.472597880311 -0.629557481175 -0.616691674520
0.504263360030 -0.129264527496 0.853820323992
0.149335563174 0.688742589708 -0.709459325610
-0.694717409994 0.344411933027 0.631465074759";
	
	my $set;
	my @f = split ' ', $string;
	my $setcount = 1;
	my $vcount = 0;
	my $i = 0;
	
	foreach (@f) {
	
		if ($_ =~ 'x') {
		
			++$setcount;
			$vcount = 0;
			next;
		}
	
		$set->[$setcount][$vcount][$i] =  $_;
		
		++$i;
		if ($i == 3) {
			$i = 0;
			++$vcount;
		}
	}
	
	return $set->[$num];;
}

sub calc_area_sphere_with_holes {

	#<
	#? Calculate the area of a sphere with holes 
	#. Does not consider hole overlap
	#; Requires: Radius, number of holes, hole radius
	#; Returns: Radius
	#>

	my $rad = $_[0]; # Sphere radius (a)
	my $nholes = $_[1] || 0; # Num holes
	my $hole_rad = $_[2] || 0; # Hole radius

	silico_msg('d', "Hole radius ($hole_rad) is greater than sphere radius ($rad)\n") if $hole_rad > $rad;
	
	my $asphere = 4*$Silico::Pi*$rad**2;
	
	# Height of cap
	my $h = $rad - sqrt($rad**2-$hole_rad**2);
	
	# Area of each cap
	my $cap_area = 2 * $Silico::Pi * $rad * $h;
	
	# Total area of caps
	my $cap_total = $cap_area*$nholes;
	
	# Remaining surface ares
	my $asphere_with_holes = $asphere-$cap_total;
	
	print heading("Sphere with $nholes holes\n");
	printf "Sphere radius: %.1f\n", $rad;
	printf "Sphere area:   %.1f\n", $asphere;
	printf "Number holes:  %d\n",$nholes;
	printf "Hole radius:   %.1f\n",$hole_rad;
	printf "Cap height:    %.1f\n",$h;
	printf "Cap area:      %.1f\n",$cap_area;
	printf "Cap total:     %.1f\n",$cap_total;
	printf "Surface area:  %.1f\n\n", $asphere_with_holes;
	
	return $asphere_with_holes;

}

return 1;
