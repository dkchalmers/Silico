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
#! silico_fp.pm
#? Routines to generate molecular fragments and make Tanimoto comparisons etc
#. $Revision: 1.46.2.4.2.1 $
#>

use strict;
package Silico;

sub setup_fp_flags {

	#<
	#? Make flags required for fingerprints
	#. Requires: max fragment size (default 8), min fragment size (default 3)
	#. Returns: nothing
	#>
	
	my $max = $_[0] || 10;
	my $min = $_[1] || 1;
	my $typer = $_[2] || 'simple';
	
	my $h;
	my $val;
	my $s = 'polar';
	
	$Silico::fp_h ||= 1;
	if ($Silico::fp_h  == 0) {
		$s = 'none';
	} elsif ($Silico::fp_h  == 2) {
		$s = 'all';
	}
	
        $Silico::fp_min = make_flag('fp_min', 'fp-min-fragment-size', 'Min fragment size', 1, $min);
        $Silico::fp_max = make_flag('fp_max', 'fp-max-fragment-size', 'Max fragment size', 1, $max);
	$Silico::fp_fragment_typer = make_flag('fp_t', 'fp-fragment-typers', 'Fragment typing method(s) [simple, steric]', 1, $typer);
	$h         =  make_flag('fp_h',  'fp-hydrogen-treatment', 'Consider hydrogens in fingerprint generation [all, polar, none]', 1, $s);
	
	$Silico::fp_h = 1;
	$Silico::fp_h = 0 if $h eq 'none';
	$Silico::fp_h = 2 if $h eq 'all';
	
	return undef;
}



sub tanimoto_matrix {

	#<
	#? Calculate a matrix of Tanimoto coefficients for two ensembles
	#. Requires: ensemble1, ensemble2 (optional), (alpha, beta coefficiencts for Tversky, optional)
	#. Returns: array
	#>

	my $mols1 = $_[0];
	my $mols2 = $_[1];
	my $alpha = $_[2];
	my $beta = $_[3];
	
	my $c;
	my $i;
	my $j;
	my $mol1;
	my $mol2;
	my $tan;
	
	# Find if mol1 and mol2 are identical 
	$c = 1 if !defined $mols2 || $mols1 == $mols2;
	
	$mols2 = $mols1 if !defined $mols2;

	$i = 0;
	foreach $mol1 (@$mols1) {

		$j = 0;
		foreach $mol2 (@$mols2) {
	
			last if !$c && $j > $i;
			
			#fragment_printout($mol1);
			#fragment_printout($mol2);
		
			my $t = tanimoto($mol1->{FRAGMENT_HASH}, $mol2->{FRAGMENT_HASH}, $alpha, $beta);
			$tan->[$i][$j] = $t;
			$tan->[$j][$i] = $t if !$c;
		
			++$j;
		}
	
		++$i;
	}
	
	return $tan
}


sub mol_tanimoto {

	#<
	#? Calculate the Tanimoto coefficient given two molecules
	#. Tanimoto coeff = Num common fragments / Num fragments in mol1 + Num fragments in mol2 - Num common fragments
	#. Tversky coeff  = Num common fragments / alpha (Num fragments in mol1 - Num common fragments) 
	#                   + beta (Num fragments in mol2 - Num common fragments) + Num common fragments
	#. alpha = 0 and beta = 1 => Molecule 2 as a substructure of Molecule 1
	#. alpha = 1 and beta = 2 => Molecule 2 as a superstructure of Molecule 1
	#. Requires: hash1, hash2, (alpha, beta coefficiencts for Tversky, optional)
	#. Returns: Tanimoto coefficient
	#>
	
	my $mol1 = $_[0];
	my $mol2 = $_[1];
	my $alpha = $_[2];
	my $beta =  $_[3];

	return tanimoto($mol1->{FRAGMENT_HASH}, $mol2->{FRAGMENT_HASH}, $alpha, $beta);
}

sub tanimoto {

	#<
	#? Calculate the Tanimoto (or Tversky) coefficient given two hashes of fragments. 
	#  See http://www.chemaxon.com/jchem/doc/user/query_similarity.html#metrics
	#. Tanimoto coeff = Num common fragments / Num fragments in mol1 + Num fragments in mol2 - Num common fragments
	#. Tversky coeff  = Num common fragments / alpha (Num fragments in mol1 - Num common fragments) 
	#                   + beta (Num fragments in mol2 - Num common fragments) + Num common fragments
	#. alpha = 0 and beta = 1 => Molecule 2 as a substructure of Molecule 1
	#. alpha = 1 and beta = 2 => Molecule 2 as a superstructure of Molecule 1
	#. Requires: hash1, hash2, (alpha, beta coefficiencts for Tversky, optional)
	#. Returns: Tanimoto coefficient
	#>

	my $hash1 = $_[0];
	my $hash2 = $_[1];
	my $alpha = $_[2];
	my $beta =  $_[3];
	
	$alpha = 1 if !defined $alpha;
	$beta = 1  if !defined $beta;
	
	my $denominator;
	my $keyboth = 0;
	my $key;
	my $numkeys1; # Use of keys in scalar context gives number of elements (not indexed from zero)
	my $numkeys2; # Ditto
	my $tan;

	$numkeys1 = keys (%$hash1);
	$numkeys2 = keys (%$hash2);

	silico_msg('d', "Molecule 1 has no fragment keys\n") if $numkeys1 == 0;
	silico_msg('d', "Molecule 2 has no fragment keys\n") if $numkeys2 == 0;
	
	foreach $key (keys %$hash2) {
	
		if ($hash1->{$key}) {
			++$keyboth;
		}
	}
	
	# Avoid divide by zero
	$denominator =  $alpha * ($numkeys1 - $keyboth) + $beta * ($numkeys2 - $keyboth) + $keyboth;
	#$denominator =  $numkeys1 + $numkeys2 - $keyboth;
	
	if ($denominator != 0 ) {
		$tan = $keyboth / $denominator;
	} else {
		silico_msg('d', "Divide by zero error\n");
	}
	
		
	return $tan;
}



sub tanimoto_freq {

	#<
	#? Calculate the Tanimoto coefficient given two hashes of fragments including fragment frequencies
	#. Requires: hash of fragments 1, hash of fragments 2
	#. Returns: Tanimoto coefficient
	#>

	my $hash1 = $_[0];
	my $hash2 = $_[1];
	
	my $common;
	my $keyboth = 0;
	my $key;
	my $numkeys1; # Use of keys in scalar context gives number of elements (not indexed from zero)
	my $numkeys2; # Ditto
	my $tan;
	
	foreach $key (keys %$hash1) {
	
		$numkeys1 += $hash1->{$key};
	}
	foreach $key (keys %$hash2) {
	
		$numkeys2 += $hash2->{$key};
	}
	
	foreach $key (keys %$hash1) {
		
		if ($hash2->{$key}) {
			if ($hash1->{$key} < $hash2->{$key}) {
				$common = $hash1->{$key};
			} else {
				$common = $hash2->{$key};
			}
			$keyboth += $common;
		}
	}

	$tan = $keyboth / ($numkeys1 + $numkeys2 - $keyboth);
		
	return $tan;
}

sub tan_substructure_frac {

	#<
	#? Calculate fraction of substructure keys present in a parent molecule
	#. frac = Num common fragments / Num fragments in mol1
	#. Requires: hash of fragments 1, hash of fragments 2
	#. Returns: Tanimoto coefficient
	#>

	my $hash1 = $_[0]; # Substructure
	my $hash2 = $_[1]; # Superstructure
	
	my $keyboth = 0;
	my $key;
	my $numkeys1; # Use of keys in scalar context gives number of elements (not indexed from zero)
	
	$numkeys1 = keys (%$hash1);
	
	my $count = 0;
	foreach $key (keys %$hash1) {
	
		++$count;
	
		if ($hash2->{$key}) {
			++$keyboth;
		}
	}
		
	return $keyboth/$numkeys1;
}

sub euclidian {

	#<
	#? Calculate the Euclidian coefficient given two hashes of fragments
	#. Euclidian coeff = 1 - (Num unshared fragments mol1 + Num unshared fragements mols)/ (Total fragments mol1, mol2)
	#. Requires: hash of fragments 1, hash of fragments 2
	#. Returns: Euclidian coefficient
	#>

	my $hash1 = $_[0];
	my $hash2 = $_[1];
	
	my $euclidian;
	my $key;
	my $notshared;
	my $numkeys1; # Use of keys in scalar context gives number of elements (not indexed from zero)
	my $numkeys2; # Ditto
	my $total;
	
	(defined $hash1 && defined $hash2) || silico_msg('d', "Fragment hash is undefined");
	
	$numkeys1 = keys (%$hash1);
	$numkeys2 = keys (%$hash2);
	
	$total = 0;
	$notshared = 0;

	foreach $key (keys %$hash1) {

		++$total;
		if (!$hash2->{$key}) {
			++$notshared;
		}
	}
	foreach $key (keys %$hash2) {

		++$total;
		if (!$hash1->{$key}) {
			++$notshared;
		}
	}

	# Stop divide by zero
	return 0 if $total == 0;

	$euclidian = 1 - $notshared/$total;
}


sub ensemble_calc_fragments {

	#<
	#? Calculate molecular fragments for each molecule in an ensemble
	#. See mol_fragment for additional information. Sets $mol->{FORMULA}, 
	#  $mol->{NUMHEAVY} and  $mol->{SDF_DATA}{NUMHEAVY}
	#; Requires: Ensemble, typer (optional)
	#>

	my $mols = $_[0];
	my $typer = $_[1] || $Silico::fp_fragment_typer;

	my $i;
	my $j;
	my $k;
	my $numfrag;
	my $mol;
	my $store;
	
	my $n = $#{$mols}+1;
 
	silico_msg('c', "\n");
	silico_msg('c', "Fragmenting $n compounds\n");
	silico_msg('c', "Minimum fragment size: $Silico::fp_min\n");
	silico_msg('c', "Maximum fragment size: $Silico::fp_max\n");
	silico_msg('c', "Using fragment typer: $typer\n");
	silico_msg('c', "Ignoring hydrogens in fragments\n") if $Silico::fp_h == 0;

	$store = $|;
	$| = 1;
	
	if ($n >=  80) {
		print "\n";
		for ($i = 0; $i <=100; $i += 10) {
			print "$i";
			print " "x(8-length($i+10));
		}

		print "\n|-------|".((("-"x7)."|")x9)."\n";
		print " ";
	}
	
	$j = 1;
	$i = 0;
	$k = $n/80;
	my $starttime = (times)[0];
	foreach $mol (@$mols) {

		if (!defined $mol->{ATOMS}) {
			next;
		}
		
		++$i;
		
		mol_fragment($mol, $typer);

		# Calculate MW and number of heavy atoms for each molecule
		$mol->{MW} = molecule_mw ($mol, undef, 1);
		$mol->{FORMULA} = molecule_formula($mol, 1);
		$mol->{NUMHEAVY} = $mol->{NUMATOMS} - ($mol->{FORMULA}{H} || 0);
		$mol->{SDF_DATA}{NUMHEAVY} = $mol->{NUMHEAVY};
		
		$numfrag = $#{$mol->{FRAGMENT_ARRAY}}+1;
		silico_msg('g', "Molecule name: $mol->{NAME}\n",
				"Molecular weight: $mol->{MW}\n",
				"Number of Fragments: $numfrag\n");
			
		if ($n >= 80) {
		
			if ($i > int($j*$k) && $j < 80) {
				print ".";
				++$j;
			}
		}
	}

	print "\n\n";

	$| = $store;

	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
}
	

sub mol_fragment {

	#<
	#? Generate molecule fragments and save in both an array and a hash
	#. By default all unique fragments are returned.  ALL fragments (including duplicates)
	#  are returned if the 'return_paths' flag  is set.  $mol->{PATH_ARRAY} is also set 
	#  in this case.
	#; Requires: molecule, fragment typing method, flag to add hydrogens, flag to delete preexisting fragments, 
	#  flag to return paths, flag to be quiet
	#  
	#; Sets: $mol->{FRAGMENT_HASH}, $mol->{FRAGMENT_ARRAY}, $mol->{FP_HASH} and $mol->{FP_ARRAY}
	#  and optionally PATH_ARRAY
	#; Returns: Array of fragments
	#>

	my $mol = $_[0];
	my $typer ||= $_[1] || $Silico::fp_fragment_typer;
	my $addh = $_[2]; # Flag to override global $Silico::fp_h
	my $reset = $_[3]; # Delete preexisting fragments.  Otherwise new fragments are added to the existing ones
	my $return_paths = $_[4]; # Flag to also return
	my $quiet = $_[5] || 1;
	
	my $farray = ();
	my $frag;
	my $fragments;
	my $fhash;
	my $fp;
	my $fparray;
	my $fphash;
	my $num1;
	my $path;
	my $paths;

	# Reset stuff if required
	if ($reset) {
		delete $mol->{FRAGMENT_HASH};
		delete $mol->{FRAGMENT_ARRAY};
		delete $mol->{FP_ARRAY};
		delete $mol->{FP_HASH};
	}
	
	# Make sure bond orders are sane
	molecule_check_and_fix_connectivity($mol);
	
	# Make sure there are hydrogens if they are required
	if ($Silico::fp_h == 1 || $addh) {
		(undef, undef, $num1) = mol_add_hydrogens($mol, 1, 1);
		print "Added $num1 hydrogen(s)" if !$quiet;
	}
	if ($Silico::fp_h == 2 || $addh) {
		(undef, undef, $num1) = mol_add_hydrogens($mol, 1, 0);
		print "Added $num1 hydrogen(s)" if !$quiet;
	}
	
	@$fragments = ();
	$paths = molecule_path($mol);

	foreach $path (@$paths) {
	
		my $frags = mol_path_to_fragments($mol, $path, $typer);
		
		foreach $frag (@$frags) {
		
			++$fhash->{$frag};
			push @$farray, $frag;
		}
	}
	
	$mol->{FRAGMENT_HASH} = $fhash;
	$mol->{FRAGMENT_ARRAY} = $farray;
	$mol->{FP_ARRAY} = $fparray;
	$mol->{FP_HASH} = $fphash;
	$mol->{FP_TYPER} = $typer;

	if ($return_paths) {
		$mol->{PATH_ARRAY} = $paths;
		return $farray, $paths;	
	}
	
	#print "eee $mol->NAME\n";
	#molecule_printout($mol);
	
	return $farray;
}

sub mol_fragment_and_path {

	#<
	#? Generate ALL molecule fragments and corresponding atom paths
	#. Paths and fragment arrays are stored as mol->{FRAGMENT_ARRAY} and {PATH_ARRAY}
	#; Requires: molecule, optional: atomlist, fragment typer, min fragment size
	#  max fragment size, flag to stop sorting fragments by size, flag to 
	#  force generation of fragments rather than use cached ones
	#; Returns: Array of fragments
	#>

	my $mol = $_[0];
	my $atomlist = $_[1] || $mol->{ATOMS};
	my $typer = $_[2] ||  $Silico::fp_fragment_typer;
	my $min = $_[3] || $Silico::fp_min;
	my $max = $_[4] || $Silico::fp_max;
	my $nosort = $_[5]; 	# Flag to prevent sorting fragments by size
	my $force = $_[6];      # Flag to force generation of fragments rather than use cached ones

	my $frag;
	my $frag_list;
	my $hash;
	my $path;
	my $path_list;
	my $path_list_new;
		
	return ($mol->{FRAGMENT_ARRAY}, $mol->{PATH_ARRAY}) if !$force && defined $mol->{FRAGMENT_ARRAY} && defined $mol->{PATH_ARRAY};

	$path_list = molecule_path($mol, $atomlist, $min, $max, $nosort);
	
	foreach $path (@$path_list) {
		my $frags = mol_path_to_fragments($mol, $path);
		foreach $frag (@$frags) {
			push @$frag_list, $frag;
			push @$path_list_new, $path;
		}
	}

	$mol->{FRAGMENT_ARRAY} = $frag_list;
	$mol->{PATH_ARRAY} = $path_list_new;
	return ($frag_list, $path_list_new);
}


sub mol_path_to_fragments {

	#<
	#? Convert an atom path into fragments
	#. Uses $Silico::fp_fragment_typer.
	#; Atom types; element (atomic element), p_np (polar [p] or nonpolar [n]) or none (all set to X)
	#; Bond types; simple (set to single, double/aromatic, triple), none (all set to single)
	#; Requires: molecule, path, typer (optional otherwise uses  $Silico::fp_fragment_typer);
	#; Returns: fragment
	#>

	my $mol = $_[0];
	my $path = $_[1];
	my $typer = lc($_[2] || $Silico::fp_fragment_typer);
	
	silico_msg('d', "No typer provided") if !$typer;
	
	my $atom;
	my $atom_typer;
	my @atom_typers;
	my $b;
	my $bo;
	my $bond_typer;
	my @bond_typers;
	my $frag;
	my $frags;
	my $i;
	my $md5_hashes;
	my $tok;
	
	if ($typer eq 'simple') {
	
		@atom_typers = ('element');
		@bond_typers = ('simple');
	
	} elsif ($typer eq 'steric') {

		@atom_typers = ('element', 'subst');
		@bond_typers = ('simple', 'simple');

	} elsif ($typer eq 'valence') {

		@atom_typers = ('subst');
		@bond_typers = ('simple');

	} else {
		silico_msg('d', "No valid typer '$typer'\n");
	}
	
	while ($atom_typer = shift @atom_typers) {
	
		$frag = '';
	
		$bond_typer = shift @bond_typers;
	
		$i = -1;
		foreach $atom (@$path) {
	
			my $el = $atom->{ELEMENT};

			++$i;
			
			# Atom typing
			next if !$Silico::fp_h && $el eq 'H';

			if ($atom_typer eq 'element') {

				$tok = "$atom->{ELEMENT}";

			} elsif ($atom_typer eq 'p_np') {
			
				if ($el eq 'N' || $el eq 'O' || $el eq 'P') {
			
					$tok = 'p';
				} else {
					$tok = 'n';
				}
			
			} elsif ($atom_typer eq 'subst') {

				# Count substituents if it hasn't been done already
				if (!defined $atom->{SUBST}) {
					my $c = 0;
					foreach (@{$atom->{CONNECT}}) {
						next if $mol->{ATOMS}[$_]{ELEMENT} eq 'H';
						++$c;
					}
					$atom->{SUBST} = $c;
				}

				$tok= $atom->{SUBST};
					
			} elsif ($atom_typer eq 'none') {
				$tok = 'X';
			} else {
				silico_msg('d', "Atom typer has not been set\n");
			}
	
			# Bond typing
			if ($bond_typer eq 'simple') {
			
				# Make sure 6-membered rings have aromatic bond orders
				# Note: this is somewhat problematic because it changes the actual bonding in the molecule
				# which must be UNDONE to write to SDF format.  
				# It is also a problem that 5-membered rings are ignored
				
				#make_aromatic_bonds($mol) if !$mol->{AROMATIC_BONDS};
				
				# Workaround
				molecule_find_rings($mol, 6);
				
				if ($i > 0) {
					$bo = find_bondorder2($path->[$i-1], $atom, $mol);
					
					if ($path->[$i-1]{PLANAR_RING} && $atom->{PLANAR_RING}) {
						# Note this will create aromatic bonds BETWEEN aromatic rings
						# but this may not be too much of a problem (because it is somewhat true)
						$b = ":";
					} elsif ($bo == 1) {
						$b = '';
					} elsif ($bo == 2) {
						$b = ':';
					} elsif ($bo == 3) {
						$b = '%';
					} elsif ($bo == 4) {
						$b = ':' 
					} else {
						$b = '!';
					}
				} else {
					$b = '';
				}
			} elsif ($bond_typer eq 'none') {
		
				$b = '';
				# single (ie all bonds are the same);
			} else {
				silico_msg('d', "Bond typer has not been set\n");
			}

			#if (!defined $tok) {
			#	atom_printout($atom);
			#	print "b $b tok $tok typer $typer atom_typer $atom_typer el $el i $i\n";
			#}
			
			$frag .= $b.$tok;
		
		}
		
		push @$frags, $frag;
	}

	return ($frags);
}

sub frag_aray_to_hashes {

	#<
	#? Generate md5 hashes from fragments
	#; Requres: array of fragments
	#; Returns: array of md5 hashes
	#>

	my $fragments = $_[0];
	
	my $frag;
	my $hashes;
	
	use Digest::MD5  qw(md5_hex);
	
	foreach $frag (@$fragments) {
	
		push @$hashes, frag_to_hash($frag);
	}
	
	return $hashes;
}

sub frag_to_hash {

	use Digest::MD5  qw(md5_hex);
	return md5_hex($_[0]);

}


sub molecule_path {

	#<
	#? Call recursive routine to find all paths in a molecule
	#; Requires: molecule, optional pointer to array of atom numbers
	#  Max fragment size, min fragment size hydrogen treatment
	#  are optionally set using globals $Silico::fp_max, $Silico::fp_min 
	#. Note will return a previously calculated, cached path unless the 'force' flag is set
	#; Returns: pointer to array of atom number arrays (paths)
	#>
	
	my $mol = $_[0];
	my $atomlist = $_[1] || $mol->{ATOMS};
	my $min = $_[2] || $Silico::fp_min;
	my $max = $_[3] || $Silico::fp_max;
	my $nosort = $_[4]; 	# Flag to prevent sorting fragments by size

	my $atom;
	my $atom_paths;
	my $path;

	# Make sure that we have bonds
	molecule_check_and_fix_connectivity($mol);
	
	silico_msg('g', "Fp_max $Silico::fp_max Fp_min $Silico::fp_min\n");

	@$atom_paths = ();
	
	foreach $atom (@$atomlist) {
		@$path = ($atom);
		atom_path_traverse($mol, $path, $atom_paths, $max, $min);
	}	

	if (!$nosort) {
		@$atom_paths = sort { $#{$b} <=> $#{$a}} @$atom_paths;
	}

	return $atom_paths;
}


sub atom_path_traverse {

	#<
	#? Recursive subroutine to enumerate LINEAR atom paths starting from a particular atom
	#; Requires: molecule, pointer to array containing starting atom number, path list
	# (pointer to empty array to start)
	#; Returns: nothing
	#; Sets: path_list to an array of atom paths
	#>
	
	my $mol = $_[0];	# Molecule
	my $path = $_[1];	# Current path (pointer to array of atoms) 
	my $path_list = $_[2];	# Pointer to array of paths collected so far (pointer to array of paths)
	my $max = $_[3];
	my $min = $_[4];
	
	my $el;	
	my $i;
	my $natom;
	
	my $atom = $path->[0];
	my $parent = $path->[1] || -1; # Parent = -1 if not defined

	if ($atom->{ELEMENT} eq 'H') {
		
		# Fragments contain no hydrogens
		return if $Silico::fp_h == 0;
		
		# Fragments contain polar hydrogens
		$el = $mol->{ATOMS}[$atom->{CONNECT}[0]]{ELEMENT};
		return if $Silico::fp_h == 1 && !($el eq 'N' || $el eq 'O' || $el eq 'S');
	}

	# Exit. Path has looped back on itself
	foreach $i (1..$#{$path}) {
		return if ($path->[$i] == $atom);
	}

	# Save path if .ge. fp_min
	if ($#{$path}+1 >= $min) {
		my $new;
		@$new = @$path;
		push @$path_list, $new;
	}

	# Exit if we have reached maximum length
	return if ($#{$path}+1 == $max);
	
        # Add new atom to path and continue
        foreach (@{$atom->{CONNECT}}) {

		$natom = $mol->{ATOMS}[$_];
                next if ($natom == $parent); # Skip parent
		my $npath;
		@$npath = ($natom, @$path);
                atom_path_traverse ($mol, $npath, $path_list, $max, $min) ;
        }
}

sub atom_path_traverse_branched {

	#<
	#? Recursive subroutine to enumerate BRANCHED atom paths starting from a 
	# particular atom
	#; Requires: molecule, pointer to array containing starting atom number, 
	#  path list
	# (pointer to empty array to start)
	#; Returns: nothing
	#; Sets: path_list to an array of atom paths
	#>

	my $mol = shift;	# Molecule
	my $path_list = shift;	# Pointer to array of paths collected so far 
			        # (arrays of atom numbers)
	my @path = @_;		# Current path (array of atom numbers)

	my $con;
	my $connum;
	my $el;	
	my $c = -1;
	my $natom;
	my $patom;
	my $pconnum;
	
	my $length = $#path+1;
	
	print "\nstart: path '@path' ";
	foreach (@path) {
		print "$mol->{ATOMS}[$_]{ELEMENT} ";
	}
	print "\n";
	
	# Loop over atoms in path
	foreach $pconnum (@path) {
	
		++$c;
		next if $c == 0;	# Do not branch at first atom of path
		next if $c == $#path;   # Do not branch at last atom of path
		
		$patom = $mol->{ATOMS}[$pconnum];
		
		print "\tpatom $pconnum $patom->{ELEMENT}\n";
	
		 # Loop over all atoms attached to parent atom
      		LOOP: foreach $connum (@{$patom->{CONNECT}}) {
		
			$con = $mol->{ATOMS}[$connum];
			
			# Deal with hydrogens
			if ($con->{ELEMENT} eq 'H') {
		
				# Fragments contain no hydrogens
				next if $Silico::fp_h == 0;
		
				# Fragments contain polar hydrogens
				$el = $mol->{ATOMS}[$con->{CONNECT}[0]]{ELEMENT};
				next if $Silico::fp_h == 1 && !($el eq 'N' || $el eq 'O' || $el eq 'S');
			}
			
			print "\t\tatom $connum $con->{ELEMENT} : ";
		
			# Path has looped back on itself
			foreach (@path) {
				if ($connum == $_) {
					print "loop\n";
					next LOOP;
				}
			}
			
			my $new;
			@$new = (@path[0..$c-1], '(', $c, $connum, ')', @path[($c+1)..$#path]);
			
			# Save path if .ge. fp_min
			if ($length >= $Silico::fp_min) {
				
				print "saved '@$new' ";
				foreach (@$new) {
					print "$mol->{ATOMS}[$_]{ELEMENT} ";
				}
				print "\n";
				
				push @$path_list, $new;
			}
		}
	}
}


sub rms_fp {

	#<
	#? Calculate the fragment-based RMS for two molecules
	#; Requires: molecule1, molecule2
	#; Returns: best fragment rms or -1 if no fit is found
	#; Sets: {SDF_DATA}{RMS_FP_AV} {SDF_DATA}{MATCH_NAME} ;
	#>

	my $mol1 = $_[0];
	my $mol2 = $_[1];
	my $thresh = $_[2] || 2;

	my $av;
	my $atom;
	my $bestf;
	my $besti;
	my $bestj;
	my $bestsumsq;
	my $d;
	my $frag1;
	my $frag2;
	my $fraglist1;
	my $fraglist2;
	my $frag_count = 0;
	my $frag_count_thresh = 0;
	my $frag_rms;
	my $i;
	my $j;
	my $len1;
	my $len2;
	my $maxsize;
	my $numfit;
	my $numfit_total = 0;
	my $numfrag1;
	my $numfrag2;
	my $path1;
	my $path2;
	my $path_list1;
	my $path_list2;
	my $sumsq;
	my $sumsq_total;
	my $t1 = 0;
	my $t2 = 0;
	
	clear_fmatch($mol1);
	clear_fmatch($mol2);

	($fraglist1, $path_list1) = mol_fragment_and_path($mol1);	
	($fraglist2, $path_list2) = mol_fragment_and_path($mol2);	

	$numfrag1 = $#{$fraglist1}+1;
	$numfrag2 = $#{$fraglist2}+1;
	
	fragmentlist_printout($fraglist1) if $Silico::debug;
	fragmentlist_printout($fraglist2) if $Silico::debug;
	
	# Get maximium size fragment from molecule 1
	$maxsize = $#{$path_list1->[0]};
	#print "Using maximum fragment size $maxsize for molecule 2\n";

	print "\nMol1: $mol1->{NAME}  Fragments: $numfrag1\nMol2: $mol2->{NAME}  Fragments: $numfrag2\n";
		
	$i = -1;
		
	LOOP: foreach $path1 (@$path_list1) {

		++$i;
		$j = -1;
		$len1 = $#{$path1}+1;
		$frag1 = $fraglist1->[$i];

		$bestsumsq = 999999;
		foreach $path2 (@$path_list2) {

			++$j;
			$frag2 = $fraglist2->[$j];
			next if $frag1 ne $frag2;

			$len2 = $#{$path2}+1;
			last if $len2 < $len1;
			next if $len1 != $len2;

			# Calculate sum of squares for pair of fragments
			$sumsq = rms_atomlist($path1, $path2, 1);

			# Retain only the result with lowest sumsq
			if ($sumsq < $bestsumsq) {
			
				$bestsumsq = $sumsq;
				$besti = $i;
				$bestj = $j;
				$bestf = $frag2;
				$numfit = $len1;
			}
		}

		# There was no matching fragment
		next if $bestsumsq == 999999;

		$frag_rms = sqrt($sumsq/$numfit); 		# RMS for fragment pair
		++$frag_count;					# Total number of matching fragments
		++$frag_count_thresh if $frag_rms < $thresh;	# Matching fragments less than a given threshold
		$numfit_total += $numfit;			# Total number of fitted atom pairs
		$sumsq_total += $bestsumsq;

		# Mark matched atoms
		foreach $atom (@{$path_list1->[$besti]}, @{$path_list2->[$bestj]}) {
			$atom->{FMATCH} = 1;
		}

		if ($Silico::debug) {
			print "Fragment: $bestf\n";
			print "Atoms mol1: ";
			foreach (@{$path_list1->[$besti]}) {
				print "$_->{NUM}$_->{ELEMENT} ";
			}
			print "\nAtoms mol2: ";
			foreach (@{$path_list2->[$bestj]}) {
				print "$_->{NUM}$_->{ELEMENT} ";
			}

			print "\nFragment RMS: $frag_rms using $numfit atoms\n\n";
		}
	}

	foreach $atom (atoms($mol1)) {
		++$t1 if $atom->{FMATCH};
	}
	foreach $atom (atoms($mol2)) {
		++$t2 if $atom->{FMATCH};
	}

	if ($numfit_total == 0) {
		print "\nNo fit found!\n\n";
		return -1;
	} 
		
	$av = sqrt($sumsq_total/$numfit_total);
	print "Matched $frag_count fragments\n";
	print "Matched $t1 of ".mol_count_heavy_atoms($mol1)." in Mol1\n";
	print "Matched $t2 of ".mol_count_heavy_atoms($mol2)." in Mol2\n";
	printf "Fragment RMS %.3f\n", $av;
	
 	$mol2->{SDF_DATA}{RMS_FP_NUMATOMS} = $numfit_total;  
	$mol2->{SDF_DATA}{MATCH_NAME} = $mol2->{NAME} || '';
	$mol2->{SDF_DATA}{RMS_FP} = $av;
		
	return $av;
}

sub clear_fmatch {

	my $mol= $_[0];

	my $atom;

	foreach $atom (@{$mol->{ATOMS}}) {
		delete $atom->{FMATCH};
	}
}

sub fragment_printout {

	#<
	#? Print out list of fragments
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];
	
	my $num = keys %{$mol->{FRAGMENT_HASH}};
	
	print heading("Fragments $mol->{NAME}\t$num\n");
	
	sub sortfunc {
		return 1 if length($a) > length($b);
		return -1 if length($b) > length($a);
		return $a cmp $b
	}
	
	foreach (sort sortfunc keys %{$mol->{FRAGMENT_HASH}}) {
	
		printf "%30s %d\n", $_, $mol->{FRAGMENT_HASH}{$_};
	
	}

}

sub fragmentlist_printout {

	#<
	#? Print out list of fragments
	#; Requires: pointer to fragments
	#; Returns: nothing
	#>
	
	my $f = $_[0];
	
	my $i = 0;
	my $frag;
	
	print "Fragments\n\n";
	
	foreach $frag (@$f) {
		++$i;	
		print "$i\t$frag\n";
	}
	
	print "\n";
}

return 1;
