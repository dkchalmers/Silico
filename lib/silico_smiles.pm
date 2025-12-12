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
#! silico_smiles.pm
#? Silico library dealing with SMILES strings
#. $Revision: 1.20.2.9.2.20 $
#>

use strict;

package Silico;

#
#  Smiles reading routines
#

sub read_smiles {

	#<
	#? Simple routine to read smiles strings (one per line) from file and return as an ensemble
	#. Reads multiple structures from a single file. Information
	#, QUIET - do not print 'Reading' line
	#>

	my $infile = $_[0];
	my $start = $_[1] || get_flag('ss', 's') || 0; 
        my $max_molecules = $_[2] || get_sflag('ms');
	my $options = uc($_[3] || '');

	my $mols;

	my $fr = open_molfile(undef, $infile, undef, undef, $options);
	return undef if !$fr;

	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Reading smiles file: $infile");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}
	if ($options =~ /\bDTAB\b/) {
		$fr->{DELIMITER} = "\t";
	}
	if ($options =~ /\bDCOMMA\b/) {
		$fr->{DELIMITER} = ',';
	}
	if ($options =~ /\bDSPACE\b/) {
		$fr->{DELIMITER} = ' ';
	}

	my $i = -1;
	my $molcount = 0;
	
	while (my $mol = read_smiles_single($fr, $options)) {
		++$i;
		next if $start and ($i < $start);
		++$molcount;
		push @$mols, $mol;
		last if (defined $max_molecules) && $molcount == $max_molecules;
	}
	
	silico_msg('c', "Read $molcount molecules\n");

	return $mols;
}

sub read_smiles_single {

	#<
	#? Read in a single smiles record from file.
	#; Requires: file record, options
	#. See read_pdb for general description
	#. Options:
	#; Returns: molecule or undef
	#>

	my $fr = $_[0];
	my $options = uc($_[1] || '');
	
	my $line;
	while (1) {

		$line = silico_readline($fr);

		return if !defined $line; # We have come to the end of the file

		chomp $line;
		++$fr->{LINECOUNT};

		# Read header if it is present ( starts with '#' or contains the word 'smiles' )
		if ($fr->{LINECOUNT} == 1) {

			if ($line =~ /^#/ || $line =~ /smiles/i) {

				$line =~ s/^#//;

				silico_msg('c', "Found header line\n");

				if ($line =~ /SMARTS/i) {
					silico_msg('e', "Smarts strings are not supported\n");
					return undef;
				} 

				if (! defined $fr->{DELIMITER}) {
					if ($line =~ "\t") {
						silico_msg ('c', "Guessing tab delimited. Setting delimiter\n");
						$fr->{DELIMITER} = "\t";
					
					} elsif ($line =~ ',') {
						
						silico_msg ('c', "Guessing comma delimited. Setting delimiter\n");
						$fr->{DELIMITER} = ",";
						
					} else {
						silico_msg ('c', ". Guessing space delimited. Setting delimiter\n");
						$fr->{DELIMITER} = " ";
					}
				}
				
				@{$fr->{HEADERS}} = split "$fr->{DELIMITER}", $line;
				
				print "Headers: @{$fr->{HEADERS}}\n";
					
				next;

			} else {

				silico_msg('c', "No header line present\n");
				if (! defined $fr->{DELIMITER}) {
					silico_msg ('c', ". Assuming space delimited. Setting delimiter\n");
					$fr->{DELIMITER} = " ";
				}
			}
		}

		next if $line =~ /^\s*$/; # Ignore blank lines
		next if $line =~ /^\s*#/; # Ignore comments
		last;
	}

	my $sdf_data;
	my $smiles;
	my $name = '';
	
	#print "d '$fr->{DELIMITER}'\n";
	
	if (defined $fr->{DELIMITER}) {
	
		# We have a delimited file (with headers)
		my @vals = split /$fr->{DELIMITER}/, $line;
		
		#print "heads ".(join "x",@{$fr->{HEADERS}})."\n";
		#print "vals ".(join "x",@vals)."\n";
		
		my $i = 0;
		$smiles = $vals[0];	
		foreach my $val (@vals) {
			if ($i > 0) {
				$fr->{HEADERS}[$i] = "header_$i" if !defined $fr->{HEADERS}[$i];
				$sdf_data->{$fr->{HEADERS}[$i]} = ($val || '');	
			}
			++$i;
		}
		foreach (qw(Name name NAME)) {
			$name = $sdf_data->{$_} if $sdf_data->{$_};
		}
		
	} else {

		# We have a non-delimited file
		($smiles, $name) = split(' ', $line, 2);
	}

	#print "smiles '$smiles'\n";
	
	silico_msg('w', "Smiles string '$smiles' contains tab\n") if $smiles =~ "\t";
	silico_msg('w', "Smiles string '$smiles' contains comma\n") if $smiles =~ ",";
	
	$smiles ||= '';
	
	my $mol = parse_smiles($smiles);
	
	return undef if !defined $mol;
	
	$mol->{SDF_DATA} = $sdf_data;

	# Process errors in here
	
	++$fr->{MOL_READ_COUNT};
	$mol->{SOURCE_FILE_TYPE} = 'smi';
	$mol->{SOURCE_SMILES} = $smiles;
	$mol->{SMILES} = $smiles;
	$mol->{SOURCE_FILE_NAME} = $fr->{FILENAME};
	
	return $mol;
}

sub parse_smiles {

	#<
	#? Parse a smiles string and return a molecule
	#. Note. Can not handle chirality.
	#; Requires: string
	#; Returns: molecule
	#>

	my $smiles = $_[0];

	if (!defined $smiles) {
		silico_msg('w', "Undefined smiles string passed to subroutine\n");
		#carp();
		return undef;
	}


	my $acount;
	my $bond = 0;
	my $bracketlevel;
	my $fchg;
	my $hash;
	my $mol = {};
	my $prevatom;
	my @stack;
	my $warn_chirality;

	chomp $smiles;

	$mol->{HAS_NO_COORDS} = 1;

	# Split smiles up into tokens
	my @f = $smiles =~ m/(\[.*?\]|Cl\d*|Br\d*|[A-Z]\d*|[a-z]\d*|\W)/g;

	#print "Smiles: '$smiles' Tokens: @f\n";

	my $x = 0;
	my $y = 0;
	my $z = 0;
	my $a = 0.4;
	
	foreach my $tok (@f) {

		#print "tok: '$tok'\n";
		
		my ($isotope, $el, $at, $h, $sign, $sign2, $num);

		my $atomflag = 0;
		if ($tok =~ /^\[/) {

			#($isotope, $el, $at, $h, $sign, $num) = $tok =~ /^\[(\d*)([a-z]+)(\@{0,2})(H?\d*)([+-])(\d*)/;

			my @f = $tok =~ /^\[(\d*)([A-Za-z][a-z]{0,2})(@*)(H*\d*)([+-]*)(\d*)](\d*)/;
			($isotope, $el, $at, $h, $sign, $sign2, $num) = @f;

			if (!defined $el) {
				silico_msg('e', "Unable to interpret smiles token '$tok'\n");
				return undef;
			}

			#print "i $isotope" if $isotope;
			#print " el $el";
			#print " a $at"     if $at;
			#print " h $h"      if $h;
			#print " s $sign"   if $sign;
			#print " s2 $sign2" if $sign2;
			#print "\n";

			if ($sign) {

				if ($sign eq '+') {
					$fchg = 1;
					$fchg = $fchg * ($sign2 || 1);
				} elsif ($sign eq '-') {
					$fchg = -1;
					$fchg = $fchg * ($sign2 || 1);
				} elsif ($sign eq '++') {
					$fchg = +2;
				} elsif ($sign eq '--') {
					$fchg = -2;
				} elsif ($sign eq '+++') {
					$fchg = +3;
				} elsif ($sign eq '---') {
					$fchg = -3;
				} else {
					silico_msg('e', "Unable to interpret sign in token $tok smiles $smiles\n");
					return undef;
				}
			}

			if ($h) {

				$h =~ s/^H//;
				$h ||= 1;
			}

			if ($at) {
				++$warn_chirality;
			}

			$atomflag = 1;
		}

		if ($tok =~ /^[a-z]/i) {

			($el, $num) = $tok =~ /^([a-z]+)(\d*)/i;
			$el       = ucfirst(lc($el));
			$atomflag = 1;
			$fchg     = 0;
			$h        = 0;
		}

		if ($atomflag) {

			my $atom;

			if (!defined $el) {
				silico_msg('d', "Unable to determine element in token '$tok' in smiles @f\n");
				return undef;
			}

			# Generate atom name
			++$acount->{$el};
			my $name = $el . $acount->{$el};

			#print "	atom: el $el name $name ";
			#print "num $num" if $num;
			#print "\n";

			$x = $x + 1.1;
			$y = $y + $a + (($num && 4 * $a) || 0);
			$a = -$a;

			$atom = mol_add_atom($mol, $name, $el, $x, $y, $z, 'UNK', 1);

			# Atom is aromatic (lower case)
			if ($tok =~ /^[a-z]/) {
				$atom->{AROM} = 1;
			}

			# Create bonds to previous atom
			# Note that any reference to $prevatom->{XXX} before here will
			# cause it to be defined and cause problems
			if (defined $prevatom) {

				$bond = 4 if $atom->{AROM} && $prevatom->{AROM};

				bond_create_atom($mol, $atom, $prevatom, $bond);
			}

			# Create bonds rings
			# Can only create single or aromatic bonds
			if ($num) {
			
				my @nums = split '', $num;
				
				#print "nums @nums\n";
				
				foreach my $n (@nums) {

					if ($hash->{$n}) {

						$bond = 1;
						$bond = 4 if $atom->{AROM} && $hash->{$n}{AROM};

						bond_create_atom($mol, $atom, $hash->{$n}, $bond);
						delete $hash->{$num};

					} else {

						$hash->{$n} = $atom;
					}
				}
			}

			$prevatom = $atom;
			$bond     = 1;

			# Add hydrogens
			if ($h) {
				atom_add_h($atom, $mol, 1.0, $h);
			}

			if ($isotope) {
				$atom->{ISOTOPE} = $isotope;
			}

			next;
		}

		if ($tok eq '(') {

			++$bracketlevel;
			push @stack, $prevatom;
			$z += 1.1;

			#print "	+brack: $bracketlevel\n";
		}

		elsif ($tok eq ')') {
			--$bracketlevel;
			if ($bracketlevel < 0) {
				silico_msg('e', "Too many close brackets ')' in smiles '$smiles'\n");
				return undef;
			}
			$prevatom = pop @stack;
			$z -= 1.1;

			#print "	+brack: $bracketlevel\n";
		}

		elsif ($tok eq '-') {

			$bond = 1;
		}

		elsif ($tok eq '=') {

			$bond = 2;
			#print "	bond: 2\n";
		}

		elsif ($tok eq '#') {

			$bond = 3;
			#print "	bond: 3\n";
		}
		
		elsif ($tok eq '~') {

			# Unspecified bond
			$bond = 8;
			#print "	bond: 8\n";
		}

		elsif ($tok eq '.') {

			#undef $prevatom;
		}

		elsif ($tok eq '/' || $tok eq '\\') {
		
			silico_msg("w", "Ignoring stereochemistry '$tok'\n");
		}

		else {

			silico_msg('e', "Unidentified token '$tok' in smiles string @f");
			return undef;
		}
	}

	#print "-----------------\n\n";
	#molecule_printout($mol, 'all');
	$mol->{WARN_SMILES_CHIRALITY} = 1 if $warn_chirality;

	return $mol;
}

#
# Smiles creation routines
#

sub clear_smiles_sets {

	#<
	#? Clear all markers on atoms identifying them as belonging
	#  to a SMILES set.
	#; Requires: molecule, list of atoms (default all atoms)
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $list = $_[1] || $mol->{ATOMS};

	foreach my $atom (@$list) {

		delete $atom->{SMILES_AROMATIC};
		delete $atom->{SMILES_INVARIANT};
		delete $atom->{SMILES_LABEL};
		delete $atom->{SMILES_ORDER};
		delete $atom->{SMILES_POS};
		delete $atom->{SMILES_RANK};
		delete $atom->{SMILES_RING_JOINS};
		delete $atom->{SMILES_RING_PRIORITY};
		delete $atom->{SMILES_RJ2};
	}
	delete $mol->{SMILES};
}

sub mol_label_smiles_invariants {

	#<
	#? Label all atoms in a molecule subset with their SMILES invariants
	#. This is based on the CANGEN algorithm described in
	#  Weininger, Weininger and Weininger,
	#  J. Chem. Inf. Comput. Sci. 1989, 29, (2), 97-101.
	#; Requires: molecule, optional subset of atoms to use
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $list = $_[1] || $mol->{ATOMS};

	mol_find_neighbouring_atoms($mol, $list);
	foreach my $atom (@$list) { 
		smiles_invariant($atom, $mol) 
	}
}

sub smiles_invariant {

	#<
	#? Generates a SMILES invariant for each atom, so as
	#  to assist in forming Unique SMILES strings.
	#; Requires: atom, molecule, list of atoms in SMILES set
	#; Returns: invariant string
	#>

	my $atom = $_[0];
	my $mol  = $_[1];
	
	my $connected = sprintf("%01d", ch($atom->{CONNECTED}, 'HEAVY'));
	my $numbonds = sprintf("%02d", valence($atom) - ch($atom->{CONNECTED}, 'H'));
	my $elnum = sprintf("%02d", element_symbol2number($atom->{ELEMENT}));
	my $formal_charge = $atom->{FORMAL_CHARGE} || 0;

	my $charge_sign;
	if (!$formal_charge) {
		$charge_sign = 0;
	} else {
		$charge_sign = ($formal_charge > 0) ? 1 : 2;
	}

	my $charge_sign2   = sprintf("%01d", $charge_sign);
	my $formal_charge2 = sprintf("%01d", abs($formal_charge));

	my $connected_h = sprintf("%01d", ch($atom->{CONNECTED}, 'H'));

	my $smiles_invariant = "$connected";
	$smiles_invariant .= "$numbonds";
	$smiles_invariant .= "$elnum";
	$smiles_invariant .= "$charge_sign";
	$smiles_invariant .= "$formal_charge2";
	$smiles_invariant .= "$connected_h";

	$atom->{SMILES_INVARIANT} = $smiles_invariant;
	return $smiles_invariant;
}

sub rank_smiles_set {

	#<
	#? Rank a set of atoms on the basis of their SMILES invariants.
	#. This is based on the CANGEN algorithm described in
	#  Weininger, Weininger and Weininger,
	#  J. Chem. Inf. Comput. Sci. 1989, 29, (2), 97-101.
	#; Requires: molecule, list of atoms that together form SMILES set
	#; Returns: nothing; atoms are marked in $atom->{SMILES_RANK}
	#>

	my $mol  = $_[0];
	my $list = $_[1] || $mol->{ATOMS};
	my $h    = $_[2];                    # Use explicit hydrogens

	mol_label_smiles_invariants($mol, $list);

	my $list2;
	foreach my $atom (@$list) {

		next if ($atom->{ELEMENT} eq 'H' && !$h);

		$atom->{SMILES_RANK} = $atom->{SMILES_INVARIANT};
		@{ $atom->{SMILES_RING_JOINS} } = ();
		$atom->{SMILES_RING_PRIORITY} = 1;

		push @$list2, $atom;
	}

	my @ranks = condense_smiles_ranks($list2);
	my @newranks;

	# Depending on symmetry, we may still have less ranks used than
	# number of atoms. These now have to be resolved via another method,
	# that of doubling all ranks and subtracting 1 from one of the atoms
	# with the lowest tied rank. This is followed by a resolution by
	# neighbouring primes, as above.
	while (1) {

		# While the number of ranks actually used is less than the
		# total number of atoms in the SMILES set, attempt to resolve
		# ties using the Corresponding Primes Method.
		while ($#ranks < $#{$list2}) {
			silico_msg('g', "Ranking SMILES atoms by neighbouring primes...\n");

			# Eventually, if two or more atoms really are equivalent by
			# symmetry, the Corresponding Primes Method will fail. We can
			# tell when this is. On the final iteration, the number of
			# unique ranks will be unchanged.
			@newranks = rank_by_neighbouring_primes($list2, $mol);

			last if $#newranks == $#ranks;

			# If the new ranks are an improvement over the existing ranks,
			# then update them.
			@ranks = @newranks;
		}

		# It is possible that the last attempt to resolve ties using neighbouring
		# primes was completely successful. If that is so, we don't want to do a
		# Double And Tie Break. So we exit here.
		last if $#ranks == $#{$list2};

		# Two or more atoms may still be equivalent by symmetry. If this is the case,
		# we need to resolve ties using Double And Tie Break.
		silico_msg('g', "Ranking SMILES atoms by double and tie break...\n");
		@ranks = double_and_tie_break($list2);

		# It is likely that the Double And Tie Break introduced new differences
		# that need to be resolved using Corresponding Primes. So we want to
		# go back to the start of the WHILE loop.
	}
}

sub double_and_tie_break {

	#<
	#? Attempt to resolve ties in SMILES ranks using the "double
	#  and tie break" method. This multiplies all ranks by 2, and
	#  then reduces the rank of any one of the atoms which share
	#  the lowest tied rank by 1.
	#. This routine has been modified so that the tie is broken in
	#  a consistent way that depends on the atomic coordinates
	#  (if any are present).  This is designed so that similar conformers
	#  will produce a similar smiles order, simplifying superimposition
	#; Requires: list of atoms
	#; Returns: list of ranks
	#>

	my $list = $_[0];

	foreach my $atom (@$list) {

		next if !defined $atom->{SMILES_RANK};

		# Double all relevant SMILES ranks.
		$atom->{SMILES_RANK} *= 2;
	}
	
	# Sort the list of SMILES ranks.
	# Sort on rank and coordinates if they exist to give
	# a consistent representation for similar structures
	@$list = sort {
		return -1 if $a->{SMILES_RANK} < $b->{SMILES_RANK};
		return 1  if $a->{SMILES_RANK} > $b->{SMILES_RANK};
		if (defined $a->{X} && defined $b->{X}) {
			return -1 if int($a->{X} - $b->{X}) > 0.1;
			return 1  if int($a->{X} - $b->{X}) < 0.1;
		}
		if (defined $a->{Y} && defined $b->{Y}) {
			return -1 if int($a->{Y} - $b->{Y}) > 0.1;
			return 1  if int($a->{Y} - $b->{Y}) < 0.1;
		}
		if (defined $a->{Z} && defined $b->{Z}) {
			return -1 if int($a->{Z} - $b->{Z}) > 0;
			return 1  if int($a->{Z} - $b->{Z}) < 0;
		}
		return 0;
	} @$list;

	# Determine the Lowest Tied SMILES Rank.
	# Deduct 1 from this atom's SMILES rank.
	for (my $i = 1 ; $i <= $#{$list} ; ++$i) {

		if ($list->[$i]{SMILES_RANK} == $list->[ $i - 1 ]{SMILES_RANK}) {
			--$list->[ $i - 1 ]{SMILES_RANK};
			last;
		}
	}

	# Redo the relative rankings.
	return condense_smiles_ranks($list);
}

sub rank_by_neighbouring_primes {

	#<
	#? Attempt to resolve ties in SMILES ranks using the product
	#  of the prime numbers corresponding to the ranks of the several
	#  neighbours of each of the atoms in question.
	#. Only atoms which have otherwise equivalent ranks will have their
	#  relative ranks altered by this procedure (though the absolute rank
	#  of any given atom may change as new ranks are created).
	#; Requires: list of atoms, molecule
	#; Returns: list of ranks
	#>

	my $list = $_[0];
	my $mol  = $_[1];

	my $hpr = 0;

	primes() if !$Silico::primes->[0];

	foreach my $atom (@$list) {

		next if !defined $atom->{SMILES_RANK};

		$atom->{PRIMERANK} = 1;

		foreach my $con (connected($atom, $mol)) {
			next if !defined $con->{SMILES_RANK};
			die "We will need to add more prime numbers" if $con->{SMILES_RANK} > 1000;
			$atom->{PRIMERANK} *= $Silico::primes->[ $con->{SMILES_RANK} ];
		}

		$hpr = $atom->{PRIMERANK} if $atom->{PRIMERANK} > $hpr;
	}

	foreach my $atom (@$list) {

		next if !defined $atom->{SMILES_RANK};

		# Make sure no ranking information is lost: add the primerank
		# to the SMILES rank. But before doing this, divide the primerank
		# by the highest primerank plus one (so the amount added will
		# always be less than one, and thus, existing ranks won't change).
		$atom->{PRIMERANK} /= ($hpr + 1);
		$atom->{SMILES_RANK} += $atom->{PRIMERANK};

		# Clean up
		delete $atom->{PRIMERANK};
	}

	# Return ranks to integral values
	return condense_smiles_ranks($list);
}

sub corresponding_prime_old {

	#<
	#? Determine the prime number corresponding to an integer.
	#  That is, for some number N, determine the Nth prime number.
	#; Requires: integer
	#; Returns: integer
	#>

	my $N = $_[0];

	my $i = 1;
	my $n = 0;
	my @p = ();

      	DO: while ($n < $N) {

		++$i;

		foreach (@p) {
			next DO if (sprintf("%.12f", $i / $_) == int($i / $_));
		}

		$p[$n] = $i;
		++$n;
	}

	#print "$N $i\n";

	return $i;
}

sub condense_smiles_ranks {

	#<
	#? Given an input set of SMILES ranks, remove the "air gaps"
	#  and non-integral values from the ranks so atoms are ranked
	#  1, 2, ...
	#; Requires: list of atoms
	#; Returns: list of used ranks
	#>

	my $list = $_[0];

	# Prepare a list of input ranks from the existing SMILES ranks
	# (these may still be initial invariants, but that is not a problem)
	my @input_ranks;
	foreach my $atom (@$list) {

		next if !defined $atom->{SMILES_RANK};
		push @input_ranks, $atom->{SMILES_RANK};
	}

	#print "inp @input_ranks\n";

	# Sort ranks numerically
	my @ranks1 = sort { $a <=> $b } @input_ranks;

	# Remove duplicate entries from the list of ordered ranks
	# Make a hash to relate the new rank to the previous one
	my $i = 0;
	my @ranks;
	my %new_rank_hash;
	foreach my $rank (@ranks1) {
		if (!defined $ranks[-1] || $rank != $ranks[-1]) {
			push @ranks, $rank;
			$new_rank_hash{$rank} = $i + 1;
			++$i;
		}
	}

	# Assign new ranks based on position in the list of ordered ranks
	foreach my $atom (@$list) {

		next if !defined $atom->{SMILES_RANK};
		$atom->{SMILES_RANK} = $new_rank_hash{ $atom->{SMILES_RANK} };
	}

	#print "out @ranks\n\n";

	return @ranks;
}

sub mol_aromatic_smiles {

	#<
	#? Of a list of atoms, mark each one as aromatic if it is.
	#; Requires: molecule, optional list
	#; Returns: nothing
	#>

	my $mol  = $_[0];
	my $list = $_[1];

	foreach my $atom (@$list) {

		my $con_a = 0;

		next if !$atom->{AROMATIC_RING};

		foreach my $con (connected($atom, $mol)) {

			if (find_bondorder2($atom, $con, $mol) == 4 || find_bondorder2($atom, $con, $mol) == 7) {
				++$con_a;
			}
		}

		if ($con_a >= 2) {
			$atom->{SMILES_AROMATIC} = 1;
		}
	}

	#molecule_printout($mol, 'all');
}

sub ens_smiles_string {

	#<
	#? Calculate smiles string
	#; Requires: molecule, flag to use ChemAxon molconvert to generate smiles
	#; Returns: nothing
	#; Sets: mol->{SDF_DATA}{SMILES}, mol->{SDF_DATA}{SMILES_SOURCE}
	#>

	my $mols = $_[0];
	my $c    = $_[1];

	if ($c) {
		ens_smiles_string_chemaxon($mols);
		return;
	}

	foreach my $mol (@$mols) {

		mol_smiles_string($mol);
		$mol->{SDF_DATA}{SMILES}        = $mol->{SMILES};
		$mol->{SDF_DATA}{SMILES_SOURCE} = 'Silico';
	}
}

sub ens_smiles_string_chemaxon {

	#<
	#? Calculate canonical smiles string using the Chemaxon program 'molconvert'
	#; Requires: molecule
	#; Returns: nothing
	#; Sets: mol->{SDF_DATA}{SMILES},  mol->{SDF_DATA}{SMILES_SOURCE}
	#>

	my $mols = $_[0];

	my $tempname = get_tempname();

	write_sdf($mols, "$tempname.sdf");

	my $molconvert = `which molconvert`;
	chomp $molconvert;

	if (!-x $molconvert) {
		silico_msg('d', "Program molconvert not found: 'which molconvert' returns $molconvert\n");
	}

	my @smileses = `$molconvert smiles:u $tempname.sdf`;

	foreach my $mol (@$mols) {

		$mol->{SDF_DATA}{SMILES} = shift @smileses;
		chomp $mol->{SDF_DATA}{SMILES};
		$mol->{SDF_DATA}{SMILES_SOURCE} = 'ChemAxon';
	}

	system("rm -f $tempname") if !$Silico::debug;
}

sub mol_smiles_string {

	#<
	#? Create a SMILES string from a molecule.
	#; Requires: molecule, atom list, flags to generate explicit bondorders, explicit hydrogens 
	#  (1 = all atoms, 2 = only where present), Kekule bonds, to mark missing hydrogens in molecule
	#. Sets $mol->{SMILES}
	#; Returns: string
	#>
	
	use vars qw($scount $smiles_b $smiles_h $smiles_k %smiles_organic_element_hash);

	my $mol = $_[0] || croak();
	my $list = $_[1] || $mol->{ATOMS};
	$smiles_b = $_[2];    # Globals
	$smiles_h = $_[3];
	$smiles_k = $_[4];
	my $a = $_[5];

	my $smiles;

	$scount = 0;          # Note use of global variable
	
	# List of 'organic' elements including the unspecified elements Q and A and normal valence
	%smiles_organic_element_hash = qw (B 3 C 4 N 3 O 2 P 5 S 2 F 1 Cl 1 Br 1 I 1 A 4 Q 4); # Note use of global variable

	if ($#{$list} == -1) {
		silico_msg('w', "Molecule " . ($mol->{NAME} || 'Unk') . " has no atoms\n");
		return '';
	}

	# Reinitialise if smiles has previously been calculated
	if ($mol->{SMILES}) {
		clear_smiles_sets($mol, $list);
	}

	# find rings
	molecule_find_rings($mol, $#{$list} + 1, 1);

	if ($smiles_k) {
		convert_aromatic_bonds_kekule($mol);
	} else {
		make_aromatic_bonds($mol);
	}

	# Calculate formal charges
	molecule_formal_charge($mol, $list) if !$mol->{HAS_FORMAL_CHARGES};

	# Determine whether each atom is to be aromatic or standard.
	mol_aromatic_smiles($mol, $list);

	rank_smiles_set($mol, $list);

	my $end_rank = 0;
	my $string   = '';

	while (1) {

		my $start_atom;
		my $start_rank = 99999999999;

		# Determine the starting atom: That with the first SMILES
		# rank (that has not already been used).
		foreach my $atom (@$list) {

			# Skip atoms without a rank
			next if !defined $atom->{SMILES_RANK};
			next if $atom->{SMILES_RANK} <= $end_rank;
			next if $atom->{SMILES_RANK} > $start_rank;

			$start_rank = $atom->{SMILES_RANK};
			$start_atom = $atom;
		}

		last if !defined $start_atom;

		$scount = 0;    # Global var

		@$smiles = smiles_list_recursive($start_atom, undef, $mol);

		$string .= "." if $string;
		$string .= make_smiles_string_from_list(@$smiles);

		# Establish relationship between each atom and position in smiles string
		my $i = 0;
		foreach my $atom (@$smiles) {

			next if (ref $atom ne "HASH");

			# Find the highest ranked atom in current smiles (end_rank)
			$end_rank = $atom->{SMILES_RANK} if $atom->{SMILES_RANK} > $end_rank;

			# Label each atom with their position in the smiles string: SMILES_POS
			# Hydrogens have to be added later following their heavy atoms
			$atom->{SMILES_POS} = $i;
			++$i;
			foreach my $con (connected($atom, $mol)) {
				next if $con->{ELEMENT} ne 'H';
				$con->{SMILES_POS} = $i;
				++$i;
			}
		}
	}

	if ($string eq '' && $mol->{NUMATOMS} != 0) {
		silico_msg(
			'e',
			"Unable to determine starting atom for molecule: '" . ($mol->{NAME} || "") . "' \n",
			"Perhaps SMILES ranks could not be generated.\n"
		);
	}

	$mol->{SMILES} = $string;

	return $string;
}

sub make_smiles_string_from_list {

	#<
	#? Generate, from a list of atoms and tokens, a SMILES string.
	#; Requires: List
	#; Returns: string
	#>

	my @smiles = @_;

	my $count  = 0;
	my $string = '';
	
	#print "smiles @smiles\n";

	foreach my $token (@smiles) {

		# If the token is a pointer to a hash, it is an atom
		# (at least, it better be one...)

		if (ref($token) eq 'HASH') {

			my $i;

			# Add the atom's SMILES label to the SMILES string.
			$string .= $token->{SMILES_LABEL};

			#print "t $token->{SMILES_LABEL}\n";

			# Have a look at ring joins. For each offset into the ring
			# joins array for this atom...
			for ($i = 0 ; $i <= $#{ $token->{SMILES_RING_JOINS} } ; ++$i) {

				my $rj;

				# load the atom at this offset
				$rj = $token->{SMILES_RING_JOINS}[$i][0];

				# if the join point occurs later in the SMILES string
				# than the current token...
				if ($rj->{SMILES_ORDER} > $token->{SMILES_ORDER}) {

					my $j;

					# add a new join point to the array of join point numbers
					# associated with the current token
					$token->{SMILES_RJ2}[$i] = ++$count;

					# for each join point on the connected atom...
					for ($j = 0 ; $j <= $#{ $rj->{SMILES_RING_JOINS} } ; ++$j) {

						my $rj2;

						$rj2 = $rj->{SMILES_RING_JOINS}[$j][0];
						next if $rj2 != $token;

						my $b = $rj->{SMILES_RING_JOINS}[$j][1];

						# add the current count to the appropriate offset
						# for the current token.
						$rj->{SMILES_RJ2}[$j] = $b;
						$rj->{SMILES_RJ2}[$j] .= '%' if $count > 9;
						$rj->{SMILES_RJ2}[$j] .= $count;
					}

				} elsif ($rj->{SMILES_ORDER} == $token->{SMILES_ORDER}) {
					silico_msg('d', "An atom supposedly completes a ring with a bond to itself!\n");
				}

				# If the join point occurs earlier in the SMILES string than
				# the current token, we don't have to do anything. It should
				# be taken care of.
			}

			# Sort the list of join points
			if (defined $token->{SMILES_RJ2}[0]) {
				@{ $token->{SMILES_RJ2} } = sort rj2sort @{ $token->{SMILES_RJ2} };
				foreach (@{ $token->{SMILES_RJ2} }) {
					$string .= $_;
				}
			}

		} else {
			$string .= $token;
		}
	}

	#print "s $string\n";
	return $string;
}

sub rj2sort {

	# Sort only on numeric values present in string
	my $a1 = $a =~ m/\d*/;
	my $b1 = $b =~ m/\d*/;

	return $a1 <=> $b1
}

sub smiles_list_recursive {

	#<
	#? Determine the SMILES string representation for a given atom.
	#; Requires: atom, molecule, ...
	#; Returns: string
	#>

	use vars qw($smiles_b $smiles_h $smiles_k);

	my $atom   = $_[0];
	my $parent = $_[1] || 0;
	my $mol    = $_[2];

	#my $smiles_b = $_[3];	# Use explicit bonds
	#my $smiles_h = $_[4];	# Use explicit hydrogens
	#my $smiles_k = $_[5];	# Use Kekule bonds and non-aromatic atom symbols

	my @array;
	my $bracket = 0;
	my $con_h = 0;
	my @paths;

	use vars qw($scount);

	if (!defined $atom) {
		silico_msg('e', "Starting atom is not defined!\n");
		return undef;
	}

	# Mark this atom as having been travelled to.
	push @array, $atom;
	$atom->{SMILES_ORDER} = ++$scount;

	atom_smiles_label($atom, $parent, $mol);

	# Provide information about connected atoms: whether they share
	# a ring with the current atom, and if so, by what order bond.
	# Also, push connected atoms (apart from parents and those not
	# forming part of the SMILES set) into a list.
	my @connected;
	foreach my $con (connected($atom, $mol)) {

		next if !defined $con->{SMILES_RANK};
		next if $con == $parent;

		my $conring;
		my $ringflag;
		my $flag;

		foreach my $ring (@{ $atom->{RINGS} }) {
			foreach $conring (@{ $con->{RINGS} }) {
				$ringflag = 1 if $ring == $conring;
			}
		}

		if ($ringflag) {
			my $border = find_bondorder2($atom, $con, $mol);
			if ($border == 2 || $border == 3) {
				$flag = 1;
			}
		}
		$con->{SMILES_RING_PRIORITY} = 2 if !$flag;

		push @connected, $con;
	}

	# Sort connected atoms by their SMILES ranks
	@connected = sort branchsort @connected;

        DO: for (my $i = 0 ; $i <= $#connected ; ++$i) {

		my $bstring;
		my $con;
		my $connum2;
		my $path;
		my $rjc;

		$con = $connected[$i];

		# Skip connected hydrogens unless the present atom
		# is also hydrogen (i.e., dihydrogen) or the connected
		# H functions as a bridge (has more than 1 connection
		# of its own)
		if ($con->{ELEMENT} eq 'H') {
			next DO unless ($atom->{ELEMENT} eq 'H' || $#{ $con->{CONNECT} } > 0);
		}

		# Add the next atom out

		my $border = find_bondorder2($atom, $con, $mol);

		if ($border == 1) {
			$bstring = "-";
		} elsif ($border == 2) {
			$bstring = "=";
		} elsif ($border == 3) {
			$bstring = "#";
		} elsif ($border == 4 || $border == 7) {
			$bstring = ":";
		} elsif ($border == 8) {
			$bstring = "~";
		} else {
			$bstring = "-";
		}

		# Add a ring join point (array of atom and bondorder) and move on if
		# $con already has an atom string.
		if (defined $con->{SMILES_ORDER}) {
			if ($con->{SMILES_ORDER} < $atom->{SMILES_ORDER}) {

				# Ring join bond.  Only required if:
				# 1. Explicit bond orders
				# 2. Ring join is not single or aromatic
				# Otherwise squash it to nothing
				$bstring = '' if ($border == 1 || $border == 4);

				my ($j1, $j2);
				@$j1 = ($atom, $bstring);
				@$j2 = ($con,  $bstring);
				push @{ $con->{SMILES_RING_JOINS} },  $j1;
				push @{ $atom->{SMILES_RING_JOINS} }, $j2;
			}
			next DO;
		}

		# Under what circumstances do I want to definitely give the bond
		# order?
		# 1. If explicit bond orders have been requested
		# 2. If bonded atoms are aromatic and are connected by
		#    other than an aromatic bond
		# 3. If bonded atoms are not aromatic and are connected by
		#    other than a single bond

		if (
			$smiles_b
			|| (       $atom->{SMILES_AROMATIC}
				&& $con->{SMILES_AROMATIC}
				&& !$smiles_k
				&& $border != 4
				&& $border != 7)
			|| ((!$atom->{SMILES_AROMATIC} || !$con->{SMILES_AROMATIC}) && $border != 1)
		  )
		{
			push @$path, $bstring;
		}
		push @$path, smiles_list_recursive($con, $atom, $mol);

		push @paths, $path;
	}

	for (my $i = 0 ; $i < $#paths ; ++$i) {
		unshift @{ $paths[$i] }, "(";
		push @{ $paths[$i] }, ")";
	}

	foreach my $a (@paths) {
		foreach (@$a) {
			push @array, $_;
		}
	}

	@{ $atom->{SMILES_RING_JOINS} } = sort rjsort @{ $atom->{SMILES_RING_JOINS} };

	return @array;
}

sub atom_smiles_label {

	#<
	#? Provide an atom with additional data in the SMILES string
	#  (square brackets, attached hydrogen atoms, formal charge, etc.)
	#; Requires: atom
	#; Returns: nothing
	#>

	use vars qw($smiles_h %smiles_k %smiles_organid_element_hash);

	my $atom   = $_[0];
	my $parent = $_[1];
	my $mol    = $_[2];

	#my $k = $_[4];

	my $con_h     = 0;
	my $con_heavy = 0;
	my $h;
	my $open_square_bracket;

	foreach my $con (connected($atom, $mol)) {
		next if $con == $parent;
		++$con_h if $con->{ELEMENT} eq 'H';
		++$con_heavy if $con->{ELEMENT} ne 'H' && !defined $con->{SMILES_RANK};
	}

	# Atom flag to add '!H' to atoms
	$con_heavy = $atom->{CONNECTED_NON_H} if $atom->{CONNECTED_NON_H};

	# We have specified all hydrogens to be included or we have elected to label explicit hydrogens
	if ($smiles_h) {
		$h = 1 if $smiles_h == 1;
		$h = 1 if $smiles_h == 2 && $con_h > 0;
	}

	# If we have connected non-SMILES heavy atoms, use hydrogens for greater certainty.
	$h = 1 if ($con_heavy);

	# If we have a formal charge, explicit hydrogens must be used.
	$h = 1 if ($atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE} != 0);

	# If the atom is a nitrogen, in an aromatic ring and has
	# hydrogens attached (e.g., pyrrole, imidazole), explicit
	# hydrogens must be used.
	$h = 1 if ($atom->{SMILES_AROMATIC} && $atom->{ELEMENT} eq 'N' && $con_h > 0);
	
	if (my $v = $smiles_organic_element_hash{$atom->{ELEMENT}}) {
		
		# If element is an 'organic' element with larger than normal valence
		$h = 1 if valence($atom) > $v;
		
	} else {
	
		# If element is not an 'organic' element (CHONSP-FClBrI + Q and A)
		$h = 1;
	}

	if ($h || $atom->{ELEMENT} eq 'H') {
		$atom->{SMILES_LABEL} .= "[";
		$open_square_bracket = 1;
	}

	# Add the atom element to the SMILES string. Use upper-case letters
	# unless the atom is within an aromatic ring. In that case, lower-case.
	# Also use upper-case letters if Kekule bonds have been requested.
	if ($atom->{SMILES_AROMATIC} && !$smiles_k) {
		$atom->{SMILES_LABEL} .= lc $atom->{ELEMENT};
	} else {
		$atom->{SMILES_LABEL} .= $atom->{ELEMENT};

	}

	# If the atom in question is anything but a hydrogen, we want
	# to describe hydrogens as a suffix inside square brackets.
	# If, however, the molecule in question is dihydrogen, it is
	# instead described as [H][H], so we need to ignore that.
	if ($h && $atom->{ELEMENT} ne 'H' && $con_h > 0) {
		$atom->{SMILES_LABEL} .= "H";
		$atom->{SMILES_LABEL} .= "$con_h" if $con_h > 1;
	}

	# Add symbols for non-SMILES heavy atoms.
	if ($con_heavy) {
		$atom->{SMILES_LABEL} .= "!H";
		$atom->{SMILES_LABEL} .= $con_heavy if $con_heavy > 1;
	}

	# Add the formal charge string
	if ($atom->{FORMAL_CHARGE}) {
		
		if ($atom->{FORMAL_CHARGE} > 0) {
			$atom->{SMILES_LABEL} .= ("+" x abs($atom->{FORMAL_CHARGE}));
		} elsif ($atom->{FORMAL_CHARGE} < 0) {
			$atom->{SMILES_LABEL} .= ("-" x abs($atom->{FORMAL_CHARGE}));
		}
	}

	# Add a closing square bracket
	$atom->{SMILES_LABEL} .= "]" if ($open_square_bracket);
}

sub branchsort {

	#<
	#? Sort connected atoms by their ring priorities then by their SMILES ranks	
	#>

	return 1  if $a->{SMILES_RING_PRIORITY} > $b->{SMILES_RING_PRIORITY};
	return -1 if $a->{SMILES_RING_PRIORITY} < $b->{SMILES_RING_PRIORITY};
	return 1  if $a->{SMILES_RANK} > $b->{SMILES_RANK};
	return -1 if $a->{SMILES_RANK} < $b->{SMILES_RANK};

	return 0;
}

#sub rjsort {
#
#	#<
#	#? Sort connected atoms by SMILES_ORDER
#	#>
#
#	return $a->[0]{SMILES_ORDER} <=> $b->[0]{SMILES_ORDER};
#}
sub rjsort {

        return 1 if $a->[0]{SMILES_ORDER} > $b->[0]{SMILES_ORDER};
        return -1 if $a->[0]{SMILES_ORDER} < $b->[0]{SMILES_ORDER};
        return 0;
}

sub mol_sort_smiles {

	#<
	#? Sort the atoms in a molecule into a canonical order based on the smiles representation
	#. Runs mol_smiles_string if required;
	#; Requires: molecule
	#; Returns: canonical smiles string
	#; Modifies: sorted molecule
	#>

	my $mol = $_[0];

	my $smiles = mol_smiles_string($mol) if !defined $mol->{SMILES};
	molecule_reorder($mol, "SMILES_POS");
	
	# Rings are now bad
	$mol->{BAD_RINGS} = 1;

	return $smiles;
}

sub write_smiles {

	#<
	#? Write an ensemble as a smiles file
	#. Basic implementation. Default is tab delimited
	#; Requires: ensemble or molecule, output filename, options, $b, $b, $k - smiles string options
	#; Options: QUIET - be quiet, 
	#. NOHEADER - dont write a header, 
	#. NODATA - don't write out molecule data
	#. DTAB - tab delimited (default)
	#. DCOMMA - comma delimited
	#. DSPACR - space delimieted
	#>

	my $mols    = ens($_[0]);
	my $outfile = $_[1];
	my $options = uc($_[2] || '');
	my $b       = $_[3];
	my $h       = $_[4];
	my $k       = $_[5];

	silico_msg('d', "Argument 0 is not an ensemble of molecules\n") if ref($_[0] ne 'ARRAY');
	silico_msg('d', "Incorrect arguments\n") if $_[3];
	
	my $warn = 0;

	$outfile = get_ofilebase($mols) . ".smi" if !$outfile;

	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	my $quiet = 1 if $options =~ /\bQUIET\b/;
	my $noheader = 1 if $options =~ /\bNOHEADER\b/;
	my $nodata = 1 if $options =~ /\bNODATA\b/;
	
	my $delimiter = "\t";
	$delimiter = "\t" if $options =~ /\bDTAB\b/;
	$delimiter = " " if $options =~ /\bDSPACE\b/;
	$delimiter = "," if $options =~ /\bDCOMMA\b/;

	my $fr = open_molfile(undef, $outfile, undef, 'write', $options);
	return undef if !$fr;
	return undef if $fr->{ERROR};
	my $FH = $fr->{FH};

	silico_msg('c', "Writing smiles file $outfile" . ($options ? " using $options" : '') . "\n") if !$quiet;

	# Find SDF data fields in all molecules
	my $headers;
	foreach my $mol (@$mols) {

		foreach my $key (keys %{ $mol->{SDF_DATA} }) {

			my $type = ref($mol->{SDF_DATA}{$key}) || "SCALAR";
			if ($type ne 'SCALAR') {
				
				++$warn;
			} else {
				++$headers->{$key};
			}
		}
	}
	my @headerlist = sort(keys %$headers);
	
	silico_msg('w', "Array and hash data is not supported. $warn fields were dropped\n") if $warn;
	
	# Print header
	if (!$noheader) {
	
		print $FH "Smiles".$delimiter."Name";
		
		if (!$nodata) {
			foreach my $key (@headerlist) {
			
				foreach ('smiles', 'SMILES', 'Smiles') {
					next if $key eq $_;
				}
			
				print $FH $delimiter.$key;
			}
			print $FH "\n";
		}
	}

	silico_msg('c', "\n") if !$quiet;
	
	# Print smiles and SDF_DATA values
	my $i = 0;
	foreach my $mol (@$mols) {

		++$i;

		my $smiles_string;
		if (!$mol->{SMILES}) {
			$smiles_string = mol_smiles_string($mol, undef, $b, $h, $k);
			$mol->{SDF_DATA}{SMILES_SILICO} = $smiles_string;
			$mol->{SDF_DATA}{SMILES}        = $smiles_string;
		} else {
			$smiles_string = $mol->{SMILES};
		}

		if (!$nodata) {
		
			my $name = $mol->{NAME} || '""';
			print $FH $smiles_string.$delimiter.$name;
			foreach my $key (@headerlist) {
			
				foreach ('smiles', 'SMILES', 'Smiles') {
					next if $key eq $_;
				}

				my $val = csv_fix($mol->{SDF_DATA}{$key});
				print $FH $delimiter.$val;
			}
			print $FH "\n";
		}
	}
}

sub csv_fix {
	
	#<
	#? Fix delimiters in value taken from CSV file
	#. See http://tools.ietf.org/html/rfc4180
	#>

	my $val = $_[0];

	my $quote;
	my $newval = $val;

	if (!defined $val) {
		$val = '';
	}

	if ($val eq '') {
		$newval = '';
		$quote  = 1;
	}

	if ($val =~ /"/) {
		$val =~ s/"/""/g;
		$quote = 1;
	}

	if ($val =~ /\n/) {
		$quote = 1;
	}

	if ($val =~ /,/) {
		$quote = 1;
	}

	if ($quote) {
		$newval = '"' . $newval . '"';
	}

	return $newval;
}

return 1;

