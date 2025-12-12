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
#! silico_asl.pm
#? Silico atom selection routines (atom specifier)
#. $Revision: 1.15-39-g4f8dfe8 $
#>

use strict;

package Silico;

###################
# ATOM SPECIFIERS #
###################

sub atom_specifier_setup {

	#<
	#? Set up flags for simple atom specifiers
	#; Requires: nothing
	#; Returns: list of atom specifier names, list of specifiers;
	#>

	die "atom_specifier_setup does not take any arguments\n" if $_[0];

	my $as;
	@$as = ();
	my $names;
	@$names = ();

	my $ca       = make_flag('ca',    'use-CA-atoms',         "Select protein backbone atoms (CA)");
	my $cacb     = make_flag('cacb',  'use-CACB-atoms',       "Select prot;ein backbone atoms (CA, CB)");
	my $backbone = make_flag('back',  'use-backbone-atoms',   "Select protein backbone atoms (C, CA, N, O)");
	my $heavy    = make_flag('heavy', 'heavy-atoms',          "Select heavy atoms");
	my $sas      = make_flag('sas',   'short-atom-specifier', "Short specifier [ca, cacb, back, back_disulfide, heavy]", 1);
	my $file     = make_flag('f',     'use-as-file',          "Read atom specifiers from a file", 1, undef);
	my $as_input = make_flag('as',    'atom-specifier',       "Atom specifier string", 1, undef, undef, undef, 1); # Double definition allowed - returns pointer to array

	$sas ||= '';

	if ($file) {

		($as, $names) = read_atom_specifier_file($file);

	} elsif ($as_input) {

		@$as = (@$as, @$as_input);
		@$names = (@$names,  @$as_input);

	} elsif ($backbone || ($sas eq 'back')) {
		push @$as,    'ANAME:C,CA,N,O';
		push @$names, "back";

	} elsif ($ca || $sas eq 'ca') {
		push @$as,    'ANAME:CA';
		push @$names, "ca";

	} elsif ($cacb || $sas eq 'cacb') {
		push @$as,    'ANAME:CA,CB';
		push @$names, "cacb";

	} elsif ($heavy || $sas eq 'heavy') {
		push @$as,    'ELEMENT:!H';
		push @$names, "heavy";

	} elsif ($sas eq 'back_disulfide') {
		push @$as,    'ANAME:C,CA,N,O|SUBNAME:CYS,ANAME:CB,SG';
		push @$names, "back_disulfide";
	}

	#silico_msg('c', "Using atom specifier(s): @$as\n") if defined $as->[0];

	return ($names, $as);
}

sub read_atom_specifier_file {

	#<
	#? Read an atom specifier file
	#; Requires: filename
	#; Returns: array of atom specifiers
	#! Simple atom specifier examples:
	#; ANAME:CA Returns all atoms called 'CA'.
	#; ANAME:CA,CB,CG Returns all atoms called 'CA' or 'CB' or 'CD'.
	#; ANAME:CA,RESNAME:TRP Returns all atoms called 'CA' in all residues called TRP
	#; ANAME:CA,SUBID:4 All atoms called 'CA' in residue number 4
	#. In the examples above, atoms are returned in the order found in the file.
	#  Atom specifiers are case sensitive.
	#!3 Ordered atom specifier examples:
	#. Successive atom specifications can be made. Each is separated by a '|'
	#; ANAME:CB|ANAME:CA|ANAME:CD Returns all atoms called 'CB', 'CA', 'CD' in that order.
	#; ANAME:CA,RESNUM:1|ANAME:CA,RESNUM:4 Returns 'CA' atoms from residue 1 and 4 in that order.
	#!3 Specifier file:
	#. Contains two sets of atom specifiers.  Each set starts with a 'Molecule' line
	#  which can contain a specifier name (space delimited)
	#. eg:
	#; Molecule Backbone
	#; ANAME:C 
	#; ANAME:CA
	#; ANAME: N
	#; Molecule Sidechain
	#; ANAME:CB
	#; ANAME:CG
	#; ANAME:CD
	#>

	my $infile = $_[0];

	my $specifiers;
	my $specnames;

	# Use open_file to also handle compressed files
	if (!open_file(*INFILE, $infile)) {
		file_read_error($infile);
		return undef;
	}

	my $molcount = -1;
	while (my $line = <INFILE>) {

		next if !($line =~ /\w/);
		next if $line =~ /^\s*#/;
		$line =~ s/\s//g;

		chomp $line;

		if ($line =~ /^Molecule/) {

			++$molcount;
			$specifiers->[$molcount] = '';
			my ($v) = $line =~ /^Molecule\s*(\w*)/;
			$specnames->[$molcount] = $v;
			print "a $v\n";
			next;
		}

		$specifiers->[$molcount] .= $line . "|" if $molcount >= 0;
	}

	my $i  = 0;
	foreach (@$specifiers) {
		s/\|$//;
		$specnames->[$i] = $_ if !defined $specnames->[$i];
		++$i;
	}

	silico_msg('d', "No 'Molecule' line found in specifier\n") if $molcount < 0;

	close INFILE;

	return $specifiers, $specnames;
}

sub get_atoms_using_specifier {

	#<
	#? Parse an asl string and return a list of atoms
	#; Requires: molecule, string
	#; Returns: pointer to array of atoms
	#>

	my $mol    = $_[0];
	my $string = $_[1];

	my $array = get_atomarray_using_specifier($mol, $string);

	my $atoms;

	my $i = 0;
	foreach my $atom (atoms($mol)) {

		push @$atoms, $atom if $array->[$i];
		++$i;
	}

	return $atoms;
}

sub get_atomnums_using_specifier {

	#<
	#? Parse an asl string and return a list of atom numnbers
	#; Requires: molecule, string
	#; Returns: pointer to array of atoms
	#>

	my $mol    = $_[0];
	my $string = $_[1];

	my $array = get_atomarray_using_specifier($mol, $string);

	my $atomnums;
	my $i = 0;
	foreach (@$array) {

		push (@$atomnums, $i) if $_;
		++$i;
	}

	return $atomnums;
}

sub get_atomarray_using_specifier {

	#<
	#? Parse an asl string and return an array indicating selected atoms
	#  and key values
	#. Allowed specifier keys: ANUM: ANAME: CHARGE: ELEMENT: NAME:
	#  RESNAME: SEGID: SUBID: TYPE: X: Y: Z:
	#; Requires: string
	#; Returns: array of atom numbers (indexed from zero)
	#>

	my $mol    = $_[0];
	my $string = $_[1];

	# Translate some types
	$string =~ s/\bANAME\b/NAME/g;
	$string =~ s/\bRESNAME\b/SUBNAME/g;
	$string =~ s/\bANUM\b/NUM/g;
	$string =~ s/\bRESNUM\b/SUBID/g;

	# Split string by '&' symbols

	my @andlist = split '&', $string;

	my $atomlist1;
	foreach my $i (0 .. $mol->{NUMATOMS} - 1) {
		$atomlist1->[$i] = 1;
	}

	foreach my $s1 (@andlist) {

		# Split string by '|' symbols

		my @orlist = split '\|', $s1;
		
		#print "s1 $s1\n";

		my $atomlist2;
		@$atomlist2 = ();
		foreach my $s2 (@orlist) {
		
			#print "s2 $s2\n";

			my $specifiers = asl_substring_to_specifiers($s2);
			foreach my $specifier (@$specifiers) {
				$atomlist2 = get_atomnums_spechash($mol, $specifier, $atomlist2);
			}
		}

		# AND lists together
		
		foreach my $i (0 .. $mol->{NUMATOMS} - 1) {
			$atomlist1->[$i] = $atomlist1->[$i] && $atomlist2->[$i];
		}
	}

	return $atomlist1;
}

sub asl_substring_to_specifiers {

	#<
	#? Parse an atom specifier string and return a hashes containing keys and key values
	#; Requires: string
	#; Returns: array of hashes
	#>

	my $string = $_[0];

	silico_msg('d', "No specifier string provided") if !defined $string;

	my @spectypes = qw(ANUM ANAME CHARGE CHAIN ELEMENT NAME NUM RESNAME RESNUM
	  SEGID SUBID SUBNAME TYPE X Y Z);

	# Separate tokens by spaces
	$string =~ s/:/: /g;
	$string =~ s/,/ /g;

	silico_msg('g', "Processing string: '$string'\n");

	foreach my $spectype (@spectypes) {
		silico_msg('d', "Syntax error: Missing colon from specifier type '$spectype:' in ASL '$string'\n") if $string =~ /\b$spectype[^:]/;
	}

	my @tokens = split " ", $string;

	my $linehash;
	my $specifiers;

      LINE: while (1) {

		last if !defined $tokens[0];

	      TYPE: while (1) {

			silico_msg('g', "tokens1: @tokens\n");

			my $tok = shift @tokens;
			last if !defined $tok;

			silico_msg('g', "tok1: $tok \n");
			last if $tok eq '|';

			my $flag = 0;
			foreach my $spectype (@spectypes) {

				next if $tok ne "$spectype:";
				$flag = 1;
				last;
			}

			# Does not match any spectype
			if (!$flag) {
				silico_msg('d', "Token '$tok' in '$string' does not match any specifier type\n");
				return undef;
			}

			my $type = $tok;
			$type =~ s/:$//;

			silico_msg('g', "type: $type\n");

			# Put remaining tokens into hash
			while (1) {

				silico_msg('g', "tokens2: @tokens\n");

				$tok = shift @tokens;
				last if !defined $tok;

				# Find if we have reached a new specifier
				if ($tok =~ /:$/) {

					# Put token back
					unshift @tokens, $tok;
					next TYPE;
				}

				silico_msg('g', "tok2: $tok\n");
				last TYPE if $tok eq '|';

				# Split up ranges of numbers
				if (($type eq 'ANUM' || $type eq 'SUBID') && $tok =~ /^[0-9]+-[0-9]+$/) {

					my @toks = split '-', $tok;

					if ($toks[1] < $toks[0]) {
						silico_msg('d', "Error in range '$tok'\n");
					} elsif ($toks[1] == $toks[0]) {
						push @{ $linehash->{$type} }, $toks[0];
					} else {
						my $i = $toks[0];
						while ($i <= $toks[1]) {
							push @{ $linehash->{$type} }, $i;
							++$i;
						}
					}
				} else {

					# Single value
					push @{ $linehash->{$type} }, $tok;
				}
			}
		}

		# Enumerate specifiers
		my $specs;
		my $typecount = 0;

		foreach my $type (keys %$linehash) {

			if ($typecount == 0) {

				foreach my $val (@{ $linehash->{$type} }) {

					my $spec;
					$spec->{$type} = $val;
					push @$specs, $spec;
				}

			} else {

				my $newspecs;
				foreach my $val (@{ $linehash->{$type} }) {

					foreach my $spec (@$specs) {

						my $newspec = deep_copy($spec);
						$newspec->{$type} = $val;
						push @$newspecs, $newspec;
					}
				}
				$specs = $newspecs;
			}
			++$typecount;
		}

		# Add specs to final list of specifiers
		foreach my $spec (@$specs) {

			push @$specifiers, $spec;

			#foreach (keys %$spec) {
			#	print "spec: $spec key: $_ val: $spec->{$_}\n";
			#}
		}
	}

	return $specifiers;
}

sub get_atomnums_spechash {

	#<
	#? Use a Silico atom specifier string to select sets of atoms from a molecule
	#; Requires: molecule, specifier array
	#; Returns: pointer to array of atom numbers (indexed from zero)
	#>

	my $mol      = $_[0];
	my $spechash = $_[1];
	my $atomlist = $_[2] || silico_msg('d', "Subroutine requres atomlist\n");

	my $i = -1;
      ATOM: foreach my $atom (atoms($mol)) {

		++$i;

		# Atom must match ALL keys
	      KEY: foreach my $key (keys %$spechash) {
	      
	      		#print "key $key val $spechash->{$key}\n";

			my $data = $atom->{$key};
			my $spec = $spechash->{$key};

			if (!defined $data) {
				silico_msg('e', "Unable to get data for key \"$key\"\n");
				return undef;
			}

			$data =~ s/ //g;    # Remove spaces

			# Deal with negation
			if ($spec =~ /^!/) {
				my $specmod = $spec;
				$specmod =~ s/^!//;
				next ATOM if $data =~ /^$specmod$/;

			} else {

				next ATOM if !($data =~ /^$spec$/);
			}
		}

		++$atomlist->[$i];
	}

	return $atomlist;
}

sub write_asl_atoms {
        
        #< 
        #? Write out atoms selected in index groups to a file
        #; Requires: molecule, list of index groups, filebase
        #; Returns: nothing
        #>      
                
        my $mol = $_[0];
        my $list = $_[1];
        my $filebase = $_[2] || get_filebase_short($mol->{SOURCE_FILE_NAME});
        
        my $ens;

	print "list: $list\n";
	carp();
        
        foreach my $key (sort keys %$list) {

                next if $key eq 'All';

                my $hash;

                foreach my $val (@{$list->{$key}}) {
                        ++$hash->{$val};
                }

                my $newmol = molecule_copy_noatoms($mol);

                foreach my $atom (atoms($mol)) {
                        if ($hash->{$atom->{NUM}}) {

                                my $at =  deep_copy($atom);

                                # Remove all bonds so that we dont have bonds to nonexistant atoms
                                delete $at->{CONNECT};
                                delete $at->{BORDERS};

                                push @{$newmol->{ATOMS}}, $at;
                        } else {
                                push @{$newmol->{ATOMS}}, undef;
                        }
                }

                molecule_pack($newmol, 1);
                $newmol->{NAME}.=" Index group: $key";
                push @$ens, $newmol;
        }

        if (defined $ens) {
                write_mol_any($ens, "$filebase\_ndx", 'pdb', 'NOCONECT');
        } else {
                silico_msg('w', "Nothing written for -write option because no groups were selected\n");
        }
}

return 1;

