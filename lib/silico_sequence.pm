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
#! silico_sequence.pm
#? Routines to deal with protein/dna sequences
#. Sequences are stored as hashes in the same manner as Silico molecules.
#
#. Sequence details
#, $seq->{NAME}		sequence name
#, $seq->{NUMRESGAP}	number of residues in sequence including gaps
#, $seq->{NUMRES}	number of residues in sequence excluding gaps
#, $seq->{SEQ}[]	array of residues
#, $seq->{ENS}		associated 3D molecules (pointer to ensemble)
#, $seq->{TYPE}		type of information stored in sequence record
#, $seq->{SOURCE}	format of data source
#
#. Standard residue details ($seq->{TYPE} = 'STD')
#, $res->{AA1}		One-letter code
#, $res->{AA3}		Three-letter code
#, $res->{NUM}		Residue number
#, $res->{ATOMS}	3D residue (array of atoms)
#
#. Sequence with associated structure ($seq->{TYPE} = 'STRUCT') contains
#  standard residue details plus
#, $res->{ATOMS}	Atoms associated with residue
#
#. Consensus residue details (seq->{TYPE} = 'CONSENSUS')
#, $res->{AA1}		One-letter code (Consensus sequence)
#, $res->{AA3}		Three-letter code  (Consensus sequence)
#, $res->{NUM}		Residue number
#, $res->{ATOMS}	3D residue (array of atoms)
#, $res->{AAS}[]	Array of residues found at this position
#, $res->{AA3S}[]	Array of residues found at this position
#, $res->{FREQS}[]	Frequency of each residue
#
#. Clustal conservation (seq->{TYPE} = 'CLUSTAL_CONS')
#, $res->{AA1}		Conservation symbyl (*, : , ., ' ' )
#, $res->{AA3}		Three-letter code
#, $res->{NUM}		Residue number
#. $Revision: 1.44.2.1.2.2 $
#>

use strict;
package Silico;

return 1;

# Amino acid routines #

sub aa {

	#<
	#? Return one letter code of the first amino acid in a sequence record
	#; Requires: Silico residue record
	#; Returns: one-letter code
	#>

	my $res = $_[0];

	return $res->{AA1};
}

sub is_gap {

	#<
	#? Determine if a residue is a gap
	#. Gaps can be marked with -,+ or maybe .
	#; Requires: Silico residue record
	#; Returns: 1 or 0
	#>

	my $res = $_[0];

	my $aa;

	$aa = $res->{AA1};

	if (!defined $aa) {
		return 1;
	}
	return 1 if $aa eq '-';
	return 1 if $aa eq '+';
	return 1 if $aa eq '.';

	return 0;
}

sub aa_three_letter {

	#<
	#? Return the three letter AA code from one letter code
	#. Returns 'UNK' if unknown
	#; Requires: string
	#; Returns: string
	#>
	
	return $Silico::Amino_Nucleic_Acids1{uc $_[0]} || 'UNK';
}

sub aa_one_letter {

	#<
	#? Return the one letter AA code from three letter
	#. Returns 'X' if unknown
	#; Requires: string
	#; Returns: string
	#>
	
	return $Silico::Amino_Nucleic_Acids{uc ($_[0])} || 'X';
}

# Sequence routines #

sub residues {

	#<
	#? Return array of residues from a sequence record
	#; Requires: Silico sequence record
	#; Returns: array of residues
	#>

	my $seq = $_[0];

	silico_msg('d', "Call to subroutine 'residues' with an undefined value!\n") if !defined $seq;

	if (!defined $seq->{SEQ}) {

		seq_printout($seq);
		silico_msg('d', "Sequence residues are undefined!\n");
	}

	return @{$seq->{SEQ}};
}

sub seq_fix_residue_numbers {

	#<
	#? Make sure residue numbers and the sequence are consistent.
	#. Renumbers residues from 1 if any residue numbers are undefined.
	#. Adds gaps to sequences if any residue numbers are missing
	#; Requires: sequence
	#; Returns: nothing
	#>
	
	my $seq = $_[0];
	
	my $diff;
	my $i;
	my $j;
	my $newresidues;
	my $nextres;
	my $res;
	my $res1;
	
	# Check that all residue numbers are defined
	# Renumber all residues from 1 if they are not.
	foreach $res (@{$seq->{SEQ}}) {
	
		next if (defined $res->{NUM});
		
		silico_msg('n', "Renumbering residues from 1\n");
		
		$i = 1;
		foreach $res1 (@{$seq->{SEQ}}) {
		
			$res->{NUM} = $i;
			++$i;
		}
		
		# Exit
		return;
	}
	
	# Add blank residues to make residue numbers contiguous
	$i = -1;
	foreach $res (@{$seq->{SEQ}}) {
	
		++$i;
	
		$nextres = $seq->{SEQ}[$i+1];
		
		last if !defined $nextres;

		push @$newresidues, $res;
		
		# If residue numbers are not contiguous
		if ($nextres->{NUM} > $res->{NUM}+1) {
			
			# Calculate number of missing residues
			$diff = $nextres->{NUM} - $res->{NUM} - 1;
			
			silico_msg('n', "Adding $diff gaps from residue $res->{NUM}\n");
			
			# Insert missing residues as blanks
			foreach $j (1..$diff) {
			
				my $newres;
				
				$newres->{NUM} = $res->{NUM} + $j;
				$newres->{AA1} = '-';
				$newres->{AA3} = '---';
				++$seq->{NUMRESGAP};
			
				push @$newresidues, $newres;
			}
		}
	}

	# Add last residue
	push @$newresidues, $seq->{SEQ}[-1];
	
	# Replace existing residue records
	$seq->{SEQ} = $newresidues;
}

sub seq_string {

	#<
	#? Return a sequence as a string of one-letter sequence
	#; Requires: sequence, flag to delete gaps
	#; Returns: string
	#>
	
	my $seq = $_[0];
	my $delgap = $_[1];
	
	my $res;
	my $string = '';

	foreach $res (residues($seq)) {

		next if ($delgap && is_gap($res));
	
		$string .= $res->{AA1} || '-';
	}
	
	return $string;
}

sub seq_merge {

	#<
	#? Merge fields from sequence1 into sequence2
	#; Requires; sequence1, sequence2, fields to merge
	#; Returns; nothing
	#>

	my $seq1 = $_[0];
	my $seq2 = $_[1];
	my $fields = $_[2];
	
	my $count = 1;
	my @fields;
	my $field;
	my $i;
	my $j;
	my $res1;
	my $res2;

	@fields = split " ", $fields;
	
	foreach $field (@fields) {
		
		if ($field eq 'ATOMS') {
			$seq2->{ENS} = $seq1->{ENS};
			last;
		}
	}
	
	silico_msg('c', heading("Merging fields (@fields) from $seq1->{NAME} to $seq2->{NAME}\n"));

	$i = -1; # Position in sequence1
	$j = -1; # Position in sequence2
	LOOP: foreach (@{$seq1->{SEQ}}) {
	
		# Handle gaps in sequence
		while (1) {
		
			++$i;
			
			$res1 = $seq1->{SEQ}[$i];
					
			# We have reached the end of sequence
			if ($i > $seq1->{NUMRESGAP}) {

				last LOOP;
			}

			# Skip gaps
			next if is_gap($res1);

			last;
		}
		
		while (1) {
		
			++$j;
			
			$res2 = $seq2->{SEQ}[$j];
			
			if (!defined $seq2->{NUMRESGAP}) {
				silico_msg('d', "\$seq2->{NUMRESGAP} is not defined!\n");
			}

			if ($j > $seq2->{NUMRESGAP}) {

				last LOOP;
			}

			next if (is_gap($res2));

			$res1->{AA3} ||= '---';
			$res2->{AA3} ||= '---';
			
			silico_msg('c', "($res1->{NUM} $res1->{AA3}:$res2->{NUM} $res2->{AA3}");
			++$count;
			
			if ($res1->{AA3} ne $res2->{AA3}) {
				silico_msg('c', "* Mismatch *");
			}
			
			silico_msg('c4', ") ");
			
			if ($count == 5) {
				$count = 0;
				silico_msg('c', "\n");
			}
			
			last;
		}
		
		
		last if !defined $res1 || !defined $res2;
		
		foreach $field (@fields) {
			
			$res1->{$field} = $res2->{$field} || '';
		}
	}
	
	silico_msg('c', "\n") if $count != 0;
	silico_msg('c', "\n");
}


sub seq_label_by_homology {

	#<
	#? Label residues by homology
	#. Requires that sequences_consensus & map sequence have been run first
	#. Percentage conservation is marked with the flags: H_CONS, H_HIGH
	#  H_MED, H_LOW and H_VLOW.  The occupancy of PDB files
	#  set to the fraction conservation (fraction of sequences with
	#  the most popular amino acid). The value $atom->{CONSERV}
	#  is also set to the fraction conservation
	#; Requires: Sequence molecule, Consensus sequence
	#; Returns: nothing
	#>
	
	my $seq = $_[0];
	my $conseq = $_[1];
	
	my $atom;
	my $conres;
	my $frac;
	my $i;
	my $mol;
	my $molres;
	my $numtypes;
	my $set;
	
	silico_msg('c', "Making atom sets by conservation Ftion\n");
	
	#sequences_printout($seq);

	$i = -1;
	foreach $conres (@{$conseq->{SEQ}}) {
	
		++$i;

		$molres = $seq->{SEQ}[$i];
	
		# Fraction of conservation from consensus sequence
		$frac = $conres->{FRAC} || 0;
		
		$numtypes = $conres->{NUMTYPES} || 0;
		
		# Assign conservation levels to sets
		if ($frac < 0.3) {
		
			$set = "H_VLOW";
			
		} elsif ($frac < 0.5){
		
			$set = "H_LOW";
		
		} elsif ($frac < 0.7){
		
			$set = "H_MED";
			
		} elsif ($frac < 0.9) {
			
			$set = "H_HIGH";
		} else {
		
			$set = "H_CONS";
		}
		
		foreach $atom (atoms($molres)) {
				
			$atom->{OCC} = $numtypes;
			
			$atom->{NUMTYPES} = $numtypes;
			$atom->{CONSERV} = $frac;
			$atom->{$set} = 1;
		}
	}

	# Make sets of marked atoms
	foreach $mol (@{$seq->{ENS}}) {
	
		make_atomset($mol, "H_VLOW", "H_VLOW");
		make_atomset($mol, "H_LOW", "H_LOW");
		make_atomset($mol, "H_MED", "H_MED");
		make_atomset($mol, "H_VLOW", "H_VLOW");
		make_atomset($mol, "H_HIGH", "H_HIGH");
		make_atomset($mol, "H_CONS", "H_CONS");
	}
}


sub seq_mark_start_end {

	#<
	#? Mark gaps at the start and end of as sequence with '+'
	#. Operates on non-consensus sequences
	#>
	
	my $sequences = $_[0];
	
	my $i;
	my $res;
	my $seq;
	
	foreach $seq (@$sequences) {
	
		# Start
		foreach $res (@{$seq->{SEQ}}) {
	
			last if ($res->{AA1} ne '-' && $res->{AA1} ne '+');
		
			$res->{AA1} = '+';
			$res->{AA3} = '+++';
		}
	
		# End
		for ($i = $#{$seq->{SEQ}}; $i >= 0; --$i) {
	
			$res = $seq->{SEQ}[$i];
	
			last if ($res->{AA1} ne '-' && $res->{AA1} ne '+');
			
			$res->{AA1} = '+';
			$res->{AA3} = '+++';
		}
	}
}


#######################
#                     #
#  Printout routines  #
#                     #
#######################

sub sequences_printout {

	#<
	#? Print out a set of sequences using one letter code
	#; Requires: sequences, options string
	#; Returns: nothing
	#>

	my $seqs = $_[0];
	
	my $seq;
	my $i;

	print "\n";
	print "Sequences Printout\n";

	$i = 0;
	foreach $seq (@$seqs) {
		
		print "\n";
		print "##### $i #####\n";
		print "\n";
				
		seq_printout_header($seq);
		seq_printout_one($seq, 100);
		++$i;
	}
	
	print "##################\n",
	print "\n";
}

sub seq_printout {

	#<
	#? Print out a single sequence record and sequence information
	#; Requires: Single sequence record, option string
	#; Returns: Nothing
	#>

	my $seq = $_[0];

	print "\n";
	print "Sequence Printout\n";
	print "#################\n";
	print "\n";
	
	if (ref $seq !~ 'HASH') {
	
		silico_msg('d', "The supplied sequence is not a hash!\n");
	}
	
	if (!defined $seq) {
	
		print "SEQUENCE NOT DEFINED\n";
		return;
	}
	
	seq_printout_header($seq);
	
	# Print out sequence in one-letter code
	print "\n";
	print "One-letter code\n";
	print "---------------\n";
	seq_printout_one($seq);
	
	# Print out sequence in three-letter code
	print "\n";
	print "Three-letter code\n";
	print "-----------------\n";
	seq_printout_three($seq);
}

sub seq_printout_header {

	#<
	#? Print out data from sequence hash except the
	#  actual SEQ records
	#; Requires: sequence
	#>
	
	my $seq = $_[0];
	my $i;
	my $key;
	print "Sequence location: $seq\n";
	print "\n";

	foreach $key (sort keys(%$seq)) {
	
		# Skip AA sequence record
		next if $key eq 'SEQ';

		# Scalars
		printf "\t%-15s\t",$key;
		
		if (defined $seq->{$key}) {
			print "\t'$seq->{$key}'";
		} else {
			print  "\tNOT DEFINED";
		}
		
		# Hashes
		if (ref $seq->{$key} eq 'HASH') {
			
			foreach (keys %{$seq->{$key}}) {
				print " <$_ ";
				if (defined $seq->{$key}{$_}) {
					print "'$seq->{$key}{$_}'";
				} else {
					print "\tNOT DEFINED";
				}
				print ">";
			}
			print "\n";
			next;
		}
		
		if (ref $seq->{$key} eq 'ARRAY') {
			print " <";
			for ($i=0; $i >= 0 && $i < 10 && $i <=$#{$seq->{$key}}; ++$i) {
				if (!defined $seq->{$key}[$i]) {
					print "\tNOT DEFINED";
					next;
				}
				print " '$seq->{$key}[$i]'";
			}
			print ">";
		}
		print "\n";
	}
}

sub seq_printout_one {

	#<
	#? Print out an AA sequence in one letter code
	#; Requires: Single sequence record, number of AAs per line (opt),
	#  number of AAs per group (opt)
	#; Returns: Nothing
	#>
	
	my $seq = $_[0];
	my $maxline = $_[1] || 50;
	my $maxgroup = $_[2] || 10;
	
	my $field;
	my $fieldcount = 0;
	my $i;
	my $j;
	my $linebase;
	my $res;
	my $rescount;
	
	if (ref($seq) !~ 'HASH') {
		silico_msg('d', "Sequence is not a HASH ($seq)!\n");
	}
		
	# Loop over fields
	foreach $field (seq_residue_fields($seq)) {
	
		# Skip fields that do not end in '1';
		next if $field !~ /1$/;
	
		$i = 1; # Group counter
		$j = 1; # Line counter
		$linebase = -1;
		$rescount = 0;
		
		print "\n" if $fieldcount > 0;
		print "$field:\n\n";
	
		while ($rescount <= $#{$seq->{SEQ}}) {
	
			# Print sequence in blocks of '$maxgroup', usually 10
			$i = 1;
			foreach $j (1..$maxline) {
		
				$rescount = $linebase+$j;
			
				last if $rescount > $#{$seq->{SEQ}};
	
				$res = $seq->{SEQ}[$rescount];
				
				# Print residue number at start of line
				if ($j == 1) {
					printf "%-5d  ", $res->{NUM} || 0;
				}
		
				# Print '_' if AA is undefined
				print ($res->{$field}||'_');
		
				++$i;
				if ($i > $maxgroup) {
					print " ";
				$i = 1;
				}
			}
		
			# Print residue number at end of line
			if ($rescount <= $#{$seq->{SEQ}}) {
				printf "  %-5d", $res->{NUM}||0;
			}
		
			print "\n";
		
			$linebase += $maxline;
		}
		
		print "\n" if $j != 1;
		
		++$fieldcount;
	}
}

sub seq_printout_three {

	#<
	#? Print out an AA sequence in three letter code
	#; Requires: Single sequence record
	#; Returns: Nothing
	#>
	
	my $seq = $_[0];
	
	my $i;
	my $j;
	my $res;
	
	$i = 1; # Group counter
	$j = 1; # Line counter
	
	foreach $res (@{$seq->{SEQ}}) {
		
		# Print out residue number at start of line
		if ($j == 1) {
			printf "%-5d  ", $res->{NUM}|| 0;
		}
		
		print ($res->{AA3}||'+++');
		
		# Print . for insertions
		if ($res->{NUM}) {
			print " ";
		} else {
			print ".";
		}
		
		++$i;
		if ($i == 11) {
			print " ";
			$i = 1;
		}
		
		# Print out residue number at end of line
		++$j;
		if ($j == 21) {
			printf "  %-5d\n", $res->{NUM}||0;
			$j = 1;
		}
	}
	
	print "\n" if $j != 1
}
	
sub residue_printout {
	
	#<
	#? Print out all seq_residue information
	#; Requires; residue
	#>

	my $res = $_[0];
	
	my $key;
	my $ref;
	my $val;
	
	if (!defined $res) {
		print "res: NOT DEFINED\n";
		print "---\n";
		return
	}

	print "Residue location: $res\n";
	
	KEY: foreach $key (sort (keys %$res)) {

		print "$key:\t";

		$ref = ref ($res->{$key});
		$val = $res->{$key};
	
		if (!defined $ref) {
			print "NOT DEFINED\n";
			next KEY;
		}
	
		# Scalars
		if ($ref eq '') {
			print "'$val'\n" ;
			next KEY;
		}
	
		# Arrays
		if ($ref eq 'ARRAY') {
			foreach (@$val) {
				if (!defined $_) {
					print "NOT DEFINED ";
					next;
				}
				print "'$_' ";
			}
			print "\n";
			next KEY;
		}
	
		# Hashes
		if ($ref eq 'HASH') {
			foreach (keys %$val) {
				if (!defined $_) {
					print "NOT DEFINED ";
					next;
				}
				print "$_ $val->{$_}' ";
			}
			print "\n";
			next KEY;
		}
	
		print "none of the above\n";
	}
	
	print "-----------------------------\n";
}


# Sequences routines

sub sequences_length {

	#<
	#? Find longest sequence in a set of sequences
	#; Requires: sequences
	#; Returns: length of longest array indexed from 1
	#>

	my $sequences = $_[0];

	my $maxlen;
	my $seq;

	$maxlen = 0;
	foreach $seq (@$sequences) {
	
		$maxlen = $seq->{NUMRESGAP} if $seq->{NUMRESGAP} > $maxlen;
	}

	return $maxlen;
}
	

sub sequences_consensus {

	#<
	#? Calculate a consensus sequences and other stats
	#. Nonstandard sequences (ie seq->{TYPE} ne 'STD) and sequences with names starting with
	#  an underscore are ignored
	#; Requires: Sequences
	#; Returns: consensus sequence record
	#>

	my $sequences =$_[0];
	
	my $aa;
	my $aas;
	my $conseq;		# Consensus sequence
	my $freqarray;
	my $frachash;
	my $i;
	my $key;
	my $maxlen;
	my $numaas;
	my $numtypes;
	my $value;

	$conseq = {};
	$conseq->{TYPE} = 'CONSENSUS';
	$conseq->{NAME} = '_consensus';

	$maxlen = sequences_length($sequences);

	# Loop over all positions in sequences
	foreach $i (0..$maxlen-1) {

		my $newres = {};
		$conseq->{SEQ}[$i] = $newres;
		$newres->{NUM} = $i+1;
	
		$aas = sequences_get_column($sequences, $i);
		$freqarray = get_aa_frequency($aas);

		# Total number of residue in this position
		$numaas = $#{$aas} + 1;
		
		# Total number of residue types in this position
		$numtypes = $#{$freqarray} + 1;
		$newres->{NT} = $numtypes;
		$newres->{NT} = 9 if $newres->{NT} > 9;
		
		# Consensus AA and count
		$newres->{FREQ} = $freqarray->[0][1];

		# Fraction of AAs due to the main AA
		$newres->{AA1}  = $freqarray->[0][0];
		$newres->{FRAC} = $newres->{FREQ}/$numaas;
		$newres->{FRAC_1} = sprintf "%1d", $newres->{FRAC} * 9;
		
		# Square of that
		$newres->{FRAC_SQ_1} = sprintf "%1d", $newres->{FRAC}*$newres->{FRAC} * 9;
		
		silico_msg('g', "freq $newres->{FREQ} naa $numaas frac $newres->{FRAC} frac_1 $newres->{FRAC_1} frac_sq_1 $newres->{FRAC_SQ_1}\n");

		# Fraction of AAs due to the secondary AA
		$newres->{AA2_1}  = $freqarray->[1][0] || '.';
		$newres->{FRAC2} = ($freqarray->[1][1]||0)/$numaas;
		$newres->{FRAC2_1} = sprintf "%1d", $newres->{FRAC2} * 9;

		# Fraction of AAs due to the tertiary AA
		$newres->{AA3_1}  = $freqarray->[2][0] || '.';
		$newres->{FRAC3} = ($freqarray->[2][1]||0)/$numaas;
		$newres->{FRAC3_1} = sprintf "%1d", $newres->{FRAC3} * 9;

		# Total number of different AAs found at this position
		$newres->{RES_TOTAL} = $numaas;
	}

	$conseq->{NUMRES} = $maxlen;
	$conseq->{NUMRESGAP} = $maxlen;
	@{$conseq->{PRINT_FIELDS}} = qw(AA1 NT FRAC_1 FRAC_SQ_1 AA2_1 FRAC2_1 AA3_1 FRAC3_1);
	
	return $conseq;
}

sub get_aa_frequency {

	#<
	#? Calculate the frequency of each AA type from an array of AAs
	#. AAs are sorted by frequency
	#; Requires: pointer to array of AAs
	#; Retruns: Pointer to array of arrays.  $array->[$i][0] = aa, $array->[$i][1] = count
	#>

	my $aas = $_[0];

	my $aa;
	my $array;
	my $count;
	my $counthash;
	my $i;

	foreach $aa (@$aas) {

		++$counthash->{$aa};
	}

	$i = 0;
	while (($aa, $count) = each (%$counthash)) {

		$array->[$i][0] = $aa;
		$array->[$i][1] = $count;
		++$i;
	}

	# Sort by frequency, then alphabetically
	@$array = sort {
		return 1 if $a->[1] < $b->[1];
		return -1 if $a->[1] > $b->[1];
		$a->[0] cmp $a->[0];
	} @$array;

	return $array;
}

sub sequences_get_column {

	#<
	#? Return one column from a set of sequences
	#; Requires: sequence, column number
	#; Returns: pointer to array of aas
	#>

	my $sequences = $_[0];
	my $position = $_[1];

	my $aa;
	my $aas;
	my $seq;
	my $res;

	@$aas = ();

	foreach $seq (@$sequences) {
	
		# Skip non-standard residue types
		if ($seq->{TYPE} ne 'STD' || $seq->{NAME} =~ /^_/) {
			next;
		}
		
		$res = $seq->{SEQ}[$position];
		
		if (is_gap($res)) {
			$aa = '-';
		} else {
			$aa = $res->{AA1};
		}

		push @$aas, $aa;
	}

	return $aas;
}

sub sequences_clustal_homology {

	#<
	#? Calculate the homology of a set of sequences according to the clustal scheme
	#. Sets the 'CONS' field to values of ~, ., : or * for each AA
	#. Nonstandard sequences (ie seq->{TYPE} ne 'STD') and sequences with names starting with
	#. an underscore are ignored
	#; Requires: Sequences
	#; Returns: a sequence containing clustal markers of homology
	#>

	my $sequences = $_[0];
	
	my $aas;
	my $newseq;	# Sequence containing clustal conservation markers
	my $i;
	my $maxlen;
	my $n1;

	$newseq = {};
	$newseq->{TYPE} = 'CLUSTAL';
	$newseq->{NAME} = '_clustal';

	$maxlen = sequences_length($sequences);
	
	print heading("seq0\tnumaas\tidentical_1\tidentical4\tstrong1\tstrong4\n");
	
	# Loop over all positions in sequences
	foreach $i (0..$maxlen-1) {
	
		$aas = sequences_get_column($sequences, $i);

		my $newres = {};

		$newseq->{SEQ}[$i] = $newres;
		$newres->{NUM} = $i+1;
			
		# Calculate clustal homology
		$newres->{SCORE}  = calc_conservation_group_membership($aas);

		#print "i $i  First res: $aas->[0]  ";
		
		# Calculate number of identical residues
		$n1 = calc_conservation_group_percentage($aas, 'identical');
		$newres->{IDENTICAL_1} = sprintf "%1d", $n1 * 9;
		$newres->{IDENTICAL_4} = sprintf "%4.2f", $n1 * 100;
		
		# Calculate fraction of strongly conserved residues
		$n1 = calc_conservation_group_percentage($aas, 'strong');
		$newres->{STRONG_1} = sprintf "%1d", $n1 * 9;
		$newres->{STRONG_4} = sprintf "%4.2f", $n1 * 100;
		
		# Calculate number of residue types
		$newres->{NUMRESTYPES_1} = calc_numrestypes($aas);
		
		print "$aas->[0]\t$newres->{IDENTICAL_1}\t$newres->{IDENTICAL_4}\t$newres->{STRONG_1}\t$newres->{STRONG_4}\n\n";
		
		++$newseq->{NUMRESGAP};
	}

	# Make the fields be printed
	@{$newseq->{PRINT_FIELDS}} = qw(SCORE STRONG_1 IDENTICAL_1 NUMRESTYPES_1);

	return $newseq;
}

sub sequences_STD {

	#<
	#! Select only amino acid sequences from a set of sequences
	#; Requires: sequences
	#; Returns: sequences
	#>

	my $sequences = $_[0];
	my $std;

	($std, undef) = sequences_STD_nonSTD($sequences);

	return $std;
}

sub sequences_nonSTD {

	#<
	#! Select only non-amino acid sequences from a set of sequences
	#; Requires: sequences
	#; Returns: sequences
	#>

	my $sequences = $_[0];
	my $nstd;

	(undef, $nstd) = sequences_STD_nonSTD($sequences);

	return $nstd;
}

sub sequences_STD_nonSTD {

	#<
	#! Select both aa and non-amino acid sequences from a set of sequences
	#; Requires: sequences
	#; Returns: standard sequences, non aa sequences
	#>

	my $sequences = $_[0];

	my $seq;
	my $std = ();
	my $nonstd = ();

	foreach $seq (@$sequences) {

		# Sequence type is standard and does not start with an underscore or 'file://'
		if  ($seq->{TYPE} eq 'STD' &&  $seq->{NAME} !~ /^_/ && $seq->{NAME} !~ /^file:\/\//i) {

			push @$std, $seq;

		} else {

			push @$nonstd, $seq;
		}

	}

	return $std, $nonstd;
}


# Reading routines #

sub read_clustal {

	#<
	#? Read in Clustal format NA or Protein sequence files
	#; Requires: Filename
	#; Returns: Sequences record
	#>

	my $infile = $_[0];
	
	my $aa;
	my $count;
	my $name;
	my $line;
	my $res;
	my $seq;
	my $sequences;
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Can not open file $infile for reading!\n");
		return undef;
	}

	silico_msg('c', "Reading aln file: $infile\n");

	# Get headers
	$line = <INFILE>;
	return undef if !defined $line;

	if ($line !~ /^CLUSTAL/) {
		silico_msg('e', "Header line does not begin with the text \"CLUSTAL\"!\n");
		close INFILE;
		return undef;
	}

	$count = 0;
	while (<INFILE>) {

		chomp;
	
		# Skip lines which don't have a 'word' character
		if (! /\w|-/) {

			$count = 0;
			next;
		}
		
		# Initialize if required
		if (!defined $sequences->[$count]) {
			my $newseq = {};
			$newseq->{TYPE} = 'STD';
			$newseq->{SOURCE} = 'aln';
			$newseq->{SOURCE_FILE_NAME} = $infile;
			$newseq->{NUMRESGAP} = 0;
			$newseq->{NUMRES} = 0;

			$sequences->[$count] = $newseq;
		}

		$seq = $sequences->[$count];

		# Match the first characters before a blank space as the
		# sequence name ($1)
		# Match everything up to the next space as the sequence ($2)
		
		/^(\S*)\s+(\S*)/;
		$name = $1;
		
		# Check that name matches expected name
		if (defined $seq->{NAME} && $name ne $seq->{NAME}) {
			silico_msg('d', "Sequence name $name does not match expected name $sequences->[$count]{NAME}!\n");
		}
		
		$seq->{NAME} = $name;
		
		# Add partial sequence to whole sequence
		foreach $aa (split '', $2) {
		
			# Number of residues in sequence including gaps
			++$seq->{NUMRESGAP};
			# Number of residues in sequence excluding gaps
			++$seq->{NUMRES} if ($aa ne '-');

			my $res = {} ;
			
			# Three and one-letter codes
			$res->{AA1} = uc $aa;
			$res->{AA3} =  aa_three_letter($aa);
			
			# Residue number (set to zero for spaces)
			if ($aa eq '-') {
				$res->{NUM} = 0;
			} else {
				$res->{NUM} = $seq->{NUMRES};
			}
			
			# Add to sequence
			push @{$seq->{SEQ}}, $res;
		}
		
		++$count;
	}
	
	return $sequences;
}

sub read_tcoffee_score {

	#<
	#? Read in tcoffee score file
	#; Requires: Filename
	#; Returns: Sequences record
	#>

	my $infile = $_[0];
	
	my $aa;
	my $count;
	my $line;
	my $name;
	my $seq;
	my $sequences;
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Can not open file $infile for reading!\n");
		return undef;
	}

	silico_msg('c', "Reading tcoffee file: $infile\n");

	# Get header
	$line = <INFILE>;
	return undef if !defined $line;

	if ($line !~ /^T-COFFEE/) {

		silico_msg('e', "File does not begin with \"T-COFFEE\"!\n");
		return undef;
	}

	foreach (0..5) {
		return undef if !defined <INFILE>; # Author, cpu time, score, *, BAD AVG GOOD, *
	}

	# Skip overall scores
	while (<INFILE>) {

		last if $_ !~ /\w/;
	}

	$count = 0;
	while (<INFILE>) {

		chomp;
	
		# Skip lines which don't have a 'word' character
		next if (! /\w|-/);

		# Reset to start of sequences at 'Cons' line
		if (/^Cons/) {
			$count = 0;
			next;
		}
		
		# Initialize if required
		if (!defined $sequences->[$count]) {
			my $newseq = {};
			$newseq->{SOURCE} = 'TCOFFEE';
			$newseq->{SOURC_FILE_NAME} = $infile;
			$newseq->{TYPE} = 'STD';
			$newseq->{NUMRESGAP} = 0;
			$newseq->{NUMRES} = 0;

			$sequences->[$count] = $newseq;
		}

		$seq = $sequences->[$count];

		# Match the first characters before a blank space as the
		# sequence name ($1)
		# Match everything up to the next space as the sequence ($2)
		
		/^(\S*)\s+(\S*)/;
		$name = $1;
		
		# Check that name matches expected name
		if (defined $seq->{NAME} && $name ne $seq->{NAME}) {
			silico_msg('d', "Sequence name $name does not match expected name $sequences->[$count]{NAME}\n");
		}
		
		$seq->{NAME} = $name;

		# Add partial sequence to whole sequence
		foreach $aa (split '', $2) {
		
			# Number of residues in sequence including gaps
			++$seq->{NUMRESGAP};
			# Number of residues in sequence excluding gaps
			++$seq->{NUMRES} if ($aa ne '-');

			my $res = {} ;
			
			# Three and one-letter codes
			$res->{AA1} = uc $aa;
			$res->{AA3} =  aa_three_letter($aa);
			
			# Residue number (set to zero for spaces)
			if ($aa eq '-') {
				$res->{NUM} = 0;
			} else {
				$res->{NUM} = $seq->{NUMRES};
			}
			
			# Add to sequence
			push @{$seq->{SEQ}}, $res;
		}
		
		++$count;
	}
	
	return $sequences;
}

sub read_fasta {

	#<
	#? Read in fasta sequence format
	#. INCOMPLETE
	#; Requires: Filename, $options
	#; Returns: Sequences record
	#>

	my $infile = $_[0];
	my $options = $_[1];

	my $aa;
	my $count;
	my $name;
	my $comment;
	my $seq;
	my $sequences;
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Can not open file $infile for reading!\n");
		return undef;
	}

	silico_msg('c', "Reading fasta file: $infile\n");

	$count = -1;
	while (<INFILE>) {

		# Skip blank lines
		next if !/\w/;

		chomp;

		# Header
		if (/^>/) {

			++$count;

			# Initialize

			$seq = {};
			$seq->{TYPE} = 'STD';
			$seq->{SOURCE} = 'fasta';
			$seq->{SOURCE_FILE_NAME} = $infile;
			$seq->{NUMRESGAP} = 0;
			$seq->{NUMRES} = 0;

			$sequences->[$count] = $seq;

			# Remove bracket
			($name) = /^>(\S*)/;
			($comment) = /^>\S*\s+(.*)/;

			$seq->{NAME} = $name || "seq_$count";
			$seq->{COMMENT} = $comment || '';

			next;
		}

		# Add partial sequence to whole sequence
		foreach $aa (split '', $_) {
		
			# Number of residues in sequence including gaps
			++$seq->{NUMRESGAP};
			# Number of residues in sequence excluding gaps
			++$seq->{NUMRES} if ($aa ne '-');

			my $res = {} ;
			
			# Three and one-letter codes
			$res->{AA1} = uc $aa;
			$res->{AA3} =  aa_three_letter($aa);
			
			# Residue number (set to zero for spaces)
			if ($aa eq '-') {
				$res->{NUM} = 0;
			} else {
				$res->{NUM} = $seq->{NUMRES};
			}
			
			#print "$res->{NUM}\n";
			
			# Add to sequence
			push @{$seq->{SEQ}}, $res;
		}
	}

	return $sequences;
}

#
# Writing routines #
#

sub write_clustal {

	#<
	#? Write a clustal format file
	#; Requires: Sequences, filename
	#; Returns: Nothing
	#; Options:
	#, DELIMITER='x' Delimiter between fields
	#, NUMLINE=XX Print XX residues per line
	#, WIDTH=XX Print XX characters per residue
	#>

	my $sequences = $_[0];
	my $outfile = $_[1];
	my $options = $_[2] || '';
	
	my $aa;
	my $count;
	my @fields;
	my $field;
	my $i;
	my $j;
	my $seq;
	my $maxlen = 0;
	my $maxnamelength = 30; # Length of sequence name
	my $numline = 60;
	my $name;
	my $delimiter = '';
	my $warn;
	my $width = 1;

	# Parse options
	if ($options =~ /\bNUMLINE=/) {
		($numline) = $options =~/\bNUMLINE=(\d*)/;
	}
	if ($options =~ /\bWIDTH=/) {
		($width) = $options =~/\bWIDTH=(\d*)/;
	}

	($delimiter) = $options =~ /\bDELIMITER='(.*?)'/;
	$delimiter ||= '';
		
	open ( COUT, ">$outfile");
	if (! open ( COUT, ">$outfile")) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing aln file: $outfile\n");
	
	print COUT "CLUSTAL (silico) multiple sequence alignment\n";
	print COUT "\n";
	print COUT "\n";
	
	# Find longest sequence
	$count = 0;
	foreach $seq (@$sequences)  {
	
		++$count;
		
		if (!defined $seq) {
			
			#carp();
			silico_msg('d',	"Undefined sequence!\n",
					"Sequence number: $count\n");
		}
		if (!defined $seq->{NUMRESGAP})  {
		
			silico_msg('d',	"NUMRESGAP is not defined!\n",
					"Sequence number: $count   Sequence name: '$seq->{NAME}'");
		}
		
		$maxlen = $seq->{NUMRESGAP} if $seq->{NUMRESGAP} > $maxlen;
	}

	# Loop over residues
	for ($i =0; $i < $maxlen; $i +=$numline) {
		
		# Loop over sequences
		foreach $seq (@$sequences) {
			
			#seq_printout($seq);
			
			$seq->{PRINT_FIELDS}[0] = "AA1" if (!defined $seq->{PRINT_FIELDS}[0]);
			
			# Loop over fields
			foreach $field (@{$seq->{PRINT_FIELDS}}) {
			
				# Add field name if it is not a standard amino acid field
				if ($field eq 'AA1') {
					$name = $seq->{NAME} || '';
				} else {
					
					# Why do we need to get the filebase of a non-filename
					# string anyway? - BPR 19/10/2006
					if (defined $seq->{NAME}) {
						$name = get_filebase($seq->{NAME}).'_'.lc ($field);
					} else {
						$name = lc($field);
					}
				}

				# Remove spaces
				$name =~ s/ //g;
				
				# Shorten
				if (length ($name) > $maxnamelength) {

					$name = substr ($name, 0, $maxnamelength);
					$warn = 1;
				}

				# Print sequence name
				printf COUT "%-30s    ",$name;
				print COUT $delimiter;

				# Print AAs
				for ($j=$i; $j < $i+$numline; ++$j) {
				
					$aa = $seq->{SEQ}[$j]{$field};
					$aa = '-' if !defined $aa;
			
					printf COUT "%*s", $width, $aa;
					
					last if $j == $maxlen;
					
					print COUT $delimiter;
				
				}
				print COUT "\n";
			}
		}
		print COUT "\n";
		print COUT "\n";
	}

	silico_msg('w', "Some sequence names have been shortened to $maxnamelength characters\n") if $warn;

	close COUT;
}

sub write_pir {

	#<
	#? Write a pir format file (also known as NBRF format) that is compatible with Modeller
	#; Requires: Sequences, filename
	#; Returns: Nothing
	#>

	my $sequences = $_[0];
	my $outfile = $_[1];
	my $options = $_[2] || '';
	
	my $aa;
	my $i;
	my $j;
	my $field;
	my @fields;
	my $seq;
	my $name;
	my $numline = 60;
	
	# Parse options
	if ($options =~ /\bNUMLINE=/) {
		
		($numline) = $options =~/\bNUMLINE=(\d*)/;
	}
	
	open (POUT, ">$outfile");
	if (! open ( POUT, ">$outfile")) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing pir file: $outfile\n");
	
	print POUT "C; Generated by Silico\n";
	print POUT "\n";
	
	foreach $seq (@$sequences)  {
	
		$seq->{PRINT_FIELDS}[0] = "AA1" if (!defined $seq->{PRINT_FIELDS}[0]);
		
		# Loop over fields
		foreach $field (@{$seq->{PRINT_FIELDS}}) {
		
			# Add field name if it is not a standard amino acid field
			if ($field eq 'AA1') {
				$name = $seq->{NAME} || '';
			} else {
				$name = get_filebase($seq->{NAME}).'_'.lc ($field) ;
			}

			# Should be two letter code describing sequence type P1, F1, DL, DC, RL, RC, or XX
			print POUT ">P1;$name\n";
			
			# Modeller comment string
			if ($seq->{STRUCTURE_FILE}) {
		
				print POUT "structure:";
				print POUT $seq->{STRUCTURE_FILE};
			} else {
				print POUT "sequence:";
			}
			print POUT ($seq->{START} || ' ').":";
			print POUT ($seq->{START_CHAIN} || ' ').":";
			print POUT ($seq->{END} || ' ').":";
			print POUT ($seq->{END_CHAIN} || ' ').":";
			print POUT ($seq->{TITLE} || ' ').":";
			print POUT ($seq->{SOURCE} || ' ').":";
			print POUT ($seq->{RESOLUTION} || ' ').":";
			print POUT ($seq->{RFACTOR} || ' ').":";
			print POUT "\n";
		
			# Print AAs
			LOOP: for ($i =0; $i < $seq->{NUMRESGAP}; $i += $numline) {
				
				for ($j=$i; $j < $i+$numline; ++$j) {
				
					$aa = $seq->{SEQ}[$j]{$field};
					$aa = '-' if !defined $aa;
			
					print POUT $aa;
					
					if ($j == $seq->{NUMRESGAP}) {

						if ($j < $i + $numline - 1) {
							print POUT "*\n";
						} else {
							print POUT "\n";
							print POUT "*\n";
						}
						last LOOP;
					}
				}
				
				print POUT "\n";
				
			}
			print POUT "\n";
		}
	}

	close POUT;
}

sub write_fasta {

	#<
	#? Write fasta format
	#; Requires: Sequences, filename, options
	#; Returns: Nothing
	#>

	my $sequences = $_[0];
	my $outfile = $_[1];
	my $options = $_[2] || '';
	
	my $aa;
	my $i;
	my $j;
	my $field;
	my @fields;
	my $seq;
	my $name;
	my $numline = 60;
	
	# Parse options
	if ($options =~ /\bNUMLINE=/) {
		($numline) = $options =~/\bNUMLINE=(\d*)/;
	}
	
	open (POUT, ">$outfile");
	if (! open ( POUT, ">$outfile")) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing fasta file: $outfile\n");
	
	foreach $seq (@$sequences)  {

		$seq->{PRINT_FIELDS}[0] = "AA1" if (!defined $seq->{PRINT_FIELDS}[0]);
	
		# Loop over fields
		foreach $field (@{$seq->{PRINT_FIELDS}}) {

			if ($field eq 'AA1') {
				$name = $seq->{NAME} || '';
			} else {
				$name = get_filebase($seq->{NAME}).'_'.lc ($field) ;
			}
			print POUT ">$name\n";
			
			# Print AAs
			LOOP: for ($i =0; $i < $seq->{NUMRESGAP}; $i += $numline) {
				
				for ($j=$i; $j < $i+$numline; ++$j) {
				
					$aa = $seq->{SEQ}[$j]{$field};
					$aa = '-' if !defined $aa;
			
					print POUT $aa;
					
					if ($j == $seq->{NUMRESGAP}) {
						print POUT "\n";
						last LOOP;
					}
				}
				print POUT "\n";
			}
			print POUT "\n";
		}
	}
	close POUT;
}


#  Molecule routines

sub mol_get_sequence {

	#<
	#? Extract the amino or nucleic acid sequence from a molecule
	#; Requires: Molecule, optional options string
	#; Returns: Ensemble of sequences. One sequence for each CHAIN or SUBID
	#. Options: NOHETATM - ignore all pdb HETATM records
	#>
	
	my $mol = $_[0];
	my $options = uc $_[1];
	
	my $chain;
	my $chainlist;
	my $chaincount;
	my $gapped;
	my $insertcount;
	my $molres;
	my $newreslist;
	my $nohetatm = 0;
	my $reslist;
	my $seq;
	my $seqs;
	my $seqname;

	# Options
	$nohetatm = 1 if $options =~ /\bNOHETATM\b/;
	$gapped = 1 if $options =~ /\bGAPPED\b/;

	$chainlist = molecule_get_chains($mol);

	silico_msg('c', "Found ".($#{$chainlist}+1)." chains:\n",
			"\n");

	# Loop over chains in molecule
	$chaincount = 0;
	foreach $chain (@$chainlist) {

		++$chaincount;

		# Make a new empty sequence
		$seq = {};

		# Number of residues
		$seq->{NUMRES} = 0;
		$seq->{NUMRESGAP} = 0;
			
		# Set sequence name
		$seqname = $mol->{NAME} || 'Mol';
		$seqname =~ s/^\s*//;			# Remove leading spaces
		$seqname =~ s/\s*$//;			# Remove trailing spaces
		$seqname =~ s/\s+/_/g;			# Replace spaces with underscore
		$seqname .= "_$chain->[0]{CHAIN}";	# add chain letter

		$seq->{NAME} = $seqname;
		$seq->{CHAIN} = $chaincount;
		$seq->{MOL} = $mol;
		$seq->{SOURCE_FILE_NAME} = $mol->{SOURCE_FILE_NAME};
		
		# Get all residues within chain
		$reslist = molecule_get_residues ($mol, $chain, $gapped);
		
		# Remove HETATM residues and undefined residues if required
		if ($nohetatm) {
			foreach $molres (@$reslist) {
		
				# Skip if a hetatm and flag is set
				next if defined $molres->[0]{LABEL} && $molres->[0]{LABEL} =~ 'HETATM';
				push @$newreslist, $molres;
			}
			$reslist = $newreslist;
		}


		$insertcount = 0;
		# Loop over residues within chain
		foreach $molres (@$reslist) {

			my $res; # New sequence residue
			
			if ($#{$molres} == -1) {
				$res->{AA3} = 'UNK';
				$res->{AA1} = 'X';
				$res->{ATOMS} = ();
				$res->{NUM} = 0;
				++$insertcount;
			
			} else {
				
				#print "$molres->[0]{SUBNAME}\n";
			
				# Three letter code AA
				$res->{AA3} = $molres->[0]{SUBNAME} || 'UNK';
				$res->{AA3} =~ s/ //g;
			
				# One letter code AA
				$res->{AA1} = aa_one_letter($res->{AA3});
				
				# Atoms in molecule residue
				$res->{ATOMS} = $molres;

				# Residue number from file
				$res->{NUM} = $molres->[0]{SUBID} || 0;
				$res->{NUM} =~ s/ //g;
			}

			# Add residue to sequence
			push @{$seq->{SEQ}}, $res;
				
			# Count the residues
			++$seq->{NUMRES};
			++$seq->{NUMRESGAP};
		}

		# Skip if there were zero residues
		if ($seq->{NUMRESGAP} == 0) {

			silico_msg('w', "Chain $chaincount contained no residues.\n",
					"Skipping.\n");
			next;
		}
		
		#seq_printout($seq);

		# Deal with any chain breaks
		seq_fix_residue_numbers($seq) if !$gapped;

		push (@$seqs, $seq);
	}
	
	silico_msg('n', "Inserted $insertcount UNK residues\n") if $insertcount;

	return $seqs;
}

sub mol_count_residue_contacts {

	#<
	#? Count the number of non-neigbour atom-atom contacts for each residue of a protein
	#; Requires: Molecule
	#>

	my $mol = $_[0];

	my $atom;
	my $atom1;
	my $atom2;
	my $count;
	my $dist_sq;
	my $i;
	my $j;
	my $label;
	my $name;
	my $newreslist;
	my $reslist;
	my $res;
	my $res1;
	my $res2;
	my ($x1, $y1, $z1);
	
	use vars qw(%Amino_Nucleic_Acids);

	my $starttime = (times)[0];
	
	silico_msg('c', "Counting interresidue contacts less than 4 Angstroms\n");
	
	# Get all residues within chain (gapped)
	$reslist = molecule_get_residues ($mol, undef, 1);
	
	# Delete hydrogens and non-amino acids
	# from search list
	foreach $res (@$reslist) {
	
		$name = $res->[0]{SUBNAME};
		$name =~ s/ //g;
		
		# Skip non AAs
		next if (!defined $Amino_Nucleic_Acids{$name});
		
		my $newres;
	
		foreach $atom (@$res) {
		
			next if $atom->{ELEMENT} eq 'H';
			
			push @$newres, $atom;
		}
		
		push @$newreslist, $newres;
	}
		
	# Residue 1
	RESONE: foreach $i (0..$#{$newreslist}) {
	
		$res1 = $newreslist->[$i];
	
		$count = 0;
		
		# Residue 2
		RESTWO: foreach $j (0..$#{$newreslist}) {
		
			last if $j >= $i-1;
			
			$res2 = $newreslist->[$j];
			
			#next if distance_sq($res1->[0], $res2->[0]) > 300;
			
			ATOMONE: foreach $atom1 (@$res1) {
			
				$x1 = $atom1->{X};
				$y1 = $atom1->{Y};
				$z1 = $atom1->{Z};
			
				foreach $atom2 (@$res2) {
					
					$dist_sq = (($x1-$atom2->{X})**2+($y1-$atom2->{Y})**2 + ($z1-$atom2->{Z})**2);
					
					# These cutoffs are empirically determined to speed up
					# the search with minimal loss of accuracy
					next RESTWO if $dist_sq > 300;
					next ATOMONE if $dist_sq > 178;
					
					# Contact distance is defined as 4 Angstroms
					next if $dist_sq > 16;
					
					++$count->[$i];
					++$count->[$j];
				}
			}
		}
	}
	
	foreach $i (0..$#{$newreslist}) {
		$res1 = $newreslist->[$i];
		$count->[$i] ||= 0;
		
		# Label the atoms BURIED, EXPOSED or VEXPOSED
		if ($count->[$i] >= 25) {
			$label = "BURIED";
		} elsif ($count->[$i] >= 10) {
			$label = "EXPOSED";
		} else {
			$label = "VEXPOSED";
		}
		
		foreach $atom (@$res1) {
			$atom->{TEMP} = $count->[$i];
			
			$atom->{$label} = 1;
		}
	}
	
	make_atomset($mol, "BURIED", "BURIED");
	make_atomset($mol, "EXPOSED", "EXPOSED");
	make_atomset($mol, "VEXPOSED", "VEXPOSED");
	
	silico_msg('c', "\n",
			"Time: ".calc_human_readable_time((times)[0]-$starttime)."\n");
}

#
# Routines for relating sequence and structural data
#

sub sequences_read_mols {

        #<
        #? Read in any molecule files that correspond to sequence names
        #. Only considers pdb and mol2 files
        #; Requires: sequences, flag to map residues in molecule to sequence (on by default)
        #; Returns: nothing
        #>

        my $sequences = $_[0];
	my $do_mapping = $_[1] || 1;

	my $mol;
        my $mols;
        my $numseq;
        my $seq;
        my $sname;

        $numseq = 0;
        foreach $seq (@$sequences) {
                
                $sname = $seq->{NAME};
                $sname =~ s/ //g;
                chomp $sname;

                next if ($sname !~ /pdb$/ && $sname !~ /mol2$/);

                # Skip any sequences without associated pdb file
                if (-e $sname) {
                        silico_msg('c', "Reading file $sname\n");
                } else {

                        silico_msg('w', "Molecule file name '$sname' not found in current directory. Skipping\n");
                        next;
                }

                $mols = read_mol_any ($sname, undef, undef, 'NOCHAINFIX');
                if (!defined $mols) {
                        silico_msg('w', "Can not read any molecule data from file $sname. Skipping.\n");
                        next;
                }

		if ($do_mapping) {

			foreach $mol (@$mols) {
				mol_map_sequence($mol, $seq);
			}				
		}

                $seq->{TYPE} = 'STRUCT';
                $seq->{ENS} = $mols;

                ++$numseq;
        }
}

sub sequences_write_mols {

        #<
        #? Write any molecule files attached to sequences
        #>

        my $sequences = $_[0];
        my $format = $_[1];

        my $mols;
        my $seq;
        my $sname;
	my $options;

        foreach $seq (@$sequences) {

                $sname = $seq->{NAME};
                $sname =~ s/ //g;
                chomp $sname;
		
		$mols = $seq->{ENS};

                # Skip any sequences without associated ensemble
                next if !$mols->[0];

                #ensemble_printout($mols);
		
		$format ||= get_oformat($mols);
		
		$options = '';
		$options = 'NOBOND' if ($format eq 'pdb' && !get_sflag('bond'));

                $mols = write_mol_any ($mols, undef, $format, $options);
              
        }
}

sub mol_map_sequence {

	#<
	#? Map residues of molecule file onto a sequence
	#. Adds the atoms of each molecule residue to the sequence residue array @{%res->{ATOMS}}
	#; Requires; molecule, sequence
	#; Returns: Nothing
	#>

	my $mol = $_[0];
	my $seq = $_[1];

	my $atom;
	my $count;
	my $isgap;
	my $i;			# Position in aa sequence
	my $j;			# Position in molecule sequence
	my $k;
	my $mols;
	my $molres;		# Molecule residue
	my $molseq;
	my $molseqs;
	my $mr;
	my $r;
	my $res;		# Sequence residue
	my $start1;
	my $start2;
	
	# Extract amino acid sequence from the molecule file
	$molseqs = mol_get_sequence($mol, 'GAPPED');
	
	print heading("Mapping $mol->{NAME}\n");
	
	$i = -1;  # Position in sequence
	# Remove leading gaps from aa sequence
	while (1) {
		$res = $seq->{SEQ}[$i+1];
		last if !defined $res;
		last if !is_gap($res);
		++$i;
	}
	
	$count = 0;
	my $molseqcount = 0;
	# Loop over molecule sequences
	LOOP_A: foreach $molseq (@$molseqs) {

		#seq_printout($molseq);

		++$molseqcount;
		$j = -1; # Position in molecule
		LOOP_B: foreach $molres (@{$molseq->{SEQ}}) {

			++$j;
			$molres = $molseq->{SEQ}[$j];
			last LOOP_B if !defined $molres;	

			# Remove any gaps from aa sequence
			while (1) {

				++$i;
				$res = $seq->{SEQ}[$i];		

				# We have reached the end of the sequence
				last LOOP_A if !defined $res;

				#print "$i $j res: $res->{NUM} $res->{AA3}\n";

				next if is_gap($res);
				last;
			}
			
			printf "%5s %3s -> %5s %3s  ", $res->{NUM},$res->{AA3},$res->{NUM},$molres->{AA3};
			
			++$count;
			if ($count == 5) {
				$count = 0;
				print "\n";
			}
			
			# Residues match:
			#	If one is a gap
			#	If one is an 'X'
			#	If they are both the same
			
			if ($molres->{AA1} eq 'X' || $res->{AA1} eq 'X' || is_gap($molres) || is_gap($res) || $molres->{AA3} eq $res->{AA3}) {
				
				
				# Add the molecule atoms to $res->{ATOMS}.
				# Atoms are added individually to allow 
				# multiple molecules to be coloured using
				# information from the same residue
				foreach $atom (@{$molres->{ATOMS}} ) {
					push @{$res->{ATOMS}}, $atom;
				}
				
			} else {
			
				print "-------No match-------\n";
				print "End of mapping\n";
				last LOOP_B;
			}
		}
	}
	
	# Store the associated molecules as $seq->{ENS}
	$seq->{ENS} = $mols;
}

sub seq_label_by_field {

	#<
	#? Label molecule atoms in target sequence according to properties of property sequence using
	#  a definable mapping
	#; Requires: target sequene, property sequence, field, mappiing (hash containing value, label pairs
	#  Eg clustal mapping is (KEY, VALUE) *, CLUSTAL_CONS; :,CLUSTAL_STRONG; ., CLUSTAL_WEAK
	#. Target and property sequences can be identical
	#; Returns: nothing
	#>
	
	my $target = $_[0]; # Target sequence
	my $seq = $_[1]; # Property sequence
	my $field = $_[2];
	my $mapping = $_[3];
	
	my $atom;
	my $i;
	my $key;
	my $mol;
	my $molres;
	my $num;
	my $res;
	my $set;
	my %sets;
	
	silico_msg('c', "Making atom sets by $field\n");
	
	$i = -1;
	foreach (@{$seq->{SEQ}}) {
	
		++$i;

		$res = $seq->{SEQ}[$i];
		$molres = $target->{SEQ}[$i];
		
		next if !defined $res->{$field};
		
		# Assign conservation levels to sets
		foreach $key (keys (%$mapping)) {
		
			if ($res->{$field} eq $key){
		
				$set = $mapping->{$key};
				
				foreach $atom (atoms($molres)) {
			
					$atom->{$set} = 1;
					
					++$sets{$set};
				}
			}
		}
	}

	# Make sets of marked atoms
	foreach $mol (@{$target->{ENS}}) {
	
		foreach $key (keys (%sets)) {
		
			$num = make_atomset($mol, $key, $key);
			silico_msg('c', "Added $num atoms to set $key\n");
		}
	}
}

sub seq_label_indels {

	#<
	#? Label sequence insertion and deletion points
	#. Requires that sequences_consensus & map sequence have been run first
	#. Insertions and deletions is marked as the sets: INS and DEL
	#; Requires: Sequence molecule, Consensus sequence
	#; Returns: nothing
	#>
	
	my $seq = $_[0];
	my $conseq = $_[1];
	
	my $atom;
	my $conres;
	my $i;
	my $mol;
	my $res;
	my $resnext;
	my $resprev;
	
	# loop over residues in consensus sequence
	$i = -1;
	foreach (@{$conseq->{SEQ}}) {
	
		++$i;
		
		$res = $seq->{SEQ}[$i];
		$conres = $conseq->{SEQ}[$i];
		
		# Deletions
		if ($conres->{HAS_GAP}) {
			
			# Mark all atoms in the molecule residue with a DEL flag
			foreach $atom (atoms($res->{ATOMS})) {
				
				$atom->{DEL} = 1;
			}
		}
		
		$resnext = $seq->{SEQ}[$i+1];
		
		last if !defined $resnext;
		
		# Check to see if next residue is missing a residue number
		if (is_gap($resnext)) {
			foreach $atom (atoms($res->{ATOMS})) {
				
				$atom->{INS} = 1;
			}
		}
		
		# Skip the first residue
		next if $i == 0;
		
		# Check to see if previous residue is missing a residue number
		$resprev = $seq->{SEQ}[$i-1];
		
		if (is_gap($resprev)) {
			foreach $atom (atoms($res->{ATOMS})) {
				
				$atom->{INS} = 1;
			}
		}
	}
	
	# Make sets of marked atoms
	foreach $mol (@{$seq->{ENS}}) {
	
		make_atomset($mol, "INS", "INS");
		make_atomset($mol, "DEL", "DEL");
	}
}

sub seq_residue_fields {

	#<
	#? Get all residue fields of a sequence
	#; Requires: sequence
	#; Returns: fields as an array
	#>
	
	my $seq = $_[0];
	
	my $key;
	my @list;
	my $res;
	
	$res = $seq->{SEQ}[0];
	
	foreach $key (keys %$res) {
	
		push @list, $key;
	}
	
	return sort @list;
}

#
# Misc routines
#

sub calc_conservation_group_membership {

	#<
	#? Calculate strong and weak conservation groups defined by Clustal
	#. "These are all the positively scoring groups that occur in the Gonnet Pam250
	#  matrix. The strong and weak groups are defined as strong score >0.5 and weak
	#  score =<0.5 respectively." Clustalw documentation
	#: Requires: pointer to array of amino acids
	#: Returns: *, :, . or '' (Corresponding to identity, strong conservation, weak conservation, no conservation)
	#>

	my $aas = $_[0];

	my $aa;
	my $flag;
	my $group;
	my $strong;
	my $weak;

	@$strong = qw(STA NEQK NHQK NDEQ QHRK MILV MILF HY FYW);
	@$weak = qw(CSA ATV SAG STNK STPA SGND SNDEQK NDEQHK NEQHRK FVLIM HFY);

	# Check to see if all aas are identical
	$flag = 1;
	foreach $aa (@$aas) {

		# return if we have a gap
		return ' ' if $aa eq '-';

		next if $aa eq $aas->[0];
		$flag = 0;
		last;
	}

	return "*" if $flag;

	GROUP: foreach $group (@$strong) {

		foreach $aa (@$aas) {
		
			next GROUP if $group !~ /$aa/;
		}

		return ":";
	}

	GROUP: foreach $group (@$weak) {

		foreach $aa (@$aas) {

			next GROUP if $group !~ /$aa/;
		}

		return ".";
	}

	return ' ';
}


sub calc_conservation_group_percentage {

	#<
	#?
	#>

	my $aas = $_[0];
	my $set = lc $_[1];

	my $aa;
	my $best_group;
	my $flag;
	my $hash;
	my $identical;
	my $key;
	my $group;
	my $max;
	my $members;
	my $num;
	my $strong;
	my $val;
	my $weak;

	# Strong plus single AAs
	@$identical = qw(A C D E F G H I K L M N P Q R S T V W Y);
	@$strong = qw(STA NEQK NHQK NDEQ QHRK MILV MILF HY FYW);
	@$weak = qw(CSA ATV SAG STNK STPA SGND SNDEQK NDEQHK NEQHRK FVLIM HFY);

	@$members = @$identical if $set eq 'identical';
	@$members = (@$identical, @$strong) if $set eq 'strong';
	@$members = (@$identical, @$strong, @$weak) if $set eq 'weak';

	$num = $#{$aas}+1;

	GROUP: foreach $group (@$members) {

		foreach $aa (@$aas) {

			next if $group !~ /$aa/;
			++$hash->{$group};
		}
	}

	$max = 0;
	$best_group = '';
	foreach $key (keys(%$hash)) {

		$val = $hash->{$key};
		if ($val > $max) {
			$max = $val;
			$best_group = $key;
		}
	}
	
	print "best_group $best_group num: $num max: $max frac:".($max/$num)."\n";

	return 0 if $num == 0;

	return $max/$num;
}

sub calc_numrestypes {

	my $aas = $_[0];
	
	my $aa;
	my $hash;
	
	foreach $aa (@$aas) {
	
		next if $aa eq '-';
		next if $aa eq 'X';
		++$hash->{$aa};
	}
	
	# Keys in scalar context returns
	# number of keys
	my $n = keys (%$hash);
	
	# reduce n to 9 max
	$n= 9 if $n > 9;
	
	return $n;
}
	


