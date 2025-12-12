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
#! silico_residue.pm
#? Silico routines for handling residues in molecules
#. $Revision: 1.1.2.11 $
#>

use strict;
package Silico;


####################################################################################
#
#  Routines operating on existing residues
#
####################################################################################


sub molecule_fix_and_get_residues {

	#<
	#? Return an array of residues.  Existing residues are returned if they exist and seem sane 
	#  otherwise new residues are generated
	#. Calls either 'molecule_get_residues' or 'molecule_make_residues'
	#. Sets $atom->{SUBCOUNT} (numbered from 1)
	#; Requires: molecule
	#; Returns: array of residues. Each residue is an array of atoms 
	#  belonging to that residue
	#>

	my $mol = $_[0];

	my $acount = 0;;
	my $flag;
	my $maxat = 0;
	my $osubid = -1;
	my $rescount = 0;
	
	foreach my $atom (@{$mol->{ATOMS}}) {

		if (!defined $atom->{SUBID} || !defined $atom->{SUBNAME} || $atom->{SUBNAME} !~ /\w/) {

			$flag = 1;
			last;
		}

		++$acount;

		if ($atom->{SUBID} != $osubid) {

			++$rescount;
			$maxat = $acount if $maxat < $acount;
			$acount = 0;
		}

		$osubid = $atom->{SUBID};
	}

	$flag = 1 if $rescount == 1 && $mol->{NUMATOMS} > 75;
	$flag = 1 if $maxat > 75;
	
	if ($flag) {
		silico_msg('n', "Too few residues found (rescount) for total atoms $mol->{NUMATOMS}.  Generating residues.\n");
	
		# Make new residues
		return molecule_make_residues($mol);
	} else {
	
		# Return existing residues
		return molecule_get_residues($mol) if !$flag;	
	}
}

sub molecule_get_residues {

	#<
	#? Return an array of amino acid residues that already exist within a file. Each residue is a simple array of atoms
	#. Sets $atom->{SUBCOUNT} (numbered from 1)
	#. Setting the gap flag will add empty residues for any gaps in the AA sequence.
	#  The mark_first_residue flag will mark the first atom of each residue with a 'STARTRES' flag.
	#; Requires: molecule, optional atom list containing a subset of atoms, gap flag, 
	#  flag to mark the first atom of each residue
	#; Returns: array of arrays of pointers to all atoms in each residue
	#>

	my $mol = $_[0];
	my $atomlist = $_[1];
	my $gapped = $_[2];
	my $mark_first_atom_of_residue = $_[3];
	
	my $rescount = 0;
	my $reslist;

	if (!defined $atomlist) {
		$atomlist = $mol->{ATOMS};
	}
	
	my $oldsubid   = $atomlist->[0]{SUBID}   ||  1;
	my $oldsubname = $atomlist->[0]{SUBNAME} || '';
	my $oldchain   = $atomlist->[0]{CHAIN}   || '';
	my $oldsegid   = $atomlist->[0]{SEGID}   || '';
	
	if ($mark_first_atom_of_residue) {
		foreach my $atom (@$atomlist) {
			delete $atom->{STARTRES};
		}
		$mol->{ATOMS}[0]{STARTRES} = 1;
	}

	ATOM: foreach my $atom (@$atomlist) {
	
	
		my $label   = $atom->{LABEL}   || '';
		my $subid   = $atom->{SUBID}   ||  1;
		my $subname = $atom->{SUBNAME} || '';
		my $chain   = $atom->{CHAIN}   || '';
		my $segid   = $atom->{SEGID}   || '';
		my $oldlabel;

		if ($segid ne $oldsegid || $chain ne $oldchain) {

			# End of SEGID or chain
			$oldsubname = $subname;
			$oldsubid = $subid;
			$oldchain = $chain;
			$oldsegid = $segid;
			$atom->{STARTRES} = 1 if $mark_first_atom_of_residue;
			++$rescount;

		} elsif ($subid != $oldsubid || $subname ne $oldsubname) {

			if ($gapped && $subid - $oldsubid > 1) {

				# Make unknown residues
				foreach my $i (1 .. (($subid-$oldsubid)-1)) {
		
					++$rescount;
					@{$reslist->[$rescount]} = ();
					$reslist->[$rescount][0]{LABEL}  = $oldlabel;
				}
			} else {
			
				$oldlabel = $label;
				$oldsubname = $subname;
				$oldsubid = $subid;
				$oldchain = $chain;
				$oldsegid = $segid;
			}
			
			$atom->{STARTRES} = 1 if $mark_first_atom_of_residue;
			++$rescount;
		}
		
		$atom->{SUBCOUNT} = $rescount+1;
	
		push @{$reslist->[$rescount]}, $atom;
	}
	
	$mol->{RESIDUES} = $reslist;

	return $reslist;
}

sub molecule_get_residues_hash {

	#<
	#? Return an array of residues. Each residue is a hash.  
	#. Atoms are in $res->{ATOMS}
	#. Number (count) are in $res->{NUM}
	#. NAME is in $res->{NAME}
	#. SUBID is in $res->{SUBID}
	#; Requires: molecule, optional atom list containing a subset of atoms, gap flag
	#; Returns: array of hashes
	#>

	my $mol = $_[0];
	my $atomlist = $_[1];
	my $gapped = $_[2];

	my $newresidues = ();
	my $residues = molecule_get_residues($mol, $atomlist, $gapped);

	my $i = 0;
	foreach my $res (@$residues) {

		++$i;
		my $nres;
		$nres->{ATOMS} = $res;
		$nres->{NUMATOMS} = $#{$res}+1;
		$nres->{NUM}   = $i;
		$nres->{NAME}  = $res->[0]{SUBNAME} || '';
		$nres->{SUBID} = $res->[0]{SUBID} || '';
		push @$newresidues, $nres;
		
		#atom_printout($nres);
	}

	return $newresidues;
}

####################################################################################
#
#  Routines creating new residues
#
####################################################################################

sub label_amino_acids {
	
	#<
	#? Create residues  and identify them based on a SMILES residue database
	#. Uses residues as the basis on which to match an atom string. 
	#; Requires: molecule
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
		
	my $aa;
		
	# Parse the list of amino acids.
	if (get_sflag('new')) {
		$aa = parse_residue_dictionary();
	} else {
		$aa = parse_amino_acids_smiles();
	}
		
	# Set formal charges on parent molecule so that we get correct smiles strings
	molecule_formal_charge($mol);
		
	# Create new residues. Each residue is a list of atoms.
	my $residues = molecule_make_residues($mol);
	
	# Label waters
	label_waters($mol);
	
	print heading("Matching residues\n");

	my $count = 0;
	foreach my $res (@$residues) {
	
		# Skip waters. 
		next if ($#{$res} == 2 && $res->[0]{FG}{W});
		++$count;

		# Store old name in OSUBNAME
		foreach my $atom (@$res) {
			$atom->{OSUBNAME} = $atom->{SUBNAME};
		}

		if (get_sflag('new')) {
			match_amino_acid2($mol, $res, $aa);
		} else {
			match_amino_acid($mol, $res, $aa);
		}
		
		if ($#{$res} < 0) {
			silico_msg('w', "Residue $count is empty. Skipping.\n");
			next;
		}
				
		clear_smiles_sets($mol, $res);
	}
	
	if ($mol->{UNMATCHED_RESIDUES_SMILES}) {
			
		print heading("Unmatched residues\n");
		foreach my $key (keys %{$mol->{UNMATCHED_RESIDUES_SMILES}}) {
			print "$key\n";
		}
	}

	print "\n";
}


sub molecule_make_residues {

        #<
        #? Divide a molecule into a NEW set residues in a molecule based on connectivity and amino acids.
	#  Old residues are overwritten.
	#. Sets $atom->{SUBCOUNT} (numbered from 1)
	#. Splits residues based on:
	#. Amide bonds
	#. Disulfides
	#. Ethers
	#. Esters
	#. Triazoles resulting from click reactions according to the atoms from the original
	#  azide and alkyne
        #; Requires: molecule
	#; Returns: List of residues ordered N->C as an array of arrays of pointers to all atoms in each 
	#  residue or undef if there is a problem. Also sets $mol->{RESIDUES}
        #>

        require silico_split;

        my $mol = $_[0];

	molecule_check_and_fix_connectivity($mol);
	mol_label_functional_group($mol);

	print heading("Generating residues.  These may not match the original residues.\n");

	# Remove existing SUBID's 
	foreach my $atom (@{$mol->{ATOMS}}) {
		$atom->{OSUBID} = $atom->{SUBID};
		undef $atom->{SUBID};
	}
	
	my $i = -1;
        my $subid = 0;
        ATOM: foreach my $atom (@{$mol->{ATOMS}}) {

                ++$i; 

		# Skip atom if we have already
                # assigned this atom to a substructure
                next if defined $atom->{SUBID};
		
		next if $atom->{SKIP};

		++$subid;
		
		my @conlist;
		$conlist[0] = $i;

		# Loop while we still have atoms that have not been assigned
		# to a substructure
		while (defined $conlist[0]) {

			# Loop over all atoms in conlist (ie those that have not
			# yet been assigned to a substructure)
			my @newconlist;
			foreach my $atom2num (@conlist) {

				my $atom2 = $mol->{ATOMS}[$atom2num];
		
				# Skip this atom if we have already
				# assigned it to a substructure.  This can
				# happen when an atom is included in
				# conlist twice because it is in a ring
				next if defined $atom2->{SUBID};

				# Skip atoms we have marked to skip (oxygens)
				next if $atom->{SKIP};

				# Assign catom to substructure
				$atom2->{SUBID} = $subid;

				# Find all unvisited neighbours of catom and put them into @newconlist
				# unless we are crossing to a new residuea
				LOOP: foreach my $atom3num (@{$atom2->{CONNECT}}) {

					my $atom3 = $mol->{ATOMS}[$atom3num];
			 
					next if defined $atom3->{SUBID};
		
					# Skip if we are crossing to next amide
					next if ($atom3->{FG}{S_AMIDE_C} || $atom3->{FG}{T_AMIDE_C}) && ($atom2->{FG}{S_AMIDE_N} || $atom2->{FG}{T_AMIDE_N});
					next if ($atom2->{FG}{S_AMIDE_C} || $atom2->{FG}{T_AMIDE_C}) && ($atom3->{FG}{S_AMIDE_N} || $atom3->{FG}{T_AMIDE_N});
			
					# Disulfide
					next if $atom2->{FG}{DISULFIDE_S} && $atom3->{FG}{DISULFIDE_S};
					
					# Esters and ethers
					# Skip the linking oxygen for the moment.  These are dealt with later
					
					if ($atom2->{ELEMENT} eq 'C' && $atom3->{ELEMENT} eq 'O' && ch($atom3->{CONNECTED}, 'C') == 2) {
						$atom3->{SKIP} = 1;
						next;
					}
					if ($atom3->{ELEMENT} eq 'C' && $atom2->{ELEMENT} eq 'O' && ch($atom2->{CONNECTED}, 'C') == 2){
						$atom2->{SKIP} = 1;
						next;
					}
					
					# Click linkages
					if ($atom2->{ELEMENT} eq 'C' && $atom2->{PLANAR_RING} && $atom3->{FG}{HET_AMINO_N} && ch($atom3->{CONNECTED},'N') >= 1) {
						#atom_printout($atom3);
						next;
					}
					if ($atom3->{ELEMENT} eq 'C' && $atom3->{PLANAR_RING} && $atom2->{FG}{HET_AMINO_N} && ch($atom2->{CONNECTED},'N') >= 1) {
						next;
					}
					
					# Phosphate esters
					# Phosphate only belongs to the same residue if it is connected via a methylene
					if ($atom3->{ELEMENT} eq 'P' && $atom2->{ELEMENT} eq 'O') {
					
						foreach my $atom4 (connected($atom2, $mol)) {
						
							next if $atom4 == $atom3;
							next LOOP if ch($atom4->{CONNECTED}, 'C') >= 2;
							next LOOP if ch($atom4->{CONNECTED}, 'C') == 0;
						}
					}
					
					push @newconlist, $atom3num;
				}
			}
	 
			# Assign newconlist to conlist
			@conlist = @newconlist;
		} 
		
		die if !defined $atom->{SUBID} && !$atom->{SKIP};
	}
	
	silico_msg('c', "Made $subid residues\n");
	
	# Fix any atoms we have privously marked to skip
	# For now, the linking atom is connected to the atom with (approximately) the lowest oxidation state using a rough heuristic
	# that could certainly be improved
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		next if !$atom->{SKIP};
			
		my $best = 100;
		my $winner;
		foreach my $c (connected($atom, $mol)) {
		
			my $s = 2*ch($c->{CONNECTED}, 'O') + ch($c->{CONNECTED}, 'C') - ch($c->{CONNECTED}, 'H');
			
			if ($s < $best) {
				$best = $s;
				$winner = $c;
			}
			
			$atom->{SUBID} = $winner->{SUBID};
		}
	}
	
	#Make sure that every hydrogen has the same residue name, number, segid as the parent atom
	sanitise_hydrogens($mol);
	
	# Make residue list
	my $reslist;
	@$reslist = ();
	foreach my $atom (atoms($mol)) {
		push @{$reslist->[$atom->{SUBID}-1]}, $atom;
	}
	
	# Make sure we don't have any empty residues
	# Remove empty residues and renumber
	my $t;
	my $count = 0;
	foreach my $res (@$reslist) {
		next if $#{$res} == -1;
		push @$t, $res;
		++$count;
		foreach my $atom (@$res) {
			$atom->{SUBD} = $count;
		}
	}
	$reslist = $t;
	
	# Note. At this point residues have been assigned - but not necessarily in the direction N->C
	# and possibly in some random order.
	# Work out the N->C sequence of amino acids by finding which residue is adjacent to N
	
	# Find all C->N linkages between residues
	my $nexthash;
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		# Find peptide linkages connected via amid C
		# (Could go wrong if linked through AA sidechain)
		if ($atom->{FG}{S_AMIDE_C} || $atom->{FG}{T_AMIDE_C}) {
			
			my $count = 0;
			foreach my $con (connected($atom, $mol)) {
			
				next if $con->{SUBID} == $atom->{SUBID};
				$nexthash->{$atom->{SUBID}} = $con->{SUBID};
				++$count;
			}
			silico_msg('w', "Found amino acid with $count C->N linkages\n") if $count > 1;
		}
	}
	
	# Sort into order
	my $counter;
	my %order;
	my $rnum = 1;
	
	while (1) {
	
		my @f = sort {$a<=>$b} keys (%$nexthash);
		my $start = $f[0];
		my $s = $start;
		
		# Find if we have any unassigned residues left	
		last if !defined $start;
		
		# Step through the residue list (nexthash) to find the residue chain
		while (1) {
		
			last if !$nexthash->{$s}; # There are no more residues in the chain
			$order{$s} = $rnum;
			++$rnum;
			my $olds = $s;
			$s = $nexthash->{$s};
			delete $nexthash->{$olds};
			last if $s == $start; # We have come in a loop;
		}
		
		silico_msg('d',"We have looped too many times. Dying.\n") if ++$counter == 3000; 
	}
	
	# Sort chain residues
	my $newreslist;
	@$newreslist = ();
	foreach my $atom (@{$mol->{ATOMS}}) {
		
		next if !defined $order{$atom->{SUBID}};
		
		my $resnum = $order{$atom->{SUBID}};
		push @{$newreslist->[$resnum-1]}, $atom;
	}
	
	# Add non-chain residues (e.g. non-amino acids, water)
	# Add these to the end of the chain residue list without sorting
	my $newreslist2;
	@$newreslist2 = ();
	my $rcount = -1;
	foreach my $res (@$reslist) {
	
		next if defined $order{$res->[0]{SUBID}};
		++$rcount;
		foreach my $atom (@$res) {
			push @{$newreslist2->[$rcount]}, $atom;	
		}	
	}
	
	@$newreslist = (@$newreslist, @$newreslist2);
	
	# Renumber residues
	if (1) {
	
		print "\n";
		print heading("Renumbering residues\n");
		
		my $rescount = 0;
		my ($s1, $s2) = '';
		
		foreach my $res (@$newreslist) {
	
			++$rescount;
			$s1 .= sprintf "%5d ", $res->[0]{OSUBID} || -1;
			$s2 .= sprintf "%5d ", $rescount;
			if ($rescount % 20 == 0) {
				print "Old SUBID: $s1\n";
				print "New SUBID: $s2\n\n";
				$s1 = $s2 = '';
			}
			foreach my $atom (@$res) {
				$atom->{SUBID} = $rescount;
				$atom->{SUBCOUNT} = $rescount;
			}
		}
		print "Old SUBID: $s1\n" if $s1;
		print "New SUBID: $s2\n" if $s1;
		print "\n";
	}
	
	$mol->{RESIDUES} = $newreslist;
	
        return $newreslist;
}

sub match_amino_acid {
	
	#<
	#? Match a set of atoms against an amino acid.
	#; Requires: molecule, list of atoms, hash of amino acids
	#; Returns:
	#>
	
	my $mol = $_[0];
	my $list = $_[1];
	my $aa = $_[2];
	
	my %binfo;
	my %sinfo;
	my $string;
	
	use silico_smiles;
	
	# Create smiles string from list of atoms
	my $smiles_string = mol_smiles_string($mol, $list);
	
	# Subname can sometimes be undefined
	foreach my $atom (atoms($mol)) {
		$atom->{SUBNAME} = 'UNK' if !defined $atom->{SUBNAME};
	}
	
	my $oldname = $list->[0]{SUBNAME};
	
	if (!defined $smiles_string) {
		silico_msg('e', "Could not generate a SMILES string for this residue!\n");
		return undef;
	}
	
	silico_msg('g', "SMILES string: $smiles_string\n");
	
	# List of amino acids by SMILES string. A hash.
	# The following conventions are used for suffixes:
	# - No suffix: in standard form, in a protein.
	# Side chains
	# - "H": Protonated side chain (applicable to ASP, GLU)
	# - "X": neutral side chain (applicable to ARG, LYS)
	# - "D": if HIS, protonated on ND (replaces "X")
	# - "E": if HIS, protonated on NE (replaces "X")
	# - "S": Cysteine substituted on terminal S (e.g., disulfides)
	# Backbones of N-terminal amino acids
	# - "N": N-terminal amino acid, protonated on terminal n
	# - "T": n-terminal amino acid, neuTral on terminal n
	# Backbones of C-terminal amino acids
	# - "C": C-terminal amino acid, deprotonated on terminal carboxylate
	# - "O": c-terminal amino acid, neutral
	# Backbones of stand-alone amino acids
	# - "A": Acidic isolated amino acid (fully protonated)
	# - "B": Basic isolated amino acid (fully deprotonated)
	# - "U": isolated amino acid with neUtral backbone
	# - "Z": isolated amino acid with Zwitterionic backbone
	
	$sinfo{0} = "";
	$sinfo{D} = "protonated on ND but not NE";
	$sinfo{E} = "protonated on NE but not ND";
	$sinfo{F} = "deprotonated on NE";
	$sinfo{H} = "protonated on side-chain carboxylate";
	$sinfo{J} = "deprotonated on NH";
	$sinfo{S} = "substituted on S";
	$sinfo{X} = "deprotonated on NZ";
	
	$binfo{0} = "";
	$binfo{N} = "N-terminal (protonated on backbone)";
	$binfo{T} = "N-terminal (neutral on backbone)";
	$binfo{C} = "C-terminal (deprotonated on backbone)";
	$binfo{O} = "C-terminal (neutral on backbone)";
	$binfo{A} = "free amino acid (fully acidic backbone)";
	$binfo{B} = "free amino acid (fully basic backbone)";
	$binfo{U} = "free amino acid (fully neutral backbone)";
	$binfo{Z} = "free amino acid (zwitterionic backbone)";
		
	foreach (keys %$aa) {

		next if ($smiles_string ne $_);
		
		my $array = $aa->{$_};
		my $name = $array->[0];
		
		# Additional info - sidechain and backbone
		my $s = $array->[1];
		my $b = $array->[2];
								
		if (!defined $name) {
			silico_msg('w', "SMILES string '$smiles_string' corresponds to a residue with no defined name!\n",
					"Please check SMILES data file.\n");
			next;
		}

		if (!defined $binfo{$b} || !defined $sinfo{$s}) {
			$string = "Found non-AA residue: $name";
		} else {
			$string = "Found amino acid: $name";
					
			if (defined $binfo{$b} && $binfo{$b} ne "") {
				$string .= ", $binfo{$b}";
			}
			
			if (defined $sinfo{$s} && $sinfo{$s} ne "") {
				$string .= ", $sinfo{$s}";
			}
		}

		$string .= " (smiles: $smiles_string)";
	
		# Rename the atoms in this residue. The SUBNAME and NAME fields will be altered. In addition,
		# the BACKBONE and SIDECHAIN fields will be created, with info about the backbone protonation state.
		foreach my $atom (@$list) {
			
			my $con_h = 0;
			my $i = 0;
			my $o = 0;
			
			# Determine the offset by the SMILES order.
			# This is for naming atoms.
			next if $atom->{ELEMENT} eq 'H';
			next if $atom->{ELEMENT} eq 'Du';

			if (defined $atom->{SMILES_ORDER}) {
				$o = $atom->{SMILES_ORDER} + 2;
			} else {
				my $s = atominfo($atom);
				silico_msg('w', "Atom $s is not in a SMILES string!\n");
			}

			$atom->{SUBNAME} = $name;
			$atom->{AA_BACKBONE} = $b if $b ne '*';
			$atom->{AA_SIDECHAIN} = $s if $s ne '*';

			if (defined $array->[$o]) {
				$atom->{NAME} = $array->[$o];
			} else {
				my $s = atominfo($atom);
				print "Atom $s does not exist in smiles string ($smiles_string). Renaming to 'X'\n";
				$atom->{NAME} = 'X';
			}

			my $aname = $atom->{NAME};
			my $el = $atom->{ELEMENT};
			
			$aname =~ s/^$el// unless $aname eq 'N';
			
			foreach my $con (connected($atom, $mol)) {
				++$con_h if $con->{ELEMENT} eq 'H';
			}
			
			foreach my $con (connected($atom, $mol)) {
				next if $con->{ELEMENT} ne 'H';
				++$i;
				$con->{SUBNAME} = $atom->{SUBNAME};
				$con->{AA_BACKBONE} = $atom->{AA_BACKBONE} if defined $atom->{AA_BACKBONE};
				$con->{AA_SIDECHAIN} = $atom->{AA_SIDECHAIN} if defined $atom->{AA_SIDECHAIN};
				
				# If an amide nitrogen, the name is HN2.
				if ($aname eq 'N' && $con_h == 1) {
					$con->{NAME} = 'HN2';
				} elsif ($con_h == 1) {
					$con->{NAME} = 'H'.$aname;
				} else {
					$con->{NAME} = 'H'.$aname.$i;
				}
			}
		}
		
		# Find chirality
		my ($n, $ca, $cb, $c);
		my $cb_count = 0;
		foreach my $atom (@$list) {
			
			$n = $atom if $atom->{NAME} eq 'N';
			$c = $atom if $atom->{NAME} eq 'C';
			$ca = $atom if $atom->{NAME} eq 'CA';
			$cb = $atom if $atom->{NAME} eq 'CB';
			++$cb_count if $atom->{NAME} eq 'CB';
		}
		
		if ($n && $ca && $cb && $c && $cb_count == 1) {
			my $c = aa_chirality($n, $ca, $cb, $c);
			
			$ca->{CHIRALITY} = $c;
			if ($c eq 'D') {
				foreach my $atom (@$list) {
					$atom->{SUBNAME} = "D".$atom->{SUBNAME};
				}
			}
			if ($c eq 'U') {
				sililco_msg('w', "Found residue ($ca->{SUBNAME}, $ca->{SUBID}) with uncertain chirality")
			}
		}
	}

	if (!$string) {
		
		my $chain = $list->[0]{CHAIN} || '';
		my $s = atominfo($list->[0]);
		my $string1  = "Residue  $s did not match a known residue. ";
		my $string2  = sprintf "%-50s %-4s 0 0 ", $smiles_string || '-', $list->[0]{SUBNAME} || '-';
		
		my @a = sort {($a->{SMILES_ORDER} || -1) <=> ($b->{SMILES_ORDER} || -1)} @$list;
		
		foreach (@a) {
			
			next if $_->{ELEMENT} eq 'H';
			$string2 .= sprintf "%-4s ", $_->{NAME};
		}
		
		++$mol->{UNMATCHED_RESIDUES_SMILES}{$string2};
	}
	
	my $resname = $list->[0]{SUBNAME};
	
	printf "%-5s %-30s", $resname || '-', $smiles_string;
	print " (oldname '$oldname')" if $oldname ne $resname;
	print "\n";
}



sub parse_amino_acids_smiles {
	
	#<
	#? Parse a file containing SMILES strings of amino acids,
	#  putting the entries into a hash.
	#; Requires: nothing, but takes an alternative file
	#  as an optional argument
	#; Returns: hash
	#>
	
	my $file = $_[0] || ($Silico::data_dir.'amino_acids_smiles.dat');
	
	my $aa;
	
	open_file('FILE', $file) || file_read_error($file, 1);
	
	my $i = 0;
	while (<FILE>) {
		
		my $f;
		
		++$i;
		
		my $line = $_;
		
		# Skip lines starting with hashes (comments).
		# Note that we don't want to skip everything after a hash:
		# somewhat inconveniently, a hash is SMILES code for a
		# triple bond.
		$line =~ s/^\s*#.*$//g;
		
		next if $line =~ /^\s*$/;
		
		@$f = split ' ', $line;
		
		if (!defined $f->[1]) {
			silico_msg('w', "Error in file $file: Too few elements on line $i!\n",
					"At least two elements should be present.\n",
					"This line has been skipped.\n");
			next;
		}
		
		my $f0 = shift @$f;
		
		if (defined $aa->{$f0}) {
			silico_msg('w', "Error in file $file: New definition for SMILES string on line $i!\n",
					"This line has been skipped.\n");
			next;
		} else {
			$aa->{$f0} = $f;
		}
	}
	
	close(FILE);
	
	return $aa;	
}


sub parse_residue_dictionary {
	
	#<
	#? Parse a mol2 file containing residue definitions
	#. Revised version that uses a simpler residue database and generate common protonation states N and C termini etc.
	#; Requires: Nothing
	#; Returns: Hash of molecules - key is SMILES string
	#>
	
	my $file = $Silico::data_dir.'dictionary_amino_acids.mol2';
	silico_msg('d', "Could not open file '$file'\n") if !-e  $file;
	
	my $hash;

	# Initialise hash - othewise passing the hash to a subroutine
	# does not seem to work
	%$hash = (); 

	my $mols = read_mol_any($file, undef, undef, "quiet");
	
	my @t = @$mols;
	foreach my $mol (@t) {
		my $newmols = make_alternate_aa_forms($mol);
		push (@$mols, @$newmols); # Combine arrays
	}

	$file = $Silico::data_dir.'dictionary_ligands_ions.mol2';
	silico_msg('d', "Could not open file '$file'\n") if !$file;
	
	my $ions = read_mol_any($file, undef, undef, "quiet");
	push (@$mols, @$ions); # Combine arrays

	# Make hash 
	foreach my $mol (@$mols) {

		mol_add_hydrogens($mol, 1);

		# Make list excluding added atoms from smiles generation
		my $list;
		foreach (atoms($mol)) {
			next if $_->{EXCLUDE};
			push @$list, $_;
		}

		$mol->{SMILES} = mol_smiles_string($mol, $list);
		$hash->{$mol->{SMILES}} = $mol;
	
		printf "%-60s %-60s\n", $mol->{NAME}, $mol->{SMILES};
	}

	return $hash;	
}

####################################################################################
#
#  Working towards an improved residue database
#
####################################################################################

sub match_amino_acid2 {
	
	#<
	#? Match a set of atoms against an amino acid.
	#  Revised version
	#; Requires: molecule, list of atoms, hash of amino acids
	#; Returns:
	#>
	
	my $mol = $_[0];
	my $list = $_[1];
	my $aa = $_[2];

	my %binfo;
	my $f;
	my @f;
	my $flag;
	my $key;
	my $smiles_string;
	my $string;
	
	use silico_smiles;
	
	# Create smiles string from list of atoms
	clear_smiles_sets($mol);
	$smiles_string = mol_smiles_string($mol, $list);
	
	if (!defined $smiles_string) {
		silico_msg('e', "Could not generate a SMILES string for this residue!\n");
		return undef;
	}

	#
	# Remember to remove sort later!!
	#
	foreach my $key (sort(keys %$aa)) {

		# Res is matching amino acid from template file
		my $res = $aa->{$key};
		my $subname = $res->{ATOMS}[0]{SUBNAME};

		$string .= sprintf "k %-40s s %-40s $res->{NAME}\n", "'$smiles_string'", "'$key'";

		next if ($smiles_string ne $key);
								
		$flag = 1;	
	
		# Rename the atoms in this residue. The SUBNAME and NAME fields will be altered. 
		foreach my $atom (@$list) {
			
			my $con;
			my $con_h = 0;
			my $el;
			my $o = 0;
			
			# Determine the offset by the SMILES order.
			# This is for naming atoms.
			next if $atom->{ELEMENT} eq 'H';
			next if $atom->{ELEMENT} eq 'Du';

			if (!defined $atom->{SMILES_ORDER}) {
				silico_msg('w', "Atom $atom->{NUM} (\"$atom->{NAME}\") is not in a SMILES string!\n");
			} else {
				$o = $atom->{SMILES_ORDER}-1;
			}

			$atom->{SUBNAME} = $subname;

			if (defined $res->{ATOMS}[$o]) {
				#print "o $o n $atom->{NAME} r $res->{ATOMS}[$o]{NAME}\n";
				$atom->{NAME} = $res->{ATOMS}[$o]{NAME};
			} else {
				my $s = atominfo($atom);
				print "Atom $s does not exist in smiles string ($smiles_string). Renaming to 'X'\n";
				$atom->{NAME} = 'X';
			}

			my $pos = $atom->{NAME};
			$pos =~ s/ //g;
			($pos) = $pos =~ /^.(.*)/;

			#print "b $pos\n";

			my $i = 0;
			foreach $con (connected($atom, $mol)) {
				if ($atom->{NAME} eq 'N') {
					$atom->{NAME} = 'H';
				}
				next if $con->{ELEMENT} ne 'H';
				++$i;
				$con->{SUBNAME} = $atom->{SUBNAME};
				$con->{NAME} = "H$pos$i";
			}
		}
	}

	if ($flag) {

		my $s = atominfo($list->[0]);
		print "Match: $s $smiles_string \n";

	} else {

		#foreach (@$list) {
		#	atom_printout($_);
		#}

		#print $string;
		my $s = atominfo($list->[0]);
		print "No match: $s $smiles_string \n";
	}
}

sub make_alternate_aa_forms {

	#<
	#? Generate alternate amino acid protonation states from parent structure
	#; Requires: Molecule, hash
	#; Returns: Nothing
	#>

	my $mol = $_[0];

	my $newmols;
	my $newmol;
	@$newmols = ();

	print "\n* $mol->{NAME} *\n";

	foreach my $atom (atoms($mol)) {
		$atom->{NAME} =~ s/ //g;
	}

	# Make silico_smiles take notice of formal charges
	$mol->{HAS_FORMAL_CHARGES} = 1;

	foreach my $atom (atoms($mol)) {
		$atom->{FORMAL_CHARGE} = 0;
	}
	
	#
	# Backbone states
	#

	# Neutral standard biopolymer
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	foreach my $atom (atoms($newmol)) {
		add_cn($atom, $newmol);
		add_nc($atom, $newmol);
	}

	# C-terminal deprotonated
	$newmol = deep_copy($mol);
	$newmol->{NAME} .= ' C-terminal deprotonated';
	foreach my $atom (atoms($newmol)) {
		if ($atom->{NAME} eq 'C') {
			my $oxt = mol_add_atom($newmol, 'OXT', 'O', undef, undef, undef, $atom->{SUBNAME}, $atom->{SUBID});
			bond_create_atom($newmol, $atom, $oxt, 1);
			$oxt->{FORMAL_CHARGE} = -1;
		}
		add_cn($atom, $newmol);
	}

	# C-terminal neutral
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' C-terminal neutral';
	foreach my $atom (atoms($newmol)) {
		if ($atom->{NAME} eq 'C') {
			my $oxt = mol_add_atom($newmol, 'OXT', 'O', undef, undef, undef, $atom->{SUBNAME}, $atom->{SUBID});
			bond_create_atom($newmol, $atom, $oxt, 1);
		}
		add_cn($atom, $newmol);
	}

	# C-terminal neutral missing OXT
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' C-terminal missing OXT';
	foreach my $atom (atoms($newmol)) {
		add_cn($atom, $newmol);
	}

	# N-terminal neutral
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' N-terminal neutral';
	foreach my $atom (atoms($newmol)) {
		add_nc($atom, $newmol);
		if ($atom->{NAME} eq 'N') {
			$atom->{FORMAL_CHARGE} = 0;
		}
	}

	# N-terminal protonated
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' N-terminal protonated';
	foreach my $atom (atoms($newmol)) {
		add_nc($atom, $newmol);
		if ($atom->{NAME} eq 'N') {
			$atom->{FORMAL_CHARGE} = 1;
		}
	}

	# Isolated aa zwitterion
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' isolated backbone zwitterion';
	foreach my $atom (atoms($newmol)) {
		if ($atom->{NAME} eq 'C') {
			my $oxt = mol_add_atom($newmol, 'OXT', 'O', undef, undef, undef, $atom->{SUBNAME}, $atom->{SUBID});
			bond_create_atom($newmol, $atom, $oxt, 1);
			$oxt->{FORMAL_CHARGE} = -1;
		}
		if ($atom->{NAME} eq 'N') {
			$atom->{FORMAL_CHARGE} = 1;
		}
	}

	# Isolated aa neutral
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' isolated backbone neutral';
	foreach my $atom (atoms($newmol)) {
		if ($atom->{NAME} eq 'C') {
			my $oxt = mol_add_atom($newmol, 'OXT', 'O', undef, undef, undef, $atom->{SUBNAME}, $atom->{SUBID});
			bond_create_atom($newmol, $atom, $oxt, 1);
			$oxt->{FORMAL_CHARGE} = 0;
		}
		if ($atom->{NAME} eq 'N') {
			$atom->{FORMAL_CHARGE} = 0;
		}
	}

	# Isolated aa anion
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' isolated backbone anion';
	foreach my $atom (atoms($newmol)) {
		if ($atom->{NAME} eq 'C') {
			my $oxt = mol_add_atom($newmol, 'OXT', 'O', undef, undef, undef, $atom->{SUBNAME}, $atom->{SUBID});
			bond_create_atom($newmol, $atom, $oxt, 1);
			$oxt->{FORMAL_CHARGE} = -1;
		}
		if ($atom->{NAME} eq 'N') {
			$atom->{FORMAL_CHARGE} = 0;
		}
	}

	# Isolated aa cation
	$newmol = deep_copy($mol);
	push @$newmols, $newmol;
	$newmol->{NAME} .= ' isolated backbone cation';
	foreach my $atom (atoms($newmol)) {
		if ($atom->{NAME} eq 'C') {
			my $oxt = mol_add_atom($newmol, 'OXT', 'O', undef, undef, undef, $atom->{SUBNAME}, $atom->{SUBID});
			bond_create_atom($newmol, $atom, $oxt, 1);
			$oxt->{FORMAL_CHARGE} = 0;
		}
		if ($atom->{NAME} eq 'N') {
			$atom->{FORMAL_CHARGE} = 1;
		}
	}

	#
	# Sidechain protonation states
	#

	my $count = 0;
	my @n = @$newmols;
	foreach my $m (@n) {

		next;
		mol_label_functional_group($m);

		my $i = -1;
		foreach my $atom (atoms($m)) {

			++$i;
			my $el = $atom->{ELEMENT};
			my $n = $atom->{NAME};
			$n =~ s/ //g;

			next if $el eq 'C';
			next if $el eq 'H';

			# Should modify this to skip backbone in a better way....
			next if $n eq 'N'; # Skip amide N
			next if $n eq 'O';
			next if $n eq 'CA';
			next if $n eq 'OXT';

			my $v = valence($atom);

			if ($el eq 'N' && $v < 4 && 
				($atom->{FG}{AMINO_N} || $atom->{FG}{GUANIDINIUM_N})
				) {

				my $newmol = deep_copy($m);
				$newmol->{NAME} .= " sidechain protonated $atom->{NAME}";
				$newmol->{ATOMS}[$i]{FORMAL_CHARGE} = 1;
				push @$newmols, $newmol;
			}

			#atom_printout($atom) if $atom->{ELEMENT} eq 'O';
			if ($el eq 'O' && $atom->{FG}{CARBOXYLATE_O} && !$atom->{FG}{CARBONYL_O}) {

				my $newmol = deep_copy($m);
				$newmol->{NAME} .= " sidechain deprotonated $atom->{NAME}";
				$newmol->{ATOMS}[$i]{FORMAL_CHARGE} = -1;
				push @$newmols, $newmol;
			}
		}
	}

	return $newmols;

}


sub add_cn {

	#
	# Subroutine for add_alternate_amino_acid_forms
	#

	my $atom = $_[0];
	my $mol =  $_[1];

	if ($atom->{NAME} eq 'N') {
		my $c = mol_add_atom($mol, 'CN', 'C', $atom->{X}-1, $atom->{Y}-1, $atom->{Z});
		bond_create($mol, $atom->{NUM}-1, $c->{NUM}-1, 1);
		$c->{EXCLUDE} = 1;
	}
}

sub add_nc {

	#
	# Subroutine for add_alternate_amino_acid_forms
	#

	my $atom = $_[0];
	my $mol =  $_[1];

	if ($atom->{NAME} eq 'C') {
		my $n = mol_add_atom($mol, 'NC', 'N', $atom->{X}-1, $atom->{Y}-1, $atom->{Z});
		bond_create($mol, $atom->{NUM}-1, $n->{NUM}-1, 1);
		$n->{EXCLUDE} = 1;
	}
}

####################################################################################
#
#  Misc
#
####################################################################################


sub find_substructure {
	
	#<
	#? Find the substructure around an atom.
	#. Assumes that a change in residue name or ID terminates the substructure.
	#. Assumes that bond orders are correct and that aromatic bonds and rings
	#  have been correctly identified.
	#; Requires: molecule, starting atom, maximum search depth
	#; Returns: munted SMILES string
	#>
	
	my $mol = $_[0];
	my $atom = $_[1];
	my $depth = $_[2] || 20;
	
	my $atom_tree;
	my $family_tree;
	my $i;
	my $relationships;
	
	if (!defined $atom->{TMPSUB}) {
		silico_msg('w', "A temporary substructure ID has not been defined for the atom ",
				"$atom->{NUM} (\"$atom->{NAME}\"). 1 will be used.\n");
		$atom->{TMPSUB} = 1;
	}
		
	find_atom_paths($mol, $atom, $depth);
	$relationships = construct_atom_relationships($mol, 1);
	$family_tree = construct_family_tree($relationships);
	
	$atom_tree = add_bondorders_to_atom_tree($mol, $family_tree);
	
	# Now, I need to find some way to merge branches so I don't
	# get parallel SMILES strings.
}

sub same_substructure {
	
	#<
	#? Check if two atoms are in the same substructure.
	#  Looks at SUBNAME, SUBID, SEGID and CHAIN.
	#; Requires: Atom 1, atom 2
	#; Returns: 1 if they are, 0 otherwise.
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $property;

	if (!defined $atom1 || !defined $atom2) {
		silico_msg('e', "Atoms not defined!\n");
		return undef;
	}
	
	foreach $property ("SUBNAME", "SUBID", "SEGID", "CHAIN") {
		
		if (defined $atom1->{$property} && defined $atom2->{$property}) {
			if ($atom1->{$property} ne $atom2->{$property}) {
				return 0;
			}
		} elsif (defined $atom1->{$property} && !defined $atom2->{$property}) {
			return 0;
		} elsif (defined $atom2->{$property} && !defined $atom1->{$property}) {
			return 0;
		}
	}
	
	return 1;
}




return 1;
