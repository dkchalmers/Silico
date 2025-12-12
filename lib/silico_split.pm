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
#! silico_split.pm
#. Routines to split a molecule into an ensemble of separate molecules.
#. $Revision: 1.52.2.5.2.24 $
#>

use strict;
package Silico;

sub mol_copy_rec {

	#<
	#? Recursively split atoms into a new ensemble
	#. Called from mol_split
	#. Requires: molecule, atomnumber, ensemble, submol
	#. Returns: nothing
	#>
	
	my $mol = $_[0];
	my $atom_num = $_[1];
	my $newmolecules = $_[2];
	my $submol = $_[3];
	
	my $atom = $mol->{ATOMS}[$atom_num];
	
	$atom->{SUBMOL} = $submol;
	
	# Make a copy of the atom and put it in the new ensemble
	my $newatom =  copy_atom($atom);
	$newmolecules->[$submol]{ATOMS}[$atom_num] = $newatom;
	
	foreach my $con (@{$atom->{CONNECT}}) {
		next if (defined $mol->{ATOMS}[$con]{SUBMOL});
		
		mol_copy_rec($mol, $con, $newmolecules, $submol);
	}

	return;
}

sub ensemble_split {

	#<
	#? Split an ensemble into seperate molecules by connectivity
	#. Calls molecule_split
	#; Requires: ensemble, minimum number of atoms (default 1);
	#; Returns: new ensemble containing unconnected molecules.
	#>

	my $ens = $_[0];
	my $minatoms = $_[1] ||= 1;

	my $newens;
	@$newens  = ();
	
	foreach my $mol (@$ens) {
		@$newens = (@$newens, @{molecule_split($mol, $minatoms)});
	}

	# Count molecules
	my $i = 0;
	my $smallmolcount = 0;
	foreach my $mol (@$newens) {

		++$i;
		if ($mol->{NUMATOMS} < 4) {
			++$smallmolcount;
			next;
		}
		silico_msg('g', "Molecule $i (\"$mol->{NAME}\") has $mol->{NUMATOMS} atoms and $mol->{NUMBONDS} bonds\n");
	}

	silico_msg('c', "Found a total of $i molecules\n",
			"Found $smallmolcount molecules with less than four atoms\n\n",
			"");

	return $newens;
}

sub ensemble_split_keeping_largest_mols {

	#<
	#? Retain only the largest molecule (by connectivity)
	#  from each Silico "molecule" in an ensemble
	#; Requires: ensemble, minimum number of atoms (default 1);
	#>
	
	my $mols = $_[0];

	my $i = 0;
	foreach my $mol (@$mols) {
	
		++$i;
		
		silico_msg('g', "Molecule $i of ".($#{$mols}+1)."\n");
		
		$mol = molecule_split_keeping_largest_mol($mol);
	}
}

sub molecule_split {

	#<
	#? A non-recursive subroutine to split a 'molecule' into a series of
	#  actual molecules
	#  21 October 2004: Improved to reduce memory use and produce
	#  atoms sorted in the order in which they were present in
	#  the original molecule file.
	#; Requires: Molecule, minimum molecule size to keep (default 1), flag to create 'minmal' molecules 
	#  (i.e. Dont use molecule_copy_noatoms)
	#. Note. Using the molecule_copy_noatoms is slow and memory-consuming
	#; Returns: Ensemble
	#>

	my $mol = $_[0];
	my $minatoms = $_[1] || 1;
	my $minimal = $_[2]; 

	my @newconlist;
	my $newmols;

	my $i = 0;
	foreach my $atom (atoms($mol)) {
	
		# Retain original atom number
		$atom->{ONUM} = $i;
		
		# Remove any old labels
		undef $atom->{SUBMOL};
		++$i;
	}
	
	my $skel;
	if (!$minimal) {
		$skel = molecule_copy_noatoms($mol);
	}

	$i = -1;
	my $submol = -1;
	ATOM: foreach my $atom (atoms($mol)) {
	
		++$i;
		
		# Skip atom if we have already
		# assigned this atom to a submolecule
		next if defined $atom->{SUBMOL};
		
		++$submol;

		# Copy all molecule fields except ATOMS
		if ($minimal) {
			$newmols->[$submol]=create_molecule();
		} else {
			$newmols->[$submol] = deep_copy($skel);
		}
		
		# Set new molecule name
		$newmols->[$submol]{NAME} = ($mol->{NAME} || 'MOL') ."_". (sprintf "%03d",$submol+1);
		$newmols->[$submol]{NAME} =~ s/ //g;
		
		my @conlist;
		$conlist[0] = $i;
		
		# Loop while we still have atoms that have not been assigned to a submolecule
		while (defined $conlist[0]) {
		
			# Loop over all atoms in conlist (ie those that have not
			# yet been assigned to a submolecule)
			my @newconlist;
			foreach my $catomnum (@conlist) {
			
				my $catom = $mol->{ATOMS}[$catomnum];
				
				# Skip this atom if we have already assigned it to a submol.  This can
				# happen when an atom is included in conlist twice because it is in a ring
				next if defined $catom->{SUBMOL};
				
				# Assign catom to submol
				$catom->{SUBMOL} = $submol;
				
				# Assign atom to newmolecule
				push @{$newmols->[$submol]{ATOMS}}, $catom;
				
				# Find all unvisited neighbours of catom
				# and put them into new conlist
				foreach my $conatom (@{$catom->{CONNECT}}) {

					my $con = $mol->{ATOMS}[$conatom];
					next if defined $con->{SUBMOL};
					push @newconlist, $conatom;
				}
			}
			
			# Assign newconlist to conlist
			@conlist = @newconlist;
		}
	}

	# Discard molecules with more or less than the required number of atoms
	my $newmols2;
	my $maxatoms = get_lflag('maximum-atoms');
	foreach my $nmol (@$newmols) {

	 	if ($minatoms) {
			next if $#{$nmol->{ATOMS}}+1 < $minatoms;
		} 
		if ($maxatoms) {
			next if $#{$nmol->{ATOMS}}+1 > $maxatoms;
		}
	
		push @$newmols2, $nmol;
	}

	clean_up_split_submols($newmols2, $mol);
	
	return $newmols2;
}

sub clean_up_split_submols {

	#<
	#? Clean up ensemble that results from splitting a molecule into submolecules
	#; Requires:Ensemble, original molecule
	#; Returns: Ensemble
	#>
	
	my $mols = $_[0] || croak();
	my $mol = $_[1] || croak();

	my @impropers;
	
	foreach my $nmol (@$mols) {

		my $atoms = $nmol->{ATOMS};

		# Sort atoms back into original order within each molecule
		@$atoms = sort {$a->{NUM} <=> $b->{NUM}} @$atoms;
	
		# Make translation table
		my $i = 0;
		my $ctable;
		foreach my $atom (@$atoms) {
		
			# Relate old atom number ($atom->{ONUM}) to new atom number ($i)
			$ctable->{$atom->{ONUM}} = $i;
			
			++$i;
		}
		
		# Translate old atom numbers to new atom numbers in CONECT records
		foreach my  $atom (@$atoms) {
			foreach my $connum (@{$atom->{CONNECT}}) {
			
				if (!defined $ctable->{$connum}) {

					silico_msg('w', "Connection to undefined atom $atom->{NUM}.\n",
							"Skipping.\n");
					next;
				}

				$connum = $ctable->{$connum};
			}
		}
		
		# Translate planar impropers (if any).
		if (defined $nmol->{P_IMPROPERS}) {
			
			@impropers = ();
			
			IMPR: foreach my $improper (@{$nmol->{P_IMPROPERS}}) {
				
				# Foreach atom (offset) in this improper...
				foreach my $ia (@$improper) {
					
					# ...if there is a corresponding ctable entry,
					# then change this improper entry
					# Note that entries in the ctable are OFFSETS,
					# not atom NUMBERS.
					if (defined $ctable->{($ia)}) {
						$ia = $ctable->{($ia)};
					# otherwise, this improper references an atom
					# outside this molecule
					} else {
						next IMPR;
					}
				}
				
				# Add this improper to the list of impropers in this molecule
				push @impropers, $improper;
			}
			@{$nmol->{P_IMPROPERS}} = @impropers;	
		}
		
		molecule_pack($nmol);
		
		# Renumber atom->{NUM}
		molecule_renumber($nmol);

		# Rings are now bad. 
		$mol->{BAD_RINGS} = 1;
	
		# Set numbers of atoms and bonds
		$nmol->{NUMATOMS} = $#{$nmol->{ATOMS}}+1;
		$nmol->{NUMBONDS} = molecule_count_numbonds($nmol);
		
		# Copy other properties (add here as needed)
		molecule_copy_header($mol, $nmol);
	}
	
	return $mols;
}

sub molecule_split_keeping_largest_mol {

	#<
	#? Split a Silico molecule up on the basis of connectivity
	#  and return the largest of its component molecules
	#; Requires: molecule
	#; Returns: molecule
	#>
	
	my $mol = $_[0];
	
	my $name = $mol->{NAME};

	# Split molecule
	my $ens = molecule_split($mol);
	
	silico_msg('n', "Found ".($#{$ens}+1)." fragments.\n");
	
	if ($#{$ens} > 0) {
		
		@$ens = sort {$b->{NUMATOMS} <=> $a->{NUMATOMS}} (@$ens);
		
		$mol = $ens->[0];
		$mol->{NAME} = $name;
	}
	
	return $mol;
}

sub molecule_label_submols {

	#<
	#? A non-recursive subroutine to label submolecules based on connectivity
	#. Sets $atom->{SUBMOL}.  Creates arrays of atoms in $mol->{SUBMOLS}[submol#]
	#; Requires: Molecule, flag to also set CHAIN
	#; Returns: $mol->{SUBMOLS}
	#>

	my $mol = $_[0];
	my $set_chain = $_[1];

	my $chain = 'A';

	# Molecule must have bonds!
	molecule_check_and_fix_connectivity($mol) if ! $mol->{CONNECTIVITY_CHECKED};
	
	# Remove any old labels
	foreach my $atom (atoms($mol)) {
		undef $atom->{SUBMOL};
	}

	my $i = -1;
	my $submol = -1;
	ATOM: foreach my $atom (atoms($mol)) {
	
		++$i;
		
		# Skip atom if we have already
		# assigned this atom to a submolecule
		next if defined $atom->{SUBMOL};

		++$submol;
		increment_chain($chain) if $i > 0 && $set_chain;
		
		my @conlist;
		$conlist[0] = $i;
		
		# Loop while we still have atoms that have not been assigned
		# to a submolecule
		while (defined $conlist[0]) {
		
			# Loop over all atoms in conlist (ie those that have not
			# yet been assigned to a submolecule)
			my @newconlist;
			foreach my $catomnum (@conlist) {
			
				my $catom = $mol->{ATOMS}[$catomnum];
				
				# Skip this atom if we have already
				# assigned it to a submol.  This can
				# happen when an atom is included in
				# conlist twice because it is in a ring
				next if defined $catom->{SUBMOL};
				
				# Assign catom to submol
				$catom->{SUBMOL} = $submol;
				$catom->{CHAIN} = $chain if $set_chain;
				
				# Find all unvisited neighbours of catom
				# and put them into new conlist
				foreach my $conatom (@{$catom->{CONNECT}}) {
					my $con = $mol->{ATOMS}[$conatom];
			
					next if defined $con->{SUBMOL};
								
					push @newconlist, $conatom;
				}
			}
			
			# Assign newconlist to conlist
			@conlist = @newconlist;
		}
	}
	
	foreach my $atom (atoms($mol)) {
		push @{$mol->{SUBMOLS}[$atom->{SUBMOL}]},$atom;
	}
	
	# Check to see that number of atoms in residue matches the expected value
	my $numhash;
	my $err = 0;

	foreach my $submol (@{$mol->{SUBMOLS}}) {
		my $n = $#{$submol}+1;
		my $name = $submol->[0]{SUBNAME};
		if (!defined $numhash->{$name}) {
			$numhash->{$name} = $n;
			next;
		}
		next if  $#{$submol}+1 == $numhash->{$name} ;
		silico_msg('w', "Number of atoms in submol $name ($n) does not match parent submol ($numhash->{$name})\n");
		$err = 1;
		
		foreach my $atom (@$submol) {
			$atom->{TEMP} =  $atom->{SUBMOL};
			$atom->{SEGID} = $atom->{SUBMOL};
		}
	}
	
	#if ($err) {
		#my $ens;
		#$ens->[0] = $mol;
		#write_pdb($ens, 'submols.pdb');
		##molecule_printout($mol, 'All');
		#carp();
		#silico_msg ('e', "We have some sort of error here\n");
	#}

	return $mol->{SUBMOLS};
}

sub ensemble_consolidate_small_molecules {

	#<
	#? Combine molecules with less than a specified number of atoms
	#; Requires: ensemble, min molecule size
	#; Returns: ensemble
	#>

	my $mols =  $_[0];
	my $size = $_[1] || 4;

	my $count = 0;
	my $newmols;

	my $cmols;
	@$cmols = ();

	foreach my $mol (@$mols) {

		if ($#{$mol->{ATOMS}} < $size) {

			push @$cmols, $mol;
			++$count;

		} else {

			push @$newmols, $mol;
		}
	}

	if ($count) {
		silico_msg('c', "Combining $count small molecules\n");
		my $smallmols = ensemble_consolidate($cmols);
		push @$newmols, $smallmols->[0];
	}

	return $newmols;
}

sub sdf_fastsplit {

	#<
	#? Split an sdf file into individual molecules without parsing the file
	#; Requires: input filename, output directory (optional), output style (see subroutine format_oname), 
        # max number of structures per output file, name of sdf field to rename structures, name of field to 
	# split structures, number of files to output
	#; Returns: list of files created
	#>

	my $infile = $_[0];
	my $dir = $_[1];
	my $ostyle = $_[2];
	my $size = $_[3] || 1;
	my $rename_field = $_[4] || '';
	my $split_field = $_[5] || '';
	my $nf = $_[6];
	
	my $filenames;
	my $max_structs = get_sflag('ms');

	if ($dir eq '.' || $dir eq '') {
		$dir = '';
	} else {
		silico_msg('d', "Output directory $dir not found!\n") if !-d $dir;
		$dir .= $Silico::FS if $dir !~ /$Silico::FS$/;
	}
	
	my $count;
	my $num;
	my $remainder;
	if ($nf) {
	
		if (!open_file(*INFILE, $infile)) {
			silico_msg('e', "Can not open file $infile for reading!\n");
			return undef;
		}
	
		while (my $line = <INFILE>) {
			++$count if $line =~ /^\$\$\$\$/;
		}
		close INFILE;
		
		$num = floor($count/$nf);
		$remainder = $count - $nf*$num;

		print "Total structures: $count. Structures per file: $num or ".($num+1)."\n";
	}
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Can not open file $infile for reading!\n");
		return undef;
	}
	
	my $filebase = get_filebase_short($infile);
	
	my $data = '';
	my $filecount = 0;
	my $lines = '';
	my $molcount = 0;
	my $name = '';
	my $old_split_field_val;
	my $sizecount = 0;
	
	while(1) {
	
		++$sizecount;
		++$molcount;
		
		last if $max_structs && ($molcount > $max_structs);
	
		# Read one molecule
		my ($mol_lines, $rename_field_val, $split_field_val, $name) = sdf_fastsplit_readmol($rename_field, $split_field);

		
		last if !defined $mol_lines;
		
		# Replace first line with new molecule name
		if ($split_field) {
			$mol_lines =~ s/^.*?\n//; 			# Remove first line
			$mol_lines = "$rename_field_val\n".$mol_lines;	# Replace with new name
			$name = $rename_field_val;
		}
		
		$old_split_field_val = $split_field_val if !defined $old_split_field_val;
		
		if ($split_field && ($split_field_val ne $old_split_field_val)) {
			
			++$filecount;
			#print "filecount $filecount split_field '$split_field_val' o '$old_split_field_val'\n";
			sdf_fastsplit_writemol($filebase, $dir, $ostyle, $data, $filecount, $name, $filenames);
			$data  = '';
			$old_split_field_val = $split_field_val;
			next;
		}

		#print "filecount $filecount split_field '$split_field_val' o '$old_split_field_val' \n";
		
		$data .= $mol_lines;

		# Divide structures into each file if nf flag
		if ($nf) {
			$size = $num;
			++$size if $remainder;
		}
		
		if ($size && $sizecount == $size) {
			++$filecount;
			sdf_fastsplit_writemol($filebase, $dir, $ostyle, $data, $filecount, $name, $filenames);
			$data = '';
			$sizecount = 0;
			--$remainder if $remainder;
			next;
		}
	}

	# Write any remaining data
	++$filecount;
	if ($data) {
		sdf_fastsplit_writemol($filebase, $dir, $ostyle, $data, $filecount, $name, $filenames);
	}

	return $filenames;
}

sub sdf_fastsplit_readmol {

	#<
	#? Read an individual SDF molecule from INFILE
	#; Requires: reanme field, split_field
        #; Returns: text lines of molecule, value of 'rename field' from SDF data, value of 'split field'
	#>
	
	my $rename_field = $_[0];
	my $split_field = $_[1];
	
	my $mol_lines = '';
	my $name;
	my $rename_field_val = '';
	my $split_field_val = '';

	my $i = 0;
	while (1) {
	
		++$i;
		my $line = <INFILE>;
			
		return undef if !defined $line;
			
		$name = $line if $i == 1;
		chomp $name;
	
		$mol_lines .= $line;
			
		# Find field containing molecule 'name'
		if ($rename_field && $line =~ /^>  <$rename_field>/) {
				
			$rename_field_val = <INFILE>;
			$mol_lines .= $rename_field_val;
			chomp $rename_field_val;
		}
			
		# Find field containing splitfield
		if ($split_field && $line =~ /^>  <$split_field>/) {
				
			$split_field_val = <INFILE>;
			$mol_lines .= $split_field_val;
			chomp $split_field_val;
		}
			
		# Loop until we have reached the end of the structure
		last if ($line =~ /^\$\$\$\$/);
	}
		
	return $mol_lines, $rename_field_val, $split_field_val, $name;
}

sub sdf_fastsplit_writemol {

	#<
	#? Write an individual SDF molecule 
	#; Requires: filebase, dir, output name format, text lines, filecount, name, filenames array
        #; Returns: nothing
	#>

	my $filebase = $_[0];
	my $dir = $_[1];
	my $ostyle = $_[2];
	my $lines = $_[3];
	my $filecount = $_[4];
	my $name = $_[5];
	my $filenames = $_[6];
	
	# Create dummy molecule containing only molecule name for use in format_oname
	my $dummy_mol->{NAME} = $name || 'Mol';
		 
	my $outfile = format_oname($ostyle, $dummy_mol, $filebase, $filecount, 'sdf');
	$outfile = $dir."/".$outfile if $dir;
	
	# Save a list of filenames
	push @$filenames, $outfile;
	
	open (SFILE, ">$outfile")|| silico_msg('d', "Could not create or open file $outfile for writing!\n");
	silico_msg('c', "Writing to file: $outfile\n");
	print SFILE $lines;
	close SFILE;
}


sub pdb_fastsplit {

	#<
	#? Quickly split a pdb file into component molecules/chains
	#. Only ATOM,  HETATM and CRYST lines are printed.
	#  Files are written directly to disk named  <filebase>_00001 etc, or by chain name if requested.
	#; Requires: infile, output directory ('.' by default), output style (see subroutine format_oname),
	#  flags to split on chain, END, ENDMDL and TER records, flag to name output files by chain
	#; Returns: list of files created
	#>

	my $infile = $_[0];
	my $dir = $_[1];
	my $ostyle = $_[2] || 'fnumbered';
	my $chainsplit = $_[3] || 0;
	my $endsplit = $_[4] || 0;
	my $endmdlsplit = $_[5] || 0;
	my $tersplit = $_[6] || 0;
	my $namechain = $_[7] || 0;
	
	my $filenames;
	my $hash;
	my $max_structs = get_sflag('ms');
	
	$ostyle = 'named' if $namechain;
	
	if ($dir eq '.' || $dir eq '') {
		$dir = '';
	} else {
		silico_msg('d', "Output directory $dir not found!\n") if !-d $dir;
		$dir .= $Silico::FS if $dir !~ /$Silico::FS$/;
	}
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Could not open file $infile for reading!\n");
		return undef;
	}
	
	my $filebase = get_filebase_short($infile);
	
	my $end;
	my $number = 0;
	my $numlines = 0;
	my $chain = '-';
	my $segid = '-';
	my $oldchain = '-';
	my $oldsegid = '-';
	my $string = '';
	my $watersplit = 1 if $chainsplit;
	my $waterflag = 0;
	
	while (1) {
	
		my $close = 0;
	
		$_ = <INFILE>;
		
		if (!defined $_) {
		
			$_='';
			$close = 1;
			$end = 1;
		}
	
		next if /^MODEL/;
		next if /^CONECT/;
	
		$close = 1 if (/^TER/ && $tersplit == 1);
		$close = 1 if (/^ENDMDL/ && $endmdlsplit == 1);
		$close = 1 if (/^END/ && (not /^ENDMDL/) && $endsplit == 1);
		
		if ($chainsplit && (/^ATOM/ || /^HETATM/)) {
			
			# Hack to stop substring outside string errors
			chomp;
			$_ .= "                                                                                             \n";
			$chain = substr($_,21,1) || '';	#Chain name
			$segid = substr($_,72,4) || '';	#SegID
			$chain =~ s/ //g;
			$segid =~ s/ //g;
		
			if (($numlines > 0) && (($chain ne $oldchain) || ($segid ne $oldsegid))) {
				$close = 1;
				print "$chain\n";
			}
		}
		
		# Split at the start  and end of each block of water molecules
		if ($watersplit && (/^ATOM/ || /^HETATM/)) {
		
			my $resname = substr($_,17,3)  || '-';	#Residue name
			
			if ($resname =~ /HOH/) {
			
				$close = 1 if !$waterflag;
				$waterflag = 1;  # We are in a water block
				
			} else {
			
				$close = 1 if $waterflag;
				$waterflag = 0;  # We are not in a water block
			}
		}
		
		if ($close) {
		
			$close = 0;
			++$number;
			$numlines = 0;
			
			my $x;
			my $y = '>';
			
			if ($namechain) {
				$x = $oldchain;
			} else {
				$x = $number;
			}
			# Keep track of previously opened files so that we can append if necessary
			++$hash->{$x}; 
			$y = '>>' if $hash->{$x} > 1;

			my $outfile = format_oname($ostyle, undef, $filebase, $x, 'pdb');
			$outfile = $dir."/".$outfile if $dir;
		
			if ($string) {
				open (SFILE, "$y$outfile")|| silico_msg('d', "Could not create or open file $outfile for writing!\n");
				print SFILE $string;
				close SFILE;
			
				$string = '';
			}
			
			last if $max_structs && $number == $max_structs;
			last if $end;
		}
				
		$oldchain = $chain;
		$oldsegid = $segid;
		
		# Only print ATOM, HETATOM and CRYST lines
		next if (!(/^ATOM/ || /^HETATM/ || /^CRYST/));
		
		++$numlines;
		$string .= $_;
	}
	
	@$filenames = keys (%$hash);
	
	return $filenames;
}

sub mol2_fastsplit {

	#<
	#? Split a multimolecule mol2 file into individual molecules without parsing the file
	#; Requires: infile, output directory (optional),  output style (see subroutine format_oname)
	#  max number of structures per output file
	#; Returns: list of files created
	#>

	my $infile = $_[0];
	my $dir = $_[1];
	my $ostyle = $_[2];
	my $size = $_[3]  || 1;
	
	my $filenames;
	my $max_structs = get_sflag('ms');
	my $number = 0;
	my $print = 0;
	
	if ($dir eq '.' || $dir eq '') {
		$dir = '';
	} else {
		silico_msg('d', "Output directory $dir not found!\n") if !-d $dir;
		$dir .= $Silico::FS if $dir !~ /$Silico::FS$/;
	}
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Could not open file $infile for reading!\n");
		return undef;
	}
	
	my $filebase = get_filebase_short($infile);

	my $count = 0;
	my $comments = '';
	my $flag = 0;
	while (<INFILE>) {
	
		# Keep comments
		if (/^#/) {
		
			$comments .= $_;
			next;
		}
		
		if (/^@<TRIPOS>MOLECULE/) {
		
			last if $max_structs && $number >= $max_structs;
			$print = 1;
		
			# Get name
			my $name = <INFILE>;
			$name =~ s/ //g;
			chomp $name;
			$name ||= 'Mol';
			
			# Make pseudo molecule to pass to format_oname
			my $mol;
			$mol->{NAME} = $name;
		
			# Set filename if we are the first structure or we have the
			# required number of structures
			if ($count >= $size || $flag == 0) {

				$count = 0;
			
				++$number;
				my $outfile = format_oname($ostyle, $mol, $filebase, $number, 'mol2');
				$outfile = $dir."/".$outfile if $dir;
				
				#print "outfile $outfile  number $number\n";
				
				open (SFILE, ">$outfile")|| silico_msg('d', "Could not create or open file $outfile for writing!\n");
			
				# Save a list of filenames
				push @$filenames, $outfile;
			}
			
			print SFILE $comments;
			print SFILE $_;
			print SFILE "$name\n";
		
			$comments = '';
			++$count;
			$flag = 1;
			next;
		}
	
		print SFILE $_ if $print;
	}
	
	close SFILE;
	
	return $filenames;
}

sub makedir {

        #>
        #? Make a directory
        #; Requires: filebase
        #; Returns: 1 or 0 if failed
        #>

        my $dirname = $_[0];
        my $force = $_[1] || get_sflag('force');

        my $error;

        # Check whether a directory (or file) by that name already exists.
        if (-e $dirname) {

                if (!$force) {
                        print "\n";
                        silico_msg('e', "Directory or file $dirname already exists!\n",
                                "    Use -force to overwrite.\n");
                        return 0;
                } else {
                        $error = system("rm -rf $dirname") if (-d "$dirname");
                        silico_msg('d', "Failed to remove directory $dirname!\n") if $error;
                }
        }

        # Make a new output directory
        $error = system ("mkdir $dirname");
        silico_msg('d', "Failed to create directory '$dirname'\n") if $error;

        return 1;
}



return 1;
