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
#! Silico_gromacs.pm
#? Routines specific to Gromacs
#. $Revision: 1.83.2.10.2.50 $
#>

use strict;

package Silico;

#
#  Gromacs read routines
#

sub read_gromacs_gro {

	#<
	#? Read a gromacs .gro coordinate file
	#; Requires: filename, undef, options
	#; Returns: ensemble
	#>

	my $infile        = $_[0];
	my $start         = $_[1] || get_flag('ss', 's') || 0;
	my $max_molecules = $_[2] || get_sflag('ms');
	my $options       = $_[3] || '';

	my $mols;

	$options = uc $options;
	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr = open_molfile(undef, $infile);
	return undef if !$fr;

	#silico_msg('c', "Reading Gromacs .gro format file: $infile\n") if $options !~ /\bQUIET\b/;

	my $i = -1;
	my $molcount = 0;
	while (my $mol = read_gromacs_gro_molecule($fr, $options)) {
		++$i;
		next if $start and ($i < $start);
		++$molcount;
		push @$mols, $mol;
		last if (defined $max_molecules) && $molcount == $max_molecules;
	}
	
	close_molfile($fr);

	return $mols;
}

sub read_gromacs_gro_single {

	#<
	#? Read in a single molecule record from a gromacs gro file.
	#; Requires: file record, options
	#; Returns: molecule or undef
	#>

	my $fr      = $_[0];
	my $options = uc($_[1] || '');

	return read_gromacs_gro_molecule($fr, $options);
}

sub read_gromacs_gro_molecule {

	#<
	#? Read in a single molecule record from a gromacs gro file.
	#; Requires; file_record, options
	#; Returns: molecule or undef
	#>

	my $fr      = $_[0];
	my $options = uc($_[1] || '');

	my $FH = $fr->{FH};
	my $has_velocities = 0;
	my $line;
	my $linenum = 0;
	my $mol;
	my ($x, $y, $z, $vx, $vy, $vz);

	$options = uc $options;

	#Read Header line
	$line = <$FH> || return undef;
	chomp $line;

	if ($line =~ /t\s*=/) {
		($mol->{NAME}, $mol->{TIME}) = split(/\s*t\s*=\s*/, $line);
	} else {
		$mol->{NAME} = $line;
	}

	# Remove trailing spaces
	$mol->{NAME} =~ s/\s+$//;

	# Number of atoms
	$line = <$FH>;

	if (!defined $line) {
		return undef;
	}
	chomp $line;

	$mol->{NUMATOMS} = $line;

	# Remove trailing spaces
	$mol->{NUMATOMS} =~ s/\s+$//;

	if (!defined $mol->{NUMATOMS} || $mol->{NUMATOMS} eq '' || $mol->{NUMATOMS} == 0) {
		silico_msg('e', "Molecule's atom count is zero or not defined!\n", "Skipping this file.\n");
		return undef;
	}

	my $format_width = 0;
	for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

		$line = <$FH>;
		++$linenum;

		if (!defined $line) {
			silico_msg('e', "Premature end of file!\n", "End of file was reached before $mol->{NUMATOMS} atoms were read.\n");
			return undef;
		}

		# Gromacs normaly writes out atom positions with format %8.3f
		# but it can be longer.  Measure number of places by measuring the
		# distance between decimal places.
		if ($i == 0) {

			my @f = split '', $line;
			my $flag = 0;
			foreach (@f) {

				next if ($_ ne '.') && !$flag;
				last if $_ eq '.'   && $flag;
				$flag = 1;
				++$format_width;
			}

			$format_width = 8                                                               if !$format_width;       # No decimal points found...
			silico_msg("n", "This gro file uses format \%f$format_width for coordinates\n") if $format_width != 8;

			# Find if the file has velocities
			my $l     = $line;
			my $count = ($l =~ tr/\.//);
			$has_velocities = 1 if $count == 6;
		}

		$mol->{GMX_FORMAT_WIDTH} = $format_width;

		chomp $line;

		# Prevent substring outside string errors
		$line .= '                                                                   ';

		my $atom;

		$atom->{SUBID}   = substr($line, 0,                      5);
		$atom->{SUBNAME} = substr($line, 5,                      5);
		$atom->{NAME}    = substr($line, 10,                     5);
		$atom->{NUM}     = $i + 1;
		$x               = substr($line, 20,                     $format_width);    # x
		$y               = substr($line, 20 + $format_width,     $format_width);    # y
		$z               = substr($line, 20 + 2 * $format_width, $format_width);    # z

		# Convert nm to Ang
		$atom->{X} = $x * 10;                                                       # x
		$atom->{Y} = $y * 10;                                                       # y
		$atom->{Z} = $z * 10;                                                       # z

		if ($has_velocities) {

			$vx = substr($line, 20 + 3 * $format_width, $format_width);         # x velocity
			$vy = substr($line, 20 + 4 * $format_width, $format_width);         # y velocity
			$vz = substr($line, 20 + 5 * $format_width, $format_width);         # z velocity

			# Convert to A/ps from nm/ps
			$atom->{VX} = $vx * 10;                                             # x velocity
			$atom->{VY} = $vy * 10;                                             # y velocity
			$atom->{VZ} = $vz * 10;                                             # z velocity
		}

		$atom->{LINE} = $linenum;

		$mol->{ATOMS}[$i] = $atom;
	}

	# Gromacs box vector
	# Format v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
	# v1(y), v1(z) and v2(z) must be zero
	# i.e. The first vector must lie along the X axis and the second vector must be in the xy plane
	#
	#  v1x v2x v3x
	#  0   v2y v3y
	#  0   0   v3z
	#

	$line = <$FH>;

	if (defined $line) {

		@{ $mol->{GROMACS_BOX_VECTORS} } = split ' ', $line;

		# If only the first three vectors are nonzero then
		# we have a rectangular box and we can
		# define 'mol->cell'

		my $flag = 0;
		foreach (@{ $mol->{GROMACS_BOX_VECTORS} }[ 3 .. 8 ]) {
			$flag = 1 if (defined $_ && $_ != 0);
		}

		if (!$flag) {
			@{ $mol->{CELL} } = @{ $mol->{GROMACS_BOX_VECTORS} }[ 0 .. 2 ];
			$mol->{CELL}[0]    = $mol->{CELL}[0] * 10;
			$mol->{CELL}[1]    = $mol->{CELL}[1] * 10;
			$mol->{CELL}[2]    = $mol->{CELL}[2] * 10;
			$mol->{CELL_ALPHA} = 90;
			$mol->{CELL_BETA}  = 90;
			$mol->{CELL_GAMMA} = 90;

			#silico_msg('c', "molecule cell: @{$mol->{CELL}}\n");

		} else {

			if ($options !~ /\bNOCELLWARNING\b/) {
				silico_msg(
					'w',
					"Off axis or non-rectangular cell:\n",
					"\t@{$mol->{GROMACS_BOX_VECTORS}}\n",
					"\tCell dimensions X, Y and Z are retained but treatment of periodic cells will not be correct!\n"
				);
			}

			@{ $mol->{CELL} } = @{ $mol->{GROMACS_BOX_VECTORS} }[ 0 .. 2 ];
			$mol->{CELL}[0] = $mol->{CELL}[0] * 10;
			$mol->{CELL}[1] = $mol->{CELL}[1] * 10;
			$mol->{CELL}[2] = $mol->{CELL}[2] * 10;
		}
	}

	my $g = $mol->{GROMACS_BOX_VECTORS};

	# v1_X, v1_Y, v1_Z, v2_X, v2_Y, v2_Z, v3_X, v3_Y, v3_Z
	@{ $mol->{CELL_VECTORS}[0] } = ($g->[0], $g->[3], $g->[4]);    # Vector 1
	@{ $mol->{CELL_VECTORS}[1] } = ($g->[5], $g->[1], $g->[6]);    # Vector 2
	@{ $mol->{CELL_VECTORS}[2] } = ($g->[7], $g->[8], $g->[2]);    # Vector 3

	# Work out the atom elements
	foreach my $atom (atoms($mol)) {
		atom_guess_element($atom, $mol);
	}

	mol_set_source_data($mol, $fr, 'gro');
	$mol->{NUMBONDS}      = 0;                                     # No bonds in a gro file
	$mol->{HAS_BAD_BONDS} = 1;
	++$fr->{MOL_READ_COUNT};
	$mol->{NUM} = $fr->{MOL_READ_COUNT} + 1;

	return $mol;
}

#
# RTP
#

sub read_gromacs_rtp {

	#<
	#? Read a gromacs .rtp residue type file
	#; Requires: filename, undef, options
	#; Returns: ensemble or undef if failed
	#>

	my $infile  = $_[0];
	my $options = $_[3] || '';

	my $ens;
	my $fr;
	my $mol;

	$options = uc $options;
	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	$fr = open_molfile($fr, $infile);
	return undef if !$fr;

	#silico_msg('c', "Reading Gromacs .rtp format file: $infile\n") if $options !~ /\bQUIET\b/;

	$mol = read_gromacs_rtp_molecule($fr, $options);
	$ens->[0] = $mol;
	return $ens;
}

sub read_rtp_single {

	#<
	#? Read in a single record from an rtp file.

	#; Requires: file record,  options
	#; Returns: molecule or undef
	#>

	my $fr      = $_[0];
	my $options = uc($_[1] || '');

	my $FH = $fr->{FH};
	my $hash;
	my $infile = $fr->{FILENAME};
	my $mol;
	my ($x, $y, $z);

	$options = uc $options;

	# Name
	while (1) {

		my $line = grom_readline($fr);

		if (!defined $line) {
			last;
		}

		if ($line =~ /\[ bondedtypes \]/) {

			grom_readline($fr);
			next;
		}

		if ($line =~ /\[ atoms \]/) {

			print "atoms\n";

			my $i = 0;

			while (1) {

				$line = grom_readline($fr);
				return undef if !defined $line;

				if ($line =~ /\[/) {
					silico_putline($line, $fr);
					last;
				}

				my @f = split " ", $line;

				my $atom;

				$atom->{SUBID}   = 1;
				$atom->{SUBNAME} = $mol->{NAME};
				$atom->{NAME}    = $f[0];
				$atom->{TYPE}    = $f[1];
				$atom->{CHARGE}  = $f[2];
				$atom->{NUM}     = $i + 1;
				$x               = 0;              # x
				$y               = 0;              # y
				$z               = 0;              # z

				$hash->{ $atom->{NAME} } = $i;

				$mol->{ATOMS}[$i] = $atom;
				++$mol->{NUMATOMS};
				++$i;
			}

			next;
		}

		if ($line =~ /\[ bonds \]/) {

			print "bonds\n";

			while (1) {

				$line = grom_readline($fr);
				return undef if !defined $line;

				if ($line =~ /\[/) {
					silico_putline($line, $fr);
					last;
				}

				my ($a1, $a2) = split " ", $line;

				silico_msg('d', "Creating bond to undefined atom ($a1, $a2)") if !defined $hash->{$a1} || !defined $hash->{$a2};

				bond_create($mol, $hash->{$a1}, $hash->{$a2});
				++$mol->{NUMBONDS};
			}

			next;
		}

		if ($line =~ /\[ cmap \]/) {
			print "cmap\n";
		}

		if ($line =~ /\[ impropers \]/) {
			print "impropers\n";
		}

		# If we have reached here (and mol->ATOMS is defined) than we have read the atoms, bonds, cmap and impropers records
		last if defined $mol->{ATOMS};

		if ($line =~ /\[ \w* \]/) {

			print "name\n";

			$line =~ s/\s*\[\s*//;
			$line =~ s/\s*\].*//;

			$mol->{NAME} = $line;
			next;
		}
	}

	# Work out the atom elements
	foreach my $atom (atoms($mol)) {
		atom_guess_element($atom, $mol);
	}

	mol_set_source_data($mol, $fr, 'rtp');
	$mol->{HAS_BAD_BONDORDERS} = 1;
	$mol->{HAS_NO_COORDS}      = 1;

	return $mol;
}

#
# TRR
#

sub read_gromacs_trr {

	#<
	#? Read a gromacs .trr or .xtc file using a wrapper script with trjconv
	#; Requires: filename(trr), filename(pdb,gro,etc), undef, options
	#; Returns: ensemble or undef if failed
	#>

	my $infile_trr = $_[0];
	my $infile_gro = $_[1] || get_flag("S", 's') || get_flag("struct", "l");
	my $options    = $_[3];

	my $mols;
	my $run;
	my $filebase = get_filebase($infile_trr);
	my $ofile    = "$filebase\_trjconv.pdb";

	if (!-e $ofile) {
		silico_msg('d', "Input file $infile_trr not found\n") if !-e $infile_trr;

		silico_msg('c', "Reading file $infile_trr.\n");

		foreach (qw(tpr gro tpb g96 pdb)) {
			if (-e "$filebase\.$_") {
				$infile_gro = "$filebase\.$_";
				silico_msg('c', "Using file $filebase\.$_ as structure file\n");
				last;
			}
		}
		silico_msg('c', "Corresponding structure file can be specified with the -S or --struct flags\n") if !$infile_gro;

		$infile_gro = prompt("Structure file (tpr, tpb, tpa, gro, g96 or pdb format)") if !defined $infile_gro;

		$run = "gmx trjconv -f $infile_trr -s $infile_gro -o $ofile << eof\n0\neof";

		system $run;

	} else {
		silico_msg 'c', "Output file already exists.  Using it\n";
	}

	$mols = read_mol_any($ofile, undef, undef $options);

	return $mols;
}

#
# Index files
#

sub read_gromacs_indexfile {

	#<
	#? Read a gromacs indexfile
	#; Requires: filename
	#; Returns: Molecule index
	#>

	my $filename = $_[0];

	my $atom;
	my @f;
	my $group;
	my $groupname = "System";
	my $index;

	if (!defined $filename) {
		silico_msg('e', "No index file was specified!\n");
		return undef;
	}

	if (!open_file(*INDEX, $filename)) {
		silico_msg('e', "Can not open file $filename for reading!\n");
		return undef;
	}

	#silico_msg('c', "Reading Gromacs index file $filename\n");

	# For each line in the index file:
	while (<INDEX>) {

		chomp;

		# If this line contains square brackets, it forms an index group name
		if (/^\s*\[/) {

			# So remove the square brackets, and whitespace at beginning and end...
			s/^\s*\[\s*//;
			s/\s*\].*$//;

			# ...and set the index name to what's left over
			$groupname = $_;

			# Finally, move to the next line
			next;
		}

		# Split the line up according to whitespace
		@f = split;

		# For each entry on the line, add it to the list of atom numbers
		# that are present in this current group
		foreach my $atom (@f) {
			push @{ $index->{$groupname} }, $atom;
		}
	}

	# Print out a summary, showing how many atoms in each index group
	foreach my $group (sort keys %$index) {
		silico_msg('c', "Index group \"$group\": " . ($#{ $index->{$group} } + 1) . " atoms\n");
	}

	# Return the hash of index groups and the atoms they contain
	close INDEX;
	return $index;
}

sub write_gromacs_index_file {

	#<
	#? Write gromacs index file
	#; Requires: Output filename, hash of atomlists, array of hash keys, flag to index atom numbers from zero
	#; Returns: Nothing
	#>

	my $outfile       = $_[0];
	my $atomlist_hash = $_[1];
	my $namelist      = $_[2];
	my $zero          = $_[3];

	open(OUTFILE, ">$outfile") || file_write_error($outfile, 1);

	foreach my $key (@$namelist) {

		print OUTFILE "[ $key ]\n";

		silico_msg('c', "Index group: $key, " . ($#{ $atomlist_hash->{$key} } + 1) . " atoms\n");

		my $count = 0;
		my $i     = 0;
		foreach my $val (sort { $a <=> $b } @{ $atomlist_hash->{$key} }) {

			++$i;

			# Reduce index number by 1 if indexing from zero
			$val -= 1 if $zero;
			printf OUTFILE "%6d ", $val;

			if (++$count > 10) {
				$count = 0;

				# Print a return if it is not the last value in the list
				if ($i != $#{ $atomlist_hash->{$key} } + 1) {
					print OUTFILE "\n";
				}
			}
		}
		print OUTFILE "\n";
	}

	close OUTFILE;
}

sub read_gromacs_mdp {

	#<
	#? Read an mdp file
	#; Requires: Filename
	#; Returns: Hash with key value pairs
	#>

	my $file = $_[0];
	open my $FH, '<', $file || silico_msg('d', "Can't open $file: $!");

	my %data;
	while (my $line = <$FH>) {

		chomp $line;
		next if $line =~ /^\s*$/;    # Skip blank lines
		next if $line =~ /^\s*;/;    # Skip comment lines

		my ($key, $value) = split /\s*=\s*/, $line, 2;
		my @values = split /\s+/, $value;

		if (scalar(@values) > 1) {
			$data{$key} = \@values;
		} else {
			$data{$key} = $values[0];
		}
	}
	close $FH;
	return \%data;
}

#
#
#

sub write_gromacs_gro {

	#<
	#? Write a gro file
	#; Requires: Ensemble, filename, $options.
	#; Returns: Zero if file open failed otherwise returns 1.
	#; Options:
	#, SPLIT - print out an ensemble as separate Gromacs files
	#, PRECISE - increase precision of gro file to 0.001 Angstroms
	#, VERYPRECISE - increase precision of gro file to 0.0001 Angstroms
	#; Requires: Ensemble, filename, options
	#; Returns: Undef if write failed, otherwise 1
	#>

	my $molecules = ens($_[0]);
	my $outfile   = $_[1];
	my $options   = $_[2] || '';

	my $FH;
	my $fr;
	my $split;

	$options = uc $options                        if $options;
	$outfile = get_ofilebase($molecules) . ".gro" if !$outfile;

	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Writing gro file: $outfile");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}
	$split = 1 if ($options !~ /\bSPLIT\b/);

	my $i = 0;
	my $filename;

	# Write out molecules
	foreach my $mol (@$molecules) {

		die if !defined $mol;

		if ($split && $#{$molecules} > 0) {
			++$i;
			my $count = sprintf "1%03d", $i;
			$filename = get_filebase($outfile) . "_$count." . get_fileext($outfile);
		} else {
			$filename = $outfile;
		}

		if (!defined $fr) {
			$fr = open_molfile($fr, $filename, undef, 'write', $options);
			if (!defined $fr) {
				silico_msg('e', "Could not open $filename");
				return undef;
			}
			$FH = $fr->{FH};
		}

		if (!defined $FH) {
			silico_msg('e', "Attempt to write to unopened filehandle. Perhaps you don't have write permission?\n");
			return undef;
		}

		write_gromacs_gro_molecule($fr, $mol, $options);

		close_molfile($fr) if $split;
	}

	close_molfile($fr) if !$split;

	return 1;
}

sub write_gromacs_gro_molecule {

	#<
	#? Print out a single gro format molecule
	#. Called from write_gromacs_gro
	#; Requires: Molecule, optional flag to increase precision of file
	#; Returns: nothing
	#>

	my $fr      = $_[0];
	my $mol     = $_[1];
	my $options = $_[2];

	my $FH      = $fr->{FH};
	my $warn    = 0;
	my $precise = 0;

	++$fr->{MOL_WRITE_COUNT};

	# Standard .gro files have three decimal places. This is 0.001 nm or 0.01 Angstroms.
	# PRECISE increases this to 0.001 Angstroms
	# VERYPRECISE increases this to 0.0001 Angstroms
	if ($options =~ /VERYPRECISE/) {
		$precise = 2;

		#silico_msg('c', "Increasing precision of .gro file to 0.0001 Angstroms\n");
	} elsif ($options =~ /PRECISE/) {
		$precise = 1;

		#silico_msg('c', "Increasing precision of .gro file to 0.001 Angstroms\n");
	}

	# Set width of gro floating point fields
	my $w = 8;
	$w = $mol->{GMX_FORMAT_WIDTH} if defined $mol->{GMX_FORMAT_WIDTH};
	$w = 9                        if $precise == 1;
	$w = 10                       if $precise == 2;

	my $dp  = $w - 5;
	my $dp2 = $w - 4;
	my $fmt = (("%$w.$dp" . "f") x 3) . (("%$w.$dp2" . "f") x 3);

	print $FH $mol->{NAME} || '';
	print $FH ('t= ' . $mol->{TIME}) if $mol->{TIME};
	print $FH "\n";
	print $FH $mol->{NUMATOMS} . "\n";

	my $i = 0;
	foreach my $atom (@{ $mol->{ATOMS} }) {

		++$i;
		$i = 1 if $i > 99999;    # Overflow!

		# Truncate atom or substructure names if they are too long
		my $name    = $atom->{NAME}    || 'X';
		my $subname = $atom->{SUBNAME} || 'UNK';

		if (length $name > 5) {
			$name = substr($name, 0, 5);
			$warn = 1;
		}
		if (length $subname > 5) {
			$subname = substr($subname, 0, 5);
			$warn    = 1;
		}

		printf $FH "%5d%-5s%5s%5d$fmt\n",

		  ($atom->{SUBID} || 1), $subname, $name, $i, ($atom->{X} || 0) / 10,    # x (convert to nm from Ang)
		  ($atom->{Y}  || 0) / 10,                                               # y
		  ($atom->{Z}  || 0) / 10,                                               # z
		  ($atom->{VX} || 0) / 10,                                               # x velocity  (convert to nm/ps from Ang/ps)
		  ($atom->{VY} || 0) / 10,                                               # y velocity
		  ($atom->{VZ} || 0) / 10;                                               # z velocity
	}

	silico_msg('w', "Atom or substructure names have been truncated in .gro file\n") if $warn;

	if ($mol->{GROMACS_BOX_VECTORS}[0]) {

		#silico_msg('c', "Using existing gromacs box vectors\n");
		for (3 .. 8) {
			$mol->{GROMACS_BOX_VECTORS}[$_] ||= 0;
		}

	} elsif ($mol->{CELL}) {

		if (!defined $mol->{CELL_ALPHA} || !defined $mol->{CELL_BETA} || !defined $mol->{CELL_GAMMA}) {

			silico_msg('c', "Write_gro: Cell dimensions are present, but angles are missing.  Assuming cuboid cell\n");
			for (0 .. 2) {
				$mol->{GROMACS_BOX_VECTORS}[$_] ||= $mol->{CELL}[$_] / 10;
			}
			for (3 .. 8) {
				$mol->{GROMACS_BOX_VECTORS}[$_] = 0;
			}

		} elsif ($mol->{CELL_ALPHA} == 90 && $mol->{CELL_BETA} == 90 && $mol->{CELL_GAMMA} == 90) {

			silico_msg('c', "Write_gro: Using speified cuboid cell\n");
			for (0 .. 2) {
				$mol->{GROMACS_BOX_VECTORS}[$_] ||= $mol->{CELL}[$_] / 10;
			}
			for (3 .. 8) {
				$mol->{GROMACS_BOX_VECTORS}[$_] = 0;
			}

		} else {

			foreach (@{ $mol->{CELL} }, $mol->{CELL_ALPHA}, $mol->{CELL_BETA}, $mol->{CELL_GAMMA}) {
				$_ =~ s/ //g;
			}

			silico_msg('c',
				    "Write_gro: Calculating gromacs vectors using cell: "
				  . "a=$mol->{CELL}[0], b=$mol->{CELL}[1], c=$mol->{CELL}[2], "
				  . "alpha=$mol->{CELL_ALPHA}, beta=$mol->{CELL_BETA}, gamma=$mol->{CELL_GAMMA}\n");

			@{ $mol->{GROMACS_BOX_VECTORS} } = make_gromacs_cell_from_pdb_cell_and_angles($mol->{CELL}[0], $mol->{CELL}[1], $mol->{CELL}[2],
				$mol->{CELL_ALPHA}, $mol->{CELL_BETA}, $mol->{CELL_GAMMA});
		}

	} else {

		# If not defined then set the box size to 100 x 100 x 100 Angstroms
		silico_msg('w', "No cell set. Seeting box vectors to 10 10 10.\n");
		@{ $mol->{GROMACS_BOX_VECTORS} } = (10, 10, 10, 0, 0, 0, 0, 0, 0);
	}

	#print " box vectors $mol->{GROMACS_BOX_VECTORS}\n";
	#foreach  (@{$mol->{GROMACS_BOX_VECTORS}}) {
	#print "'$_'";
	#}
	#print "\n";

	# Gromacs box vectors v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
	printf $FH "%f %f %f %f %f %f %f %f %f\n", @{ $mol->{GROMACS_BOX_VECTORS} };
}

sub make_gromacs_cell_from_pdb_cell_and_angles {

	#<
	#? Convert PDB cell information into gromacs box vectors
	#; Requires: vector - a, b, c, alpha, beta, gamma
	#; Returns: gromacs vector (6 numbers)
	#>

	my ($a, $b, $c, $alpha, $beta, $gamma) = @_;

	my $pi = 3.14159265359;

	$a /= 10;
	$b /= 10;
	$c /= 10;

	$alpha = ($alpha / 180) * $pi;
	$beta  = ($beta / 180) * $pi;
	$gamma = ($gamma / 180) * $pi;

	my $v = sqrt(1 - cos($alpha)**2 - cos($beta)**2 - cos($gamma)**2 + 2 * cos($alpha) * cos($beta) * cos($gamma)) * $a * $b * $c;

	return $a, $b * sin($gamma), $v / ($a * $b * sin($gamma)), 0, 0, $b * cos($gamma), 0, $c * cos($beta),
	  ($c / sin($gamma)) * (cos($alpha) - cos($beta) * cos($gamma));
}

sub mol_move_into_cell_using_cell_vectors {

	#<
	#? Translate a molecule into a unit cell using cell_vecltors
	#; Requires: Molecule, cell, flag to move based on position the 'middle' atom of the molecule rather than centroid (much faster).
	#; Returns: Nothing
	#>

	my $mol = $_[0];
	my $cv  = $_[1] = $mol->{CELL_VECTORS} || (molecule_printout($mol) && silico_msg('d', "No cell vectors\n"));

	print_cell_vector($cv);

	my @centre = centroid($mol);

	foreach my $i (0 .. 2) {

		my $g = 0;

		print "i $i centre $centre[$i] cv $cv->[$i][$i]\n";

		while ($centre[$i] < 0) {
			@centre = vector_add(@centre, @{ $cv->[$i] });
			++$g;
		}
		while ($centre[$i] >= $cv->[$i][$i]) {
			@centre = vector_add(@centre, scale_vector(-1, @{ $cv->[$i] }));
			--$g;
		}

		print "i $i g $g\n";

	}
}

sub print_cell_vector {

	my $cv = $_[0];

	foreach (@$cv) {
		foreach (@{$_}) {
			print "$_ ";
		}
		print "\n";
	}
}

#
# ITP & TOP files
#

sub read_gromacs_itp {

	#<
	#? Read a Gromacs .itp or .top topology file
	#. See read_gromacs_top
	#>

	return read_gromacs_top(@_);
}

sub read_gromacs_top {

	#<
	#? Read a Gromacs .itp or .top topology file
	#  Not all block types are read (yet)
	#. Note: data read from the system-wide blocks ('defaults', 'atomtypes', etc.) are stored in the GMX_TOPDATA.
	#  This information is also linked to SDF_DATA
	#  in the ensemble
	#; Requires: inputfile, undef, undef, options
	#; Returns: molecule
	#>

	my $infile  = $_[0];
	die if $_[1];
	die if $_[2];
	my $options = uc($_[3] || '');
	
	my $fr = open_molfile(undef, $infile);
	return undef if !$fr;
	
	# List of blocks
	my @blocks = qw(defaults moleculetype settles atomtypes atoms bonds pairs angles dihedrals position_restraints dihedral_restraints molecules system exclusions cmap system molecules);    
	
	my $mols = ();
	my $globaldata = {};
	my $block     = '';
	my $atomcount = 0;
	my $molcount = -1;
	my $mol;
	
        LOOP: while (my $line = grom_readline($fr)) {

		#print "itp_mol line $line\n" if $line =~ /\[/;
		
		# Check to see if we have found a block keyword
		if ($line =~ /\[/) {

			foreach my $t (@blocks) {

				if ($line =~ /\[\s*\Q$t\E\s*/) {

					$block = $t;
					next LOOP;
				}
			}
			
			silico_msg('d', "Unrecognised ITP block type $line in file $fr->{FILENAME} line $fr->{LINECOUNT}\n");
		}

		#
		# Defaults
		#
		
		if ($block eq "defaults") {

			my @f = split(' ', $line);
			
			my $r;
			
			$r->{NB_FUNC}   = $f[0] || '';
			$r->{COMB_RULE} = $f[1] || '';
			$r->{GEN_PAIRS} = $f[2] || '';
			$r->{FUDGE_LJ}  = $f[3] || '';
			$r->{FUDGE_QQ}  = $f[4] || '';

			$globaldata->{GMX_DEFAULTS} = $r;
			next;
		}

		#
		# Atomtypes
		#
		
		if ($block eq "atomtypes") {
		
			my @f = split(' ', $line);
			
			my $atype;

			$atype->{GMX_TYPE}     = $f[0];
			$atype->{ELEMENT_NUM}  = $f[1] || 'UNK';
			$atype->{MASS}         = $f[2] || 1;
			$atype->{CHARGE}       = $f[3] || 0;
			$atype->{PTYPE}        = $f[4] || 'UNK';
			$atype->{SIGMA}        = $f[5] || 0;
			$atype->{EPSILON}      = $f[6] || 0;
			$atype->{ELEMENT} = element_number2symbol($atype->{ELEMENT_NUM});

			$globaldata->{GMX_ATOMTYPES} ||= ();
			push @{ $globaldata->{GMX_ATOMTYPES} }, $atype;
			next;
		}

		#
		# System name
		#
		
		if ($block eq "system") {

			my @f = split(' ', $line);
			
			$globaldata->{GMX_SYSTEM} = $f[0];
			next;
		}

		#
		# Moleculs = numbers of residues/molecules
		#
		
		if ($block eq "molecules") {

			my @f = split(' ', $line);
			
			$globaldata->{GMX_MOLECULES}{$f[0]} = $f[1];
			next;
		}

		#
		# Moleculetypes = residues
		#

		if ($block eq "moleculetype") {
		
			my @f = split(' ', $line);
			
			++$molcount;
			
			$mols->[$molcount] ||= {};
			
			$mol = $mols->[$molcount];
			
			$mol->{GMX_NAME} = $mol->{NAME} = $f[0] || 'UNK';
			$mol->{GMX_NREXCL} = $f[1] || ++$mol->{WARN}->{'nrxcl not set'};
			$mol->{NUMATOMS} = 0;
			$mol->{ATOMS} = ();
			$mol->{HAS_NO_COORDS}      = 1;
			$mol->{HAS_BAD_BONDORDERS} = 1;
			$mol->{IS_TOPOLOGY_FILE}   = 1;
			$mol->{SOURCE_FILE_TYPE}   = $fr->{FORMAT};
			$mol->{SOURCE_FILE_NAME}   = $fr->{FILENAME};
			++$fr->{MOL_READ_COUNT};
			
			next;
		}

		#
		# Atoms
		# 

		if ($block eq 'atoms') {

			my @f = split(' ', $line);
			my $atom;
			++$atomcount;

			$atom->{NUM}          = $f[0];
			$atom->{GMX_TYPE}     = $f[1] || 'UNK';
			$atom->{TYPE}         = $f[1] || 'UNK';
			$atom->{SUBID}        = $f[2] || 1;
			$atom->{SUBNAME}      = $f[3] || 'UNK';
			$atom->{NAME}         = $f[4] || 'UNK';
			$atom->{CHARGE_GROUP} = $f[5] || 1;
			$atom->{CHARGE}       = $f[6] || 999;
			$atom->{MASS}         = $f[7] || 0;
			$atom->{TYPE_B}       = $f[8] if defined $f[8];
			$atom->{CHARGE_B}     = $f[9] if defined $f[9];
			$atom->{MASS_B}       = $f[10] if defined $f[10];
			$atom->{X}            = 0;
			$atom->{Y}            = 0;
			$atom->{Z}            = 0;

			atom_guess_element($atom, $mol);
			push @{ $mol->{ATOMS} }, $atom;
			++$mol->{NUMATOMS};

			silico_msg('w', "Atom $atomcount does not have correct atom number. It is $f[0]\n") if $f[0] != $atomcount;
			next;
		}

		#
		# Bonds
		# 

		if ($block eq 'bonds') {

			if (!defined $mol->{ATOMS}) {
				molecule_printout($mol);
				silico_msg('d', "[ atoms ] line not found before [ bonds ]\n");
			}

			my @f = split(' ', $line);
			itp_make_bond($mol, @f);
			bond_create($mol, $f[0]-1, $f[1]-1, 1, 1);
			++$mol->{NUMBONDS};
			next;
		}

		#
		# Pairs
		# 

		if ($block eq 'pairs') {

			my @f = split(' ', $line);
			itp_make_pair($mol, @f);
			next;
		}

		if ($block eq 'angles') {

			my @f = split(' ', $line);
			itp_make_angle($mol, @f);
			next;
		}

		if ($block eq 'dihedrals') {

			my @f = split(' ', $line);
			itp_make_dihedral($mol, @f);
			next;
		}
		
		if ($block eq 'dihedral_restraints') {

			my @f = split(' ', $line);
			itp_make_dihedral_restraint($mol, @f);
			next;
		}
		
		if ($block eq 'settles') {
		
			my @f = split(' ', $line);
			itp_make_settle($mol, @f);
			next;
		}

		if ($block eq 'exclusions') {

			my @f = split(' ', $line);
			my $exclusion;

			@{ $exclusion->{ATOMS} } = ($f[0] - 1, $f[1] - 1);
			$exclusion->{GMX_FUNC} = $f[2];
			$exclusion->{TYPE}     = 'exclusion';
			push @{ $mol->{GMX_EXCLUSIONS} }, $exclusion;
			next;
		}

		if ($block eq 'cmap' || $block eq 'position_restraints' || $block eq 'system' || $block eq 'molecules') {

			++$mol->{WARN}->{"Ignoring lines in block \"$block\". Be careful. This may significantly affect your .itp file"};
		}
	}
	
	
	foreach my $mol (@$mols) {
		$mol->{GMX_TOPDATA} = $globaldata;
		mol_print_warnings($mol);
	}
	
	return $mols;
}

sub make_gromacs_itp {

	#<
	#? Make basic gromacs itp data where it does not exist
	# Note: Blocks [ position_restraints ], [ molecules ] and [ system ] are ignored
	#; Requires: molecule
	#; Returns: molecule
	#>

	my $mol = $_[0];
	
	atomic_masses();
	molecule_check_and_fix_connectivity($mol);
	mol_find_all_bonds($mol);
	mol_find_all_angles($mol);
	mol_find_all_dihedrals($mol);

	foreach my $atom (atoms($mol)) {

		$atom->{MASS}          = $Silico::Atomic_Masses[$atom->{ELEMENT_NUM}];
		$atom->{TYPE}         ||= $atom->{ELEMENT} || 'UNK';
		$atom->{CHARGE_GROUP} ||= 1;
		$atom->{CHARGE}       ||= 0;
		$atom->{TYPE_B}       ||= '';
		$atom->{CHARGE_B}     ||= '';
		$atom->{MASS_B}       ||= '';
		$atom->{X}            ||= 0;
		$atom->{Y}            ||= 0;
		$atom->{Z}            ||= 0;
	}

	foreach my $bond (@{ $mol->{BONDS} }) {
	
		# Mol, a1, a2, func, Func, Length, Fconst
		tp_make_bond($mol, $bond->[0]+1, $bond->[1]+1, 2, 999, 999);
	}

	#
	# Note GMX_PAIRS is not created
	#

	foreach my $angle (@{ $mol->{ANGLES} }) {

		itp_make_angle($mol, $angle->[0]+1, $angle->[1]+1, $angle->[2]+1, 1, 999, 999);
	}

	foreach my $angle (@{ $mol->{DIHEDRALS} }) {

		my $a = $angle->{ATOMS};
	
		itp_make_dihedral($mol, $a->[0]+1, $a->[1]+1, $a->[2]+1, $a->[3]+1, 1, 999, 999, 999);
	}
	
	#
	# Not made
	# exclusions, cmap, position_restraints, dihedral_restraints, system, molecules
	#

	my $ens;
	$ens->[0] = $mol;
	return $ens;
}

sub itp_make_bond {

	#<
	#? Make a Gromacs itp bond
	#. Note: Uses atom NUMBERS (not atom indexes).
	#; Requires molecule, atom Number 1, atom Number 2, type, etc
	#
	#>

	my $mol = shift(@_);
	my $a1  = shift(@_) - 1;    # Atom1 index
	my $a2  = shift(@_) - 1;    # Atom2 index
	my @f   = @_;

	my $bond;

	$bond->{NUM} = $#{ $mol->{GMX_BONDS} } + 1;
	@{ $bond->{ATOMS} } = ($a1, $a2);
	$bond->{GMX_FUNC} = $f[0];
	$bond->{TYPE}     = 'bond';

	if (!defined $f[1]) {

		# No parameters provided
		$bond->{FORMAT} = 'null';

	} elsif (defined $f[1] && !defined $f[2]) {

		$bond->{CODE}   = $f[1];
		$bond->{FORMAT} = 'code';

	} else {

		# Full format
		$bond->{GMX_FC}   = $f[1];
		$bond->{LENGTH}   = $f[2];
		$bond->{GMX_FC_B} = $f[3];
		$bond->{LENGTH_B} = $f[4];
		$bond->{FORMAT}   = 'full';
	}

	#print "bond format $bond->{FORMAT}\n";

	push @{ $mol->{GMX_BONDS} }, $bond;

	return 1;
}

sub itp_make_pair {

	#<
	#? Make a Gromacs itp pair
	#. Note: Uses atom numbers (not atom indexes).
	#; Requires molecule, atom Number 1, atom Number 2, type
	#
	#>

	my $mol = shift;
	my @f   = @_;
	my $pair;

	$pair->{NUM} = $#{ $mol->{GMX_PAIRS} } + 1;
	@{ $pair->{ATOMS} } = ($f[0] - 1, $f[1] - 1);
	$pair->{GMX_FUNC} = $f[2];
	$pair->{TYPE}     = 'pair';
	push @{ $mol->{GMX_PAIRS} }, $pair;

	return 1;
}

sub itp_make_angle {

	#<
	#? Make a Gromacs itp angle

	#. Note: Uses atom numbers (not atom indexes).
	#; Requires molecule, atom Number 1, atom Number 2,  atom Number 3,  angle type, format, etc
	#
	#>

	my $mol = shift;
	my @f   = @_;

	my $angle;

	if (!defined $f[3] || $f[3] !~ /^[0-9]/) {

		silico_msg('d', "Missing angle function type: @f\n");
	}

	my $func = $f[3];

	$angle->{NUM} = $#{ $mol->{GMX_ANGLES} } + 1;
	@{ $angle->{ATOMS} } = ($f[0] - 1, $f[1] - 1, $f[2] - 1);
	$angle->{GMX_FUNC} = $f[3];
	$angle->{TYPE}     = 'angle';

	if (!defined $f[4]) {

		$f[4] = '';
	}

	if ($f[4] eq '') {

		# Empty value ( = use defaults )
		# Note: we do not interpret the default values.
		$angle->{FORMAT} = 'null';

	} elsif (defined $f[4] && !defined $f[5]) {

		# Has an alphanumeric code
		# Note: we do not interpret the codes. We just parrot them back.
		$angle->{CODE}   = $f[4];
		$angle->{FORMAT} = 'code';

	} else {
		$angle->{ANGLE}    = $f[4];
		$angle->{GMX_FC}   = $f[5];
		$angle->{ANGLE_B}  = $f[6];
		$angle->{GMX_FC_B} = $f[7];
		$angle->{FORMAT}   = 'full';
	}

	#print "angle format $angle->{FORMAT}\n";

	push @{ $mol->{GMX_ANGLES} }, $angle;

	return 1;
}

sub itp_make_dihedral {

	#<
	#? Make a Gromacs itp dihedral (mol->GMX_DIHEDRALS)
	#. Note: Uses atom numbers (not atom indexes).
	#; Requires molecule, atom Number 1, atom Number 2,  atom Number 3,  atom Number 4, 
	#  dihedral type (1=proper, 2=improper), format, etc
	#. Note. We do not read individual values for default values or where codes are used to define parameters
	#>

	my $mol = shift;
	my @f   = @_;

	my $angle;

	if (!defined $f[4] || $f[4] !~ /^[0-9]/) {

		silico_msg('d', "Missing angle function type: @f\n");
	}

	my $func = $f[4];

	$angle->{NUM} = $#{ $mol->{GMX_DIHEDRALS} } + 1;
	@{ $angle->{ATOMS} } = ($f[0] - 1, $f[1] - 1, $f[2] - 1, $f[3] - 1);
	$angle->{GMX_FUNC} = $f[4];

	if (!defined $f[5]) {

		$f[5] = '';
	}

	if ($f[5] eq '') {

		# Empty value ( = use defaults )
		# Note: we do not interpret the default values.
		$angle->{FORMAT} = 'null';

	} elsif (defined $f[5] && !defined $f[6]) {

		# Has an alphanumeric code
		# Note: we do not interpret the codes. We just parrot them back.
		$angle->{CODE}   = $f[5];
		$angle->{FORMAT} = 'code';

	} elsif ($func == '1' || $func == 9) {

		# Propers
		$angle->{TYPE}   = 'proper';
		$angle->{FORMAT} = 'full';

		$angle->{ANGLE}   = $f[5];
		$angle->{PHASE}   = $f[6];
		$angle->{MULT}    = $f[7];
		$angle->{ANGLE_B} = $f[8];
		$angle->{PHASE_B} = $f[9];
		$angle->{MULT_B}  = $f[10];

	} elsif ($func == 2 || $func == 4) {

		# Impropers
		$angle->{TYPE}   = 'improper';
		$angle->{FORMAT} = 'full';

		$angle->{ANGLE}    = $f[5];
		$angle->{GMX_FC}   = $f[6];
		$angle->{MULT}     = $f[7];
		$angle->{ANGLE_B}  = $f[8];
		$angle->{GMX_FC_B} = $f[9];
		$angle->{MULT_B}   = $f[10];

	} elsif ($func == '3') {

		$angle->{TYPE}   = 'proper';
		$angle->{FORMAT} = 'full';

		foreach (0 .. 5) {
			$angle->{RB}[$_] = $f[ 5 + $_ ];
		}

	} else {

		silico_msg('w', "Found dihedral angle type '$func'. Dihedral functions other than types 1-4 and 9 are not implemented\n");
	}

	#print "dihedral format $angle->{FORMAT}\n";

	push @{ $mol->{GMX_DIHEDRALS} }, $angle;

	return 1;
}

sub itp_make_dihedral_restraint {

	#<
	#? Make a Gromacs itp dihedral restraint (mol->GMX_DIHEDRAL_RESTRAINTS)
	#. Note: Uses atom numbers (not atom indexes).
	#; Requires molecule, atom Number 1, atom Number 2,  atom Number 3,  atom Number 4, 
	#  function type ( usually 1)
	#>

	my $mol = shift;
	my @f   = @_;

	my $angle;
	my $func = $f[4];

	if (!defined $func || $func !~ /^[0-9]/) {

		silico_msg('d', "Missing or nonnumeric angle function type: @f\n");
	}

	if ($func != 1) {

		silico_msg('w', "Found dihedral angle restraint type '$func'. Dihedral functions other than type 1 are not implemented\n");
	}

	$angle->{NUM} = $#{ $mol->{GMX_DIHEDRAL_RESTRAINTS} } + 1;
	@{ $angle->{ATOMS} } = ($f[0] - 1, $f[1] - 1, $f[2] - 1, $f[3] - 1);
	$angle->{GMX_FUNC} = $f[4];
	$angle->{ANGLE} = $f[5];
	$angle->{DELTA_ANGLE} = $f[6];
	$angle->{GMX_FC} = $f[7];

	push @{ $mol->{GMX_DIHEDRAL_RESTRAINTS} }, $angle;

	return 1;
}

sub itp_restrain_dihedrals {

	#<
	#? Make a set of dihedral restraints
	#; Requires: molecule, dihedral list, force constant, angle delta
	#>
	
	my $mol = $_[0];
	my $dihedrals = $_[1];
	my $fc = $_[2] || 10;
	my $delta = $_[3] || 0;
	
	foreach my $dihedral (@$dihedrals) {
	
		print "@$dihedral\n";
	}
	
	molecule_printout($mol);
	die;
}

sub itp_make_settle {

	#<
	#? Make a Gromacs settle
	#. Not implemented
	#
	#>

	my $mol = shift;
	my @f   = @_;

	if (!defined $f[0]) {
	
		silico_msg('d', "Empty value found\n");
	}

	$mol->{GMX_SETTLES}{I} = $f[0];
	$mol->{GMX_SETTLES}{FUNC} = $f[1];
	$mol->{GMX_SETTLES}{D_OH} = $f[2];
	$mol->{GMX_SETTLES}{D_HH} = $f[3];

	return;
}


#
# ITP maninipulation
#

sub itp_delete_parameters_atom {

	#<
	#? Delete Gromacs GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS or GMX_EXCLUSIONS
	#  where the specified atom is involved
	#; Requires molecule, thing, list of atom NUMBERs (not indexes)
	#>

	my $mol     = shift;
	my $thing   = shift;
	my @numbers = @_;

	my $new;
	my $skipped = 0;

	if (!defined $mol->{$thing}) {
		silico_message('w', "$thing not defined\n");
		return;
	}

      THING: foreach my $t (@{ $mol->{$thing} }) {
		foreach my $n (@{ $t->{ATOMS} }) {
			foreach my $num (@numbers) {
				if ($n == $num - 1) {

					print "a @{$t->{ATOMS}}\n";

					++$skipped;
					next THING;
				}
			}
		}

		push @$new, $t;
	}

	$mol->{$thing} = $new;

	print "Deleted $skipped $thing from atoms: " . join(', ', @numbers) . "\n";
}

sub itp_delete_angle_parameters_atom {

	#<
	#? Delete Gromacs GMX_ANGLES where the specified atom is a CENTRAL atom
	#; Requires molecule, list of atom NUMBERS (not indexes)
	#>

	my $mol     = shift;
	my @numbers = @_;

	my $skipped = 0;
	my $new;

      ANGLE: foreach my $ang (@{ $mol->{"GMX_ANGLES"} }) {
		my $n = $ang->{ATOMS}[1];
		foreach my $num (@numbers) {

			if ($n == $num - 1) {
				++$skipped;
				next ANGLE;
			}
		}

		push @$new, $ang;
	}

	$mol->{"GMX_ANGLES"} = $new;
	print "Deleted $skipped GMX_ANGLES from atom(s): " . join(', ', @numbers) . "\n";
}

sub pack_gromacs_angles {

	#<
	#? Pack new format gromacs bonds and angles, etc. following atom deletion
	#. Atoms are renumbered
	#. Any bond or angle with undefined values is deleted.
	#; Requires: Molecule, ctable (array of old and new atom numbers);
	#>

	my $mol    = $_[0];
	my $ctable = $_[1];

	my @list  = qw (GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS GMX_EXCLUSIONS);
	my $count = 0;

	#silico_msg('c', "Renumbering gromacs bonds, pairs, etc\n");

	foreach my $l (@list) {

		next if !defined $mol->{$l};
		my $new_l;

	      GMX_THING: foreach my $n (@{ $mol->{$l} }) {

			foreach (@{ $n->{ATOMS} }) {

				if (!defined $ctable->[$_]) {
					++$count;
					next GMX_THING;
				}

				$_ = $ctable->[$_];
			}

			push @$new_l, $n;
		}

		$mol->{$l} = $new_l;
	}
}

sub grom_readline {

	#<
	#? Read a gromacs .itp topology file line
	#. Filters  out blank lines and comments
	#; Requires: file_record
	#; Returns: line
	#>

	my $fr = $_[0];
	
	my $line;

	while (1) {

		$line = silico_readline($fr);

		return undef if !defined $line;

		chomp $line;

		# Remove everything after semicolon (a coment)
		$line =~ s/;.*//;

		# Remove everything after a hash (a coment)
		$line =~ s/#.*//;

		# Empty blank lines
		$line =~ s/^\s*$//;

		last if $line;
	}

	return $line;
}

sub grom_putline {

	#<
	#? Return a gromacs line
	#; Requires: line, file_record
	#; Returns: nothing
	#>
	
	my $line = $_[0];
	my $fr = $_[1] || croak();
	
	silico_putline($line, $fr);
}

sub gmx_make_heavy_h {

	#<
	#? Increase the mass of hydrogen atoms by a scale factor (default 4) and reduce the
	#  mass of the attached heavy atoms
	#; Requires: molecule, scale factor
	#; Returns: nothing
	#>

	my $mol     = $_[0];
	my $scale_h = $_[1] || 4;

	# Modify mass of hydrogen atoms
	foreach my $atom (atoms($mol)) {

		if ($atom->{MASS} < 1.1) {
			$atom->{MASS}   = $scale_h * 1.008;
			$atom->{HEAVYH} = 1;
		}
		if ($atom->{MASS_B} && $atom->{MASS_B} < 1.1) {
			$atom->{MASS_B}   = $scale_h * 1.008;
			$atom->{HEAVYH_B} = 1;
		}
	}

	# Modify mass of nonhydrogen atoms
	foreach my $atom (atoms($mol)) {

		my $hcount   = 0;
		my $hcount_b = 0;

		foreach my $con (connected($atom, $mol)) {
			++$hcount   if $con->{HEAVYH};
			++$hcount_b if $con->{HEAVYH_B};
		}

		# Subtract added mass of hydrogens, but only if parent atom has enough mass
		# this is necessary to deal with unusual topologies
		my $m = ($scale_h - 1) * 1.008 * $hcount;
		if ($hcount && $atom->{MASS} > $m + 4) {
			$atom->{MASS} -= $m;
		}

		# Subtract added mass of hydrogens, but only if parent atom has enough mass
		# this is necessary to deal with unusual topologies
		my $m_b = ($scale_h - 1) * 1.008 * $hcount_b;
		if ($hcount_b && $atom->{MASS_B} > $m_b + 4) {
			($atom->{MASS_B} -= $m_b);
		}

		if ($atom->{MASS} < 0) {
			
			silico_msg('w', "Atom with negative mass!\n");
			foreach my $con (connected($atom, $mol)) {
				print "el $con->{ELEMENT}\n";
			}
			silico_msg('d', "Attempt to create atom with mass less than 0");
		}
		if ($atom->{MASS_B} && $atom->{MASS_B} < 0) {
		
			silico_msg('w', "Atom with negative mass!\n");
			atom_printout($atom);
			silico_msg('d', "Attempt to create atom with mass_b less than 0");
		}
	}
}

#
# Write Gromacs itp/top 
#

sub write_gromacs_top {

	#<
	#? Write a gromacs full topology file
	#; Requires: Ensemble (or molecule), options
	#; Returns: Undef if write failed, otherwise 1
	#>

	my $mols    = $_[0];
	my $outfile = $_[1];
	my $options = $_[2];

	my $fr;

	$options = uc $options                   if $options;
	
	silico_msg('c', "Writing gromacs TOP file\n");
	
	$outfile = get_ofilebase($mols) . ".top" if !$outfile;
	
	$fr = open_molfile($fr, $outfile, undef, 'write', $options);
	if (!defined $fr) {
		silico_msg('e', "Could not open $outfile");
		return undef;
	}
	
	my $FH = $fr->{FH};
	
	silico_msg('c', "Writing gromacs TOP file: $outfile\n");
	
	#
	# Header
	#
	
	print $FH ";\n";
	print $FH "; Gromacs topology file.\n";
	print $FH "; Generated by SILICO version $silico_version\n";
	print $FH ";\n\n";
	
	my $globaldata = $mols->[0]{GMX_TOPDATA};

	#
	# Defaults
	#
	
	if (my $d = $globaldata->{GMX_DEFAULTS}) {
	
		printf $FH "[ defaults ]\n";
		printf $FH ";%9s %10s %10s %10s %10s\n", qw(nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ);
		printf $FH "%10d %10d %10s %10.1f %10.1f\n\n", $d->{NB_FUNC} ,$d->{COMB_RULE}, $d->{GEN_PAIRS}, $d->{FUDGE_LJ}, $d->{FUDGE_QQ};
	
	} else {
	
		++$mols->[0]{WARN}{"Missing \"defaults\" block"};
	}

	#
	# Atomtypes
	#
	
	if (my $d = $globaldata->{GMX_ATOMTYPES}) {
	
		printf $FH "[ atomtypes ]\n";
		printf $FH ";%19s %10s %10s %10s %10s %12s %12s\n", qw(type at_num mass charge ptype sigma epsilon);
		
		foreach my $t (@$d) {

			printf $FH "%-20s %10d %10.6f %10.6f %10s %12.9f %12.9f\n", 
				$t->{GMX_TYPE}, $t->{ELEMENT_NUM}, $t->{MASS}, $t->{CHARGE}, $t->{PTYPE}, $t->{SIGMA}, $t->{EPSILON};      
		}
		
	} else {
	
		++$mols->[0]{WARN}{"Missing \"atomtypes\" block"};
	}

	#
        # Molecular topologies (.itp)
        #
	
	foreach my $mol (@$mols) {
		
		if (!defined $fr) {
			silico_msg('e', "Could not open $outfile");
			return undef;
		}

		write_gromacs_itp_mol($fr, $mol, $options);
	}
	
	#
	# System
	#
	
	printf $FH "\n";
	printf $FH "[ system ]\n";
	
	if (my $d = $globaldata->{GMX_SYSTEM}) {
	
		printf $FH "$globaldata->{GMX_SYSTEM}\n";
		
		
	} else {
	
		++$mols->[0]{WARN}{"Missing \"system\" block"};
	}
	
	#
	# Molecules
	#
	
	printf $FH "\n";
	printf $FH "[ molecules ]\n";
	
	if (my $d = $globaldata->{GMX_MOLECULES}) {
	
		foreach my $key (keys %$d) {
		
			printf $FH "%-20s %10d",$key, $globaldata->{GMX_MOLECULES}{$key};
		}
		
		
	} else {
	
		++$mols->[0]{WARN}{"Missing \"system\" block"};
	}
	
	mol_print_warnings($mols->[0]);
	
	close_molfile($fr);
	
	return 1;
}

sub write_gromacs_itp {

	#<
	#? Write a gromacs included topology file
	#. Only writes first molecule
	#; Requires: Ensemble (or molecule), options
	#; Returns: Undef if write failed, otherwise 1
	#>

	my $mols    = ens($_[0]);
	my $outfile = $_[1];
	my $options = $_[2];

	my $fr;

	$options = uc $options                   if $options;
	$outfile = get_ofilebase($mols) . ".itp" if !$outfile;

	$fr = open_molfile($fr, $outfile, undef, 'write', $options);
	if (!defined $fr) {
		silico_msg('e', "Could not open $outfile");
		return undef;
	}
	
	my $FH = $fr->{FH};

	silico_msg('c', "Writing gromacs ITP file: $outfile\n");

	foreach my $mol (@$mols) {

		print $FH ";\n";
		print $FH "; Gromacs included topology file.\n";
		print $FH "; Generated by SILICO version $silico_version\n";
		print $FH ";\n\n";
		
		write_gromacs_itp_mol($fr, $mol, $options);
	}
	
	return 1;
}

sub write_gromacs_itp_mol {

	#<
	#? Write a gromacs included topology file
	#. Assumes that bonds are present and that bond orders are correct
	#. Data is stored in atoms (several fields), GMX_BONDS, GMX_PAIRS, GMX_ANGLES, GMX_DIHEDRALS and GMX_EXCLUSIONS
	#; Requires: file_record, molecule, options
	#; Returns: Undef if write failed, otherwise 1
	#; Note that these subroutines use a new BONDS and ANGLES format currently called GMX_BONDS etc
	#>

	my $fr      = $_[0];
	my $mol     = $_[1];
	my $options = $_[2] || '';


	my $FH = $fr->{FH};
	++$fr->{MOL_WRITE_COUNT};

	molecule_find_rings($mol);

	my $atom;
	
	
	if ($mol->{COMMENTS}) {
	
		print $FH "\n";
		my $c = $mol->{COMMENTS};
		chomp $c;
		$c = "$c\n";
		my @f = split "\n", $c;
		
		foreach (@f) {
		
			s/^;//;
			s/^ //;
			print $FH "; $_\n";
		}
	}

	#
	# Molecule
	#
	print $FH "\n[ moleculetype ]\n";
	print $FH ";\n";
	print $FH "; Name";
	print $FH " " x (length(($mol->{NAME} || '')) + 4);
	print $FH "nrexcl\n";
	printf $FH "  %s        %d\n", ($mol->{NAME} || ''), ($mol->{NREXCL} || 3);
	print $FH "\n";

	print $FH "[ atoms ]\n";
	print $FH "; AtNum  AtType ResNum ResNme AtName  ChgGrp   Charge     Mass ; TotChg\n";

	#
	# Atoms
	#
	my $i = 0;
	my $totchg;
	my $warn;

	foreach my $atom (atoms($mol)) {

		# Skip water atoms
		next if $atom->{FG}{W};

		++$i;

		# Checks
		my $charge;

		my $group = $atom->{CHARGE_GROUP} || 999;

		if (defined $atom->{CHARGE}) {
			$charge = sprintf "%8.5f", $atom->{CHARGE};
			$totchg += $atom->{CHARGE};
			$totchg = 999 if $totchg > 999;

		} else {

			#$charge = "UNDEF";
			$charge = 0;
			++$warn->{ATOM_WITH_MISSING_CHARGE};
			$totchg = 999;
		}

		if (!defined $atom->{MASS}) {
			$atom->{MASS} = 999;
			++$warn->{ATOM_WITH_MISSING_MASS};
		}

		$charge = 0 if $options =~ "ZEROCHARGE";

		printf $FH "  %-6d %-6s %-6d %-7s %-6s %-6d %8.5f %8.4f ", $i,
		  ($atom->{TYPE}    || 'UNK'),
		  ($atom->{SUBID}   || 1),
		  ($atom->{SUBNAME} || 'UNK'),
		  ($atom->{NAME}    || 'UNK'),
		  $group,
		  $charge,
		  $atom->{MASS};
		  
		if ($atom->{TYPE_B}) {
			printf $FH " %-6s %8.5f %8.4f", ($atom->{TYPE_B} || ''), ($atom->{CHARGE_B} || ''), ($atom->{MASS_B} || ''),;
		}

		printf $FH " ; %6.3f\n", $totchg;
	}

	#
	# Pairs
	#
	
	print $FH "\n";
	print $FH "[ pairs ]\n";
	print $FH "; Atom1 Atom2 Func ; Type1 Type2    Name1 Name2\n";

	foreach my $pair (@{ $mol->{GMX_PAIRS} }) {

		#atom_printout($pair);

		my $a1 = $mol->{ATOMS}[ $pair->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $pair->{ATOMS}[1] ];

		printf $FH "  %-5d %-5d %-4s ", $pair->{ATOMS}[0] + 1, $pair->{ATOMS}[1] + 1, $pair->{GMX_FUNC};

		printf $FH ";  %-5s %-5s   %-5s %-5s\n", ($a1->{TYPE} || 'UNK'), ($a2->{TYPE} || 'UNK'), ($a1->{NAME} || 'UNK'), ($a2->{NAME} || 'UNK');
	}
	
	#
	# Bonds
	#
	
	print $FH "\n";
	print $FH "[ bonds ]\n";
	print $FH "; Atom1 Atom2 Func  Length  FConst    ; Type1 Type2    Name1 Name2 BOrd\n";

	foreach my $bond (@{ $mol->{GMX_BONDS} }) {
	
		#atom_printout($bond);

		my $a1 = $mol->{ATOMS}[ $bond->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $bond->{ATOMS}[1] ];

		my $bo = find_bondorder2($a1, $a2, $mol);

		printf $FH "  %-5d %-5d %-4s ", $bond->{ATOMS}[0] + 1, $bond->{ATOMS}[1] + 1, $bond->{GMX_FUNC};

		if ($bond->{FORMAT} eq 'null') {
			print $FH '                   ';
		} elsif ($bond->{FORMAT} eq 'full') {
			printf $FH "%8.4f %10.4g ", $bond->{GMX_FC}, $bond->{LENGTH};
			if (defined $bond->{LENGTH_B} && $bond->{LENGTH_B} ne '') {
				printf $FH "%8.4f %10.4g ", $bond->{GMX_FC_B}, $bond->{LENGTH_B};
			}

		} elsif ($bond->{FORMAT} eq 'code') {
			printf $FH "%-6s             ", $bond->{CODE};
		} else {
			die("Should not reach here. bond->{FORMANT} = '$bond->{FORMAT}'\n");
		}

		# Atom type and name comments
		printf $FH ";  %-5s %-5s   %-5s %-5s %-4d\n",
		  ($a1->{TYPE} || 'UNK'), ($a2->{TYPE} || 'UNK'),
		  ($a1->{NAME} || 'UNK'), ($a2->{NAME} || 'UNK'),
		  $bo;
	}

	#
	# Angles
	#
	
	print $FH "\n";
	print $FH "[ angles ]\n";
	print $FH "; Atom1 Atom2 Atom3 Func    Angle   FConst  ;   Res Subid Type1 Type2 Type3   Name1 Name2 Name3\n";
	
	foreach my $angle (@{ $mol->{GMX_ANGLES} }) {

		#atom_printout($angle);

		my $a1 = $mol->{ATOMS}[ $angle->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $angle->{ATOMS}[1] ];
		my $a3 = $mol->{ATOMS}[ $angle->{ATOMS}[2] ];

		printf $FH "  %-5d %-5d %-5d %-4d ", $angle->{ATOMS}[0] + 1, $angle->{ATOMS}[1] + 1, $angle->{ATOMS}[2] + 1, $angle->{GMX_FUNC};

		if ($angle->{FORMAT} eq 'null') {
			print $FH '                   ';
		} elsif ($angle->{FORMAT} eq 'full') {
			printf $FH " %8.2f %8.2f ", $angle->{ANGLE}, $angle->{GMX_FC};
			if (defined $angle->{ANGLE_B} && $angle->{ANGLE_B} ne '') {
				printf $FH " %8.2f %8.2f ", $angle->{ANGLE_B}, $angle->{GMX_FC_B};
			}
		} else {
			printf $FH " %-6s            ", $angle->{CODE};
		}

		printf $FH ";  %4s %-5s %-5s %-5s %-5s   %-5s %-5s %-5s\n",
		  ($a1->{SUBNAME} || 'UNK'), ($a1->{SUBID} || 'UNK'),
		  ($a1->{TYPE}    || 'UNK'), ($a2->{TYPE}  || 'UNK'), ($a3->{TYPE} || 'UNK'),
		  ($a1->{NAME}    || 'UNK'), ($a2->{NAME}  || 'UNK'), ($a3->{NAME} || 'UNK');
	}

	#
	# Dihedrals
	#
	
	print $FH "\n";
	print $FH "[ dihedrals ]\n";

	my $prev = '';
	foreach my $angle (@{ $mol->{GMX_DIHEDRALS} }) {

		#atom_printout($angle);

		my $a1 = $mol->{ATOMS}[ $angle->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $angle->{ATOMS}[1] ];
		my $a3 = $mol->{ATOMS}[ $angle->{ATOMS}[2] ];
		my $a4 = $mol->{ATOMS}[ $angle->{ATOMS}[3] ];

 		if ($angle->{GMX_FUNC} == 1 && $prev ne 'proper') {
			
			print $FH "; Proper dihedrals\n";
			print $FH "; Atom1 Atom2 Atom3 Atom4 Func  Phase     cp  Mult ;  Res Subid Type1 Type2 Type3 Type4  Name1 Name2 Name3 Name4\n";
		}
		 
		if ($angle->{GMX_FUNC} == 2 && $prev ne 'improper') {
			
			print $FH "; Improper dihedrals\n";
			print $FH "; Atom1 Atom2 Atom3 Atom4 Func Angle FConst        ;  Res Subid Type1 Type2 Type3 Type4  Name1 Name2 Name3 Name4\n";
		}
		
		$prev = 'proper'   if ($angle->{GMX_FUNC} == 1);
		$prev = 'improper' if ($angle->{GMX_FUNC} == 2);

		printf $FH "  %5d %5d %5d %5d %-4s ",
		  $angle->{ATOMS}[0] + 1,
		  $angle->{ATOMS}[1] + 1,
		  $angle->{ATOMS}[2] + 1,
		  $angle->{ATOMS}[3] + 1,
		  $angle->{GMX_FUNC};

		if ($angle->{FORMAT} eq 'null') {
			print $FH '                    ';
		} elsif ($angle->{FORMAT} eq 'full') {

			if ($angle->{GMX_FUNC} == 1 || $angle->{GMX_FUNC} == 9) {

				printf $FH "%6.2f %6.2f %5d ", $angle->{ANGLE}, $angle->{PHASE}, $angle->{MULT};

				if (defined $angle->{ANGLE_B} && $angle->{ANGLE_B} ne '') {
					printf $FH "%6.2f %6.2f %5d ", $angle->{ANGLE_B}, $angle->{PHASE_B}, $angle->{MULT_B};
				}

			} elsif ($angle->{GMX_FUNC} == 2) {

				printf $FH "%6.2f %6.2f      ", $angle->{ANGLE}, $angle->{GMX_FC};

			} elsif ($angle->{GMX_FUNC} == 4) {

				if ($angle->{ANGLE} ne '') {
					printf $FH "%6.2f %6.2f %6d      ", $angle->{ANGLE}, $angle->{GMX_FC}, $angle->{MULT};

					if (defined $angle->{ANGLE_B} && $angle->{ANGLE_B} ne '') {
						printf $FH "%6.2f %6.2f %6d      ", $angle->{ANGLE_B}, $angle->{GMX_FC_B}, $angle->{MULT_B};
					}
				}
			} elsif ($angle->{GMX_FUNC} == 3) {

				foreach (0 .. 5) {
					printf $FH "%6.2f ", $angle->{RB}[$_];
				}

			} else {
				silico_msg("d", "Angle function $angle->{GMX_FUNC} has not been implemented\n");
			}

		} else {
			printf $FH "%-6s              ", $angle->{CODE},;
		}

		printf $FH "; %4s %-5s %-5s %-5s %-5s %-5s  %-5s %-5s %-5s %-5s\n",
		  ($a1->{SUBNAME} || 'UNK'), ($a1->{SUBID} || 'UNK'),
		  ($a1->{TYPE}    || 'UNK'), ($a2->{TYPE}  || 'UNK'), ($a3->{TYPE} || 'UNK'), ($a4->{TYPE} || 'UNK'),
		  ($a1->{NAME}    || 'UNK'), ($a2->{NAME}  || 'UNK'), ($a3->{NAME} || 'UNK'), ($a4->{NAME} || 'UNK');
	}
	
	if ( $mol->{GMX_DIHEDRAL_RESTRAINTS} ) {
	
		print $FH "\n";
		print $FH "[ dihedral_restraints ]\n";
		print $FH ";";
		print $FH "; Atom1 Atom2 Atom3 Atom4 Func Angle FConst        ;  Res Subid Type1 Type2 Type3 Type4  Name1 Name2 Name3 Name4\n";
		print $FH "\n";
	
		foreach my $angle (@{ $mol->{GMX_DIHEDRAL_RESTRAINTS} }) {
		
			my $a1 = $mol->{ATOMS}[ $angle->{ATOMS}[0] ];
			my $a2 = $mol->{ATOMS}[ $angle->{ATOMS}[1] ];
			my $a3 = $mol->{ATOMS}[ $angle->{ATOMS}[2] ];
			my $a4 = $mol->{ATOMS}[ $angle->{ATOMS}[3] ];
		
			printf $FH "  %5d %5d %5d %5d %-4s %5d %5d %5d ",
				$angle->{ATOMS}[0] + 1,
				$angle->{ATOMS}[1] + 1,
				$angle->{ATOMS}[2] + 1,
				$angle->{ATOMS}[3] + 1,
				$angle->{GMX_FUNC},
				$angle->{ANGLE},
				$angle->{DELTA_ANGLE},
				$angle->{GMX_FC};
				
				
			printf $FH "; %4s %-5s %-5s %-5s %-5s %-5s  %-5s %-5s %-5s %-5s\n",
		  		($a1->{SUBNAME} || 'UNK'), ($a1->{SUBID} || 'UNK'),
		  		($a1->{TYPE}    || 'UNK'), ($a2->{TYPE}  || 'UNK'), ($a3->{TYPE} || 'UNK'), ($a4->{TYPE} || 'UNK'),
		  		($a1->{NAME}    || 'UNK'), ($a2->{NAME}  || 'UNK'), ($a3->{NAME} || 'UNK'), ($a4->{NAME} || 'UNK');
		}
		
	}

	print $FH "\n";
	print $FH "[ exclusions ]\n";
	print $FH "; Atom1 Atom2 Type\n";
	print $FH "\n";

	if (defined $mol->{GMX_EXCLUSIONS}[0]) {
		foreach my $ex (@{ $mol->{GMX_EXCLUSIONS} }) {
			my $a1 = $mol->{ATOMS}[ $ex->{ATOMS}[0] ];
			my $a2 = $mol->{ATOMS}[ $ex->{ATOMS}[1] ];
			printf $FH "  %-5d %-5d ", $ex->{ATOMS}[0] + 1, $ex->{ATOMS}[1] + 1;
			printf $FH "%6s ", $ex->{GMX_FUNC} || '';
			printf $FH ";  %-5s %-5s %-5s %-5s\n", ($a1->{TYPE} || 'UNK'), ($a2->{TYPE} || 'UNK'), ($a1->{NAME} || 'UNK'), ($a2->{NAME} || 'UNK');
		}
	}
	
	
	#
	# Settles
	#
	
	if (defined $mol->{GMX_SETTLES}) {
		print $FH "\n";
		print $FH "[ settles ]\n";
		print $FH ";i  funct   dOH  dHH\n";
		
		my $s = $mol->{GMX_SETTLES};
		
		print $FH "$s->{I} $s->{FUNC} $s->{D_OH} $s->{D_HH}\n";
	}
	
}


sub gmx_sort_angles {

	my $mol = $_[0];

	foreach (@{ $mol->{GMX_PAIRS} }) {
		@{ $_->{ATOMS} } = reverse(@$_->{ATOMS}) if $_->{ATOMS}[0] > $_->{ATOMS}[1];
	}
	@{ $mol->{GMX_PAIRS} } = sort {
		return 1  if $a->{ATOMS}[0] > $b->{ATOMS}[0];
		return -1 if $a->{ATOMS}[0] < $b->{ATOMS}[0];
		return $a->{ATOMS}[1] <=> $b->{ATOMS}[1]
	} @{ $mol->{GMX_PAIRS} };

	foreach (@{ $mol->{GMX_BONDS} }) {
		@{ $_->{ATOMS} } = reverse(@{ $_->{ATOMS} }) if $_->{ATOMS}[0] > $_->{ATOMS}[1];
	}
	@{ $mol->{GMX_BONDS} } = sort {
		return 1  if $a->{ATOMS}[0] > $b->{ATOMS}[0];
		return -1 if $a->{ATOMS}[0] < $b->{ATOMS}[0];
		return $a->{ATOMS}[1] <=> $b->{ATOMS}[1]
	} @{ $mol->{GMX_BONDS} };

	foreach (@{ $mol->{GMX_ANGLES} }) {
		@{ $_->{ATOMS} } = reverse(@{ $_->{ATOMS} }) if $_->{ATOMS}[0] > $_->{ATOMS}[2];
	}
	@{ $mol->{GMX_ANGLES} } = sort {
		return 1  if $a->{ATOMS}[0] > $b->{ATOMS}[0];
		return -1 if $a->{ATOMS}[0] < $b->{ATOMS}[0];
		return $a->{ATOMS}[2] <=> $b->{ATOMS}[2]
	} @{ $mol->{GMX_ANGLES} };

	@{ $mol->{GMX_DIHEDRALS} } = sort {
		return 1  if $a->{GMX_FUNC} > $b->{GMX_FUNC};
		return -1 if $a->{GMX_FUNC} < $b->{GMX_FUNC};
		return 1  if $a->{ATOMS}[0] > $b->{ATOMS}[0];
		return -1 if $a->{ATOMS}[0] < $b->{ATOMS}[0];
		return $a->{ATOMS}[1] <=> $b->{ATOMS}[1]
	} @{ $mol->{GMX_DIHEDRALS} };
}

sub gmx_merge_coords_into_top {

	#<
	#? Merge coordinates into ITP file
	#; Requires: Topology file, coordinate file
	#>

	my $imol = deep_copy($_[0]);    # Itp (topology) molecule
	my $cmol = $_[1];               # Coordinate molecule

	# Generate and print atom keys from coordinate file
	my $coord_hash = atom_keys($cmol, 'coordinate');

	# Print out atom keys from topology file
	atom_keys($imol, 'topology');
	print "\n";

	my $death = 0;
	foreach my $iatom (atoms($imol)) {

		# Find matching atom from coordinate molecule
		my $catom = $coord_hash->{ $iatom->{KEY} };

		if ($catom) {

			# Copy coordinates to topology file
			$iatom->{X}     = $catom->{X};
			$iatom->{Y}     = $catom->{Y};
			$iatom->{Z}     = $catom->{Z};
			$iatom->{MATCH} = $catom;

			# copy bonds from coordinate file
			$iatom->{BORDERS} = $catom->{BORDERS};
			$iatom->{CONNECT} = $catom->{CONNECT};

			# Save the matching atom in $iatom->{MATCH}
			++$catom->{MATCH_COUNT};

			if ($catom->{MATCH_COUNT} > 1) {
				silico_msg('w', "Multiple matches ($catom->{MATCH_COUNT}) for atom: $iatom->{KEY}\n");
				++$death;
			}
		} else {

			silico_msg('w', "No matching coordinate atom found for topology atom: $iatom->{KEY}\n");
			++$death;
		}
	}

	delete $imol->{HAS_NO_COORDS};

	if ($death) {

		silico_msg('c', "Writing modified coordinate file to file err.mol2");
		write_mol_any($cmol, "err.mol2", 'mol2');
		silico_msg('d', "Problems encountered\n");
	}

	return $imol;
}

sub atom_keys {

	#>
	#? Create hash of atom name and subid
	#; Requires: molecule, identifier string
	#>

	my $mol    = $_[0];
	my $string = $_[1] || '';

	my $hash;

	foreach my $atom (atoms($mol)) {

		my $n = $atom->{NAME};
		my $s = $atom->{SUBID};
		$n =~ s/ //g;
		$s =~ s/ //g;
		my $key = "$n\_$s";

		$hash->{$key} = $atom;
		$atom->{KEY} = $key;
	}

	print heading("Atom keys from '$string' file\n");
	my @f = keys(%$hash);
	my $i = 0;
	foreach (sort @f) {
		++$i;
		print "$_ ";
		if ($i == 20) {
			$i = 0;
			print "\n";
		}
	}
	print "\n" if $i != 0;

	return $hash;
}

sub gromacs_check_itp {

	#<
	#? Check a gromacs included topology data and make sure that it is consistent with the bonded molecular structure
	#; Requires: molecule
	#; Returns: Undef if write failed, otherwise 1
	#>

	my $mol = $_[0];

	#
	# Check atoms
	#
	my $i = 0;
	my $totchg;

	foreach my $atom (atoms($mol)) {

		# Skip water atoms
		next if $atom->{FG}{W};

		++$i;

		# Checks
		my $charge;

		my $group = $atom->{CHARGE_GROUP} || 999;

		if (defined $atom->{CHARGE}) {
			$charge = sprintf "%8.5f", $atom->{CHARGE};
			$totchg += $atom->{CHARGE};
			if (abs($totchg) > 999) {
				$mol->{WARN}{"MOLECULE HAS A CHARGE PROBLEM"} = 1;
			}
		} else {
			++$mol->{WARN}{"ATOMS WITH UNDEFINED CHARGES"};
			$charge = "UNDEF";
			$totchg = 999;
		}
	}

	#
	# Check GMX BONDS, ANGLES and DIHEDRALS are consistent with molecular structure
	#

	foreach my $angle (@{ $mol->{GMX_BONDS} }, @{ $mol->{GMX_ANGLES} }, @{ $mol->{GMX_DIHEDRALS} }) {
		check_itp_is_bonded($angle, $mol);
	}

	#
	# Check the molecular structure is consistent with GMX BONDS, ANGLES and DIHEDRALS
	#

	foreach my $atom (atoms($mol)) {
		check_atom_has_itp_data($atom, $mol);
	}

	foreach (keys(%{ $mol->{WARN} })) {
		silico_msg('w', "Error '$_' found $mol->{WARN}{$_} times\n") if $mol->{WARN}{$_};
	}
}

sub check_atom_has_itp_data {

	my $atom = $_[0];
	my $mol  = $_[1];

	my @connected = connected($atom, $mol);

	my $match;
	foreach my $con (@connected) {

		$match = 0;
		my $matches;
		foreach my $bond (@{ $mol->{GMX_BONDS} }) {

			if (       ($bond->{ATOMS}[0] == $atom->{NUM} - 1 && $bond->{ATOMS}[1] == $con->{NUM} - 1)
				|| ($bond->{ATOMS}[0] == $con->{NUM} - 1 && $bond->{ATOMS}[1] == $atom->{NUM} - 1))
			{

				++$match;
				push @$matches, $bond;
			}
		}

		if ($match == 0) {
			++$mol->{WARN}{"Bond missing from GMX_BONDS"};

			print heading("Bond $atom->{NUM} $con->{NUM} missing from GMX_BONDS\n");
			atom_printout($atom);
			print "\n----\n\n";
			atom_printout($con);

			print heading("GMX_BONDS\n");

			foreach (sort { $a->{ATOMS}[0] <=> $b->{ATOMS}[0] } @{ $mol->{GMX_BONDS} }) {
				print "@{$_->{ATOMS}}, ";
			}
			print "\n";
			carp();
			die;
		}

		if ($match > 1) {

			++$mol->{WARN}{"Bond found multiple times in GMX_BONDS"};

			print heading("Bond $atom->{NUM} $con->{NUM} found multiple times ($match) in GMX_BONDS\n");
			atom_printout($atom);
			print "\n----\n\n";
			atom_printout($con);

			print heading("GMX_BONDS:\n");
			print "@{$mol->{GMX_BONDS}}\n";

			print heading("Matching atoms:\n");
			foreach (@$matches) {
				atom_printout($_);
				print "\n----\n\n";
			}
		}
	}
}

sub check_itp_is_bonded {

	#<
	#? Check that atoms in a GMX_BOND, GMX_ANGLE, or GMX_DIHEDRAL ($angle) are actually bonded
	#>

	my $angle = $_[0];
	my $mol   = $_[1];

	my $a1;
	my $a2;
	my @atomnums = @{ $angle->{ATOMS} };
	my $count    = 0;

	# Propers

	$a1 = shift @atomnums;

	while ($a2 = shift @atomnums) {

		my $bo = find_bondorder($a1, $a2, $mol);

		if ($bo < 1) {
			++$count;
		}

		$a1 = $a2;
	}

	# Chiral impropers
	# Note: improper atoms may  bonded I-J-K-L in the case of
	# ring impropers or I(-J)(-K)-L or L(-I)(-J)-K where the improper
	# is a constraint on a tetrahedral or trigonal centre.

	# Check to see if this is a sane chiral improper I(-J)(-K)-L if we failed the Propers test
	if ($count && $angle->{TYPE} eq 'improper') {

		$count = 0;
		foreach (1 .. 3) {
			my $bo = find_bondorder($angle->{ATOMS}[0], $angle->{ATOMS}[$_], $mol);
			++$count if $bo < 1;
		}
	}

	# Check to see if this is a sane chiral improper L(-I)(-J)-K if we failed the Propers test
	if ($count && $angle->{TYPE} eq 'improper') {

		$count = 0;
		foreach (0 .. 2) {
			my $bo = find_bondorder($angle->{ATOMS}[3], $angle->{ATOMS}[$_], $mol);
			++$count if $bo < 1;
		}
	}

	if ($angle->{WARN}{"Not bonded"}) {

		print heading("GMX_ANGLE contains unbonded atoms\n");
		atom_printout($angle);

		print heading("Atoms referenced by GMX_ANGLE\n");
		foreach (@{ $angle->{ATOMS} }) {
			print "---\n";
			atom_printout($mol->{ATOMS}[$_]);
		}
		die;

	}

	$angle->{WARN}{"Not bonded"}                        = $count;
	$mol->{WARN}{"ITP $angle->{TYPE} atoms not bonded"} = $count;
}

sub gromacs_regenerate_pairs {

	#<
	#? Regenerate 1-4 atom pairs list
	#; Requires: molecule
	#>

	my $mol = $_[0];
	my $angle;
	my $atom1;
	my $atom2;

	my $dihedrals = mol_find_all_dihedrals($mol);
	my $oldpairs  = $mol->{GMX_PAIRS};
	undef $mol->{GMX_PAIRS};

	my $pairs14;
      DIH: foreach my $dihedral (@$dihedrals) {

		my $pair;

		# Skip if the dihedral loops back to the same atom (eg. in 4 membered rings?)
		next DIH if ($dihedral->[0] == $dihedral->[3]);

		# Skip if the dihedral loops back to a directly bonded atom (thus a 1-2 relationship)
		foreach my $bond (@{ $mol->{BONDS} }) {
			next DIH if ($dihedral->[0] == $bond->[0] && $dihedral->[3] == $bond->[1]);
			next DIH if ($dihedral->[3] == $bond->[0] && $dihedral->[0] == $bond->[1]);
		}

		# Skip if the dihedral loops back to an atom in an angle (thus a 1-3 relationship)
		foreach my $angle (@{ $mol->{ANGLES} }) {
			next DIH if ($dihedral->[0] == $angle->[0] && $dihedral->[3] == $angle->[2]);
			next DIH if ($dihedral->[3] == $angle->[0] && $dihedral->[0] == $angle->[2]);
		}

		# Skip if the dihedral has already been listed. This has been a problem in certain
		# aromatic rings, e.g. benzene, where, for example, the pairs arising from the
		# dihedral C1-C2-C3-C4 and the dihedral C4-C5-C6-C1 were both getting a listing (thus
		# two instances of C1-C4).
		foreach my $pair14 (@$pairs14) {
			next DIH if ($dihedral->[0] == $pair14->[0] && $dihedral->[3] == $pair14->[1]);
			next DIH if ($dihedral->[3] == $pair14->[0] && $dihedral->[0] == $pair14->[1]);
		}

		# Make sure the lower atom number comes first
		@$pair = sort { $a <=> $b } ($dihedral->[3], $dihedral->[0]);

		# Add to the list of pairs
		push @$pairs14, $pair;
	}

	# Sort the pair list into a sensible order
	@$pairs14 = sort {
		if    ($a->[0] > $b->[0]) { return 1; }
		elsif ($a->[0] < $b->[0]) { return -1; }
		else {
			if    ($a->[1] > $b->[1]) { return 1; }
			elsif ($a->[1] < $b->[1]) { return -1; }
			else                      { return 0; }
		}
	} @$pairs14;

	foreach my $pair14 (@$pairs14) {

		$atom1 = $mol->{ATOMS}[ $pair14->[0] ];
		$atom2 = $mol->{ATOMS}[ $pair14->[1] ];

		# Commented out -not sure of thie function???? DKC
		# Do not generate exclusions for HC-X-X-X (??)
		#next DIH if ($atom1->{TYPE} eq 'HC' || $atom2->{TYPE} eq 'HC');

		itp_make_pair($mol, $pair14->[0] + 1, $pair14->[1] + 1, 1);
	}

	# Print things out
	if (0) {
		foreach ($oldpairs, $mol->{GMX_PAIRS}) {
			my $count = 0;
			foreach my $pair14 (sort { $a->{ATOMS}[0] <=> $b->{ATOMS}[0] } @$_) {
				print "$pair14->{ATOMS}[0] $pair14->{ATOMS}[1], ";
				++$count;
			}
			print "\n$count pairs\n\n";
		}
	}
}

sub dump_itp_data {

	#<
	#? Write out debugging file with all itp bonds, pairs, angles etc
	#; Requires: molecule
	#; Writes out sdf file for each bond, angle etc to a director
	#>

	my $mol = deep_copy($_[0]);

	my $mols;

	my $filebase = get_filebase($mol->{SOURCE_FILE_NAME});

	# Delete existing bonds
	foreach my $atom (atoms($mol)) {
		delete $atom->{CONNECT};
		delete $atom->{BORDERS};
		delete $atom->{MATCH};
	}
	delete $mol->{RINGS};

	# Bonds
	my $i = 0;
	foreach my $bond (@{ $mol->{GMX_BONDS} }) {

		$i = sprintf "%03d", ++$i;
		my $a1 = $mol->{ATOMS}[ $bond->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $bond->{ATOMS}[1] ];

		my $nm;

		push @{ $nm->{ATOMS} }, deep_copy($a1), deep_copy($a2);
		bond_create($nm, 0, 1);

		$nm->{NUMATOMS} = 2;
		$nm->{NUMBONDS} = 1;
		$nm->{NAME}     = "Bond_$i\_$a1->{SUBID}_$a2->{SUBID}_$a1->{NAME}_$a2->{NAME}";
		push @$mols, $nm;
	}

	# Pairs
	$i = 0;
	foreach my $pair (@{ $mol->{GMX_PAIRS} }) {

		$i = sprintf "%03d", ++$i;
		my $a1 = $mol->{ATOMS}[ $pair->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $pair->{ATOMS}[1] ];
		my $nm;

		push @{ $nm->{ATOMS} }, deep_copy($a1), deep_copy($a2);
		bond_create($nm, 0, 1);

		$nm->{NUMATOMS} = 2;
		$nm->{NUMBONDS} = 1;
		$nm->{NAME}     = "Pair_$i\_$a1->{SUBID}_$a2->{SUBID}-$a1->{NAME}_$a2->{NAME}";
		push @$mols, $nm;
	}

	# Angles
	$i = 0;
	foreach my $angle (@{ $mol->{GMX_ANGLES} }) {

		$i = sprintf "%03d", ++$i;
		my $a1 = $mol->{ATOMS}[ $angle->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $angle->{ATOMS}[1] ];
		my $a3 = $mol->{ATOMS}[ $angle->{ATOMS}[2] ];

		my $nm;

		push @{ $nm->{ATOMS} }, deep_copy($a1), deep_copy($a2), deep_copy($a3);
		bond_create($nm, 0, 1);
		bond_create($nm, 1, 2);

		$nm->{NUMATOMS} = 3;
		$nm->{NUMBONDS} = 2;
		my $s = $angle->{CODE} || "$angle->{ANGLE}_$angle->{GMX_FC}";
		$nm->{NAME} = "Angle_$i\_$a1->{SUBID}_$a2->{SUBID}_$a3->{SUBID}-$a1->{NAME}_$a2->{NAME}_$a3->{NAME}_$s";
		push @$mols, $nm;
	}

	# Dihedrals
	$i = 0;
	foreach my $angle (@{ $mol->{GMX_DIHEDRALS} }) {

		$i = sprintf "%03d", ++$i;
		my $a1 = $mol->{ATOMS}[ $angle->{ATOMS}[0] ];
		my $a2 = $mol->{ATOMS}[ $angle->{ATOMS}[1] ];
		my $a3 = $mol->{ATOMS}[ $angle->{ATOMS}[2] ];
		my $a4 = $mol->{ATOMS}[ $angle->{ATOMS}[3] ];

		my $nm;

		push @{ $nm->{ATOMS} }, deep_copy($a1), deep_copy($a2), deep_copy($a3), deep_copy($a4);
		bond_create($nm, 0, 1);
		bond_create($nm, 1, 2);
		bond_create($nm, 2, 3);

		$nm->{NUMATOMS} = 4;
		$nm->{NUMBONDS} = 3;
		if ($angle->{GMX_FUNC} == 1) {
			$nm->{NAME} = "Dihedral";
		} elsif ($angle->{GMX_FUNC} == 2) {
			$nm->{NAME} = "Improper";
		} elsif ($angle->{GMX_FUNC} == 3) {
			$nm->{NAME} = "RB";
		} else {
			$nm->{NAME} = "UNK";
		}

		my $s;
		if ($angle->{GMX_FUNC} == 1) {
			if ($angle->{CODE}) {
				$s = $angle->{CODE};
			} else {
				$s = "$angle->{ANGLE}_$angle->{PHASE}_$angle->_{MULT}";
			}
		} elsif ($angle->{GMX_FUNC} == 2) {
			if ($angle->{CODE}) {
				$s = $angle->{CODE};
			} else {
				$s = "$angle->{ANGLE}_$angle->{PHASE}_$angle->_{MULT}";
			}
		} elsif ($angle->{GMX_FUNC} == 3) {
			$s = join "_", @{ $angle->{RB} };
		} else {
			silico_msg('w', "Angle type $angle->{GMX_FUNC} not implemented\n");
		}

		$nm->{NAME} .= "_$i\_$a1->{SUBID}_$a2->{SUBID}_$a3->{SUBID}_$a4->{SUBID}-$a1->{NAME}_$a2->{NAME}_$a3->{NAME}_$a4->{NAME}_$s";
		push @$mols, $nm;
	}

	# Exclusions
	if (defined $mol->{GMX_EXCLUSIONS}[0]) {
		foreach my $ex (@{ $mol->{GMX_EXCLUSIONS} }) {

			$i = sprintf "%03d", ++$i;
			my $a1 = $mol->{ATOMS}[ $ex->{ATOMS}[0] ];
			my $a2 = $mol->{ATOMS}[ $ex->{ATOMS}[1] ];

			my $nm;

			push @{ $nm->{ATOMS} }, deep_copy($a1), deep_copy($a2);
			bond_create($nm, 0, 1);

			$nm->{NUMATOMS} = 2;
			$nm->{NUMBONDS} = 1;
			$nm->{NAME}     = "Excl_$i\_$a1->{SUBID}_$a2->{SUBID}-$a1->{NAME}_$a2->{NAME}";
			push @$mols, $nm;
		}
	}

	my $dir = "$filebase\_itp-data.dir";
	system("rm -rf $dir $dir\.tar.bz2");
	system("mkdir $dir");

	my $ens;
	my $count = 0;
	my $files;
	$i = 0;
	foreach (@$mols) {
		++$count;
		$count = sprintf "%05d", $count;
		$ens->[0] = $_;
		my $name = "$_->{NAME}";
		$name =~ s/ /_/g;
		write_sdf($ens, "$dir/$name.sdf");
		push @$files, "$name.sdf";
	}

	open(OUTFILE, ">$dir/0load.pml") || die;
	foreach (@$files) {
		print OUTFILE "load $_\n";
	}
	close OUTFILE;
	system "tar -cf $dir.tar $dir; bzip2 $dir.tar";

}


###############
# Subroutines #
###############



sub merge_itp_coord {

	#<
	#? Merge coordinates into ITP molecule
	#; Requires: ITP molecule, coordinate molecule
	#; Returns: molecule
	#>
	
	#
	# Note there seems to be a simliar routine in silico_gromacs already
	#

	my $mol  = $_[0];    # ITP molecule
	my $cmol = $_[1];    # Coordinate molecule

	die ('Error');
	# Note: subroutine is using a global hash. This needs to be fixed.
	my $hash;

	my $warn = 0;

	# Delete old MATCH and MATCH_COUNT values
	foreach my $atom ($cmol) {
		delete $atom->{MATCH_COUNT};
		delete $atom->{MATCH};
	}

	my $death = 0;
	foreach my $atom (atoms($mol)) {
		my $n = $atom->{NAME};
		my $s = $atom->{SUBID};
		$n =~ s/ //g;
		$s =~ s/ //g;
		my $key   = "$n\_$s";
		my $tatom = $hash->{$key};    # Template atom

		if ($tatom) {

			$atom->{X}     = $tatom->{X};
			$atom->{Y}     = $tatom->{Y};
			$atom->{Z}     = $tatom->{Z};
			$atom->{MATCH} = $tatom;
			++$tatom->{MATCH_COUNT};
			if ($tatom->{MATCH_COUNT} > 1) {
				silico_msg('w', "Multiple matches ($tatom->{MATCH_COUNT}) for atom: $key\n");
				++$death;
			}
		} else {
			silico_msg('w', "No matching itp atom found for atom: $key\n");
			++$death;
		}
	}
	
	silico_msg('d', "Problems encountered\n") if $death;

	foreach my $atom (atoms($mol)) {
		$atom->{SMILES_POS} = $atom->{MATCH}{SMILES_POS};
	}

	delete $mol->{HAS_NO_COORDS};

	return $mol;
}

sub read_coord_mol {

	#<
	#? Read in coordinate molcule, check connectivity, fix ATB problems
	#; Requires: filename
	#; Returns: molecule
	#>

	my $cmol = $_[0];

	die('Error');
	# Subroutine is using a global hash. This needs to be fixed
	my $hash;

	my $ens = read_mol_any($cmol);
	die "Could not read file $cmol\n" if !defined $ens->[0];

	$cmol = $ens->[0];
	molecule_check_and_fix_connectivity($cmol);

	my $flag;
	foreach my $atom (atoms($cmol)) {

		# Fix ATB pdb file residue numbering which starts from 0
		$flag = 1 if ($atom->{SUBID} == 0);
		++$atom->{SUBID} if $flag;

		my $n = $atom->{NAME};
		my $s = $atom->{SUBID};
		$n =~ s/ //g;
		$s =~ s/ //g;
		my $key = "$n\_$s";
		$hash->{$key} = $atom;
	}

	print heading("Atom keys (<atom name>_<substructure id>)");
	my @f = keys(%$hash);
	my $i = 0;
	foreach (@f) {
		++$i;
		print "$_ ";
		if ($i == 20) {
			$i = 0;
			print "\n";
		}
	}
	print "\n" if $i != 0;
	print "\n";

	return ($cmol);
}

sub sort_itp_by_smiles {

	#<
	#? Sort itp file by canonical smiles
	#  ascending order by atom bumber
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];

	mol_smiles_string($mol) if !defined $mol->{SMILES};

	# Smiles pos calculated based on template molecule
	my $ctable = molecule_reorder($mol, "SMILES_POS");    #

	my @list = qw(GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS GMX_EXCLUSIONS);

	foreach my $l (@list) {
		foreach my $thing (@{ $mol->{$l} }) {
			foreach (@{ $thing->{ATOMS} }) {
				$_ = $ctable->[$_];
			}
		}
	}

	sort_itp_by_atom_number($mol);
}

sub sort_itp_by_atom_number {

	#<
	#? Sort GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS GMX_EXCLUSIONS in
	#  ascending order by atom bumber
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];

	my @list = qw(GMX_BONDS GMX_PAIRS GMX_ANGLES GMX_DIHEDRALS GMX_EXCLUSIONS);

	# Sort GMX_* into numerical order by atom
	foreach my $l (@list) {
		@{ $mol->{$l} } = sort atom_number_sort_function @{ $mol->{$l} };
	}
}

sub atom_number_sort_function {

	#<
	#? Sort function for sort_itp_by_atom_number()
	#>

	foreach my $i (0 .. 3) {

		return 0  if !defined $a->{ATOMS}[$i];
		return 1  if $a->{ATOMS}[$i] > $b->{ATOMS}[$i];
		return -1 if $a->{ATOMS}[$i] < $b->{ATOMS}[$i];
	}

	return 0;

}

sub gmx_replace_peptide_dihedrals {

	#<
	#? Replace the dihedral parameters of peptide backbone atoms
	#; Requires: molecule
	#; Returns: number of dihedrals replaced
	#>

	my $mol = $_[0];

	my $count = 0;

	my $dihedral_hash;
	%$dihedral_hash = qw (N_CH1_C_N gd_45  CH1_N_C_CH1 gd_14  C_CH1_N_C gd_44);

      DIHEDRAL: foreach $a (@{ $mol->{DIHEDRALS} }) {

		next if $a->{GMX_FUNC} == 2;
		my $a1;
		my $a2 = $mol->{ATOMS}[ $a->{ATOMS}[1] ];
		my $a3 = $mol->{ATOMS}[ $a->{ATOMS}[2] ];
		my $a4;

		foreach $a1 (connected($a2, $mol)) {
			next if $a1 == $a3;
			foreach $a4 (connected($a3, $mol)) {
				next if $a4 == $a2;

				my @types;
				foreach ($a1, $a2, $a3, $a4) {
					push @types, $_->{TYPE};
				}

				my @types2 = reverse @types;
				my $key1   = join "_", @types;
				my $key2   = join "_", @types2;

				#print "key: $key1\n";
				#print "key: $key2\n";

				my $key;
				foreach $key ($key1, $key2) {
					if ($dihedral_hash->{$key}) {

						$a->{CODE}   = $dihedral_hash->{$key};
						$a->{FORMAT} = 'code';
						@{ $a->{ATOMS} } = ($a1->{NUM} - 1, $a2->{NUM} - 1, $a3->{NUM} - 1, $a4->{NUM} - 1);
						++$count;
						next DIHEDRAL;
					}
				}
			}
		}
	}
	return $count;
}

return 1;
