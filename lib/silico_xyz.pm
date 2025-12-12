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
#! silico_xyz.pm
#? Silico routines specific to Jay Ponder's TINKER and SLEUTH software
#. $Revision: 1.1.2.3 $
#>

use strict;

package Silico;

sub read_xyz {

	#<
	#? XYZ read routine.
	#; Requires: filename, (optional) start, max molecules, options string
	#; Returns: ensemble or undef if failed
	#>

	my $infile        = $_[0];
	my $start         = $_[1] || get_flag('ss', 's') || 0;
	my $max_molecules = $_[2] || get_flag('ms', 's');
	my $options       = uc($_[3] || '');
	my $mols;

	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr = open_molfile($fr, $infile, undef, undef, $options);
	if (!defined $fr) {
		silico_msg('e', "Can not open file $infile for reading\n");
		return undef;
	}

	my $FH = $fr->{FH};

	if ($options !~ /\bQUIET\b/) {

		silico_msg('c', "Reading xyz file: $infile\n");
	}

	my $i = -1;
	my $molcount = 0;
	while (my $mol = read_xyz_single($fr, $options)) {
		++$i;
		next if $start and ($i < $start);
		++$molcount;
		push @$mols, $mol;
		last if (defined $max_molecules) && $molcount == $max_molecules;
	}

	close($FH);

	return $mols;
}

sub read_xyz_single {

	#<
	#? Read in a single record from an xyz file
	#. Tinker format has no comment line (comments are on first line after number of atoms)
	#; Requires: file record, options
	#; Returns: molecule or undef
	#>

	my $fr      = $_[0];
	my $options = $_[1] || '';

	#croak();

	my $FH = $fr->{FH};

	my $mol;

	# Read Header
	my $line = <$FH>;

	return if !defined $line;

	++$fr->{MOL_READ_COUNT};

	chomp $line;

	$line =~ s/^\s*//;

	my @f = split(" ", $line);
	$mol->{NUMATOMS} = $f[0];

	$mol->{NAME} = "Mol_" . ($fr->{MOL_READ_COUNT} + 1);

	$line = <$FH>;
	chomp $line;
	$mol->{SDF_DATA}{COMMENTS} = $mol->{COMMENTS} = $line;

	for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

		$line = <$FH>;

		if (!defined $line) {
			silico_msg('e', "Premature end of file!\n");
			return undef;
		}

		chomp $line;
		my @f = split(" ", $line);

		my $atom;

		$atom->{NAME} = $f[0];    #Atom name
		$atom->{X}    = $f[1];    #x
		$atom->{Y}    = $f[2];    #y
		$atom->{Z}    = $f[3];    #z
		$atom->{NUM}  = $i + 1;

		$atom->{SUBNAME} = 'SUB';
		$atom->{SUBID}   = 1;

		push @{ $mol->{ATOMS} }, $atom;
	}

	# Work out the atom elements
	foreach my $atom (atoms($mol)) {

		atom_guess_element($atom, $mol);
	}

	mol_set_source_data($mol, $fr, 'xyz');

	# No bondorder information in a .xyz file
	$mol->{HAS_BAD_BONDORDERS} = 1;

	#molecule_printout($mol);
	#die if ++$Silico::count == 10;

	return $mol;
}

sub write_xyz {
	
	#<
	#? Write an xyz file
	#; Requires: Ensemble or molecule, filename, $options.
	#; Returns: Zero if file open failed otherwise returns 1.
	#>

	my $molecules = ens($_[0]);
	my $outfile = $_[1];
	my $options = $_[2] || '';


	$options = uc $options if $options;
	$outfile = get_ofilebase($molecules).".xyz" if !$outfile;
	
	my $fr = open_molfile($fr, $outfile, undef, 'write', $options);
	if (!defined $fr) {
		silico_msg('e', "Could not open $outfile");
		return undef;
	}
	my $FH = $fr->{FH};
	
	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Writing xyz file: $outfile");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}
	
	# Write out molecules
	my $molcount = 0;
	foreach my $mol (@$molecules) {
	
		die if !defined $mol;
	
		++$molcount;
		if (!defined $FH) {
			silico_msg ('e',"Attempt to write to unopened filehandle. Perhaps you don't have write permission?\n");
			return undef;
		}
		
		write_xyz_molecule($fr, $mol, $options);
	}

	close_molfile($fr);
	
	return 1;
}


sub write_xyz_molecule {

	#<
	#? Write xyz molecule
	#; Requires: ensemble (or molecule), filename, forcefield file (optional) to
	#  type atoms.
	#; Returns: undef if file open failed otherwise returns 1.
	#>

	my $fr  = $_[0] || silico_msg('d', "File record not defined\n");
	my $mol = $_[1] || silico_msg('d', "Molecule not defined\n");
	my $options = $_[2] || '';

	if (!defined $mol->{NUMATOMS}) {

		silico_msg('e', "\$mol->NUMATOMS is not defined!\n");
		return undef;
	}
	
	my $FH = $fr->{FH};
	++$fr->{MOL_WRITE_COUNT};

	# Print header
	print $FH "$mol->{NUMATOMS}\n";
	printf $FH "%s\n", ($mol->{NAME} || 'Mol');

	# Print atoms
	for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

		my $atom = $mol->{ATOMS}[$i];

		printf $FH "%6d  %4s %12.6f %12.6f %12.6f\n", $i + 1, $atom->{ELEMENT}, $atom->{X}, $atom->{Y}, $atom->{Z};
	}

	return 1;
}

return 1;
