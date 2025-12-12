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
#! silico_mmod.pm
#? Silico routines specific to macromodel format files
#. $Revision: 1.15-73-g3bc62ee $
#>

use strict;

package Silico;

##################################################################
#
#	Macromodel format read routines
#
##################################################################

sub read_mmod {

	#<
	#? Macromodel read routine.
	#. Can read Macromodel 'short' text format.
	#; Requires: filename, start molecule (currently ignored), max molecules, option string
	#; Returns: ensemble or undef if failed
	#>

	my $infile        = $_[0];
	my $start         = $_[1] || get_flag('ss', 's');    # Ignored!!!!
	my $max_molecules = $_[2] || get_flag('ms', 's');
	my $options       = uc($_[3] || '');

	my $ensemble;
	my %hash;

	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr = open_molfile(undef, $infile, undef, undef, $options);
	if (!defined $fr) {
		silico_msg('e', "Can not open file $infile for reading\n");
		return undef;
	}
	my $FH = $fr->{FH};

	silico_msg('c', "Reading mmod file: $_[0]\n") if $options !~ /\bQUIET\b/;

	# Macromodel atom type to element conversion
	%Silico::silico_mmod_atypes = qw(

	  1	C	2	C	3	C	4	C	5	C
	  6	C	7	C	8	C	9	C	10	C
	  11	C	12	C	13	C	14	C	15	O
	  16	O	17	O	18	O	19	O	20	O
	  21	O	22	O	23	O	24	N	25	N
	  26	N	27	N	28	N	29	N	30	N
	  31	N	32	N	33	N	34	N	35	N
	  36	N	37	N	38	N	39	N	40	N
	  41	H	42	H	43	H	44	H	45	H
	  46	H	47	H	48	H	49	S	50	S
	  51	S	52	S	53	P	54	B	55	B
	  56	F	57	Cl	58	Br	59	I	60	Si
	  61	Du	62	Z0	63	Lp	65	Li	66	Na
	  67	K	68	Rb	69	Cs	70	Ca	71	Ba
	);

	my $molcount = 0;
	while (1) {

		my $mol;
		++$molcount;

		#
		# Needs to be modified to use standart max_molecules test subroutine
		#

		# Exit if we have reached the maximum number of molecules
		if ($max_molecules && ($molcount > $max_molecules)) {
			silico_msg('c', "Read maximum number of molecules $max_molecules. Finishing\n") if $options !~ /\bQUIET\b/;
			last;
		}

		undef %hash;

		# Read Header
		my $line = <$FH>;

		last if !defined $line;

		# Make sure that there is a space after the '=' sign
		# In some cases this is not present
		$line =~ s/=/= /;

		my @f = split(' ', $line);

		my $compressed;
		if ($f[0] >= 0) {
			$compressed = 0;
			$mol->{NUMATOMS} = $f[0];
		} else {
			$compressed = 1;
			$mol->{NUMATOMS} = -$f[0];
		}

		$mol->{NAME} = "";

		# Get molecule name and energy (value following E =)
		for (my $i = 1 ; $i <= $#f - 2 ; ++$i) {

			if ($f[$i] eq "E" && $f[ $i + 1 ] eq "=") {
				$mol->{ENERGY}       = $f[ $i + 2 ];
				$mol->{ENERGY_UNITS} = "kJ/mol";
				$mol->{ENERGY_TYPE}  = "Macromodel";
				last;
			}

			$mol->{NAME} = join " ", $mol->{NAME}, $f[$i];
		}

		# Create all the atoms first because we will need to
		# create bonds to atoms that have not been read in yet.

		for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {
			my $at;
			$at->{NUM} = $i + 1;       #Atom number
			$mol->{ATOMS}[$i] = $at;
		}

		# Read macromodel uncompressed and compressed formats
		if ($compressed == 0) {

			# Uncompressed format
			# -------------------

			for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

				my $atom = $mol->{ATOMS}[$i];

				my $line = <$FH>;

				defined $line || silico_msg('d', "Premature end of file $infile!\n");

				# Pad out line with spacs
				chomp $line;
				$line .=
'                                                                                                                              ';

#                0     1        2  3        4  5        6  7        8  9        10 11       12  13            14            15           16     17  18  19   20         21          22     23
#                type  con      bo con      bo con      bo con      bo con      bo con      bo  x             y             z            subid  res chn colr occ        temp        subn   name
				my @f = $line =~
/^(....)(......).(.)(......).(.)(......).(.)(......).(.)(......).(.)(......).(.).(...........).(...........).(...........)(......)(.)(.)(....)(.........)(.........).(.....)(....)/;

				$atom->{NUM}       = $i + 1;
				$atom->{MMOD_TYPE} = $f[0];     #Atom type
				$atom->{SUBNAME}   = $f[20];    #Residue name
				$atom->{X}         = $f[13];    #x
				$atom->{Y}         = $f[14];    #y
				$atom->{Z}         = $f[15];    #z

				$atom->{SUBID}   = ($f[16]);            #Substructure number
				$atom->{CHAIN}   = ($f[18]);            #chain
				$atom->{OCC}     = ($f[20] || 0);
				$atom->{TEMP}    = ($f[21] || 0);
				$atom->{SUBNAME} = ($f[22] || 'UNK');
				$atom->{NAME}    = $f[23];

				$atom->{SUBID} = 1 if $atom->{SUBID} eq ' ';

				# Remove spaces
				$atom->{MMOD_TYPE} =~ s/ //g;
				$atom->{SUBID}     =~ s/ //g;
				$atom->{SUBNAME}   =~ s/ //g;

				$atom->{ELEMENT}     = $Silico::silico_mmod_atypes{ $atom->{MMOD_TYPE} };
				$atom->{ELEMENT}     = "Du" if (!defined $atom->{ELEMENT});
				$atom->{ELEMENT_NUM} = element_symbol2number($atom->{ELEMENT});

				# Create bonds
				for (my $j = 1 ; $j <= 12 ; $j = $j + 2) {

					last if ($f[$j] == 0);

					# The internal connectivity matrix is indexed
					# from zero. The connectivities in the input
					# file are indexed from 1

					my $atnum = $f[$j] - 1;

					# Next if these atoms already bonded
					next if ($hash{"$atnum $i"});

					++$hash{"$i $atnum"};

					# Translate bond order here!!
					# Not done at present
					my $bo = $f[ $j + 1 ];

					# Create bond.  Remember that the connectivity
					# matrix is indexed from zero
					bond_create($mol, $i, $atnum, $bo);
				}
			}
		}

		if ($compressed) {

			# Compressed format
			# -------------------

			$mol->{NUMBONDS} = $ensemble->[0]{NUMBONDS};

			for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

				my $atom = $mol->{ATOMS}[$i];

				my $line = <$FH>;

				defined $line || silico_msg('d', "Premature end of file $infile!\n");

				@f = split(' ', $line);

				$atom->{NUM} = $i + 1;
				$atom->{X}   = $f[1];    #x
				$atom->{Y}   = $f[2];    #y
				$atom->{Z}   = $f[3];    #z

				# Call the identical atom in the first
				# molecule of the ensemble the 'parent'
				# Then get values from there

				my $parent = $ensemble->[0]{ATOMS}[$i];

				$atom->{NAME}        = $parent->{NAME};
				$atom->{MMOD_TYPE}   = $parent->{MMOD_TYPE};
				$atom->{ELEMENT}     = $parent->{ELEMENT};
				$atom->{ELEMENT_NUM} = $parent->{ELEMENT_NUM};
				$atom->{SUBID}       = $parent->{SUBID};
				$atom->{SUBNAME}     = $parent->{SUBNAME};
				@{ $atom->{CONNECT} } = @{ $parent->{CONNECT} };
				@{ $atom->{BORDERS} } = @{ $parent->{BORDERS} };
			}
		}

		mol_rename_atoms($mol);    # Macromodel produces stupid atom names
		mol_set_source_data($mol, $fr);
		push @$ensemble, $mol;
	}

	close_molfile($fr);

	return $ensemble;
}

##################################################################
#
#	Macromodel format write routines
#
##################################################################

sub write_mmod {

	#<
	#? Macromodel write routine.
	#; Requires: ensemble (or molecule), filename.
	#; Returns: undef if file open failed otherwise returns 1.
	#>

	my $ensemble = ens($_[0]);
	my $outfile  = $_[1];

	my $warn_atomname = 0;
	my $warn_subname  = 0;
	my $warn_kekule   = 0;

	silico_msg('d', "Argument 0 is not an ensemble of molecule\n") if ref($_[0] ne 'ARRAY');

	$outfile = get_ofilebase($ensemble) . ".dat" if !$outfile;

	# Macromodel colours of elements.  Many are missing
	my %colourtable = qw(C 2 N 4 O 16);

	my $success = open(MMODOUT, ">$outfile");
	if (!$success) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing Macromodel file: $outfile\n");

	foreach my $mol (@$ensemble) {

		# Fix bonds if we have bad bondorders (ie we came from a pdb file)
		molecule_check_and_fix_connectivity($mol);

		make_atom_names_unique($mol) if get_flag('unique-atom-names', 'l');

		# Convert bonds to kekule if we have any aromatic bondorders
		ATOM: foreach my $atom (atoms($mol)) {
			foreach my $bo (@{ $atom->{BORDERS} }) {
				if ($bo > 3) {
					++$warn_kekule;
					convert_aromatic_bonds_kekule($mol);
					last ATOM;
				}
			}
		}

		if (!defined $mol->{NUMATOMS}) {
			silico_msg('e', "\$mol->NUMATOMS is not defined!\n", "Aborting write.\n");
			return undef;
		}

		# Generate macromodel atom types
		mol_type_mmod($mol) if !defined $mol->{ATOMS}[0]{MMOD_TYPE};

		# Print header
		printf MMODOUT " %5d", $mol->{NUMATOMS};
		print MMODOUT "  $mol->{NAME}"      if defined $mol->{NAME};
		print MMODOUT " E = $mol->{ENERGY}" if defined $mol->{ENERGY};
		print MMODOUT "\n";

		for (my $i = 0 ; $i < $mol->{NUMATOMS} ; ++$i) {

			my $atom = $mol->{ATOMS}[$i];

			printf MMODOUT " %3s", $atom->{MMOD_TYPE};

			# Bonds and orders
			for (my $j = 0 ; $j <= 5 ; ++$j) {

				if (defined $atom->{CONNECT}[$j]) {

					# Translate bond order here (NOT DONE YET!)

					my $bo = $atom->{BORDERS}[$j];

					printf MMODOUT " %5d %1d", $atom->{CONNECT}[$j] + 1, $bo;

				} else {

					print MMODOUT "     0 0";
				}
			}

			# Atom colours
			my $colour = $colourtable{ $atom->{ELEMENT} } || 2;

			# Truncate residue name if it is too long
			my $subname = $atom->{SUBNAME} || 'UNK';
			if (length $subname > 3) {
				$subname = substr($subname, 0, 3);
				++$warn_subname;
			}

			# Truncate atom name if it is too long
			my $atomname = $atom->{NAME} || 'X';
			if (length $atomname > 4) {
				$atomname = substr($atomname, 0, 4);
				++$warn_atomname;
			}

			printf MMODOUT " %11.6f %11.6f %11.6f %5dX%1s %3d %8.5f %8.5f %3s %4s",
			  $atom->{X},
			  $atom->{Y},
			  $atom->{Z},
			  $atom->{SUBID},
			  ($atom->{CHAIN} || ' '),
			  $colour,
			  ($atom->{OCC}  || 0),
			  ($atom->{TEMP} || 0),
			  $subname,
			  $atomname;

			# Missing some optional fields in here!!

			print MMODOUT "\n";
		}
	}

	silico_msg('n', "Molecules ($warn_kekule) contained aromatic bondorders. Converting to Kekule.\n") if $warn_kekule;
	silico_msg('n', "Substructure names have been truncated to three characters\n")                    if $warn_subname;
	silico_msg('n', "Atom names have been truncated to four characters\n")                             if $warn_atomname;

	close(MMODOUT);

	return 1;
}

sub mol_type_mmod {

	#<
	#? Macromodel atom type routine
	#. Assumes that structures are three-dimensional (ie not planar 2D representations)
	#  and have complete bond types.  Atom geometry is used to identify aromatic rings
	#  and planar nitrogen atoms.  Does not assume that hydrogens are present.
	#  Makes primary amines N.4 and makes O-C=O into carboxylate groups if explicit
	#  hydrogens are not present.
	#; Requires: molecule.
	#; Returns: Macromodel atom type.
	#>

	my $mol = $_[0];

	my $error = 0;

	# Tag planar atoms
	molecule_find_planar_atoms($mol);

	# Returns: 0 if no hydrogens, 1 if polar hydrogens, 2 nonpolar hydrogens
	# (ie connected to carbon) or -1 if hydrogens are present but bonds are not.
	my $has_hydrogens = mol_has_h($mol);

	#print "Molecule has hydrogens = $has_hydrogens\n";

	foreach my $atom (atoms($mol)) {

		my $numaromatic     = 0;
		my $numbonds        = 0;
		my $numsingle       = 0;
		my $numdouble       = 0;
		my $numtriple       = 0;
		my $con_planar_atom = 0;
		my $con_n           = 0;
		my $con_o           = 0;
		my $con_h           = 0;
		my $con_c           = 0;
		my $con_s           = 0;

		my $el = $atom->{ELEMENT};

		if (!$el) {
			silico_msg('w', "Atom number $atom->{NUM} has no defined element!\n", "Setting this atom's Macromodel type to 61.\n");
			$atom->{MMOD_TYPE} = '61';
			next;
		}

		# Get types of bonds and connected atoms
		my $i = 0;
		foreach my $connum (@{ $atom->{CONNECT} }) {

			my $bo = $atom->{BORDERS}[$i];

			# Skip if zero order or no bonds defined
			next if (!defined $bo);
			next if ($bo == 0);

			my $con = $mol->{ATOMS}[$connum];

			if ($con->{ELEMENT} eq 'O') {
				++$con_o;
			}
			if ($con->{ELEMENT} eq 'H') {
				++$con_h;
			}
			if ($con->{ELEMENT} eq 'C') {
				++$con_c;
			}
			if ($con->{ELEMENT} eq 'S') {
				++$con_s;
			}
			if ($con->{ELEMENT} eq 'N') {
				++$con_n;
			}
			if ($con->{PLANAR_ATOM}) {
				++$con_planar_atom;
			}

			++$numsingle   if ($bo == 1);
			++$numdouble   if ($bo == 2);
			++$numtriple   if ($bo == 3);
			++$numaromatic if ($bo == 4);

			++$numbonds;

			++$i;
		}

		if ($el eq 'H') {

			if ($con_c || $con_s) {
				$atom->{MMOD_TYPE} = '41';
				next;
			}
			if ($con_o) {
				$atom->{MMOD_TYPE} = '42';
				next;
			}
			if ($con_n) {
				$atom->{MMOD_TYPE} = '42';
				next;
			}

			warn_atom_type_error($mol, $atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{MMOD_TYPE} = 48;
		}

		# Carbon
		if ($el eq 'C') {

			# United atom
			if ($has_hydrogens == 0 || $has_hydrogens == 1) {

				# UA CH sp3
				if ($numsingle == 3 && $numbonds == 3) {
					$atom->{MMOD_TYPE} = '4';
					next;
				}

				# UA CH2 sp3
				if ($numsingle == 2 && $numbonds == 2) {
					$atom->{MMOD_TYPE} = 5;
					next;
				}

				# UA CH3 sp3
				if ($numsingle == 1 && $numbonds == 1) {
					$atom->{MMOD_TYPE} = 6;
					next;
				}

				# UA CH sp2
				if ($numdouble == 1 && $numbonds == 2) {
					$atom->{MMOD_TYPE} = 7;
					next;
				}

				# UA CH2 sp2
				if ($numdouble == 1 && $numbonds == 1) {
					$atom->{MMOD_TYPE} = 8;
					next;
				}

				# UA CH2 sp2
				if ($numtriple == 1) {
					$atom->{MMOD_TYPE} = 9;
					next;
				}
			}

			if ($numdouble == 0 && $numtriple == 0) {
				$atom->{MMOD_TYPE} = '3';
				next;
			}
			if ($numdouble >= 1) {
				$atom->{MMOD_TYPE} = 2;
				next;
			}
			if ($numtriple == 1) {
				$atom->{MMOD_TYPE} = 1;
				next;
			}

			warn_atom_type_error($mol, $atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{MMOD_TYPE} = 61;
			$error = 1;
			next;
		}

		# Nitrogen
		if ($el eq 'N') {

			# N+ sp3
			if ($numsingle == 4) {
				$atom->{MMOD_TYPE} = 32;
				next;
			}

			# N sp
			if ($numtriple == 1) {
				$atom->{MMOD_TYPE} = 24;

			}

			if ($numdouble == 0 && $numtriple == 0) {

				#  N sp2
				if ($atom->{PLANAR_ATOM}) {
					$atom->{MMOD_TYPE} = 25;
					next;
				}

				#  N sp2 attached to aromatic
				if ($con_planar_atom) {
					$atom->{MMOD_TYPE} = 25;
					next;
				}

				# N+ sp3 but no hydrogens
				# Promote nitrogens in molecules without hydrogens to N.4
				if (!$has_hydrogens) {
					$atom->{MMOD_TYPE} = 32;
					next;
				}

				# N sp3
				$atom->{MMOD_TYPE} = 26;
				next;
			}

			# N sp2
			if ($numdouble >= 1) {
				$atom->{MMOD_TYPE} = 25;
				next;
			}

			# N sp
			if ($numtriple == 1) {
				$atom->{MMOD_TYPE} = 24;
				next;
			}

			warn_atom_type_error($mol, $atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{MMOD_TYPE} = 61;
			$error = 1;
			next;
		}

		# Oxygen
		if ($el eq 'O') {

			# O sp2
			if ($numdouble == 1) {
				$atom->{MMOD_TYPE} = 15;
				next;
			}

			# O sp3
			if ($numsingle == 2) {
				$atom->{MMOD_TYPE} = 16;
				next;
			}

			# Alkoxide or carboxylate
			if ($numsingle == 1 && $numdouble == 0) {
				$atom->{MMOD_TYPE} = 18;
				next;
			}

			warn_atom_type_error($mol, $atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
			$atom->{MMOD_TYPE} = 61;
			$error = 1;
			next;

		}

		# Sulfur
		if ($el eq 'S') {

			# Thiolate
			if ($numdouble == 0 && $numsingle == 1) {
				$atom->{MMOD_TYPE} = 51;
				next;
			}

			# All other S
			$atom->{MMOD_TYPE} = 49;
			next;
		}

		# Phosphorous
		if ($el eq 'P') {
			$atom->{MMOD_TYPE} = 53;
			next;
		}

		# Fluorine
		if ($el eq 'F') {
			if ($numbonds == 0) {
				$atom->{MMOD_TYPE} = 56;
			} else {

				# Fluoride
				$atom->{MMOD_TYPE} = 104;
			}
			next;
		}

		# Chlorine
		if ($el eq 'Cl') {
			if ($numbonds == 0) {
				$atom->{MMOD_TYPE} = 57;
			} else {

				# Chloride
				$atom->{MMOD_TYPE} = 102;
			}
			next;
		}

		# Bromine
		if ($el eq 'Br') {
			if ($numbonds == 0) {
				$atom->{MMOD_TYPE} = 105;
			} else {
				$atom->{MMOD_TYPE} = 58;
			}
			next;
		}

		# Iodine
		if ($el eq 'I') {
			if ($numbonds == 0) {
				$atom->{MMOD_TYPE} = 106;
			} else {
				$atom->{MMOD_TYPE} = 59;
			}
			next;
		}

		# Silicon
		if ($el eq 'Si') {
			$atom->{MMOD_TYPE} = 60;
			next;
		}

		# Lithium
		if ($el eq 'Li') {
			$atom->{MMOD_TYPE} = 65;
			next;
		}

		# Sodium +
		if ($el eq 'Na') {
			$atom->{MMOD_TYPE} = 66;
			next;
		}

		# Potassium +
		if ($el eq 'K') {
			$atom->{MMOD_TYPE} = 67;
			next;
		}

		# Rubidium +
		if ($el eq 'Rb') {
			$atom->{MMOD_TYPE} = 68;
			next;
		}

		# Caesium
		if ($el eq 'Cs') {
			$atom->{MMOD_TYPE} = 69;
			next;
		}

		# Calcium 2+
		if ($el eq 'Ca') {
			$atom->{MMOD_TYPE} = 70;
			next;
		}

		# Barium 2+
		if ($el eq 'Ba') {
			$atom->{MMOD_TYPE} = 71;
			next;
		}

		# Magnesium 2+
		if ($el eq 'Mg') {
			$atom->{MMOD_TYPE} = 72;
			next;
		}

		warn_atom_type_error($mol, $atom, $numsingle, $numdouble, $numtriple, $numaromatic, $numbonds);
	}
}

##################################################################
#
#	Maestro format read routines
#
##################################################################

sub read_maestro {

	#<
	#? Schroedinger maestro read routine.
	#; Requires: filename, start molecule (currently ignored), max molecules, option string
	#; Returns: ensemble or undef if failed
	#. Options:
	#  Returns an empty array when there are no more structures.
	#, Will respond to a -me (maximum energy) and max molecule (-ms) flags
	#>

	my $infile        = $_[0];
	my $start         = $_[1] || get_flag('ss', 's');
	my $max_molecules = $_[2] || get_flag('ms', 's');
	my $options       = uc($_[3] || '');

	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	my $energy_max = get_flag('me', 's');

	my $fr = open_molfile(undef, $infile, undef, undef, $options);
	return undef if !$fr;

	if ($options !~ /\bQUIET\b/ && $fr->{MOL_READ_COUNT} < 1) {
		silico_msg('c', "Reading maestro file: $_[0]");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}

	my $i        = -1;
	my $molcount = 0;
	my $mols;
	while (my $mol = read_maestro_single($fr)) {
		++$i;

		# Needs to be replaced with a routine similar to "test_energy_molcount"
		next if $start and ($i < $start);
		++$molcount;
		push @$mols, $mol;
	}

	silico_msg('w', "Renamed $fr->{MAE_WARN_RENAME} atoms with missing or empty atom names\n") if $fr->{MAE_WARN_RENAME};

	return $mols;
}

sub read_maestro_single {

	#<
	#? Schroedinger maestro read routine.
	#; Requires: filehandle, options
	#; Returns: molecule
	#. Options:
	#>

	my $fr      = $_[0];
	my $options = uc($_[1] || '');

	my $infile = $fr->{FILENAME};
	my $mol;    # New molecule
	my $molcount = 0;

	++$fr->{MOL_READ_COUNT};

      LOOP: while (1) {

		my $tok = mae_get_tok($fr);

		# Exit at legitimate end of file
		if (!defined $tok) {
			last LOOP;
		}

		# Read header for first molecule
		if ($fr->{MOL_READ_COUNT} == 0) {
			if ($tok eq '{') {
				$tok = mae_get_tok($fr);
				if ($tok !~ 'm2io_version') {
					silico_msg('e', "File $infile does not appear to be a maestro file!\n");
					close_file();
					return undef;
				}
				mae_ignore_block($fr);    # Ignore remainder of block
			}
		}

		#
		# Read ct data blocks
		#

		next if $tok !~ /_ct$/;

		if ($tok eq 'f_m_ct') {

			++$molcount;
			if ($molcount > 1) {
				mae_put_tok($fr, $tok);
				last;
			}

			# Pass an empty molecule to mae_parse_ct_data
			@{ $mol->{ATOMS} } = ();
			$mol->{NUMATOMS} = 0;
			$mol = mae_parse_ct_data($fr, $mol);

			# Retain pointer to previous molecule so that we can use
			# it as a template if required
			$fr->{MAE_PREV_MOLECULE} = $mol;

		} elsif ($tok eq 'p_m_ct') {

			++$molcount;
			if ($molcount > 1) {
				mae_put_tok($fr, $tok);
				last;
			}

			# Make a copy of the previous molecule
			# to use as a template for for the
			# partial molecule

			silico_msg('d', "Malformed file.  No f_m_ct record\n") if !defined $fr->{MAE_PREV_MOLECULE};
			my $pmol = deep_copy($fr->{MAE_PREV_MOLECULE});
			$mol = mae_parse_ct_data($fr, $pmol);
		}
	}

	# End of file or read error
	if (!defined $mol) {
		return;
	}

	my $sd = $mol->{SDF_DATA};

	if (defined $sd->{PDB_CRYST1_a}) {

		my @f = qw(a b c alpha beta gamma Space_Group);

		$mol->{CELL}[0]     = $sd->{PDB_CRYST1_a};
		$mol->{CELL}[1]     = $sd->{PDB_CRYST1_b};
		$mol->{CELL}[2]     = $sd->{PDB_CRYST1_c};
		$mol->{CELL_ALPHA}  = $sd->{PDB_CRYST1_alpha};
		$mol->{CELL_BETA}   = $sd->{PDB_CRYST1_beta};
		$mol->{CELL_GAMMA}  = $sd->{PDB_CRYST1_gamma};
		$mol->{SPACE_GROUP} = space_group_pdb_to_silico($sd->{PDB_CRYST1_Space_Group});
	}

	mol_set_source_data($mol, $fr, 'mae');

	# Rename atoms if any atoms have empty or undefined names
	my $count = 0;
	foreach my $atom (atoms($mol)) {

		next if $atom->{NAME} && $atom->{NAME} =~ /\w/;
		++$count;
		$atom->{NAME} = $atom->{ELEMENT} . $count;
		++$fr->{WARN}{MAE_RENAME};
	}

	# Test that we have not exeeded maximum molecule energy or count
	# set by flags
	return undef if test_energy_molcount($fr, $mol);

	return $mol;
}

sub test_energy_molcount {

	#<
	#? Test to see if we have exceeded molecule maximum energy or maximum molecule count
	#  set by flags
	#; Requires: filehandle, molecule, max_energy (optional or from flag), max molcount (optional or from flag)
	#; Returns; true if we have exeeded counts/energy
	#>

	my $fr           = $_[0];
	my $mol          = $_[1];
	my $max_energy   = $_[2];
	my $max_molcount = $_[3];

	if (!$fr->{MES_FLAG}) {

		$fr->{ME} = $max_energy || get_sflag('me');
		if ($fr->{ME}) {
			if (!defined $mol->{ENERGY}) {
				molecule_printout($mol);
				silico_msg('d', "\nFlag -me has ben set, but energy value not found in first structure\n");
			}
			$fr->{ENERGY_INIT} = $mol->{ENERGY};
		}
		$fr->{MS}       = $max_molcount || get_sflag('ms');
		$fr->{MES_FLAG} = 1;
	}

	if (defined $fr->{MS} && $fr->{MOL_READ_COUNT} >= $fr->{MS}) {
		silico_msg('c', "\nRead maximum number of molecules $fr->{MS}. Finishing\n");
		return 1;
	}

	if (defined $fr->{ME}) {
		my $en = $mol->{ENERGY} - $fr->{ENERGY_INIT};
		if ($en > $fr->{ME}) {
			silico_msg('c',
				"\nEnergy $en exeeeded maximum energy $fr->{ME} $mol->{ENERGY_UNITS} after $fr->{MOL_READ_COUNT} molecules. Finishing\n");
			return 1;
		}
	}

	return 0;

}

sub mae_get_tok {

	#<
	#? Read mae file, process tokens and return one at a time
	#; Requires: filehandle
	#; Returns: token (string)
	#>

	my $fr = $_[0];

	my $lines;
	my $FH = $fr->{FH};

	# Read in more of the input file if there are no more tokens
	# in the data array
	if ($#{ $fr->{MAE_INFILE_TOKENS} } < 0) {

		while (<$FH>) {

			# Ignore comments
			next if /^\s*#/;

			$lines .= $_;

			if (/_m_ct/) {

				mae_split_data($fr, $lines);
				$lines = '';
				last;
			}
		}

		# Process any remaining data
		mae_split_data($fr, $lines) if defined $lines;
	}

	my $tok = shift @{ $fr->{MAE_INFILE_TOKENS} };

	return $tok;
}

sub mae_put_tok {

	my $fr  = $_[0];
	my $tok = $_[1];

	unshift @{ $fr->{MAE_INFILE_TOKENS} }, $tok;
}

sub mae_split_data {

	#<
	#? Split data into tokens
	#; Requires: lines from input file (string)
	#; Returns: array
	#>

	my $fr    = $_[0];
	my $lines = $_[1];

	# Match on quoted strings OR delimited by run of nonwhitespace characters
	my @data = $lines =~ m/(".*?"|\S+)/gs;

	# Remove double quotes
	foreach (@data) {

		s/^"(.*)"$/$1/;
	}

	foreach (@data) {

		push @{ $fr->{MAE_INFILE_TOKENS} }, $_;
	}
}

sub mae_parse_ct_data {

	#<
	#? Parse a maestro connection table data block (a molecule)
	#; Requires: filehandle, molecule (can be empty);
	#; Returns: molecule
	#>

	my $fr  = $_[0];
	my $mol = $_[1];

	my $hash;
	my $hash_type;

	mae_get_tok($fr);    # Discard bracket

	my ($headers, $types, $progs) = mae_parse_block_header($fr);

	my $i = -1;
	foreach my $head (@$headers) {

		++$i;
		my $tok = mae_get_tok($fr);

		#silico_msg('g', "ct head: $head\n", "tok2: (ct data) $tok\n");

		if (!defined $tok) {
			silico_msg('e', "Premature end of file!\n");
			return undef;
		}

		$hash->{$head}      = $tok;
		$hash_type->{$head} = $types->[$i];
	}

	# Maestro seems to have several 'name' entries.  Different names can be used under different
	# circumstances
	$mol->{NAME} = $hash->{title} || $hash->{entry_name} || $hash->{NAME} || 'Mol';

	# Put molecule data into SDF data fields and get an energy value
	foreach my $key (keys %$hash) {

		my $val = $hash->{$key};
		$val                        = '' if !defined $val;
		$mol->{SDF_DATA}{$key}      = $val;
		$mol->{SDF_DATA_TYPE}{$key} = $hash_type->{$key};

		#print "key '$key'\n";

		#
		# Units need to be fixed. Macromodel produces kJ/mol other Schrodinger progs produce kCal/mol
		#

		# Obtain potential energy value and type
		if ($key =~ /^potential_energy/i) {
			$mol->{ENERGY}       = $val;
			$mol->{ENERGY_TYPE}  = $key;
			$mol->{ENERGY_UNITS} = 'kcal/mol or kJ/mol';
		}

		# Obtain glide energy
		if ($key =~ /^glide_energy/i) {
			$mol->{ENERGY}       = $val;
			$mol->{ENERGY_TYPE}  = $key;
			$mol->{ENERGY_UNITS} = 'kcal/mol or kJ/mol';
		}

		# MINTA free energy
		if ($key =~ /^MINTA_Free_Energy-/i) {
			$mol->{ENERGY}       = $val;
			$mol->{ENERGY_TYPE}  = $key;
			$mol->{ENERGY_UNITS} = 'kcal/mol';
		}
		if ($key =~ /^MINTA_Free_Energy_Error/i) {
			$mol->{ENERGY_ERROR} = $val;
			$mol->{ENERGY_UNITS} = 'kcal/mol or kJ/mol';
		}
	}

	if (defined $mol->{ENERGY}) {
		$mol->{SDF_DATA}{ENERGY}       = $mol->{ENERGY};
		$mol->{SDF_DATA}{ENERGY_TYPE}  = $mol->{ENERGY_TYPE};
		$mol->{SDF_DATA}{ENERGY_UNITS} = $mol->{ENERGY_UNITS};
	}

	# Read additional blocks
	my $bracketlevel = 0;
	while (1) {

		my $tok = mae_get_tok($fr);

		if (!defined $tok) {

			#silico_msg('g', "End of file\n");
			last;
		}

		#silico_msg('g', "tok2a: (ct data) $tok\n");

		if ($tok eq '{') {
			++$bracketlevel;
			next;
		}
		if ($tok eq '}') {
			--$bracketlevel;
			next if $bracketlevel > 0;
			last;
		}

		if ($tok =~ /_atom\[(.*)\]/) {
			mae_parse_atoms($fr, $mol, $1);
			next;
		}
		if ($tok =~ /_bond\[(.*)\]/) {
			mae_parse_bonds($fr, $mol, $1);
			next;
		}
		mae_ignore_block($fr);
	}

	return $mol;
}

sub mae_parse_block_header {

	#<
	#? Parse the header information in a maestro data block
	#; Requires: filehandle
	#; Returns: array of headers, array of header types (r: real, i: integer, etc),
	#  array of header source programs
	#>

	my $fr = $_[0];

	my ($headers, $types, $progs);

	# Read headers
	while (1) {

		my $tok = mae_get_tok($fr);

		#silico_msg('g', "tok3: (block header) $tok\n");

		if (!defined $tok) {
			silico_msg('e', "Premature end of file!\n");
			return undef;
		}

		last if $tok eq ":::";

		# Remove field type and owner from header field
		my ($type, $prog) = $tok =~ m/^(.+?)_(.+?)_/;
		$tok =~ s/^(.+?)_(.+?)_//;

		# Replace '\_' with '_'
		$tok =~ s/\\_/_/g;

		push @$headers, $tok;
		push @$types,   $type;
		push @$progs,   $tok;
	}

	return ($headers, $types, $progs);
}

sub mae_ignore_block {

	#<
	#? Read a maestro data block but ignore the contents
	#; Requires: filehandle
	#; Returns: nothing
	#>

	my $fr = $_[0];

	my $bracketlevel = 0;

	while (1) {

		my $tok = mae_get_tok($fr);

		last if !defined $tok;

		#silico_msg('g', "tok4: (ignore block) $tok\n");

		++$bracketlevel if $tok eq '{';
		--$bracketlevel if $tok eq '}';

		last if $bracketlevel <= 0;
	}
}

sub mae_parse_atoms {

	#<
	#? Read a maestro atom block
	#. Can handle full atom records (f_m_ct) and partial atom records (p_m_ct)
	#; Requires: filehandle, molecule, number of atoms
	#; Returns: nothing
	#>

	my $fr  = $_[0];
	my $mol = $_[1];
	my $num = $_[2];    # Total number of atoms in molecule (there may be fewer atoms in block)

	my $hash;

	# Discard bracket
	mae_get_tok($fr);

	my ($headers, $types, $progs) = mae_parse_block_header($fr);

	# Add index column not declared in fields(!)
	unshift @$headers, '_idx';

	# Loop over atoms
      LOOP: while (1) {

		undef $hash;

		# Loop over atom fields
		foreach my $head (@$headers) {

			my $tok = mae_get_tok($fr);

			if (!defined $tok) {

				silico_msg('e', "Premature end to atom data!\n");
				return undef;
			}

			#silico_msg('g', ("$head : '$tok'  "));

			last LOOP if $tok eq ":::";

			$hash->{$head} = $tok;
		}
		
		
		#print "------\n";
		#foreach my $key (sort keys(%$hash)) {
		#
		#	print "$key $hash->{$key}\n";
		#}

		#silico_msg('g',"\n");

		my $atom;

		my $idx = $hash->{_idx};

		if (defined $mol->{ATOMS}[$idx-1]) {
		
			# p_m_ct record

			$atom = $mol->{ATOMS}[$idx-1];
			
			$atom->{X} = $hash->{x_coord};    # x
			$atom->{Y} = $hash->{y_coord};    # y
			$atom->{Z} = $hash->{z_coord};    # z

		} else {

			# f_m_ct record

			$atom->{NUM}           = $idx;                                            # Atom number
			$atom->{NAME}          = $hash->{pdb_atom_name} || $hash->{atom_name};    # Atom name
			$atom->{SUBNAME}       = $hash->{pdb_residue_name} || 'UNK';              # Residue name
			$atom->{SUBID}         = $hash->{residue_number};                         # Residue number
			$atom->{X}             = $hash->{x_coord};                                # x
			$atom->{Y}             = $hash->{y_coord};                                # y
			$atom->{Z}             = $hash->{z_coord};                                # z
			$atom->{CHAIN}         = $hash->{chain_name} || ' ';                      # Chain
			$atom->{ALT}           = $hash->{insertion_code} || '';                   # Alternate atom indicator?
			$atom->{ELEMENT_NUM}   = $hash->{atomic_number} || 0;                     # Atom atomic number
			$atom->{ELEMENT}       = element_number2symbol($atom->{ELEMENT_NUM});     # Element symbol
			$atom->{MMOD_TYPE}     = $hash->{mmod_type} || 0;                         # Atom mmod type
			$atom->{CHARGE}        = $hash->{charge1} || 0;                           # Atom charge
			$atom->{FORMAL_CHARGE} = $hash->{formal_charge} || 0;                     # Atom formal charge
			$atom->{MMOD_COLOUR}   = $hash->{color} || 0;                             # Macromodel atom colour code

			if ($atom->{SUBNAME} =~ ' ') {
				$atom->{SUBNAME} =~ s/ //g;
				++$mol->{WARN}{RENAME_SUBNAME};
			}
			
			$atom->{SUBID} = 1 if !defined $atom->{SUBID};

			push @{ $mol->{ATOMS} }, $atom;

			++$mol->{NUMATOMS};
		}
	}

	# Discard close bracket
	mae_get_tok($fr);
}

sub mae_parse_bonds {

	#<
	#? Read a maestro bond data block
	#; Requires: filehandle, molecule, number of bonds
	#; Returns: nothing
	#>

	my $fr  = $_[0];
	my $mol = $_[1];
	my $num = $_[2];    # Number of bonds

	my $hash;

	# Discard bracket
	mae_get_tok($fr);

	my ($headers, $types, $progs) = mae_parse_block_header($fr);

	# Add index column not declared in fields(!)
	# Named _idx so that it does not clash with i_rdk_index header
	unshift @$headers, '_idx';

	# Loop over bonds
      LOOP: while (1) {

		undef $hash;

		foreach my $head (@$headers) {
			my $tok = mae_get_tok($fr);

			#silico_msg('g', "tok6  (bonds): $tok\n");

			if (!defined $tok) {

				silico_msg('e', "Premature end to bond data!\n");
				return undef;
			}

			last LOOP if ($tok eq ':::');

			$hash->{$head} = $tok;
		}
		
		#foreach my $key (sort keys(%$hash)) {
		#
		#	print "$key $hash->{$key}\n";
		#}

		my $a1 = $hash->{from} - 1;    # Atom1
		my $a2 = $hash->{to} - 1;      # Atom2
		my $bo = $hash->{order};       # Bondorder

		# MMOD bonds are included twice in file
		next if $a2 < $a1;

		# Macromodel uses kekule aromatic representation
		# Bondorder 4 is of undefined order
		$bo = 0 if $bo == 4;

		bond_create($mol, $a1, $a2, $bo);
	}

	# Discard close bracket
	mae_get_tok($fr);
}

##################################################################
#
#	Maestro write routines
#
##################################################################

sub write_mae {

	#<
	#? Basic Maestro write routine.
	#; Requires: ensemble (or molecule), filename.
	#; Returns: undef if file open failed otherwise returns 1.
	#; Options: GZ (writes maegz format)
	#>

	my $ensemble = ens($_[0]);
	my $outfile  = $_[1];
	my $options  = uc($_[2] || '');

	silico_msg('d', "Argument 0 is not an ensemble of molecules\n") if ref($_[0] ne 'ARRAY');

	my $fr = open_molfile(undef, $outfile, undef, 'write', $options);
	my $FH = $fr->{FH};
	if (!defined $fr) {
		silico_msg('e', "Could not open $outfile");
		return undef;
	}

	silico_msg('c', "Writing Maestro file: $outfile");
	silico_msg('c', " using $options") if $options;
	silico_msg('c', "\n");

	foreach my $mol (@$ensemble) {
		write_mae_molecule($fr, $mol, $options);
		mol_print_warnings($mol);
	}

	mae_write_warnings($fr);
}

sub write_maegz {

	#<
	#? Basic Maestro maegz write routine (wrapper for write_mae)
	#; Requires: ensemble (or molecule), filename.
	#; Returns: undef if file open failed otherwise returns 1.
	#>

	my $ensemble = ens($_[0]);
	my $outfile  = $_[1];
	my $options  = uc($_[2] || '');

	write_mae($ensemble, $outfile, "$options GZ");
}

sub mae_write_warnings {

	my $fr = $_[0];

	silico_msg('n', "$fr->{WARN}{KEKULE} Molecules contained aromatic bondorders. Converted to Kekule.\n") if $fr->{WARN}{KEKULE};

}

sub write_mae_molecule {

	#<
	#? Write a single record from a pdb file.
	#; Requires: file record, molecule, options
	#. See read_pdb for general description
	#; Returns: 1 or undef if failed
	#>

	my $fr  = $_[0] || silico_msg('d', "File record not defined\n");
	my $mol = $_[1] || silico_msg('d', "Molecule not defined\n");
	my $options = $_[2] || '';    # Unused

	my $key;
	my $type;
	my %colourtable;

	my $FH = $fr->{FH};
	if (!defined $FH) {
		silico_msg('e', "FH not defined\n");
		fr_printout($fr);
		croak();
	}
	++$fr->{MOL_WRITE_COUNT};     # Note. Indexed from zero

	# Macromodel colours
	my %col = qw (H 21 Li 4 B 10 C 2 N 43 O 70 F 8 Mg 4 Si 14 P 15 S 13 Cl 9 Br 22);

	# Fix bonds
	molecule_check_and_fix_connectivity($mol);

	make_atom_names_unique($mol) if get_flag('unique-atom-names', 'l');

	# Generate macromodel atomtypes if any are missing
	# (macromodel fails to read files if we do not provide atomtypes);
	foreach my $atom (atoms($mol)) {
		next if $atom->{MMOD_TYPE};
		mol_type_mmod($mol);
		last;
	}

	# Convert bonds to kekule if we have any aromatic bondorders
      ATOM: foreach my $atom (atoms($mol)) {
		foreach my $bo (@{ $atom->{BORDERS} }) {
			if ($bo > 3) {

				#silico_msg('c', "Molecule contains aromatic bondorders. Converting to Kekule representation.\n");
				++$fr->{WARN}{KEKULE};
				convert_aromatic_bonds_kekule($mol);
				last ATOM;
			}
		}
	}

	# Make sure a single file does not contain duplicate compound names which
	# creates problems when reading in Pymol
	my $n = $mol->{NAME} || 'mol';
	if (defined $fr->{NAMES}{$n}) {

		$n .= "_" . ($fr->{MOL_WRITE_COUNT} + 1);
		$mol->{NAME} = $n;
	}
	++$fr->{NAMES}{$n};

	#print "Name: $mol->{NAME}\n";

	# Write file header
	if (!$fr->{MAE_HEADER_WRITTEN}) {

		print $FH "{\n  s_m_m2io_version\n  :::\n  2.0.0\n}\n\n";
		$fr->{MAE_HEADER_WRITTEN} = 1;
	}

	print $FH "f_m_ct {\n";
	print $FH "  s_m_title\n";

	#
	# Note - we need a better treatment of the cell data for maestro files
	#
	#

	# Crystal headers
	if (defined $mol->{CELL}[0]) {
		print $FH "   s_pdb_PDB_ID\n";
		print $FH "   r_pdb_PDB_CRYST1_a\n";
		print $FH "   r_pdb_PDB_CRYST1_b\n";
		print $FH "   r_pdb_PDB_CRYST1_c\n";
		print $FH "   r_pdb_PDB_CRYST1_alpha\n";
		print $FH "   r_pdb_PDB_CRYST1_beta\n";
		print $FH "   r_pdb_PDB_CRYST1_gamma\n";
		print $FH "   s_pdb_PDB_CRYST1_Space_Group\n";
		print $FH "   i_pdb_PDB_CRYST1_z\n";
		print $FH "   s_pdb_PDB_CLASSIFICATION\n";    # ??
	}

	# SDF_DATA headers
	if ($mol->{SDF_DATA}) {
		foreach $key (sort(keys %{ $mol->{SDF_DATA} })) {
			$type = $mol->{SDF_DATA_TYPE}{$key} || 's';
			next if $key eq 'title';
			print $FH "  $type\_m_$key\n";
		}
	}
	print $FH "  :::\n";
	my $nm = $mol->{NAME} || '';
	print $FH "  \"$nm\"\n";

	# Crystal data
	if (defined $mol->{CELL}[0]) {

		my $space_group = space_group_silico_to_pdb($mol->{SPACE_GROUP});
		if (!defined $space_group) {
			silico_msg('w', "Setting space group to P 1\n");
			$space_group       = 'P 1';
			$mol->{CELL_ALPHA} = 90;
			$mol->{CELL_BETA}  = 90;
			$mol->{CELL_GAMMA} = 90;
		}

		$mol->{CHAINS_PER_CELL} ||= 1;

		print $FH "   \"" . ($mol->{SDF_DATA}{PDB_CODE} || '') . "\"\n";
		print $FH "   $mol->{CELL}[0]\n";
		print $FH "   $mol->{CELL}[1]\n";
		print $FH "   $mol->{CELL}[2]\n";
		print $FH "   $mol->{CELL_ALPHA}\n";
		print $FH "   $mol->{CELL_BETA}\n";
		print $FH "   $mol->{CELL_GAMMA}\n";
		print $FH "   \"$space_group\"\n";
		print $FH "   $mol->{CHAINS_PER_CELL}\n";
		print $FH "   1\n";    # ??
	}

	# SDF_DATA
	if ($mol->{SDF_DATA}) {
		foreach $key (sort(keys %{ $mol->{SDF_DATA} })) {
			next if $key eq 'title';
			print $FH "  \"$mol->{SDF_DATA}{$key}\"\n";
		}
	}

	print $FH "  m_atom[$mol->{NUMATOMS}] {\n";
	print $FH "    # First column is atom index #\n";
	print $FH "    i_m_mmod_type\n";
	print $FH "    r_m_x_coord\n";
	print $FH "    r_m_y_coord\n";
	print $FH "    r_m_z_coord\n";
	print $FH "    i_m_residue_number\n";
	print $FH "    s_m_chain_name\n";
	print $FH "    i_m_formal_charge\n";
	print $FH "    i_m_atomic_number\n";
	print $FH "    r_m_charge1\n";
	print $FH "    s_m_pdb_residue_name\n";
	print $FH "    s_m_pdb_atom_name\n";
	print $FH "    i_m_color\n";    # Colour is included because maestro complains if it is missing

	print $FH "    :::\n";

	my $i = 0;
	foreach my $atom (@{ $mol->{ATOMS} }) {

		++$i;
		my $atomname = $atom->{NAME};
		$atomname =~ s/ //g;
		$atomname ||= 'X';
		$atomname = "\"$atomname\"";

		my $subname = $atom->{SUBNAME} || '';
		$subname =~ s/ //g;

		#$subname = substr($subname,0,4) if length($subname)> 4;
		$subname ||= "UNK";
		$subname = "\"$subname\"";

		my $chain = $atom->{CHAIN} || '';
		$chain =~ s/ //g;
		$chain ||= " ";
		$chain = "\"$chain\"";

		# Includes some rough fixes for undefined values
		# that could be improved
		printf $FH " %6i %4i %12.6f %12.6f %12.6f %5i %3s %8.3f %4i %8.3f %-8s %-8s %s\n",
		  $i, ($atom->{MMOD_TYPE} || '64'), $atom->{X}, $atom->{Y}, $atom->{Z}, ($atom->{SUBID} || 1), $chain,
		  ($atom->{FORMAL_CHARGE} || 0), ($atom->{ELEMENT_NUM} || 0),
		  ($atom->{CHARGE} || 0), ($subname || 'UNK'), ($atomname || 'X'),
		  ($atom->{MMOD_COLOUR} || $col{ $atom->{ELEMENT} } || 20);
	}

	print $FH "    :::\n";
	print $FH "  }\n";

	# Bonds

	print $FH "  m_bond[" . (2 * $mol->{NUMBONDS}) . "] {\n";

	#print $FH "  m_bond[$mol->{NUMBONDS}] {\n";
	print $FH "    # First column is bond index #\n";
	print $FH "    i_m_from\n";
	print $FH "    i_m_to\n";
	print $FH "    i_m_order\n";
	print $FH "    :::\n";

	$i = 0;
	foreach my $atom (@{ $mol->{ATOMS} }) {

		my $j = -1;

		foreach my $con (connected($atom, $mol)) {

			++$j;
			++$i;

			printf $FH " %6i %6i %6i %2i\n", $i, $atom->{NUM}, $con->{NUM}, $atom->{BORDERS}[$j];
		}
	}
	print $FH "    :::\n";
	print $FH "  }\n";

	print $FH "}\n";
}

##################################################################
#
#	Macromodel command-file/wrapper routines
#
##################################################################

sub find_schrodinger_exe {

	#<
	#? Find the Schrodinger executableis
	#; Reqires:  Executable name
	#; Returns: Path to Macromodel executable or undef
	#>

	my $name = $_[0];

	if (!defined $ENV{SCHRODINGER}) {
		silico_msg('e', "Environment variable '\$SCHRODINGER' is not defined." . " Macromodel does not appear to be installed on this system.\n");
		return undef;
	}

	my $exe = $ENV{SCHRODINGER} . "/$name";

	if (!-x $exe) {
		silico_msg('e', "$exe is not a valid executable file!\n");
		return undef;
	}

	return $exe;
}

sub set_debg17 {

	#<
	#? Work out if macromodel debug command DEBG 17 should be set
	#; Reqires:  molecule
	#; Returns: true or false
	#>

	my $mol = $_[0];

	my $debg;

	# Set DEBG 17 if FXAT is defined and any of the force constants are less than 100
	foreach my $atom (atoms($mol)) {
		if (defined $atom->{MMOD}{FXAT}) {
			next
			  if defined $atom->{MMOD}{FXAT_FORCECONST}
			  && (     $atom->{MMOD}{FXAT_FORCECONST} < 0
				|| $atom->{MMOD}{FXAT_FORCECONST} > 100);
			$debg = 1;
			last;
		}
	}

	return $debg;
}

sub mmod_min {

	#<
	#? Minimise a molecule using macromodel
	#. If filebase is not specified, temporary files are written to $Silico::temp_dir, unless $Silico::debug
	#  is set, in which case they are written to the current working directory.
	#; Requires:  Molecule, filebase, max number of steps, final gradient
	#  forcefield, solvent, minimsation method
	#; Returns: ensemble or undef if things go wrong
	#>

	my $mol      = $_[0];
	my $filebase = $_[1] || get_tempname();
	my $steps    = $_[2];
	my $grad     = $_[3];
	my $ff       = $_[4] || "opls3";
	my $solv     = $_[5] || "water";
	my $method   = $_[6];

	my $subs = 1 if $mol->{MMOD}{SUBS} || $mol->{MMOD}{FXAT};
	my $ens;
	$ens->[0] = $mol;

	my $val = write_mmod_com_setup($mol, "$filebase\.com", "$filebase\.maegz", "$filebase\_out.maegz", $ff, $solv);
	silico_msg('e', "Error writing mmod com file $filebase.com\n") if !defined $val;
	return undef                                                   if !defined $val;

	write_mmod_com_mini_stage("$filebase\.com", $steps, $grad, $method);

	write_maegz($ens, "$filebase\.maegz");

	return mmod_run($filebase);

}

sub mmod_min_h {

	#<
	#? Minimise hydrogens using macromodel.  All other atoms are constrained
	#; Requires:  Molecule, filebase, max number of steps, final gradient
	#  forcefield, solvent, minimsation method
	#; Returns: ensemble or undef if things go wrong
	#>

	foreach my $atom (atoms($_[0])) {
		next if $atom->{ELEMENT} eq 'H';
		$atom->{MMOD}{FXAT} = 1;
	}

	$_[0]->{MMOD}{FXAT} = 1;

	return mmod_min(@_);
}

sub mmod_mult_min {

	#<
	#? Use Batchmin to do a multiconformer minimisation
	#; Requires: molecule, output filebase, number of steps, solvent
	#; Returns: ensemble or undef if run failed
	#>

	my $mols     = $_[0];
	my $filebase = $_[1] || get_tempname();
	my $steps    = $_[2];
	my $grad     = $_[3];
	my $ff       = lc $_[4];
	my $solv     = lc $_[5];

	my $comfile = "$filebase.com";

	my $val = write_mmod_com_setup($mols->[0], $comfile, undef, undef, $ff, $solv);
	silico_msg('e', "Error writing mmod com file $filebase.com\n") if !defined $val;
	return undef                                                   if !defined $val;

	# Conjugate gradient minimisation
	write_mmod_com_mini_stage("$filebase\.com", $steps, $grad, 9, 1);
	open(OUTFILE, ">>$comfile") || die;

	write_maegz($mols, "$filebase\.maegz");

	my $nmols = mmod_run($filebase);

	my $i = 0;
	foreach my $mol (@$mols) {
		mol_copy_mmod_atom_data($mol, $nmols->[$i]);
		++$i;
	}

	silico_msg('e', "Error running mmod com file $filebase.com\n") if !defined $nmols;
	return $nmols;
}

sub mmod_anneal {

	#<
	#? Write Macromodel command file for simulated annealing run
	#. The first output structure produces is the initial, minimised structures
	#  one structure is then written after each annealing cycle
	#; Requires: molecule, output filebase, minimum temp, maximum temp,
	#  total cycle time (ps), number of annealing cycles, forcefield,
	#  solvent, molecular dynamics timestep (fs), number of structures to be saved from the
	#  cooling stage of the annealing run
	#; Returns: ensemble or undef if run failed
	#>

	my $mol       = $_[0];
	my $filebase  = $_[1] || get_tempname();
	my $mintemp   = $_[2];
	my $maxtemp   = $_[3];
	my $cycletime = $_[4] || 50;
	my $numcycles = $_[5] || 1;
	my $ff        = lc $_[6];
	my $solv      = lc $_[7];
	my $ts        = lc $_[8];
	my $numstruct = $_[9] || 0;

	my $steps    = 50;
	my $comfile  = "$filebase.com";
	my $uptime   = $cycletime / 5;
	my $downtime = 3 * $uptime;

	my $val = write_mmod_com_setup($mol, $comfile, undef, undef, $ff, $solv);
	silico_msg('e', "Error writing mmod com file $filebase.com\n") if !defined $val;
	return undef                                                   if !defined $val;

	write_mmod_com_mini_stage($comfile, $steps);

	for (1 .. $numcycles) {

		# Heat up using stochastic dynamics
		write_mmod_com_mdyn_stage($comfile, $uptime, undef, $ts, $mintemp, $maxtemp, undef, 1);

		# Run at max temp
		write_mmod_com_mdyn_stage($comfile, $uptime, undef, $ts, $maxtemp, $maxtemp, undef, 1);

		# Cool down
		write_mmod_com_mdyn_stage($comfile, $downtime, undef, $ts, $maxtemp, $mintemp, undef, 1, $numstruct);
		write_mmod_com_write_stage($comfile) if !$numstruct;    # Make sure that the final structure is written
	}

	my $ens;
	$ens->[0] = $mol;
	write_maegz($ens, "$filebase\.maegz");

	my $newmols = mmod_run($filebase);
	foreach my $nmol (@$newmols) {
		mol_copy_mmod_atom_data($mol, $nmol);
	}

	silico_msg('e', "Error running mmod com file $filebase.com\n") if !defined $newmols;
	return $newmols;
}

sub mmod_calc_energy {

	#<
	#? Use Batchmin to to a 'single point' energy calculation
	#. The output .mmo file is read and molecular energy vaules are stored as SDF_DATA fields
	#
	#; Requires: molecule, output filebase, number of steps, solvent
	#; Returns: ensemble or undef if run failed
	#>

	my $mol      = $_[0];
	my $filebase = $_[1] || get_tempname();
	my $ff       = lc $_[2];
	my $solv     = lc $_[3];

	my $comfile = "$filebase.com";

	my $val = write_mmod_com_setup($mol, $comfile, undef, undef, $ff, $solv);
	silico_msg('e', "Error writing mmod com file $filebase.com\n") if !defined $val;
	return undef                                                   if !defined $val;

	# Turn of sbc file
	$mol->{MMOD}{FXAT} = 0;
	$mol->{MMOD}{SUBS} = 0;

	open(OUTFILE, ">>$comfile") || die;
	pmmod('ELST');
	close(OUTFILE);

	my $ens;
	$ens->[0] = $mol;
	write_maegz($ens, "$filebase\.maegz");

	my $newmols = mmod_run($filebase);
	silico_msg('e', "Error running mmod com file $filebase.com\n") if !defined $newmols;

	my $hash = read_mmod_mmo("$filebase\_out.mmo");

	return undef if !defined $hash;

	foreach (keys(%$hash)) {
		$mol->{SDF_DATA}{$_} = $hash->{$_};
	}

	return $mol;
}

sub mmod_csearch {

	#<
	#? Write Macromodel command file for a conformational search run
	#; Requires: molecule or ensemble (all molecules must be identical).
	#  Optional or via flags: solvent, search_steps, filebase(otherwise temp filenames are used)
	#  Other arguments are passed via flags
	#; Returns: ensemble unless nowait flag is specified (then 1)
	#>

	my $mols         = ens($_[0]);
	my $solv         = $_[1] || get_lflag('solvent') || croak();
	my $search_steps = $_[2] || get_lflag('search-steps') || croak();
	my $filebase     = $_[3];
	my $nowait       = $_[4];

	my $noclean = 1 if $filebase;
	my $comfile = "$filebase.com";

	# Wait for license to become available
	my $used = schrodinger_wait_for_license_availability('bmin');

	# Write command files
	write_maegz($mols, "$filebase\.maegz");
	mmod_csearch_write_com($mols, $solv, $search_steps, $filebase);

	my $newmols = mmod_run($filebase, $nowait, $noclean, $used);

	return 1 if $nowait;

	silico_msg('e', "Error running mmod com file $filebase.com\n") if !defined $newmols;

	return $newmols;
}

sub mmod_csearch_write_com {

	#<
	#? Write Macromodel command file for a conformational search run
	#; Requires: molecule or ensemble (all molecules must be identical).
	#  Optional or via flags: solvent, search_steps, filebase(otherwise temp filenames are used)
	#  Other arguments are passed via flags
	#; Returns: ensemble or undef if run failed
	#>

	my $mols         = ens($_[0]);
	my $solv         = $_[1] || get_lflag('solvent') || croak();
	my $search_steps = $_[2] || get_lflag('search-steps') || croak();
	my $filebase     = $_[3] || get_tempname();

	my $ff     = get_lflag('force-field')   || croak();
	my $method = get_lflag('search-method') || croak();
	my $flap      = get_lflag('ring-flap', "Use FLAP option for ring searchers in mcmm and low-mode searches");
	my $min_steps = get_lflag('min-steps') || croak();
	my $grad      = get_lflag('min-gradient-cutoff') || croak();
	my $window    = get_lflag('energy-window') || 50;

	my $noclean = 1 if $filebase;

	silico_msg('w', "Flap is not implemented (yet)") if $flap;

	my $comfile     = "$filebase.com";
	my $mini_method = 9;                 # TNCG

	# Macromodel MCSS structure selection
	my $struct_select = 2;               #
	$method = lc $method;

	my $val = write_mmod_com_setup($mols->[0], $comfile, undef, undef, $ff, $solv, 0);
	silico_msg('e', "Error writing mmod com file $filebase.com\n") if !defined $val;
	return undef                                                   if !defined $val;

	open(OUTFILE, ">>$comfile") || die;

	pmmod('CRMS', 0, 0, 0, 0, 0, 0.5);    # Convergence RMS 0.5 Ang

	my $mcop_frac = 0;
	if ($method eq 'lmcs') {

		# LMCS low mode search
		pmmod('LMCS', $search_steps, 0, 0, 0);

	} elsif ($method eq 'lmcs-shake') {

		# LMCS shake (Do LMCS with minimisation to reduce distortion)
		pmmod('LMCS', $search_steps, 0, 0, 0, 1);

	} elsif ($method eq 'lmcs-mcmm') {

		# Mixed LMCS low mode MCMM 50% frequency of each
		pmmod('LMCS', $search_steps, 0, 0, 0, 0, 0, 3, 6);
		$mcop_frac = 0.5;

	} elsif ($method eq 'mcmm') {

		# Monte Carlo multiple minimum
		pmmod('MCMM', $search_steps, 0, 0, 0, -1);

	} elsif ($method eq 'spmc') {

		# Systematic pseudo-Monte Carlo multiple minimum
		pmmod('SPMC', $search_steps, 0, 0, 0, -1);

	} else {
		silico_msg('d', "Invalid method ($method) chosen\n");
	}

	pmmod('NANT');
	pmmod('MCNV', 1, 5);                          # Monte Carlo number of variables
	pmmod('MCSS', 2, 0, 0, 0, $window);
	pmmod('MCOP', 1, 0, 0, 0, $mcop_frac);
	pmmod('DEMX', 0, 833, 0, 0, $window, 100);    # Delta E max energy windowing
	pmmod('MSYM');
	pmmod('AUTO', 0, 2, 1, 1, 0, -1, 0, 3);       # No torsional constraints (Arg 8 = 3)
	pmmod('CONV', 2, 0, 0, 0, 0.1);               # Convergence criterea
	pmmod('MINI', $mini_method, 0, $min_steps);

	close OUTFILE;

	return 1;
}

sub mmod_run {

	#<
	#? Run macromodel
	#. Will wait for output and read it in. Setting the nowait flag will cause it to return as soon as the
	#  job has started, without returning output
	#; Requires:  filebase of macromodel com file, flag not to wait for job to complete, flag not to remove run files,
	#; number of used licenses (schrodinger_wait_for_license_availability has already been run);
	#; Returns: molecule or 1 if nowait flag is set.
	#>

	my $filebase = $_[0] || croak();
	my $nowait   = $_[1];
	my $noclean  = $_[2];
	my $used     = $_[3];

	my $bmin_exe = find_schrodinger_exe('bmin');

	if (!$bmin_exe) {
		silico_msg('d', "Macromodel executable not found.\n");
	}

	$used = schrodinger_wait_for_license_availability('bmin') if !defined $used;

	silico_msg('c', heading("Running Macromodel"));

	my $run = "$bmin_exe $filebase -NO_REDIRECT";
	$run .= " -WAIT"           if !$nowait;
	$run .= " > $filebase.log" if (!get_flag("print-macromodel", 'l') || $Silico::debug);
	$run .= " &"               if $nowait;

	print "run: $run\n";
	my $failure = system($run);

	if ($failure) {
	 	system ("mv $filebase.* .");
		silico_msg('d', "Macromodel failed to execute correctly. Moved files $filebase.* to current working directory. Please check log file\n");
		return undef;
	}

	if ($nowait) {

		# If we are not waiting for job to complete, then
		# make sure that license is checked out before we return
		# to start a new job
		schrodinger_wait_for_license_use('bmin', $used);
		return 1;
	}

	my $outfile = "$filebase\_out.maegz";
	my $mols    = read_mol_any($outfile);

	#ensemble_printout($mols);

	silico_msg('d', "Macromodel appears to have failed. Could not read any molecules from file $outfile") if !defined $mols;

	# Clean up temporary files if debug is not set and no filebase was provided
	system "rm $filebase.maegz $filebase.com $filebase\_out.maegz $filebase.log" if !($noclean || $Silico::debug);

	return $mols;
}

sub mmod_logfile_normal_termination {

	#<
	#? Find if Batchmin logfile terminated properly
	#>

	my $infile = $_[0];

	open(LOGFILE, "$infile") || return 0;

	while (<LOGFILE>) {
		return 1 if /BatchMin: normal termination/;
	}

	close(LOGFILE);

	silico_msg('c', "Found unfinished logfile $infile\n");

	return 0;
}

sub mmod_logfile_num_conf {

	#<
	#? Count number of Conf lines in conformational search logfile
	#>

	my $infile = $_[0];

	open(LOGFILE, "$infile") || return 0;

	my $num = 0;
	while (<LOGFILE>) {
		next if !/Conf/;
		($num) = $_ =~ /Conf\s*(\d*)/;
	}

	close(LOGFILE);

	return $num;
}

sub mmod_make_complex {

	#<
	#? Create a receptor/ligand complex for use with Macromodel
	#. Can be used to specify fixed and frozen atoms in the receptor
	#. Defines atom->{MMOD}{FXAT} flag for receptor atoms and
	#  atom->{MMOD}{SUBS} flag for ligand.  Also sets flags
	#  mol->{MMOD}{FXAT} and mol->{MMOD}{SUBS}
	#. DEBG 17 needs to be set so that all intramolecular interactions are not discarded
	#  by Macromodel
	#; Requires: Protein (receptor), ligand, forceconst1, forceconst2, forceconst3
	#  freeze distahce, constraint scheme, flag to constrain hydrogens
	#  freeze atoms more than this distance from the ligand
	#; Returns:  Molecule or undef if things go wrong
	#>

	my $prot        = deep_copy($_[0]);
	my $lig         = deep_copy($_[1]);
	my $forceconst  = $_[2] || 1;
	my $forceconst2 = $_[3] || 0;
	my $forceconst3 = $_[4] || 0;
	my $freezedist  = $_[5];
	my $scheme      = $_[6] || 'A';
	my $fixh        = $_[7] || 0;

	my $ignore_dist = 2;

	silico_msg('c', "Using constraint scheme: $scheme\n");

	# Chain of ligand is changed to Z to ensure that naming
	# is unique (required by Macromodel)
	foreach my $atom (atoms($prot)) {
		$atom->{SEGID} = 'PROT';
		$atom->{CHAIN} = 'z' if $atom->{CHAIN} eq 'Z';
	}
	foreach my $atom (atoms($lig)) {
		$atom->{SEGID} = 'LIG';
		$atom->{CHAIN} = 'Z';
	}

	my $ens;
	@$ens = ($prot, $lig);
	my $newens  = ensemble_consolidate($ens);
	my $complex = $newens->[0];

	label_aa_backbone($complex);
	label_waters($complex);

	foreach my $atom (atoms($complex)) {

		if ($atom->{SEGID} eq 'PROT') {

			my $dist = distance_atoms_min($atom, atoms($lig));
			$atom->{MMOD}{FXAT} = 1;
			$newens->[0]{MMOD}{FXAT} = 1;

			# Set original coordinates for atom constraints
			$atom->{MMOD}{FXAT_X} = $atom->{X};
			$atom->{MMOD}{FXAT_Y} = $atom->{Y};
			$atom->{MMOD}{FXAT_Z} = $atom->{Z};

			# Constrained atoms using force constant scheme
			if ($dist <= $freezedist) {
				if ($scheme eq 'A') {

					#
					# Single valued force constant
					#

					if (!$fixh && $atom->{ELEMENT} eq 'H') {
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst3;
					} else {
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst;
					}

				} elsif ($scheme eq 'B') {

					#
					# Distance dependent force constant
					#

					if (!$fixh && $atom->{ELEMENT} eq 'H') {

						# Hydrogens
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst3;
					} else {

						# Everything else
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst * ($dist / $freezedist)**2;
					}

				} elsif ($scheme eq 'C') {

					#
					# Constrained backbone with very weak constraints on other atoms
					#

					if (!$fixh && $atom->{ELEMENT} eq 'H') {

						# Hydrogens
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst3;
					} elsif ($atom->{FG}{AA_N} || $atom->{FG}{AA_C} || $atom->{FG}{AA_CA} || $atom->{FG}{AA_O}) {

						# Backbone atoms
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst;
					} elsif ($atom->{FG}{W} || $atom->{SUBNAME} =~ 'SO4') {

						# Sulfates and waters
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst2;
					} else {

						# Everything else
						$atom->{MMOD}{FXAT_FORCECONST} = $forceconst3;
					}

				} else {
					silico_msg('d', "Invalid constraint scheme: $scheme\n");
				}
			}

			# Frozen atoms
			if ($dist > $freezedist) {
				$atom->{MMOD}{FXAT}            = 1;
				$atom->{MMOD}{FXAT_FORCECONST} = -1;
				$newens->[0]{MMOD}{FXAT}       = 1;
			}

			# Ignored atoms
			if ($dist >= $freezedist + $ignore_dist) {
				$atom->{MMOD}{IGNORE} = 1;
			}
		}

		if ($atom->{SEGID} eq 'LIG') {
			$atom->{MMOD}{SUBS} = 1;
			$newens->[0]{MMOD}{SUBS} = 1;
		}
	}

	return $newens->[0];
}

sub mol_copy_mmod_atom_data {

	my $mol1 = $_[0];    # Source
	my $mol2 = $_[1];    # Target

	$mol2->{MMOD} = deep_copy($mol1->{MMOD});

	my $i = 0;
	foreach my $atom1 (atoms($mol1)) {

		my $atom2 = $mol2->{ATOMS}[$i];

		next if !defined $atom1->{MMOD};
		$atom2->{MMOD} = deep_copy($atom1->{MMOD});

		++$i;
	}
}

##################################################################
#
# Routines to write out macromodel command file commands
#
##################################################################

sub write_mmod_com_setup {

	#<
	#? Write initial section of Macromodel command file for a minimisation or other simulation
	#; Requires: molecule, output filename (required), structure input filename (opt),
	#  structure output filename (opt), forcefield (opls3, opls2005, opls2001, amber94, amber, mm3, mm2),
	#  solvent (water, chloroform, octanol, none), minimisation method (macromodel codes, default 20),
	#; Returns: undef if write failed
	#>

	my $mol                = $_[0];
	my $outfile            = $_[1] || silico_msg('d', "No output filename provided\n");
	my $struct_input_file  = $_[2] || (get_filebase($outfile) . ".maegz");
	my $struct_output_file = $_[3] || (get_filebase($outfile) . "_out.maegz");
	my $ff                 = lc $_[4] || 'opls3';
	my $solv               = lc $_[5] || 'water';

	my $ffcode;
	my $solvcode;

	my $sbcfile = $outfile;
	$sbcfile =~ s/.com$//;
	$sbcfile .= ".sbc";

	if ($ff eq 'opls3' || $ff eq 'opls' || !$ff) {
		$ffcode = 16;
	} elsif ($ff eq 'opls2005') {
		$ffcode = 14;
	} elsif ($ff eq 'opls2001') {
		$ffcode = 11;
	} elsif ($ff =~ /mmff/) {
		$ffcode = 10;
	} elsif ($ff eq 'amber94') {
		$ffcode = 5;
	} elsif ($ff =~ /^amber/) {
		$ffcode = 3;
	} elsif ($ff =~ /^mm3/) {
		$ffcode = 2;
	} elsif ($ff =~ /^mm2/) {
		$ffcode = 1;
	} else {
		silico_msg('d', "Error specifying forcefield; $ff\n");
	}

	if ($solv eq 'none') {
		$solvcode = 0;
	} elsif ($solv eq 'water') {
		$solvcode = 1;
	} elsif ($solv eq 'chloroform') {
		$solvcode = 5;
	} elsif ($solv eq 'octanol') {
		$solvcode = 9;
	} else {
		silico_msg('d', "Error specifying solvent; $solv\n");
	}

	my $subs = write_mmod_sbc($mol, $sbcfile);

	my $success = open(OUTFILE, ">$outfile");
	if (!$success) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing Macromodel command file: $outfile\n");

	print OUTFILE "$struct_input_file\n";
	print OUTFILE "$struct_output_file\n";

	pmmod('MMOD', 0, 1);
	pmmod('FFLD', $ffcode, 1, 0, 0, 1);
	pmmod('SOLV', 3, $solvcode) if $solvcode;
	pmmod('SUBS') if $subs;
	pmmod('EXNB');
	pmmod('BDCO', 0, 0, 0, 0, 89.4427, 99999.0);
	pmmod('READ');
	pmmod('DEBG', 17) if set_debg17($mol);

	close(OUTFILE);
	return 1;
}

sub pmmod {

	#<
	#? Write a formatted line of a macromodel com file
	#; Requires: arguments as an array
	#>

	my $n = 8 - $#_;

	silico_msg('d', "Wrong number of arguments ($#_) provided to pmmod\n") if $n < 0;

	for (1 .. $n) {
		push @_, 0;
	}

	#print "pmmod: @_\n";

	printf OUTFILE " %4s  %6d %6d %6d %6d %10.4f %10.4f %10.4f %10.4f\n", @_;
}

sub write_mmod_com_mini_stage {

	#<
	#? Write the minimisation stage of a macromodel com file
	#; Requires: number of steps, final gradient, macromodel minimisation method code.
	#  Flag to loop minimisation
	#>

	my $outfile = $_[0];
	my $steps   = $_[1] || 500;
	my $grad    = $_[2] || 0.1;
	my $method  = $_[3] || 9;     # TNCG
	my $mult    = $_[4];

	open(OUTFILE, ">>$outfile") || die;

	pmmod('BGIN') if $mult;
	pmmod('READ') if $mult;
	pmmod('CONV', 2, 0, 0, 0, $grad);
	pmmod('MINI', $method, 0, $steps);
	pmmod('END ') if $mult;

	close OUTFILE;
}

sub write_mmod_com_rewind_stage {

	#<
	#? Write a rewind step
	#; Requires: flag to rewind input (otherwise rewinds output)

	open(OUTFILE, ">>$_[0]") || die;
	pmmod('RWND', $_[1]);

	#pmmod ('READ');
	close OUTFILE;
}

sub write_mmod_com_write_stage {

	#<
	#? Write a macromodel 'WRIT' command to a macromodel com file
	#; Requires:outfile
	#>

	my $outfile = $_[0];

	open(OUTFILE, ">>$outfile") || die;
	pmmod('WRIT');
	pmmod('READ');
	close OUTFILE;
}

sub write_mmod_com_mdyn_stage {

	#<
	#? Write the molecular dynamics stage of a macromodel com file
	#; Requires: output filename, simulation time (ps),
	#  temperature (K), timestep, initial temp (opt), final temp (opt),
	#  shake flag (0 - none, 1 - hydrogens, 2 - all bonds),
	#  stochastic dynamics flag, output frequency (ps)
	#>

	my $outfile  = $_[0] || die;
	my $time     = $_[1] || 500;    # ps
	my $temp     = $_[2] || 300;
	my $timestep = $_[3] || 1.5;
	my $inittemp = $_[4];
	my $finaltemp  = $_[5];
	my $shake      = $_[6] || 0;    #  0 none, 1 hydrogens, 2 all bonds
	my $stochastic = $_[7] || 0;
	my $numstruct  = $_[8] || 0;

	$temp = $inittemp if $inittemp;

	open(OUTFILE, ">>$outfile") || die;

	# Output (sampling) options for MD
	pmmod('DEBG', 17);              # Prevent FXAT from removing all interatomic interactions;
	pmmod('MDSA', $numstruct);

	# Initial and final temperatures
	pmmod('MDIT', 0, 0, 0, 0, $inittemp)  if $inittemp;
	pmmod('MDFT', 0, 0, 0, 0, $finaltemp) if $finaltemp;

	# Dynamics
	pmmod('MDYN', 0, $shake, $stochastic, 0, $timestep, $time, $temp);

	close OUTFILE;
}

sub write_mmod_sbc {

	#<
	#? Write Macromodel constraint commands (SUBS, FXAT) to .sbc or .com files
	#. Uses flags $atom->{MMOD}{XXX} where XXX is SUBS, FXAT, etc
	#. Force constants and other MMOD parameters are stored on a per atom basis
	#  in $atom->{MMOD}{ZZZ} where ZZZ can be FXAT_WELL_HALF_WIDTH, FXAT_FORCECONST, FXAT_X, FXAT_Y, FXAT_Z
	#; Requires: molecule, filename
	#; Returns: undef if write failed
	#>

	my $mol        = $_[0];
	my $outfile    = $_[1];
	my $forceconst = $_[2] || 10;    # Default forceconst if not provided in molecule file

	my @array;
	my $success;

	return 0 if !$mol->{MMOD}{SUBS} && !$mol->{MMOD}{FXAT};

	# Overwrite any existing .sbc file but append to .com files
	# This allows insertion of SUBS and FXAT commands into either
	# sbc or com files
	if ($outfile =~ /sbc$/) {
		$success = open(SBCFILE, ">$outfile");
	} else {
		$success = open(SBCFILE, ">>$outfile");
	}
	if (!$success) {
		silico_msg('e', "Can not create or open file $outfile for writing!\n");
		return undef;
	}

	silico_msg('c', "Writing Macromodel sbc file: $outfile\n");

	# SUBS command
	my $i = 0;
	my $j = -1;
	foreach my $atom (atoms($mol)) {

		++$j;
		next if !$atom->{MMOD}{SUBS};

		++$i;
		push @array, $atom->{NUM};

		if ($i == 4 || $j == $#{ $mol->{ATOMS} }) {

			printf SBCFILE " SUBS  %6i %6i %6i %6i %10.5f %10.5f %10.5f %10.5f\n", @array, 0, 0, 0, 0, 0, 0, 0, 0;
			$i     = 0;
			@array = ();
		}
	}

	# FXAT command
	$i = 0;
	foreach my $atom (atoms($mol)) {
		next if !$atom->{MMOD}{FXAT};
		next if $atom->{MMOD}{IGNORE};
		printf SBCFILE " FXAT  %6i %6i %6i %6i %10.5f %10.5f %10.5f %10.5f\n",
		  $atom->{NUM}, 0, 0, ($atom->{MMOD}{FXAT_WELL_HALF_WIDTH} || 0),
		  ($atom->{MMOD}{FXAT_FORCECONST} || $forceconst), ($atom->{MMOD}{FXAT_X} || 0),
		  ($atom->{MMOD}{FXAT_Y} || 0), ($atom->{MMOD}{FXAT_Z} || 0);
	}

	# FXDI FXBA FXTA not implemented (yet)

	close SBCFILE;
	return 1;
}

sub read_mmod_mmo {

	#<
	#? Read Macromodel .mmo energy file
	#; Requires: filename,
	#; Returns: ensemble or undef if failed
	#>

	my $infile = $_[0];

	my $hash;

	my $fr = open_molfile(undef, $infile, undef, 'read');
	my $FH = $fr->{FH};
	if (!defined $fr) {
		silico_msg('e', "Could not open $infile");
		return undef;
	}

	silico_msg('c', "Reading mmod mmo file: $infile\n");

	while (<$FH>) {

		my $flag;
		if (/^Total energy:/) {
			$flag = 1;
		}

		next if !$flag;

		last if !/\w/;

		s/kJ.*//g;    # Discard line after 'kJ'
		s/\s+/ /g;
		my @f = split ":";    # Split on colon
		$f[0] =~ s/^\s*//;    # Removes spaces
		$f[0] =~ s/\s*$//;
		$f[0] =~ s/\s+/_/g;
		$f[0] = lc $f[0];     # Make key lower case
		$f[0] =~ s/-/_/g;     # Remove dashes from key;
		next if !$f[0];
		$hash->{"MMOD_$f[0]"} = $f[1] || '';

		#print "'$f[0]' '$f[1]'\n";
	}

	silico_msg('e', "No data was read from mmo file $_[0]\n") if !defined $hash;

	return $hash;
}

sub mmod_qikprop {

	#<
	#? Run qikprop
	#; Requires:  filebase of macromodel com file, flag not to remove run files
	#; Returns: molecule or undef if things go wrong
	#>

	my $mols     = ens($_[0]);
	my $filebase = $_[1] || get_tempname();
	my $fast     = $_[2];

	my $options = '';
	$options = '-nosim -fast -nsim 0';

	my $qikprop_exe = find_schrodinger_exe('qikprop');

	if (!-x $qikprop_exe) {
		silico_msg('w', "Qikprop executable not found.\n");
		return undef;
	}

	my $used = schrodinger_wait_for_license_availability('qikprop');

	write_maegz($mols, "$filebase.maegz");

	silico_msg('c', heading("Running Qikprop"));

	my $run = "$qikprop_exe $options $filebase.maegz -NO_REDIRECT > $filebase.log";
	print "run: $run\n";
	my $failure = system($run);

	schrodinger_wait_for_license_use($used);

	return $failure;
}

#
# Queueing Schrodinger jobs
#

sub schrodinger_wait_for_license_availability {

	#<
	#? Wait until there a licenses available
	#; Requires: program_name, quiet flag
	#; Returns: Nothing
	#>

	my $prog  = $_[0];
	my $quiet = $_[1];

	$prog = lc($prog);

	silico_msg('c', "Waiting until license is free\n") if !$quiet;

	# Maximum number of schrodinger tokens to use
	my $maxtok = maxtok_file() || get_lflag('maxtok') || 0;
	silico_msg('c', "Maximum number of Schrodinger tokens is set to $maxtok\n") if $maxtok;

	my $requested = 0;
	my $name;
	if ($prog eq 'glide') {
		$requested = 5;
		$name = "GLIDE_MAIN";
	} elsif ($prog eq 'bmin') {
		$requested = 2;
		$name = "MMOD_MACROMODEL";
	} elsif ($prog eq 'qikprop') {
		$requested = 2;
		$name = "QIKPROP_MAIN";
	} else {
		silico_msg('d', "Unknown Schrodinger program '$prog'");
	}

      LOOP: while (1) {

		my ($issued, $used) = schrodinger_check_tokens($name, $requested, $quiet);

		my $max = $issued;
		$max = $maxtok if $maxtok && $maxtok < $max && $maxtok >= 0;

		if ($max - $used - $requested >= 0) {

			silico_msg('c', "  Let's go!\n") if !$quiet;
			return $used;
		}

		# Delay between license checks
		sleep(60);
	}

	return 0;
}

sub maxtok_file {

	my @f = glob "maxtok*";

	return undef if !defined $f[0];

	my $n = $f[0];
	$n =~ s/\D//g;
	$n = int(abs($n));

	silico_msg('c', "File $f[0] sets maxtok to $n\n");

	return $n;
}

sub schrodinger_wait_for_license_use {

	#<
	#? Wait until requested license has been used
	#; Requires: program_name, quiet flag
	#; Returns: Nothing
	#>

	my $name = $_[0];
	my $oldused = $_[1];
	my $quiet   = $_[2];

	silico_msg('c', "Waiting until license is used\n") if !$quiet;

	my $count = 0;
      LOOP: while (1) {

		++$count;
		if ($count == 4) {
			silico_msg('c', "Giving up waiting\n") if !$quiet;
			last;
		}

		sleep(20);
		my (undef, $used) = schrodinger_check_tokens($name, 0, $quiet);
		last if $oldused != $used;
	}

	return 1;
}

sub schrodinger_check_tokens {

	my $name = $_[0];
	my $requested = $_[1];
	my $quiet     = $_[2];

	# Redirect STDOUT and STDERR to $string
	# lictool can print errors to stderr
	my $string = `\$SCHRODINGER/run lictool status 2>&1`;

	if (!defined $string) {
		silico_msg('e', "licadmin returned undefined value\n");
		return (undef, undef);
	}

	my @f = split "\n", $string;

	my $issued;
	my $used;

	foreach (@f) {

		next if !/Users of $name/;
		m/Total of(.*)licenses issued; Total of(.*)licenses in use/;
		$issued = $1;
		$used = $2;
		print "$name licenses issued: $issued used: $used\n";
		last;
	}

	if (!defined $issued || !defined $used) {

		print "$string\n";
		silico_msg('d', "Error reading license output\n");
	}

	$issued =~ s/ //g;
	$used   =~ s/ //g;

	silico_msg('c', "Schrodinger licenses. Requested: $requested.  In use: $used. Issued: $issued.\n") if !$quiet;

	return $issued, $used;
}

return 1;
