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
#! silico_tinker.pm
#? Silico routines specific to Jay Ponder's TINKER and SLEUTH software
#. $Revision: 1.40.2.1.2.4 $
#>

use strict;
package Silico;

sub read_tinker_xyz {

	

	#<
	#? Tinker read routine.
	#; Requires: filename, (optional) start, max molecules, options string, Tinker format flag
	#; Returns: ensemble or undef if failed
	#>

	my $infile = $_[0];
	my $start = $_[1] || get_flag('ss', 's') || 0; 
        my $max_molecules = $_[2] || get_flag('ms', 's');
	my $options = uc ($_[3] || '');
	
	my $mols;
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
		
	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
        if (!defined $fr) {
                silico_msg('e', "Can not open file $infile for reading\n");
                return undef;
        }
	
	my $FH = $fr->{FH};
	
	if ($options !~ /\bQUIET\b/) {
	
		
		silico_msg('c', "Reading tinker xyz file: $infile\n");
	}
	
	my $i = -1;
	my $molcount = 0;
	while (my $mol = read_tinker_xyz_single($fr, $options)) {
		++$i;
		next if $start and ($i < $start);
		++$molcount;
		push @$mols, $mol;
		last if (defined $max_molecules) && $molcount == $max_molecules;
	}
	
	close ($FH);
     
	return $mols;
}
	
	
sub read_tinker_xyz_single {

	#<
	#? Read in a single record from a Tinker xyz file
	#. Tinker format has no comment line (comments are on first line after number of atoms)
	#; Requires: file record, options, flag to read Tinker format xyz
	#; Returns: molecule or undef
	#>

	my $fr = $_[0];
	my $options = $_[1] || '';
	
	my $FH = $fr->{FH};
	
	my $mol;
	
	# Read Header
	my $line = <$FH>;
	
	return if !defined $line;
	
	++$fr->{MOL_READ_COUNT};
	
	chomp $line;
	
	$line =~ s/^\s*//;
	
	my @f = split (" ", $line);
	$mol->{NUMATOMS} = $f[0];
	
	$mol->{NAME} = join (":", @f[1..($#f)]);
	chomp $mol->{NAME};
	

	# Create all the atoms first because we will need to
	# create bonds to atoms that have not been read in yet.
	
	for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {
		
		my $at;
		$at->{NUM} = $i+1;		#Atom number
		$mol->{ATOMS}[$i] = $at;
	}

	for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {

		$line = <$FH>;
		
		if (!defined $line) {
			silico_msg('e', "Premature end of file!\n");
			return undef;
		}
		
		chomp $line;
		
		my @f = split (" ", $line);

		my $atom = $mol->{ATOMS}[$i];

		$atom->{NUM} = $i+1;
		$atom->{NAME} = $f[1];		#Atom name
		$atom->{X} = $f[2];		#x
		$atom->{Y} = $f[3];		#y
		$atom->{Z} = $f[4];		#z
		$atom->{TYPE} =  $f[5];
		
		$atom->{SUBNAME} = 'SUB';
		$atom->{SUBID} =1;
		
		# Create bonds;
		for (my $j = 6; $j <= 9; ++$j) {

			last if !defined $f[$j];
			last if $f[$j] eq "";

			# The internal connectivity matrix is indexed
			# from zero. The connectivities in the input
			#  file are indexed from 1

			my $atnum = $f[$j]-1;

			# Create bonds in only one direction
			next if $atnum <= $i;

			my $bo = 1;
			bond_create($mol, $i, $atnum, $bo);
		}
	}

	# Work out the atom elements
	foreach my $atom (atoms($mol)) {

		atom_guess_element($atom, $mol);
	}

	mol_set_source_data($mol, $fr, 'xyz');

	# No bondorder information in a .xyz file
	$mol->{HAS_BAD_BONDORDERS} = 1;
	
	return $mol;
}


sub write_tinker_xyz {

	#<
	#? Tinker write routine.
	#; Requires: ensemble (or molecule), filename, forcefield file (optional) to
	#  type atoms.
	#; Returns: undef if file open failed otherwise returns 1.
	#>

	my $ensemble = ens($_[0]);
	my $filename = $_[1];
	my $options = uc ($_[2] || '');

	my $atom;
	my $ext;
	my $filebase;
	my $i;
	my $j;
	my $mol;
	my $molcount;
	my $outfile;
	my $success;
	my $type;
	my $warn;

	$filename = get_outfile($ensemble, 'xyz') if !$filename;
	$filebase = get_filebase($filename);
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	ensemble_printout($ensemble, 'all') if $Silico::debug;

	$molcount = -1;

	foreach $mol (@$ensemble) {

		++$molcount;
		$warn = 0;
		
		if (!defined $mol->{NUMATOMS}) {

			silico_msg('e', "\$mol->NUMATOMS is not defined!\n");
			return undef;
		}

		if ($molcount > 0) {
			$outfile = $filebase.sprintf "%04d", $molcount.".".$ext;
		} else {
			$outfile = $filename;
		}

		$success = open (TINOUT, ">$outfile");
		
		if (!$success) {
			silico_msg('e', "Can not create or open file $outfile for writing!\n");
			return undef;
		}

		if ($options !~ /\bQUIET\b/) {
			silico_msg('c', "Writing xyz file $outfile");
			silico_msg('c', " using $options") if $options;
			silico_msg('c',  "\n");
		}

		# Print header
		print TINOUT " $mol->{NUMATOMS} $mol->{NAME}\n";
		
		# Generate tinker atom names and atom types
		for ($i=0; $i < $mol->{NUMATOMS}; ++$i) {

			$atom = $mol->{ATOMS}[$i];

			$atom->{NAME} =~ s/ //;
			
			# If atom type is not numeric then change it to zero
			$type = $atom->{TINKER_TYPE};
			if (!defined $type || $type !~ /^\d+$/) {
				$type = 0;
				++$warn;
			}
			$type =~ s/ //g;

			printf TINOUT "%6d  %4s %12.6f %12.6f %12.6f %6d",
				$i+1,$atom->{NAME}, $atom->{X}, $atom->{Y}, $atom->{Z}, $type;

			# Bonds
			for ($j=0; $j<=3; ++$j) {

				last if (!defined $atom->{CONNECT}[$j]);
				
				printf TINOUT " %6d", $atom->{CONNECT}[$j]+1;
			}
			print TINOUT "\n";
		}
		close (TINOUT);
		
		if ($warn) {
			silico_msg('w', "$warn atom types were set to \"0\". Check output file.\n");
		}
	}
	

	return 1;
}

sub write_tinker_prm {
	
	#<
	#? Write a TINKER force field parameter file
	#; Requires: Ensemble, filename
	#; Returns: 0 if any molecules failed, 1 if success
	#>
	
	my $molecules = $_[0];
	my $filename = $_[1];
	my $ff = $_[2];
	my $options = $_[3] || '';
	
	my $count;
	my $mol;
	my $failure;
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	foreach $mol (@$molecules) {
	
		my $prmfile;
		
		if ($#{$molecules} > 0) {
			++$count;
			$count = sprintf "%03d", $count;
			$prmfile = get_filebase($filename)."_$count.prm";
		} else {
			$prmfile = $filename;
		}
		
		if (!write_tinker_prm_mol($mol, $prmfile, $ff, $options)) {
			silico_msg('w', "Could not write a Tinker parameter file for molecule \"$mol->{NAME}\"\n",
					"($count of ".($#{$molecules}+1).")!\n",
					"Skipping.\n");
			++$failure;
			next;
		}
	}
	
	return 0 if $failure;
	return 1;
}

sub write_tinker_prm_mol {

	my $mol		= $_[0];
	my $filename	= $_[1];
	my $ff		= $_[2];
	my $options	= uc ($_[3] || '');
	
	my $success;
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	$success = open(PRM, ">$filename");
	
	if (!defined $success) {
		silico_msg('e', "Can not create or open file $filename for writing!\n",
				"Aborting write.\n");
		return undef;
	} else {
		silico_msg('c', "Writing Tinker format parameter file: $filename\n")
			if ($options !~ /\bQUIET\b/);
		
		print PRM "\n";
		print PRM "      ##############################\n";
		print PRM "      ##                          ##\n";
		print PRM "      ##  Force Field Definition  ##\n";
		print PRM "      ##                          ##\n";
		print PRM "      ##############################\n";
		print PRM "\n";
		
		if ($ff eq 'oplsaa') {write_tinker_prm_oplsaa()}
		
		else {
			silico_msg('e', "No instructions for force field $ff.\n",
					"Aborting write.\n");
			close PRM;
			return undef;
		}
	}

	write_tinker_prm_mol_atoms($mol, $ff, $options);
	write_tinker_prm_mol_vdw($mol, $ff);
	write_tinker_prm_mol_bonds($mol, $ff);
	write_tinker_prm_mol_angles($mol, $ff);
	write_tinker_prm_mol_impropers($mol, $ff);
	write_tinker_prm_mol_dihedrals($mol, $ff);
	write_tinker_prm_mol_charges($mol, $ff);
	
	return 1;
}

sub write_tinker_prm_mol_atoms {

	my $mol = $_[0];
	my $ff = $_[1];
	my $options = $_[2] || '';
	
	write_tinker_prm_mol_atoms_oplsaa($mol, $options) if $ff eq 'oplsaa';

}

sub write_tinker_prm_mol_atoms_oplsaa {

	my $mol = $_[0];
	my $options = $_[1] || '';
	
	my $array;
	my $c = 0;
	my %class;
	my $atom;
	my $hash;
	my $heavycount;
	my $i;
	my $list;
	my $mass;
	my $type;
	
	# Load atomic masses
	atomic_masses();
	
	print PRM "\n";
	print PRM "      #############################\n";
	print PRM "      ##                         ##\n";
	print PRM "      ##  Atom Type Definitions  ##\n";
	print PRM "      ##                         ##\n";
	print PRM "      #############################\n";
	print PRM "\n";
	print PRM "\n";
	
	foreach $atom (atoms($mol)) {
		
		my $list2;
		@$list2 = ();
		undef $mass;
		undef $type;
		
		if (!defined $atom->{OPLS_TYPE}) {
			silico_msg('w', "OPLS-AA specific type missing from atom!\n",
					"Number: $atom->{NUM}\n",
					"Specific type will be set to \"UNK\" in parameter file.\n");
		}
		
		if (defined $atom->{OPLS_TYPE} && !defined $atom->{OPLS_TYPE2}) {
			silico_msg('w', "OPLS-AA generic type missing from atom!\n",
					"Number: $atom->{NUM}  Specific type: \"$atom->{OPLS_TYPE}\"\n",
					"Generic type will be set to \"UNK\" in parameter file.\n");
		}
		
		# Prepare the TINKER_CLASS.
		if (!$class{"$atom->{OPLS_TYPE}"}) {
			++$c;
			$class{"$atom->{OPLS_TYPE}"} = $c;
		}
		$atom->{TINKER_CLASS} = $class{"$atom->{OPLS_TYPE}"};
		
		# Prepare a number for the numerical sort
		$type = $atom->{OPLS_TYPE};

		++$hash->{$type};
		next if ($hash->{$type} > 1);
		
		# Do further operations on the number: remove whitespace
		# and ensure that problem atom types appear first
		$type =~ s/ //g;
		$type = 0 if ($type =~ /\W/);
		
		# Modify atomic masses if required
		$mass = $Silico::Atomic_Masses[$atom->{ELEMENT_NUM}];
		
		if ($atom->{ELEMENT} eq 'H' && $options =~ /\bHEAVYH\b/) {
			++$heavycount;
			$mass = 2.0141;
		}
		
		if (!defined $mass) {
			silico_msg('w', "No mass defined for atom!\n",
					"Number: $atom->{NUM}  Element: $atom->{ELEMENT}\n",
					"Mass will be set to zero in parameter file.\n");
		}
		
		if (!defined $atom->{CONNECT}) {
			silico_msg('w', "Number of connections not defined for atom!\n",
					"Number: $atom->{NUM}\n",
					"Number of connections will be set to zero in parameter file.\n");
		}
		
		@{$list2} = (	$type,
				($atom->{OPLS_TYPE} || 'UNK'),
				($atom->{TINKER_CLASS} || 'UNK'),
				($atom->{OPLS_TYPE2} || 'UNK'),
				($atom->{ELEMENT_NUM} || 0),
				($mass || 0),
				(($#{$atom->{CONNECT}}+1) || 0));
				
		push @$array, $list2;
	}
	
	silico_msg('c', "Used deuterium mass for $heavycount hydrogen atoms\n")
		if ($heavycount && $options !~ /\bQUIET\b/);
	
	# Sort by the modified atom type.
	
	@$array = sort {
		$a->[0] <=> $b->[0];
		
		} @$array;
	
	for ($i = 0; $i <= $#{$array}; ++$i) {
	
		my $atom;
		
		$list = $array->[$i];
	
		#$list->[2] = $i + 1;
		shift @$list;
		
		printf PRM "atom %6s %5s %5s      \"\" %26d %10.3f %5d\n", @$list;
	}
	
	print PRM "\n";
}

sub write_tinker_prm_mol_vdw {

	my $mol = $_[0];
	my $ff = $_[1];
	
	write_tinker_prm_mol_vdw_oplsaa($mol) if ($ff eq 'oplsaa');

}

sub write_tinker_prm_mol_vdw_oplsaa {

	my $mol = $_[0];
	
	my @array;
	my $atom;
	my $hash;
	my $list;
	
	print PRM "\n";
	print PRM "      ################################\n";
	print PRM "      ##                            ##\n";
	print PRM "      ##  Van der Waals Parameters  ##\n";
	print PRM "      ##                            ##\n";
	print PRM "      ################################\n";
	print PRM "\n";
	print PRM "\n";
	
	foreach $atom (atoms($mol)) {
		
		my $list;
		my $type;
		
		if (!defined $atom->{TINKER_CLASS} || !defined $atom->{SIGMA}
			|| !defined $atom->{EPSILON}) {
			silico_msg('w', "VDW parameters missing from atom!\n",
					"Number: $atom->{NUM}  Type: ".($atom->{OPLS_TYPE} || 'UNDEFINED')."\n",
					"Sigma and Epsilon will be set to zero in parameter file.\n");
		}
		
		# Prepare a number for the numerical sort
		$type = $atom->{TINKER_CLASS};

		++$hash->{$type};
		next if ($hash->{$type} > 1);
		
		# Do further operations on the number: remove whitespace
		# and ensure that problem atom types appear first
		$type =~ s/ //g;
		$type = 0 if ($type =~ /\W/);
		
		@{$list} = (	$type,
				($atom->{TINKER_CLASS} || 'UNK'),
				($atom->{SIGMA} || 0),
				($atom->{EPSILON} || 0));
		
		push @array, $list;
	}
	
	# Sort by the modified atom type.
	
	@array = sort {
		$a->[0] <=> $b->[0];
		
		} @array;
		
	foreach $list (@array) {
	
		shift @$list;
		
		printf PRM "vdw %10s %19.4f %10.4f\n", @$list;
	}
		
	print PRM "\n";
}

sub write_tinker_prm_mol_charges {

	my $mol = $_[0];
	my $ff = $_[1];
	
	write_tinker_prm_mol_charges_oplsaa($mol) if $ff eq 'oplsaa';

}

sub write_tinker_prm_mol_charges_oplsaa {

	my $mol = $_[0];
	
	my @array;
	my $atom;
	my $hash;
	my $list;
	
	print PRM "\n";
	print PRM "      ########################################\n";
	print PRM "      ##                                    ##\n";
	print PRM "      ##  Atomic Partial Charge Parameters  ##\n";
	print PRM "      ##                                    ##\n";
	print PRM "      ########################################\n";
	print PRM "\n";
	print PRM "\n";
	
	foreach $atom (atoms($mol)) {
		
		my $list;
		my $type;
		
		if (!defined $atom->{OPLS_TYPE} || !defined $atom->{CHARGE}) {
			silico_msg('w', "Charge missing from atom!\n",
					"Number: $atom->{NUM}  Type: ".($atom->{OPLS_TYPE} || 'UNDEFINED')."\n",
					"Charge will be set to zero in parameter file.\n");
		}
		
		# Prepare a number for the numerical sort
		$type = $atom->{OPLS_TYPE};

		++$hash->{$type};
		next if ($hash->{$type} > 1);
		
		# Do further operations on the number: remove whitespace
		# and ensure that problem atoms appear first
		$type =~ s/ //g;
		$type = 0 if ($type =~ /\W/);
		
		@{$list} = (	$type,
				($atom->{OPLS_TYPE} || 'UNK'),
				($atom->{CHARGE} || 0));
		
		push @array, $list;
	}
	
	# Sort by the modified atom type.
	
	@array = sort {
		$a->[0] <=> $b->[0];
		
		} @array;
	
	foreach $list (@array) {
	
		shift @$list;
		
		printf PRM "charge %7s %13.4f\n", @$list;
	}
	
	print PRM "\n";
}

sub write_tinker_prm_mol_bonds {

	my $mol = $_[0];
	my $ff = $_[1];
	
	write_tinker_prm_mol_bonds_oplsaa($mol) if $ff eq 'oplsaa';

}

sub write_tinker_prm_mol_bonds_oplsaa {

	my $mol = $_[0];
	
	my @array;
	my $bond;
	my $list;
	my %printed;
	my $printedbond;
	my @printedbonds;
	
	print PRM "\n";
	print PRM "      ##################################\n";
	print PRM "      ##                              ##\n";
	print PRM "      ##  Bond Stretching Parameters  ##\n";
	print PRM "      ##                              ##\n";
	print PRM "      ##################################\n";
	print PRM "\n";
	print PRM "\n";
	
	foreach $bond (@{$mol->{BONDS}}) {
	
		my $atom1;
		my $atom2;
		my $atype1;
		my $atype2;
		my $k;
		my $list;
		my $r;

		# Get atoms involved in this bond, and their types
		$atom1 = $mol->{ATOMS}[$bond->[0]];
		$atype1 = $atom1->{TINKER_CLASS};
		
		$atom2 = $mol->{ATOMS}[$bond->[1]];
		$atype2 = $atom2->{TINKER_CLASS};
		
		# Sort to put lowest atom type first in each bond
		if ($atype1 gt $atype2) {
			($atype1, $atype2) = ($atype2, $atype1);
		}
	
		# Provide variables for force constant and reference distance
		$k = $bond->[2];
		$r = $bond->[3];
		
		# Skip and display a warning if there are any parameters missing for this bond
		if (!defined $k || !defined $r) {
		
			silico_msg('w', "Parameters missing for bond!\n",
					"Atoms: $atom1->{NUM}, $atom2->{NUM}\n",
					"Types: $atom1->{OPLS_TYPE} ($atom1->{OPLS_TYPE2}), $atom2->{OPLS_TYPE} ($atom2->{OPLS_TYPE2})\n",
					"No entry created in parameter file.\n");
			next;
		}
		
		# Check whether this bond is already parametrised
		# If it is not, then add it to the list of bonds to print
		# If it is then do nothing (will automatically move on to the next bond)
		foreach $printedbond (@printedbonds) {
			$printed{$printedbond} = 1
		}
		
		unless ($printed{"$atype1-$atype2"}) {
			@{$list} = ($atype1, $atype2, $k, $r);
			push @array, $list;
			push (@printedbonds, "$atype1-$atype2");
		}
	}
	
	# Sort the list of bonds to print: Atom 1, then Atom 2
	@array = sort {
		return  1 if $a->[0] gt $b->[0];
		return -1 if $a->[0] lt $b->[0];
		return  1 if $a->[1] gt $b->[1];
		return -1 if $a->[1] lt $b->[1];
		return 0;
		} @array;
	
	# Print out each bond in the list
	# Prints out the two atoms, then the force constant and reference distance
	foreach $list (@array) {
		
		my $atype1;
		my $atype2;
		my $k;
		my $r;
		
		($atype1, $atype2, $k, $r) = @$list;
		
		# Do sprintf on numeric values.  This lets us handle undefined (UNK) values
		if ($k !~ /unk/i) {
			$k = sprintf "%.1f", $k;
		}
		if ($r !~ /unk/i) {
			$r = sprintf "%.4f", $r;
		}
		
		printf PRM "bond %9s %4s %14s %10s\n", $atype1, $atype2, $k, $r;
	}
	
	print PRM "\n";
}

sub write_tinker_prm_mol_angles {

	my $mol = $_[0];
	my $ff = $_[1];
	
	write_tinker_prm_mol_angles_oplsaa($mol) if $ff eq 'oplsaa';

}

sub write_tinker_prm_mol_angles_oplsaa {

	my $mol = $_[0];
	
	my $angle;
	my @array;
	my $join;
	my $list;
	my %printed;
	my $printedangle;
	my @printedangles;
	
	print PRM "\n";
	print PRM "      ################################\n";
	print PRM "      ##                            ##\n";
	print PRM "      ##  Angle Bending Parameters  ##\n";
	print PRM "      ##                            ##\n";
	print PRM "      ################################\n";
	print PRM "\n";
	print PRM "\n";
	
	# Get OPLS-AA numeric types for each atom in every angle
	foreach $angle (@{$mol->{ANGLES}}) {
	
		my $atom1;
		my $atom2;
		my $atom3;
		my $atype1;
		my $atype2;
		my $atype3;
		my $k;
		my $r;
		
		# Get atoms involved in this angle, and their types
		$atom1 = $mol->{ATOMS}[$angle->[0]];
		$atype1 = $atom1->{TINKER_CLASS};
		
		$atom2 = $mol->{ATOMS}[$angle->[1]];
		$atype2 = $atom2->{TINKER_CLASS};
		
		$atom3 = $mol->{ATOMS}[$angle->[2]];
		$atype3 = $atom3->{TINKER_CLASS};
		
		# Sort to put lowest atom type first in each angle
		if ($atype1 gt $atype3) {
			($atype1, $atype3) = ($atype3, $atype1);
		}
		
		# Skip if there are any atom types missing for this angle
		#next if (!$atom1->{OPLS_TYPE2} || !$atom2->{OPLS_TYPE2} || !$atom3->{OPLS_TYPE2});
		
		# Skip and display a warning if there are any parameters missing for this angle
		if (!defined $angle->[3] || !defined $angle->[4]) {
		
			silico_msg('w', "Parameters missing for angle!\n",
					"Atoms: $atom1->{NUM}, $atom2->{NUM}, $atom3->{NUM}\n",
					"Types: $atom1->{OPLS_TYPE} ($atom1->{OPLS_TYPE2}), $atom2->{OPLS_TYPE} ($atom2->{OPLS_TYPE2}), $atom3->{OPLS_TYPE} ($atom3->{OPLS_TYPE2})\n",
					"No entry created in parameter file.\n");
			next;
		}
		
		# Provide variables for force constant and reference angle
		$k = $angle->[3];
		$r = $angle->[4];
		
		# Check whether this angle is already parametrised
		# If it is not, then add it to the list of angles to print, and mark it as parametrised
		# If it is, then do nothing (will automatically move on to the next angle)
		foreach $printedangle (@printedangles) {
			$printed{$printedangle} = 1
		}
		
		unless ($printed{"$atype1-$atype2-$atype3"}) {
			my $list;
			@{$list} = ($atype1, $atype2, $atype3, $k, $r);
			push @array, $list;
			push (@printedangles, "$atype1-$atype2-$atype3");
		}
	}
	
	# Sort the list of angles to print: Atom 2, then Atom 1, then Atom 3
	@array = sort {
		return  1 if $a->[1] gt $b->[1];
		return -1 if $a->[1] lt $b->[1];
		return  1 if $a->[0] gt $b->[0];
		return -1 if $a->[0] lt $b->[0];
		return  1 if $a->[2] gt $b->[2];
		return -1 if $a->[2] lt $b->[2];
		return 0;
		} @array;
	
	# Print out each angle in the list
	# Prints out the three atoms, then the force constant and reference angle
	foreach $list (@array) {
	
		my $atype1;
		my $atype2;
		my $atype3;
		my $k;
		my $r;
		
		($atype1, $atype2, $atype3, $k, $r) = @$list;
		
		# Do sprintf on numeric values.  This lets us handle undefined (UNK) values
		if ($k !~ /unk/i) {
			$k = sprintf "%.2f", $k;
		}
		
		if ($r !~ /unk/i) {
			$r = sprintf "%.2f", $r;
		}
		
		printf PRM "angle %8s %4s %4s %9s %10s\n", $atype1, $atype2, $atype3, $k, $r;
	}
	
	print PRM "\n";
}

sub write_tinker_prm_mol_dihedrals {

	my $mol = $_[0];
	my $ff = $_[1];
	
	write_tinker_prm_mol_dihedrals_oplsaa($mol) if $ff eq 'oplsaa';

}

sub write_tinker_prm_mol_dihedrals_oplsaa {
	
	#<
	#? Writes out parameters for each dihedral type in the molecule
	#; Requires: Molecule
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	
	my @array;
	my $dihedral;
	my $list;
	my %printed;
	my $printeddih;
	my @printeddihs;
	
	print PRM "\n";
	print PRM "      ############################\n";
	print PRM "      ##                        ##\n";
	print PRM "      ##  Torsional Parameters  ##\n";
	print PRM "      ##                        ##\n";
	print PRM "      ############################\n";
	print PRM "\n";
	print PRM "\n";
	
	# Get OPLS-AA numeric types for each atom in every dihedral
	foreach $dihedral (@{$mol->{DIHEDRALS}}) {
			
		# Get atoms involved in this dihedral, and their types
		my $atom1 = $mol->{ATOMS}[$dihedral->[0]];
		my $atype1 = $atom1->{TINKER_CLASS};
		
		my $atom2 = $mol->{ATOMS}[$dihedral->[1]];
		my $atype2 = $atom2->{TINKER_CLASS};
		
		my $atom3 = $mol->{ATOMS}[$dihedral->[2]];
		my $atype3 = $atom3->{TINKER_CLASS};
		
		my $atom4 = $mol->{ATOMS}[$dihedral->[3]];
		my $atype4 = $atom4->{TINKER_CLASS};
		
		# Sort this dihedral internally by OPLS type (bonded atoms, then terminal atoms)
		if (($atype2 gt $atype3) || ($atype2 eq $atype3 && $atype1 gt $atype4)) {
			($atype1, $atype2, $atype3, $atype4) = ($atype4, $atype3, $atype2, $atype1);
		}
		
		# Provide variables for the four Fourier coefficients
		my $v1 = $dihedral->[4];
		my $v2 = $dihedral->[5];
		my $v3 = $dihedral->[6];
		my $v4 = $dihedral->[7];
		my $type = $dihedral->[8];
		my $atoms = $dihedral->[9];
		my $comment = $dihedral->[10];
		
		# Check whether this dihedral is already parametrised
		# If it is not, then add it to the list of dihedrals to print
		# If it is then do nothing (will automatically move on to the next dihedral)
		foreach $printeddih (@printeddihs) {
			$printed{$printeddih} = 1
		}
		
		unless ($printed{"$atype1-$atype2-$atype3-$atype4"}) {
			my $list;
			@{$list} = ($atype1, $atype2, $atype3, $atype4, $v1, $v2, $v3, $v4, $type, $atoms, $comment);
			push @array, $list;
			push (@printeddihs, "$atype1-$atype2-$atype3-$atype4");
		}
	}
	
	# Sort the list of dihedrals to print: Atom 2, then Atom 3, then Atom 1, then Atom 4
	@array = sort {
		return  1 if $a->[1] gt $b->[1];
		return -1 if $a->[1] lt $b->[1];
		return  1 if $a->[2] gt $b->[2];
		return -1 if $a->[2] lt $b->[2];
		return  1 if $a->[0] gt $b->[0];
		return -1 if $a->[0] lt $b->[0];
		return  1 if $a->[3] gt $b->[3];
		return -1 if $a->[3] lt $b->[3];
		return 0;
		} @array;
	
	# Print each dihedral in the list: Four atom types, then the Fourier coefficient, order, and offset
	# Prints out a single line for each order
	foreach $list (@array) {
		
		my ($atype1, $atype2, $atype3, $atype4, $v1, $v2, $v3, $v4, $type, $atoms, $comment) = @$list;
		
		printf PRM "torsion %6s %4s %4s %4s", $atype1, $atype2, $atype3, $atype4;
		
		if ($v1 != 0) {
			printf PRM "%11.3f 0.0 1", $v1;
		} elsif ($v2 != 0 || $v3 != 0 || $v4 != 0) {
			print PRM  "                 ";
		}
		
		if ($v2 != 0) {
			printf PRM "%9.3f 180.0 2", $v2;
		} elsif ($v3 != 0 || $v4 != 0) {
			print PRM  "                 ";
		}
			
		if ($v3 != 0) {
			printf PRM "%9.3f 0.0 3", $v3;
		} elsif ($v4 != 0) {
			print PRM  "               ";
		}
		
		if ($v4 != 0) {
			printf PRM "%9.3f 180.0 4", $v4;
		}
		
		print PRM "\n";
	}
	
	print PRM "\n";
}

sub write_tinker_prm_mol_impropers {

	my $mol = $_[0];
	my $ff = $_[1];
	
	write_tinker_prm_mol_impropers_oplsaa($mol) if $ff eq 'oplsaa';

}

sub write_tinker_prm_mol_impropers_oplsaa {

	my $mol = $_[0];
	
	my @array;
	my $improper;
	my $list;
	my %printed;
	my $printedimpr;
	my @printedimprs;
	
	print PRM "\n";
	print PRM "      #####################################\n";
	print PRM "      ##                                 ##\n";
	print PRM "      ##  Improper Torsional Parameters  ##\n";
	print PRM "      ##                                 ##\n";
	print PRM "      #####################################\n";
	print PRM "\n";
	print PRM "\n";
	
	# Get OPLS-AA numeric types for each atom in every planar improper
	foreach $improper (@{$mol->{P_IMPROPERS}}) {
			
		
		# Get atoms involved in this dihedral, and their types
		my $atom1 = $mol->{ATOMS}[$improper->[0]];
		my $atype1 = $atom1->{TINKER_CLASS};
		
		my $atom2 = $mol->{ATOMS}[$improper->[1]];
		my $atype2 = $atom2->{TINKER_CLASS};
		
		my $atom3 = $mol->{ATOMS}[$improper->[2]];
		my $atype3 = $atom3->{TINKER_CLASS};
		
		my $atom4 = $mol->{ATOMS}[$improper->[3]];
		my $atype4 = $atom4->{TINKER_CLASS};
				
		# Sort this dihedral internally by OPLS type (bonded atoms, then terminal atoms)
		if (($atype2 gt $atype3) || ($atype2 eq $atype3 && $atype1 gt $atype4)) {
			($atype1, $atype2, $atype3, $atype4) = ($atype4, $atype3, $atype2, $atype1);
		}
		
		# Provide variables for the four Fourier coefficients
		# Because of Silico internal data handling (two impropers for
		# each planar atom), these need to be divided by 2
		my $v1 = $improper->[4]/2;
		my $v2 = $improper->[5]/2;
		my $v3 = $improper->[6]/2;
		my $v4 = $improper->[7]/2;
		my $type = $improper->[8];
		my $atoms = $improper->[9];
		my $comment = $improper->[10];
		
		# Check whether this dihedral is already parametrised
		# If it is not, then add it to the list of dihedrals to print
		# If it is then do nothing (will automatically move on to the next dihedral)
		foreach $printedimpr (@printedimprs) {
			$printed{$printedimpr} = 1
		}
		
		unless ($printed{"$atype1-$atype2-$atype3-$atype4"}) {
			my $list;
			@{$list} = ($atype1, $atype2, $atype3, $atype4, $v1, $v2, $v3, $v4, $type, $atoms, $comment);
			push @array, $list;
			push (@printedimprs, "$atype1-$atype2-$atype3-$atype4");
		}
	}
	
	# Sort the list of dihedrals to print: Atom 2, then Atom 3, then Atom 1, then Atom 4
	@array = sort {
		return  1 if $a->[1] gt $b->[1];
		return -1 if $a->[1] lt $b->[1];
		return  1 if $a->[2] gt $b->[2];
		return -1 if $a->[2] lt $b->[2];
		return  1 if $a->[0] gt $b->[0];
		return -1 if $a->[0] lt $b->[0];
		return  1 if $a->[3] gt $b->[3];
		return -1 if $a->[3] lt $b->[3];
		return 0;
		} @array;
	
	# Print each dihedral in the list: Four atom types, then the Fourier coefficient, order, and offset
	# Prints out a single line for each order
	foreach $list (@array) {
		
		
		my ($atype1, $atype2, $atype3, $atype4, $v1, $v2, $v3, $v4, $type, $atoms, $comment) = @$list;
		
		printf PRM "imptors %6s %4s %4s %4s", $atype1, $atype2, $atype3, $atype4;
		
		if ($v1 != 0) {
			printf PRM "%11.3f 0.0 1", $v1;
		} elsif ($v2 != 0 || $v3 != 0 || $v4 != 0) {
			print PRM  "                 ";
		}
		
		if ($v2 != 0) {
			printf PRM "%9.3f 180.0 2", $v2;
		} elsif ($v3 != 0 || $v4 != 0) {
			print PRM  "                 ";
		}
			
		if ($v3 != 0) {
			printf PRM "%9.3f 0.0 3", $v3;
		} elsif ($v4 != 0) {
			print PRM  "               ";
		}
		
		if ($v4 != 0) {
			printf PRM "%9.3f 180.0 4", $v4;
		}
		
		print PRM "\n";
	}
	
	print PRM "\n";
}

sub write_tinker_prm_oplsaa {

	#<
	#? Write those parts of a TINKER parameter file that
	#  pertain to the OPLS-AA force field
	#; Requires:
	#; Returns:
	#>
	
	print PRM "\n";
	print PRM "forcefield              OPLS-AA\n";
	print PRM "\n";
	print PRM "vdwtype                 LENNARD-JONES\n";
	print PRM "radiusrule              GEOMETRIC\n";
	print PRM "radiustype              SIGMA\n";
	print PRM "radiussize              DIAMETER\n";
	print PRM "epsilonrule             GEOMETRIC\n";
	print PRM "torsionunit             0.5\n";
	print PRM "vdw-14-scale            2.0\n";
	print PRM "chg-14-scale            2.0\n";
	print PRM "dielectric              1.0\n";
	print PRM "\n";
	print PRM "\n";
	print PRM "      #############################\n";
	print PRM "      ##                         ##\n";
	print PRM "      ##  Literature References  ##\n";
	print PRM "      ##                         ##\n";
	print PRM "      #############################\n";
	print PRM "\n";
	print PRM "\n";
	print PRM "W. L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, \"Development\n";
	print PRM "and Testing of the OPLS All-Atom Force Field on Conformational\n";
	print PRM "Energetics and Properties of Organic Liquids\", J. Am. Chem. Soc.,\n";
	print PRM "117, 11225-11236 (1996)\n";
	print PRM "\n";
	print PRM "D. S. Maxwell, J. Tirado-Rives and W. L. Jorgensen, \"A Comprehensive\n";
	print PRM "Study of the Rotational Energy Profiles of Organic Systems by Ab Initio\n";
	print PRM "MO Theory, Forming a Basis for Peptide Torsional Parameters\", J. Comput.\n";
	print PRM "Chem., 16, 984-1010 (1995)\n";
	print PRM "\n";
	print PRM "W. L. Jorgensen and N. A. McDonald, \"Development of an All-Atom Force\n";
	print PRM "Field for Heterocycles. Properties of Liquid Pyridine and Diazenes\",\n";
	print PRM "THEOCHEM-J. Mol. Struct., 424, 145-155 (1998)\n";
	print PRM "\n";
	print PRM "N. A. McDonald and W. L. Jorgensen, \"Development of an All-Atom\n";
	print PRM "Force Field for Heterocycles. Properties of Liquid Pyrrole, Furan,\n";
	print PRM "Diazoles, and Oxazoles\", J. Phys. Chem. B, 102, 8049-8059 (1998)\n";
	print PRM "\n";
	print PRM "R. C. Rizzo and W. L. Jorgensen, \"OPLS All-Atom Model for Amines:\n";
	print PRM "Resolution of the Amine Hydration Problem\", J. Am. Chem. Soc.,\n";
	print PRM "121, 4827-4836 (1999)\n";
	print PRM "\n";
	print PRM "M. L. P. Price, D. Ostrovsky and W. L. Jorgensen, \"Gas-Phase and\n";
	print PRM "Liquid-State Properties of Esters, Nitriles, and Nitro Compounds with\n";
	print PRM "the OPLS-AA Force Field\", J. Comput. Chem., 22, 1340-1352 (2001)\n";
	print PRM "\n";
	print PRM "Most parameters contained in this file are from \"OPLS and OPLS-AA\n";
	print PRM "Parameters for Organic Molecules, Ions, and Nucleic Acids\" as provided\n";
	print PRM "by W. L. Jorgensen, Yale University, July 2004\n";
	print PRM "\n";
}

sub read_sleuth_abstract {

	#<
	#? Read a SLEUTH abstract file
	#; Requires: input file
	#; Returns: sequence record or undef if file open failed otherwise returns 1.
	#>

	my $infile = $_[0];

	my $newseq;
	my $seq = {};

	$newseq->{TYPE} = 'SLEUTH';
	$newseq->{NUMRESGAP} = 0;
	$newseq->{NUMRES} = 0;

	if (!open_file(*IO, $infile)) {
		silico_msg('e', "Can not open file $infile for reading!\n");
		return undef;
	}

	$seq->{NAME} = get_filebase_short($infile);
	
	LOOP: while  (<IO>) {

		# Add amino acids to squence
		if (/Amino Acid Sequence :/) {

			<IO>; <IO>; <IO>;

			my $i = 0;
			while (<IO>) {

				next LOOP if $_ !~ /\w/;

				my @f = split;
				shift @f;
				foreach my $aa (@f) {

					++$i;

					my $res = {} ;

					++$seq->{NUMRESGAP};
					++$seq->{NUMRES};

					# Three and one-letter codes
					$res->{AA1} = aa_one_letter($aa);
					$res->{AA3} =  $aa;

					$res->{NUM} = $i;
				
					# Add to sequence
					push @{$seq->{SEQ}}, $res;
				}
			}
		}
		
		# Add alpha coordinates
		if (/Alpha Carbon Cartesian Coordinates :/) {

			<IO>; <IO>; <IO>;

			my $i = -1;
			while (<IO>) {

				next LOOP if $_ !~ /\w/;

				my @f = split;
				shift @f;
				shift @f;
				
				++$i;
				my $res = $seq->{SEQ}[$i];
				
				($res->{X}, $res->{Y}, $res->{Z}) = @f;
				
				
			}
		}

		# Add secondary structure assignment
		if (/Secondary Structure Assignment :/) {

			my @f = split;

			my $type = "SS_".$f[0];

			#print "type $type\n";

			<IO>; <IO>; <IO>;

			my $i = -1;
			while (<IO>) {

				next LOOP if $_ !~ /\w/;
				
				chomp;
				
				# Stop regular expression falling outside string
				$_ .= '                                                                               ';
				my @f = /^...........(...)..(...)..(...)..(...)..(...)..(...)..(...)..(...)..(...)..(...)/;


				foreach my $aa (@f) {

					++$i;
					
					last if !defined($seq->{SEQ}[$i]);
					
					my $res = $seq->{SEQ}[$i];

					$aa =~ s/ //g;
					$aa = '.' if $aa eq '';

					$res->{$type} = $aa;
				}
			}
		}

		# Add secondary structure assignment
		if (/Buried Surface Area Percentages :/) {

			<IO>; <IO>; <IO>;

			my $i = -1;
			while (<IO>) {

				next LOOP if $_ !~ /\w/;

				++$i;

				my @f = split;

				my $res = $seq->{SEQ}[$i];

				$res->{BURIED_AREA} = $f[2]." ";
				$res->{BURIED_PC} = $f[3]." ";
				$res->{POLAR_COVER_PC} = $f[4]." ";

				# Single digit values
				$res->{BURIED_SINGLE} = int($f[3]/100*9+0.5); # Scaled 0-9
				$res->{POLAR_COVER_SINGLE} = int($f[4]/100*9+0.5); # Scaled 0-9
			}
		}
	}
	
	if (!$seq->{NUMRESGAP}) {
	
		silico_msg('e', "No residues read!\n");
		return undef;
	}

	@{$seq->{PRINT_FIELDS}} = qw(AA1 SS_Kabsch-Sander BURIED_SINGLE POLAR_COVER_SINGLE);

	#seq_printout($seq);
	
	return $seq;
}

sub tinker_space_group {
	
	#<
	#? Convert a Silico space group (numeric code) to Tinker format
	#  based on Hermann-Mauguin symbols.
	#; Requires: Space group in Silico format
	#; Returns: Space group in Tinker format; undef if error
	#>
	
	my $i;
	my $sg_name;
	my @tinker_groups;
	
	space_groups();
	
	@tinker_groups = qw|	P1
				P1(-)
				P21
				Cc
				P21/a
				P21/n
				P21/c
				C2/c
				P212121
				Pna21
				Pn21a
				Cmc21
				Pccn
				Pbcn
				Pbca
				P41
				I41/a
				P4(-)21c
				P4(-)m2
				R3c
				P6(3)/mcm
				Fm3(-)m
				Im3(-)m
			|;
	
	if (defined $Silico::Space_Groups{$_[0]}) {
		$sg_name = $Silico::Space_Groups{$_[0]};
	} elsif (defined $_[0]) {
		silico_msg('e', "Space group $_[0] is not recognised!\n");
		return undef;
	} else {
		silico_msg('e', "Space group is not defined!\n");
		return undef;
	}
	
	# Replace "sub1" with "1"
	$sg_name =~ s/sub1/1/g;
	
	# Replace "sub2", "sub3" etc with "(2)", "(3)", etc
	for ($i = 2; $i <= 6; ++$i) {
		$sg_name =~ s/sub$i/($i)/g;
	}
	
	# Replace "-" with "(-)"
	$sg_name =~ s/-/(-)/g;
	
	# Check to see if there's a match
	foreach (@tinker_groups) {
		return $sg_name if $sg_name eq $_;
	}
	
	silico_msg('e', "Space group $_[0] is not an allowed Tinker space group!\n");
	return undef;
		
}


return 1;
