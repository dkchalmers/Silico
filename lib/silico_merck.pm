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
#! silico_merck.pm
#? Routines to read and write Merck format files for use with MMFF94 in Charmm
#. $Revision: 1.34.2.1.2.2 $
#>

use strict;
package Silico;

##################################################################
#
#	Merck read routines
#
##################################################################

sub read_merck {

	#<
	#? Read in Merck format file
	#; Requires: filename, unused, unused, , options string;
	#; Returns: ensemble
	#>

	my $infile = $_[0];
	my $options = $_[3] || '';
	
	my $at1;
	my $at2;
	my $atomcount;
	my $bo;
	my $bondcount;
	my $header_bonds;
	my $charge_code;
	my @ctable;
	my $i;
	my $line;
	my $molecules;
	my $molcount = 0;
	my @f;
	
	$options = uc $options if $options;

	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
        if (!defined $fr) {
                silico_msg('e', "Can not open file $infile for reading\n");
                return undef;
        }
        my $FH = $fr->{FH};

	silico_msg('c', "Reading Merck format file: $_[0]\n") if $options !~ /\bQUIET\b/;
	
	while (1) {
	
		my $mol;

		#Read Header 1 (ignoring date etc)
		$line = <$FH>;
		last if !$line;
		
		chomp $line;
		if (length($line) > 70) {
		$line = substr($line,0,70);
		}
		# Remove trailing spaces, stars etc
		$line =~ s/ +$//;
		$line =~ s/\*//g;
		$mol->{NAME} = $line;

		# Read header 2
		$line = <$FH>;
		if ($line !~ /^MOL/) {
			silico_msg('e', "File $infile is missing a MOL line! Possibly not a Merck file.\n");
			return undef;
		}
	
		# Atoms and bonds
		$line = <$FH>;
		@f = split (" ",$line);
		$mol->{NUMATOMS} = $f[0];
		$mol->{NUMBONDS} = 0;
		mol_set_source_data($mol, $fr);
		$header_bonds = $f[1];
	
		silico_msg('g', "numatoms: $f[0] numbonds: $f[1]\n");
	
		$atomcount = 0;

		for ($i=0; $i < $mol->{NUMATOMS}; ++$i) {

			$line = <$FH>;
			chomp $line;
			if (!defined $line) {
				silico_msg('e', "Premature end of file!\n",
						"This file will be skipped.\n");
				return undef;
			}
		
			if (length $line < 70) {
				silico_msg('e', "Atom record ".($i+1)." is too short!\n",
						"Line: \"$line\"\n",
						"This file will be skipped.\n");
				return undef;
			}

			my $atom;

			$atom->{X}			= substr($line,0,10);	# X
			$atom->{Y}			= substr($line,11,10);	# Y
			$atom->{Z}			= substr($line,22,10);	# Z
			$atom->{ELEMENT_NUM}		= substr($line,33,5);	# Element number
			$atom->{MMFF_TYPE}		= substr($line,39,2);	# MMFF atom type
			$charge_code			= substr($line,42,1);	# MMFF chg code
			$atom->{NUM}			= substr($line,44,5);	# Atom sequence num
			$atom->{NAME}			= substr($line,50,4);	# Atom name
			$atom->{SUBID}			= substr($line,54,4);	# Residue number
			$atom->{SUBNAME}		= substr($line,58,4);	# Residue name
			$atom->{CHARGE}			= substr($line,62,10);	# Partial charge
			$atom->{SEGID}			= substr($line,76,10);	# SEGID
		
			# Charge code
			$atom->{FORMAL_CHARGE} = merck_charge_code_to_formal($charge_code);

			# Remove leading spaces from SUBNAME
			$atom->{SUBNAME} =~ s/^ *//g;
		
			#Generate translation table between atom sequence
			# number and atom number in file
			$ctable[$atom->{NUM}] = $atomcount;

			# Renumber atoms
			$atom->{NUM} = $atomcount+1;

			silico_msg('g', "Atom number: ".($atom->{NUM} || 'UNDEFINED'),
					"   Atom name: ".($atom->{NAME} || 'UNDEFINED')."\n",
					"X coord: ".($atom->{X} || 'UNDEFINED'),
					"    Y coord: ".($atom->{Y} || 'UNDEFINED'),
					"    Z coord: ".($atom->{Z} || 'UNDEFINED')."\n",
					"Element number: ".($atom->{ELEMENT_NUM} || 'UNDEFINED'),
					"    Formal charge: ".($atom->{CHARGE} || 'UNDEFINED'),
					"    MMFF Type: ".($atom->{MMFF_TYPE} || 'UNDEFINED'),
					"    MMFF Charge Code: ".($atom->{MMFF_CHARGE_CODE} || 'UNDEFINED')."\n",
					"Residue number: ".($atom->{SUBID} || 'UNDEFINED'),
					"    Residue name: ".($atom->{SUBNAME} || 'UNDEFINED'),
					"    Segment ID: ".($atom->{SEGID} || 'UNDEFINED')."\n",
					"\n");
		
			$atom->{ELEMENT} = element_number2symbol($atom->{ELEMENT_NUM});
			
			$mol->{ATOMS}[$i] = $atom;
			++$atomcount;
		}
	
		if ($atomcount != $mol->{NUMATOMS}) {
			silico_msg('e', "Number of atom records does not match file header!\n",
					"Number of atom records: $atomcount; File header: $mol->{NUMATOMS}\n",
					"Skipping.\n");
			return undef;
		}
	
		$bondcount = 0;
		BOND: while (1) {
		
			last BOND if $bondcount == $header_bonds;
			
			$line = <$FH>;
			
			if (!defined $line) {
			
				silico_msg('e', "Premature end of file while reading bonds!\n",
						"Skipping.\n");
				return undef;
			}
	
			chomp $line;
	
			@f = split " ", $line;
			
			# Loop over 5 atom records in each line
			for ($i=0; $i<5; ++$i) {

				++$bondcount;
				last BOND if $bondcount > $header_bonds;
				
				if (!defined $f[$i*3]) {
					silico_msg('w', "Number of bond records does not match file header!\n",
							"Number of bond records: $bondcount; File header: \"$header_bonds\"\n",
							"Skipping bonds.\n");
					last BOND;
				}

				# Translate atom number to atom sequence
				$at1 = $ctable[$f[$i*3]];
				$at2 = $ctable[$f[$i*3+1]];

				# Bond orders
				$bo =  $f[$i*3+2];
			
				bond_create($mol, $at1, $at2, $bo);
			}
		}

		if ($header_bonds != $mol->{NUMBONDS}) {
			silico_msg('w', "Number of bonds does not match file header!\n",
					"Number of bonds: $mol->{NUMBONDS}; file header: \"$header_bonds\"\n");
		}

		silico_msg('c', "$mol->{NUMATOMS} atoms and $mol->{NUMBONDS} bonds have been read.\n");
	
		
		$molecules->[$molcount] = $mol;
		++$molcount;
	}
	
	return $molecules;
}



##################################################################
#
#	Merck write routines
#
##################################################################

sub write_merck {

	#<
	#? Print out an ensemble as separate Merck molecules
	#; Requires: Ensemble (or molecule), filename
	#; Returns: undef if failed
	#>

	my $molecules = ens($_[0]);
	my $file = $_[1];
	my $options = $_[2] || '';
	
	my $count;
	my $success;
	my $mol;
	my $outfile;
	
	$options = uc $options if $options;

	$outfile = get_ofilebase($molecules).".mrk" if !$outfile;
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Writing Merck file: $file");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}

	$count = 0;
	foreach $mol (@$molecules) {
	
		if ($#{$molecules} > 0) {
			++$count;
			$count = sprintf "%03d", $count;
			$outfile = get_filebase($file)."_$count.".get_fileext($file);
		} else {
			$outfile = $file;
		}
	
		$success = open (OUTFILE, ">$outfile");
		if (!defined $success) {
			silico_msg('e', "Can not create or open file $outfile for writing!\n");
			return undef;
		}

		write_merck_mol($mol);
		close OUTIFLE;
	}
}

sub write_merck_mol {

	#<
	#? Print out a single Merck format molecule
	#. Currently crudely attempts to assign Merck formal charge parameters.
	#  This makes the assumption that all hydrogens are present and all
	#  bond valences are correct
	#. Attempts to fix valence of dataive bonds to single and use alternating
	#  double and single bonds for aromatic rings
	#; Requires: Molecule, optional output options string
	#; Returns: 1
	#>
      
	my $mol = $_[0];
	my $options = $_[1] || '';
	
	my $aromflag = 0;
	my $atom;
	my $bo;
	my $bond;
	my $bondorder;
	my %bondlist;
	my $count;
	my @f;
	my $fix;
	my $hour;
	my $i;
	my $j;
	my $linecount;
	my $min;
	my $sortflag;
	my $warn;
	my $yday;
	my $year;
	
	$mol->{NUMBONDS} = molecule_count_numbonds($mol);
	
	# Correct bonding if required
	molecule_check_and_fix_connectivity($mol) if ($options !~ /NOBOND/);
	
	# Ensure that SEGID is defined for all atoms
	# This reduces problems with the sort routine (I hope)
	foreach $atom (atoms($mol)) {
		$atom->{SEGID} ||= 'MOL1';
		
	}

	# If there are any bonds with aromatic bondorder (4)
	# then convert them to kekule representation
	LOOP: foreach $atom (atoms($mol)) {
	
		foreach $bondorder (@{$atom->{BORDERS}}) {
			next if ($bondorder < 4);
			convert_aromatic_bonds_kekule($mol);
			last LOOP;
		}
	}

	# Ensure that molecule is written out with substructures
	# in numerical order.  This _is_ necessary for some
	# Charmm routines (eg parallel shake)
	if ($options !~ /\bNOSORT\b/) {
	
		$i = 0;
		foreach $atom (atoms($mol)) {
	
			if (!$sortflag && $i > 0 && compare_atom_order ($mol->{ATOMS}[$i-1], $atom) > 0 ) {
				silico_msg('n', "Substructures are not in order.");
				$| = 1;
				silico_msg('c', "Sorting substructures...");
				pdb_molecule_sort_atoms($mol);
				silico_msg('c', "done.\n");
				$| = 0;
				last;
			}

			++$i;
		}
	}
	
	# Assign S and P bond orders as required for merck files
	$fix = fix_dative_bonds ($mol, 1);
	silico_msg('c', "Modified bond order in $fix dative bonds\n");
	
	# Merck Header 1
	# --------------
	# Get last two digits of year
	(undef,$min,$hour,undef,undef,$year,undef,$yday,undef) = gmtime(time);
	$year +=1900;
	$year = substr ($year,2);
	
	$mol->{NAME} ||= 'MOL';
	# There may be a bug in perl with regards to the hour field?
	printf OUTFILE "%-70s%02d%03d%02d%02d1\n", $mol->{NAME}, $year, $yday, $hour, $min;
	
	# Merck Header 2
	print OUTFILE "MOL silico\n";
	
	# Number of atoms and bonds
	printf OUTFILE "%5d %5d\n", $mol->{NUMATOMS}, $mol->{NUMBONDS};
	
	$warn = 0;
	# Atoms
	for ($i=0; $i < $mol->{NUMATOMS}; ++$i) {

		$atom = $mol->{ATOMS}[$i];
		
		# Calculate merck formal charge code if not defined
		if (!defined $atom->{FORMAL_CHARGE}) {
			$atom->{MMFF_CHARGE_CODE} = merck_calculate_formal_charge($atom);
			return undef if !defined $atom->{MMFF_CHARGE_CODE};
		} else {
			$atom->{MMFF_CHARGE_CODE} = merck_formal_to_charge_code($atom->{FORMAL_CHARGE});
		};
		
		$f[0]	= $atom->{X};			# X
		$f[1]	= $atom->{Y};			# Y
		$f[2]	= $atom->{Z};			# Z
		$f[3]	= $atom->{ELEMENT_NUM};		# Element number (0 for lone pair)
		$f[4]	= $atom->{MMFF_TYPE} || 0;	# MMFF atom type
		$f[5]	= $atom->{MMFF_CHARGE_CODE};	# MMFF charge code
		$f[6]	= $atom->{NUM};			# Atom sequence number
		$f[7]	= $atom->{NAME};			# Atom name
		$f[8]	= $atom->{SUBID};		# Residue number
		$f[9]	= $atom->{SUBNAME};		# Residue name
		$f[10]	= $atom->{CHARGE} || 0;		# Partial charge
		$f[11]	= $atom->{SEGID} || 'MOL0';	# SEGID
		
		# Remove spaces
		$f[7] =~ s/ //g;
		$f[8] =~ s/ //g;
		$f[9] =~ s/ //g;
		
		# Truncate fields that are too long
		if (length $f[7] > 4) {
			$f[7] = substr($f[7],0,4);
			$warn = 1;
		}
		if (length $f[8] > 4) {
			$f[8] = substr($f[8],0,4);
			$warn = 1;
		}
		if (length $f[9] > 4) {
			$f[9] = substr($f[9],0,4);
			$warn = 1;
		}
		
		printf OUTFILE "%10.4f %10.4f %10.4f %5d %2d %1d %5d %4s%4s%-4s %8.4f     %-4s\n", @f;
			
	}
	
	$warn && silico_msg('w', "Some fields have been truncated. Please check your output file.\n");
	
	# Hash to keep track of bonds we have already printed
	undef %bondlist;
	
	$count = 0;	# Bond counter
	$linecount = 0;	# Counter to print out five bonds per row
	for ($i=0; $i <$mol->{NUMATOMS}; ++$i) {
			
		$atom = $mol->{ATOMS}[$i];
		$j=0;
		foreach $bond (@{$atom->{CONNECT}}) {
			
			if (!defined $bondlist{"$i $bond"}) {
				++$count;
					
				# Fix up the bond orders
				$bo = ${$atom->{BORDERS}}[$j];
					
				# Skip ir zero order or no bonds defined
				next if (!defined $bo);
				next if ($bo == 0);
			
				$bondorder = $bo;
				# Note need to fix aromatic bondorders in here!
				$bondorder = 1 if ($bo == 4);
				$bondorder = 1 if ($bo == 5 || $bo == 6 || $bo == 8);
				$bondorder = 2 if ($bo == 7);
				# Note we are ignoring single order bonds here
					
				printf OUTFILE "%5d %5d %2d  ",
					($i+1),($mol->{ATOMS}[$bond]{NUM}),$bondorder;
					
				# Split in to 5 bonds per line or last carriage return
				# Note!! Charmm requires that there are no blank lines
				# following the bond list
				++$linecount;
				if ($linecount == 5 || $count == $mol->{NUMBONDS}) {
					print OUTFILE "\n";
					$linecount = 0;
				}
				
				++$bondlist{"$bond $i"};
			}
			++$j;
		}
	}
	
	if ($count != $mol->{NUMBONDS}) {
		silico_msg('w', "Number of bonds written doesn't match \$mol->{NUMBONDS}!\n",
				"Bonds written: $count; Number of bonds: $mol->{NUMBONDS}\n");
	}
	
	return 1;
}

sub merck_charge_code_to_formal {

	#<
	#? Convert a Merck charge code to a formal charge
	#>

	my $code = $_[0];

	return 0 if $code == 0;
	return 1 if $code == 1;
	return -1 if $code == 2;
	return 0 if $code == 3; # Radical
	return 2 if $code == 4;
	return -2 if $code == 5;
	return 3 if $code == 6;
	return -3 if $code == 7;
	return 4 if $code == 8;
	return -4 if $code == 9;

	silico_msg('w', "Unable to assign charge to charge code $code\n");

	return 0;
}

sub merck_formal_to_charge_code {

	#<
	#? Convert a formal charge to a Merck charge code for .mrk files
	#>

	my $fcharge = $_[0];

	return 0 if $fcharge == 0;
	return 1 if $fcharge == 1;
	return 2 if $fcharge == -1;
	return 4 if $fcharge == 2;
	return 5 if $fcharge == -2;
	return 6 if $fcharge == 3;
	return 7 if $fcharge == -3;
	return 8 if $fcharge == 4;
	return 9 if $fcharge == -4;

	silico_msg('w', "Unable to assign Merck charge code for charge $fcharge\n");

	return 0;
}

sub merck_calculate_formal_charge {

	#<
	#? Calculate the merck charge code for .mrk files
	#. The valid charge codes are (From the charmm documentation):
	#
	#;       Code            Charge Code
	#
	#;         0              Neutral
	#;         1               +1
	#;         2               -1
	#;         3              Radical
	#;         4               +2
	#;         5               -2
	#;         6               +3
	#;         7               -3
	#;         8               +4
	#;         9               -4
	#
	#
	#; Requires: atom
	#; Returns: charge code
	#>

	my $atom = $_[0];
	
	my $valence;
	my $el;
	
	# Calculate valence
	
	$el = $atom->{ELEMENT};
	$valence = 0; # No bonds to start with
	foreach (@{$atom->{BORDERS}}) {
	
		if ($_ == 1 || $_ == 2 || $_ == 3) {
			$valence += $_;
			next;
		}
		if ($_ == 4) {
			$valence += 1.5;	# Aromatic
			next;
		}
		
		silico_msg('w', "Bond orders of atom $atom->{NUM} (name \"$atom->{NAME}\" are poorly defined.\n");
	}
	
	if ($el eq "C") {
	
		return 0 if $valence == 4; # Neutral
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 4    Actually is: $valence\n");
		return 0;
	}
	
	if ($el eq 'H') {
	
		return 0 if $valence == 1; # Neutral
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 1    Actually is: $valence\n");
		return 0;
	}
	if ($el eq "N") {
	
		return 0 if $valence == 3; # Neutral
		return 1 if $valence == 4; # +1
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 3 or 4    Actually is: $valence\n");
		return 0;
	}
	if ($el eq "O") {
	
		return 0 if $valence == 2; # Neutral
		return 2 if $valence == 1; # -1
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 2 or 1    Actually is: $valence\n");
		return 0;
	}
	if ($el eq "S") {
	
		return 0 if $valence == 2; # Neutral
		return 0 if $valence == 4; # Neutral
		return 2 if $valence == 1; # -1
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 2, 4 or 1    Actually is: $valence\n");
		return 0;
	}
	if ($el eq "P") {
	
		return 0 if $valence == 2; # Neutral
		return 0 if $valence == 4; # Neutral
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 2 or 4    Actually is: $valence\n");
		return 0;
	}
	if ($el eq "F" || $el eq "Cl" || $el eq "Br" || $el eq "I") {
		
		return 2 if $valence == 0; # Ion
		return 0 if $valence == 1; # Monovalent
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 0 or 1    Actually is: $valence\n");
		return 0;
	}
	if ($el eq "Na" || $el eq "K") {
	
		return 1 if $valence == 0; # +1 ion
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 0    Actually is: $valence\n");
		return 0;
	}
	
	# Divalent cations (some assumptions made here....)
	if ($el eq "Mg" || $el eq "Ca" || $el eq "Mn" ||$el eq "Zn") {
	
		return 4 if $valence == 0; # +2 ion
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 2    Actually is: $valence\n");
		return 0;
	}
	
	# Trivalent cations (some assumptions made here....)
	if ($el eq "Fe") {
	
		return 6 if $valence == 0; # +3 ion
		silico_msg('w', "Atom $atom->{NUM} (name \"$atom->{NAME}\", residue \"$atom->{SUBNAME}\", element \"$el\") has an incorrect valence!\n",
				"Should be: 0    Actually is: $valence\n");
		return 0;
	}
	
	silico_msg('w', "Element of atom $atom->{NUM} (name \"$atom->{NAME}\") is not known!\n",
			"Element set to: \"$atom->{ELEMENT}\".");
	return undef;
	
}

return 1;
