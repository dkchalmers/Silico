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
#! silico_pdb.pm
#? Silico routines to read and write pdb format files.
#. $Revision: 1.15-57-gb152f62 $
#>

use strict;
package Silico;

##################################################################
#
#	PDB read routines
#
##################################################################

sub read_pdb {

	#<
	#? Read in a pdb file.
	#. Reads pdb files containing multiple MODEL records.  Each is 
	#  read into a separate molecule file
	#. Similarly, reads pdb files with multiple molecules separated by END 
	#  records (as produced by VMD)
	#. Reads Autodock output files (this is now old and has not been checked
	#  recently)
	#. The $mol->{HAS_BAD_BONDORDERS} is set 
	#. Saves some information from the PDB header records.  This
	#  is  these are:
	#, HEADER
	#, COMPND
	#, HETNAM
	#, USER (Autodock energies)
	#, REMARK (PDB_VERSION, PDB_RESOLUTION, PDB_EXPT_TYPE)
	#, JRNL (PDB_JRNL_AUTH, PDB_JRNL_TITL, PDB_JRNL_REF)
	#. Reads CRYST records and saves in Silico format
	#. Assumes that atom names starting with Q are pseudo atoms.  These
	#  atoms are set to dummy (Du) type.
	#.  The PDB file is read in two passes.  First to get all non atom information
	#  (header info and CONECT records) then to read atoms.  Header and
	#  CONECT record information is read until an END record is encountered.
	#  It is then copied into all subsequent molecules.  
	#  This lets files containing MODELS be read.
	#; Requires: Filename, start structure (numbered from 0), max_molecules, option string
	#. Options:
	#,  QUIET - do not print 'Reading' line or warnings
	#,  DELALT - delete alternate atoms
	#; Returns: Ensemble containing molecule or undef if file open fails
	#>

	my $infile = $_[0];
	my $start = $_[1] || get_flag('ss', 's') || 0; 
        my $max_molecules = $_[2] || get_sflag('ms');
	my $options = $_[3] || '';
	
	my $mols;
	my $starttime = (times)[0];
	
	# Options
	$options = uc $options;
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	my $quiet = 1 if ($options =~ /\bQUIET\b/);
	
	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
	return undef if !$fr;
	
	if (!$quiet) {
		if ($options) {
			silico_msg('c', "Reading pdb file: $infile using options $options\n");
		} else {
			silico_msg('c', "Reading pdb file: $infile\n");
		}
	}
	
	my $i = -1;
	my $molcount = 0;

	while (my $mol = read_pdb_molecule($fr, $options)) {
		++$i;
		next if $start and ($i < $start);
		++$molcount;
		push @$mols, $mol;
		last if (defined $max_molecules) && $molcount == $max_molecules;
	}
	
	if ($#{$mols} > 0) {
		silico_msg('n', "Found ".($#{$mols}+1)." models.\n",
				"\tThese have been read as separate molecules.\n");
	}
	
	close_molfile($fr);

	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
	return $mols;
}

sub read_pdb_single {

	#<
	#? Read in a single record from a pdb file.
	#; Requires: file record, options
	#. See read_pdb for general description
	#. Options:
	#; Returns: molecule or undef
	#>

	my $fr = $_[0];
	my $options = $_[1] || '';

	$options = uc $options if $options;

	my $mol = read_pdb_molecule($fr, $options);

	return undef if !defined $mol;

	# Print warnings
	mol_print_warnings($mol);
	
	return $mol;
}

sub pdb_rewind {

	#<
	#? Rewind a pdb file to the start
	#; Requires: file record
	#; Requires: nothing
	#>

	my $fr = $_[0];

	my $FH = $fr->{FH};

	# Close and reopen file to rewind (Seek does not work on pipes used to open compressed files);
	$FH->close;
	$FH->open($fr->{OPEN_COMMAND});
}

sub read_pdb_molecule {

	#<
	#? Read a single PDB molecule up to and END or ENDMDL record
	#; Requires: file record, options
	#; Returns: molecule, undef if end of file
	#>
	
	my $fr = $_[0];
	my $options = $_[1];

	my $FH = $fr->{FH};
	my $hmol;
	my $mol;
	my $tflag = 1;
	my $reslen = 3;

	my $quiet = 1 if $options && $options =~ /\bQUIET\b/;
	$reslen = 4 if get_lflag('pdb-long-resname')  || $options && $options =~ /\bLONG_RESNAME\b/;
	
	my $stride = get_sflag('stride');
	if ($stride) {
		while (1) {
			last if $stride == 1;
	  		$_ = <$FH>;
               		if (!defined $_) {
			    $fr->{END} = 1;
        	   	    last;
               		}
			if (/^END\b/ || /^ENDMDL\b/) {
				--$stride;
				++$fr->{MOL_READ_COUNT};
               		}
		}
	}

	# Initialise hash variables (doesn't work if we don't do this)
	%$mol = (); # Molecule with coordinate data
	
	# Molecule with header data so we can read multi molecule records
	if (defined $fr->{PDB_HEADERS}) {
		$hmol = $fr->{PDB_HEADERS};
	} else {
		%$hmol = (); 
		$fr->{PDB_HEADERS} = $hmol;
	} 
	      
       	while (1) {

               	$_ = <$FH>;
               	if (!defined $_) {
		       $fr->{END} = 1;
        	       last;
               	}
		
		# Use this to remove Windows line breaks
		$_ =~ s/\R//;
               	#chomp;

		if (/^ATOM/ || /^HETATM/) {
        	       read_pdb_atom($_, $mol, $quiet, $reslen);
		}  elsif (/^CONECT/) {
			$fr->{PDB_CONECT} = read_pdb_conect($_, $mol);
                } elsif (/^HEADER/) {
                        read_pdb_header_line($_, $hmol);	
                } elsif (/^TITLE/) {
                        read_pdb_title($_,$hmol);
                } elsif (/^COMPND/) {
                        read_pdb_compnd($_, $hmol);
                } elsif (/^CRYST1/) {
                        read_pdb_cryst($_, $hmol);
                } elsif (/^HETNAM/) {
                        read_pdb_hetnam($_,$hmol);
                } elsif (/^REMARK/) {
                        read_pdb_remark($_, $hmol);
                } elsif (/^JRNL/) {
                        read_pdb_jrnl($_, $hmol);
                } elsif (/^EXPDTA/) {
                        read_pdb_expdta($_, $hmol);
                } elsif (/^USER/) {
                        read_pdb_user($_, $hmol);
		} elsif (/^TER\b/) {
               
        	       # Revised mechanism for TER records
        	       $mol->{ATOMS}[-1]{TER} = 1 if defined $mol->{ATOMS}[-1];

               	} elsif (/^END\b/) {
        	       	last;
               	} elsif (/^ENDMDL\b/) {
        	       	last;
               	}
	}

	# End of file
	if (!$mol->{NUMATOMS}) {
		return undef;
	}
	
	# Create bonds from CONECT records and CTABLE
	#  Pdb files have two records for each bond (ie atom 1 bonds to 3 and atom 3 bonds to 1).
	#  Some programs (eg Schrodinger) encode double bonds by entering a bond twice.

	foreach my $a (@{$mol->{PDB_CONECT}}) {
	
		# Translate atom number b to atom sequence
		my $anuma = $mol->{PDB_CTABLE}[${$a}[0]];
		my $anumb = $mol->{PDB_CTABLE}[${$a}[1]];
	
		if (!defined $anuma) {
			#print "anuma $anuma $anumb a0  ${$a}[0] a1 ${$a}[1]\n";
			silico_msg('w', "Atom number ${$a}[0] is present in CONECT records but is missing from ATOM records! Skipping.\n");
			next;
		}
		if (!defined $anumb) {
			silico_msg('w', "Atom number ${$a}[1] is present in CONECT records but is missing from ATOM records! Skipping.\n");
			next;
		}

		my $atoma = $mol->{ATOMS}[$anuma];
		
		# Make (half) bond
		my $flag  = 0;
		my $i = 0;
		
		# Check to see if bond already exists and if so increment bondorder
		foreach my $connum (@{$atoma->{CONNECT}}) {
			if ($connum == $anumb) {
				$flag = 1;
				mol_warn($mol, 'Bondorder modified');
				++$atoma->{BORDERS}[$i];
				last;
			}
			++$i;
		}
		
		# Bond does not already exist.  Make a new single bond
		if (!$flag) {
		
			# Make (half) bond
			push @{$atoma->{CONNECT}}, $anumb;
			push @{$atoma->{BORDERS}}, 1;
		}
	}
	
	delete ($mol->{PDB_CONECT});
	delete ($mol->{PDB_CTABLE});
	
	# Copy header data to molecule
	foreach my $key (keys %$hmol) {
		$mol->{$key} = $hmol->{$key};
	}
	
	$mol->{NUMBONDS} = molecule_count_numbonds($mol);
	$mol->{SOURCE_FILE_TYPE} = 'pdb';
	$mol->{SOURCE_FILE_NAME} = $fr->{FILENAME};
	$mol->{HAS_BAD_BONDORDERS} = 1;
	$mol->{HAS_BAD_BONDS} = 1 if ($mol->{NUMBONDS} < ($mol->{NUMATOMS}-1)/2);
	$mol->{NUM} = $fr->{MOL_READ_COUNT}+1;
	
	name_pdb_molecule($mol);

	return $mol;
}

sub read_pdb_atom {

	#<
	#? Read in an ATOM line from a pdb  file.
	#; Requires; line, molecule, quiet flag, residue length (default 3)
	#; Returns; new atom
	#>

	my $line = $_[0];
	my $mol = $_[1];
	my $quiet = $_[2];
	my $reslen = $_[3] || 3;

	my $atom;

	++$mol->{NUMATOMS};

	# Pad out line to prevent 'substring outside string' errors
	$line .= '                              ';

	# The original (FORTRAN) format for pdb files is:
	# ( A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2 )
	
	my ($aname, $subid, $c, $el);
	my @f = ($atom->{LABEL}, 
		$atom->{NUM}, 
		$aname, 
		$atom->{ALT}, 
		$atom->{SUBNAME}, 
		$atom->{CHAIN}, 
		$subid,
		$atom->{EXT},
		$atom->{X},
		$atom->{Y},
		$atom->{Z},
		$atom->{OCC},
		$atom->{TEMP},
		$atom->{SEGID},
		$el,
		$c
	) = $line =~ 
		#  Label   Anum    AnameAltSubnm Chn Subid Ext    X         Y         Z         Occ     Temp    Null Segid El  Chg
		/^(......)(.....).(....)(.)(....)(.)(....)(.)...(........)(........)(........)(......)(......)......(....)(..)(..)/;
		
		#print "line $line\n";
		#foreach (@f) {
		#	if (defined $_) {
		#		print "'$_'  "; 
		#	} else {
		#		print "''  ";
		#	}
		#}
		#print "\n\n";
		#die;

	#$atom->{LABEL}		= substr($line,0,6);		#Record type
	#$atom->{NUM}		= substr($line,6,5);		#Atom number 
	#my $aname		= substr($line,12,4);		#Atom name
	#$atom->{ALT}		= substr($line,16,1);		#Alt location indicator
	#$atom->{SUBNAME}	= substr($line,17,$reslen);	#Residue name
	#$atom->{CHAIN}		= substr($line,21,1);		#Chain name
	#my $subid		= substr($line,22,4);		#Residue number
	#$atom->{EXT}		= substr($line,26,1);		#Residue number extension
	#$atom->{X}		= substr($line,30,8);		#x
	#$atom->{Y}		= substr($line,38,8);		#y
	#$atom->{Z}		= substr($line,46,8);		#z
	#$atom->{OCC}		= substr($line,54,6);		#Occupancy
	#$atom->{TEMP}		= substr($line,60,6);		#Temp Factor
	#$atom->{SEGID}		= substr($line,72,4);		#SegID
	#$atom->{ELEMENT}	= substr($line,76,2);		#Element
	#my $c			= substr($line,78,2);		#Formal charge

	# Remove spaces from atom name
	#$aname =~ s/(\b) +(\b)/_/;
	$aname =~ s/ //g;
	$atom->{NAME} = $aname;
	
	# Sometimes even very important fields are empty
	$subid = 1 if $subid eq '    ';
	# Prune the atom substructure ID if it is not a number.
	# That is, replace all non-numbers with the empty string.
	$subid =~ s/\D//g;
	$atom->{SUBID} = $subid;
	
	# Formal charge (could be in the PDB format '1-' or a normal number '-1')
	if ( $c =~ /^([ +-]\d|\d[ +-])$/ ) {
		
		$c =~ s/\+//;
		if ($c =~ '-') {
			$c =~ s/\-//;
			$c = -$c;
		}
		
		#print "c $c\n";
		
		$atom->{FORMAL_CHARGE} = $c;
		mol_warn("Atom with large formal charge ($c)") if abs($c) > 3;
	}

	# Others are lacking OCC, TEMP, SEGID or ELEMENT fields
	# Undefine them if they are empty
	undef $atom->{OCC}     if $atom->{OCC} eq '      ';
	undef $atom->{TEMP}    if $atom->{TEMP} eq '      ';
	undef $atom->{SEGID}   if $atom->{SEGID} eq '    ';
	
	# VMD can use hexadecimal atom numbers in files with more than 99999 atoms (!)
	# Use our own numbers after this
	if ($mol->{NUMATOMS} > 99999) {
		$atom->{NUM} = $mol->{NUMATOMS};
	}
	
	#
	# Work out the element
	#
	
	if ($el) {
	
		$el =~ s/ //g;	
		$el = ucfirst(lc($el));
	
		my $n = $Silico::Atomic_Elements{$el};

		if ($n) {
			$atom->{ELEMENT} = $el;
			$atom->{ELEMENT_NUM} = $n;
		
		} else {
			atom_guess_element($atom, $mol, $quiet);
		}
	
		#print "a $atom->{ELEMENT} $atom->{ELEMENT_NUM}\n";
	}

	# Make translation table between atom index and pdb file atom number
	if (defined $mol->{PDB_CTABLE}[$atom->{NUM}]) {
		mol_warn($mol, "Found duplicate atom");
	}
	$mol->{PDB_CTABLE}[$atom->{NUM}] =  $mol->{NUMATOMS}-1;
	
	# Renumber atoms
	$atom->{NUM} = $mol->{NUMATOMS};

	push @{$mol->{ATOMS}}, $atom;
}

sub read_pdb_compnd {

	#<
	#? Read PDB compnd line
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
	
	$line =~ s/^COMPND\s*//;
	return $mol->{COMPND} .= $line;
}

sub read_pdb_conect {

	#<
	#? Read in a CONECT line from a pdb  file and save it in mol->{PDB_CONECT}
	#; Requires; line, molecule
	#; Returns; mol->{PDB_CONECT}
	#>
	
	my $line = $_[0];
	my $mol = $_[1];
	
	$line .= '                                              ';
	my @f = $line =~ /(......)(.....)(.....)(.....)(.....)(.....)(.....)/;
	
	shift @f;
	
	# Remove spaces
	foreach (@f) {
		next if !$_;
		s/ //g;
	}
	
	if (!$f[0]) {
		silico_msg('w', "Error in CONECT record. Skiping\n");
		return undef;
	}
	
	# Get parent atom
	my $ca = shift @f;
	
	# Save bonds in an array;
	foreach my $connum (@f) {
		next if !$connum;
		my $a;
		@$a = ($ca, $connum);
		push @{$mol->{PDB_CONECT}}, $a;
	}

	return $mol->{PDB_CONECT};
}
	
	
sub read_pdb_cryst {

	#<
	#? Read PDB CRYST line
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
	
	my $ncperuc;
	my $space_group;
	my $flag;
	
	$line .= '                                                                                            ';
	my @f = $line =~ /(......)(.........)(.........)(.........)(.......)(.......)(.......)(............)(....)/;
	
	$mol->{CELL}[0] = $f[1];
	$mol->{CELL}[1] = $f[2];
	$mol->{CELL}[2] = $f[3];
	
	# Set next three entries to unit cell angles
	$mol->{CELL_ALPHA} = $f[4];
	$mol->{CELL_BETA} = $f[5];
	$mol->{CELL_GAMMA} = $f[6];

	foreach (@{$mol->{CELL}}, $mol->{CELL_ALPHA}, $mol->{CELL_BETA}, $mol->{CELL_GAMMA}) {
		next if $_ != 0;
		delete $mol->{CELL};
		delete $mol->{CELL_ALPHA};
		delete $mol->{CELL_BETA};
		delete $mol->{CELL_GAMMA};
		silico_msg('e', "Molecule cell data contains zero value! Ignoring cell data\n");
		return;
	}
	
	foreach (@f) {
		s/^\s*//;
		s/\s*$//;
	}

	# Try and match space group
	# Space group must contain at least two elements, space separated.
	# Number of chains in unit cell provides an additional element.
	# Thus, there must be at least 10 space-separated values on this line.
	# If not, we fall back to default values.
	if ($#f < 8) {
		silico_msg('w', "Badly formed PDB CRYST1 record.\n",
			"Space group information, if present, must include a lattice centering,\n",
			"at least one symmetry element, and the number of main chains per unit cell.\n",
			"Continuing, using default values for space group and number of chains.\n");
	 
	} else {
		
		$space_group = $f[7];
		$mol->{SPACE_GROUP} = space_group_pdb_to_silico($space_group);
		
		$ncperuc = $f[8]; # Number of chains per unit cell is last entry on line.
		if ($ncperuc eq '' || !check_data_type($ncperuc, 'INT > 0')) {
		
			$mol->{CHAINS_PER_CELL} = 1;
			
		} else {
		
			$mol->{CHAINS_PER_CELL} = $ncperuc;
		}
	}

	# Set space group and chains per unit cell to default values.
	# That is, "P1" (Silico code 0011) for space group, and 1 for number of chains
	# per unit cell.
	# Only do this if they are otherwise not defined.
	$mol->{SPACE_GROUP} ||= "0011";
	$mol->{CHAINS_PER_CELL} ||= 1;
}

sub read_pdb_header_line {

	#<
	#? Read PDB HEADER line
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];

	my $date;
	my $code;
	
	$line .= "                                                                                     ";
	$mol->{HEADER} = substr($line, 7, 40); 
	$mol->{HEADER} =~ s/\s*$//;
	$mol->{HEADER} =~ s/^\s*//;
	$date = substr($line,50,9);
	$date =~ s/ //g;
	$mol->{SDF_DATA}{PDB_DATE} = $date if $date ne '';
	$code = substr($line,62,4);
	$code =~ s/ //g;
	$mol->{SDF_DATA}{PDB_CODE} = $code if $code ne '';
	
	return $mol;
}

sub read_pdb_hetnam {

	#<
	#? Read PDB HETNAM line
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
	
	clean_pdb_line($line, 'HETNAM');
	
	$mol->{HETNAM} .= " " if $mol->{HETNAM};
	$mol->{HETNAM} .= $line;
	$mol->{SDF_DATA}{PDB_HETNAM} = $mol->{HETNAM};
}

sub read_pdb_expdta {

	#<
	#? Read PDB EXPDTA line
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
	
	clean_pdb_line($line, 'EXPDTA');
	
	$mol->{EXPDTA} .= $line;
	$mol->{SDF_DATA}{PDB_EXPDTA} .= $line;
}

sub read_pdb_jrnl {

	#<
	#? Read PDB JRNL record
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
	
	 if (/^JRNL\s*AUTH/) {
	 	 $line =~ s/^JRNL\s*AUTH...//;
	 	 $line =~ s/\s*$//;
	 	 $mol->{SDF_DATA}{PDB_JRNL_AUTH} .= $line; 
	 }
	 if (/^JRNL\s*TITL/) {
	 	 $line =~ s/^JRNL\s*TITL...//;
	 	 $line =~ s/\s*$/ /;
	 	 $mol->{SDF_DATA}{PDB_JRNL_TITL} .= $line; 
	 }
	 if (/^JRNL\s*REF /) {
	 	 $line =~ s/^JRNL\s*REF....//;
	 	 $line =~ s/\s*$//;
	 	 $line =~ s/\s+/ /g;
	 	 $mol->{SDF_DATA}{PDB_JRNL_REF} .= $line; 
	 }
}

sub read_pdb_remark {

	#<
	#? Read PDB REMARK record
	#; Extracts: PDB file version, Resolution and experiment type
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];

	# PDB Format
	if ($line =~ /\s*4.*COMPLIES WITH FORMAT V\.\s*(\S*)/) {
		$mol->{SDF_DATA}{PDB_VERSION} = $1;
	}
	#  Resolution
	if ($line =~ /\s*2\s*RESOLUTION.\s*(\S*)\s*ANGSTROMS/) {
		$mol->{SDF_DATA}{PDB_RESOLUTION} = $1;
	}
	# Autodock Vina
	if ($line =~ /VINA RESULT:/) {
		my @f = split ' ', $line;
		$mol->{SDF_DATA}{ENERGY} = $f[3];
		$mol->{SDF_DATA}{VINA_RESULT} = $f[3];
	}
	# Experiment type
	if ($line =~ /EXPERIMENT TYPE\s*:\s*(.*)\s*/) {
		$mol->{SDF_DATA}{PDB_EXPT_TYPE} = $1;
		clean_pdb_line($mol->{SDF_DATA}{PDB_EXPT_TYPE});
	}
}

sub read_pdb_title {

	#<
	#? Read PDB TITLE record
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
	my $tflag = $_[2];
	
	clean_pdb_line($line);

	$mol->{TITLE} = '' if $tflag;
	$tflag = 0;

	$mol->{TITLE} .= ' ' if $mol->{TITLE};
	
	$mol->{TITLE} .= $line;
	$mol->{SDF_DATA}{PDB_TITLE} = $mol->{TITLE};
}


sub read_pdb_user {

	#<
	#? Read USER record from pdb file
	#; Requires; line, molecule
	#>

	my $line = $_[0];
	my $mol = $_[1];
					
	my @f = split (/\s+/,$line);
	
	if (/Final Docked Energy/) {
		$mol->{ADOCK_BINDING} = $f[5] if defined $f[5];
		# Also assign this value to ENERGY
		$mol->{ENERGY} = $f[5] if defined $f[5];
		$mol->{ENERGY_TYPE} = "Autodock";
	}
	if (/Free Energy of Binding/) {
		$mol->{ADOCK_FREE} = $f[7] if defined $f[7];
	}
}

sub name_pdb_molecule {

	#<
	#? Generate molecule name for pdb molecule
	#. 	1. HEADER
	# 	2. COMPND
	# 	or 3. filename
	#>
	
	my $mol = $_[0];
	
	my $name = $mol->{HEADER} if (defined $mol->{HEADER});
		
	if (!$name && defined $mol->{COMPND}) {
		if ($mol->{COMPND} =~ 'MOL_ID') {
			my @f = split "COMPND", $mol->{COMPND}; 
			foreach (@f) {
				next if ! /MOLECULE:/;
				
				s/.*MOLECULE:\s*//;
				s/;\s*$//;
				
				print "$_\n";
				$name .= "$_ ";
				print "name $name\n";
			}
		} else {
			$name = $mol->{COMPND};
		}
		$name = substr("$name                                       ",0,30);
		$name =~ s/ *$//; # Remove trailing spaces
	}
	
	if (!$name) {
		$name = $mol->{SOURCE_FILE_NAME};
		# Remove 3-letter extensions
		$name =~ s/\.pdb$|\.ent$//;
	}
	
	$mol->{NAME} = $name;
}

#
# Fast read routines
#

sub read_pdb_fast {

	#<
	#? Read in a pdb file.
	#. Is able to read files containing multiple MODEL records and Autodock
	#  output files.  Sets the $mol->{HAS_BAD_BONDORDERS} flag
	#. Also saves some information from the PDB header records these are:
	#, HEADER
	#, COMPND
	#, HETNAM
	#, USER
	#. Assumes that atom names starting with Q are pseudo atoms.  These
	#  atoms are set to dummy (Du) type.
	#; Requires: Filename, unused, unused, option string
	#; Returns: Ensemble containing molecule or undef if file open fails
	#>

	my $infile = $_[0];
	
	my $mols;
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Can not open file $infile!\n");
		return undef;
	}
	
	MOL: while (1) {

		my $newmol = read_pdb_molecule_fast($infile);
		
		last if $newmol == 0;

		push @$mols, $newmol;
	}
	
	close_file();

	if (!defined $mols->[0]{NUMATOMS}) {
		silico_msg('e', "Could not determine the number of atoms in file $infile!\n",
				"There may be an error in the ATOM and/or HETATM records.\n");
		return undef;
	}
	
	return $mols;
}

sub read_pdb_molecule_fast {

	#<
	#? Read a PDB molecule
	#; Requires: nothing
	#; Returns: molecule, 0 if no more atom lines, undef if error

	my $infile = $_[0];
	
	my $atomlines = 0;  # Number of ATOM records
	my $linecount = 0;
	my $mol;

	while (1) {

		++$linecount;
		$_ = <INFILE>;
		last if !defined $_;
		chomp;

		if (/^ATOM/ || /^HETATM/) {
			my $atom;
			$atom = read_pdb_atom_fast($_, $mol);
			++$atomlines;
			push @{$mol->{ATOMS}}, $atom;
		} elsif (/^END/) {
			last;
		} elsif (/^MODEL/ && $atomlines) {
			last;
		} elsif (/^HEADER/) {
			pdb_read_header($_, $mol);
		} elsif (/^TITLE/) {
			pdb_read_title($_, $mol);
		} elsif  (/^COMPND/) {
			pdb_read_compnd($_, $mol);
		} elsif (/^CRYST1/) {
			pdb_read_cryst($_, $mol);
		} elsif (/^HETNAM/) {
			pdb_read_hetnam($_, $mol);
		}
	}

	return 0 if (!$mol->{NUMATOMS});

	# Set name if not defined
	if (!defined $mol->{NAME}) {
		$mol->{NAME} = $infile;
	
		# Remove 3-letter extensions
		$mol->{NAME} =~ s/\.pdb$|\.ent$//;
	}

	return $mol;
}

sub read_pdb_atom_fast {

	#<
	#? Read in an ATOM line from a pdb file, getting only some fields.
	#; Requires: line, molecule
	#; Returns: An atom
	#>

	my $line = $_[0];
	my $mol = $_[1];

	my $atom;

	++$mol->{NUMATOMS};

	# Pad out line to prevent 'substring outside string' errors
	$line .= '                              ';

	# The original (FORTRAN) format for pdb files is:
	# ( A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2 )

	$atom->{NUM}		= substr($line,6,5);	#Atom number
	$atom->{NAME}		= substr($line,12,4);	#Atom name
	$atom->{SUBNAME}	= substr($line,17,3);	#Residue name
	$atom->{SUBID}		= substr($line,22,4);	#Residue number
	$atom->{X}		= substr($line,30,8);	#x
	$atom->{Y}		= substr($line,38,8);	#y
	$atom->{Z}		= substr($line,46,8);	#z
	$atom->{ELEMENT}	= substr($line,76,2);	#Element
	
	undef $atom->{ELEMENT} if $atom->{ELEMENT} eq '  ';
	
	# Work out the element
	atom_guess_element($atom, $mol);
	
	# Renumber atoms (Required for CONECT records to work)
	$atom->{NUM} = $mol->{NUMATOMS};

	# Some pathological pdb files have spaces
	# in the atom name! Get rid of them.
	$atom->{NAME} =~ s/(\b) +(\b)/_/;

	return $atom;
}

sub clean_pdb_line {

	#<
	#? Remove spaces and line name from PDB record
	#>

	my $string = $_[1] || '';
	
	$_[0] =~ s/$string//;
	$_[0] =~ s/ +/ /g;
	$_[0] =~ s/^ //;
	$_[0] =~ s/ $//;
}
	

sub molecule_connect_pdb_atoms_by_name {

	#<
	#? Bond atoms together using a pdb lookup table
	#. This routine bonds together atoms proteins using a lookup table
	#. This routine depends on the atom names being correct (See also molecule_pdb_rename_h).
	#  Terminal residues etc are not present in amino_acid.dat file
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];

	my @connected;
	my $cutoff = 2;
	my $distance_sq;
	my $starttime = (times)[0];
	
	# Read in file of pdb connectivites
	read_pdb_connected();
	
	# Delete any bonds already present
	for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {
		
		my $atom1 = $mol->{ATOMS}[$i];
		@{$atom1->{CONNECT}} = ();
		@{$atom1->{BORDERS}} = ();
	}
	
	for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {
	
		my $atom1 = $mol->{ATOMS}[$i];
		
		my $aname = $atom1->{NAME};
		$aname =~ s/ //g;
		my $resname = $atom1->{SUBNAME};
		$resname =~ s/ //g;
		
		my $hash = "$resname $aname";
		if (defined $Silico::pdb_connected->{$hash}) {
			@connected = @{$Silico::pdb_connected->{$hash}};
		} else {
			# Don't print warning for hetatms
			if ($atom1->{LABEL} =~ /ATOM/) {
				silico_msg('w', "Unknown atom type!\n",
						"Atom: $atom1->{NUM} Name: \"$aname\" Substructure: \"$resname\"\n");
			}
			undef @connected;
		}
		
		# Make cutoff longer for S or P
		my $acutoff = $cutoff;
		$acutoff += 0.5 if ($atom1->{ELEMENT} eq 'S' || $atom1->{ELEMENT} eq 'P');
		my $acutoff2 = $acutoff**2;

		# Do a distance based search to find connected atoms
		for (my $j=0; $j < $i; ++$j) {
		
			my $atom2 = $mol->{ATOMS}[$j];
			
			next if abs($atom1->{X} - $atom2->{X}) > $acutoff;
			next if abs($atom1->{Y} - $atom2->{Y}) > $acutoff;
			next if abs($atom1->{Z} - $atom2->{Z}) > $acutoff;
			
			next if ($distance_sq = distance_sq($atom1, $atom2) > $acutoff2);
					
			# Check to see if there should be a pdb
			#  bond between atom1 and atom2
			my $bo = undef;
			$aname = $atom2->{NAME};
			$aname =~ s/ //g;
							
			for (my $k=0;$k<=$#connected;++$k) {
				# Note some of the connected atom names have +, * or @ in them
				# to denote connection to a different residue
				if ($connected[$k] =~/^$aname/) {
					$bo = $Silico::pdb_bondorder->{$hash}[$k];
					last;
				}
			}
			
			next if !defined $bo;
			
			bond_create($mol, $i, $j, $bo);
		}
	}
	
	silico_msg('t', "Time in ".subname().": ".calc_human_readable_time((times)[0]-$starttime)."\n");
}




##################################################################
#
#	PDB write routines
#
##################################################################


sub write_pdb {
	
	#<
	#? Write a pdb file
	#; Requires: Ensemble or molecule, filename, $options.
	#; Returns: Zero if file open failed otherwise returns 1.
	#; Options: FAST, NOSORT, NOCONECT, BOND, MODEL
	#. Note: this routine renumbers the ensemble of pdb files
	#. Also note: the handling of pdb chains is somewhat confused
	#. A TER record will be written following an atom with atom->{TER} set
	#>

	my $molecules = ens($_[0]);
	my $outfile = $_[1];
	my $options = $_[2] || '';

	my $molcount = 0;

	$options = uc $options if $options;
	$outfile = get_ofilebase($molecules).".pdb" if !$outfile;
	
	my $fr = open_molfile(undef, $outfile, undef, 'write', $options);
	if (!defined $fr) {
		silico_msg('e', "Could not open $outfile");
		return undef;
	}
	my $FH = $fr->{FH};
	
	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Writing pdb file: $outfile");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}
	
	# Write out molecules
	foreach my $mol (@$molecules) {
	
		die if !defined $mol;
	
		++$molcount;
		if (!defined $FH) {
			silico_msg ('e',"Attempt to write to unopened filehandle. Perhaps you don't have write permission?\n");
			return undef;
		}
		
		write_pdb_molecule($fr, $mol, $options);
	}
	
	close_molfile($fr);
	
	return 1;
}

sub write_pdb_molecule {

	#<
	#? Write a single record from a pdb file.
	#; Options: FAST, NOSORT, BOND, NOCONECT, DUCONECT, MODEL
	#; Requires: file record, molecule, options
	#. See read_pdb for general description
	#; Returns: 1 or undef if failed
	#>

	my $fr = $_[0]  || silico_msg('d', "File record not defined\n");
	my $mol = $_[1] || silico_msg('d', "Molecule not defined\n");
	my $options = $_[2] || '';

	my $atomcount = 0;
	my $dupflag;
	my $reslen = 4;
	
	my $FH = $fr->{FH};
	
	++$fr->{MOL_WRITE_COUNT};
	
	# Set standard output options based on flags
	add_option("FAST", $options) 	  if get_lflag('fast');
	add_option("NOSORT", $options) 	  if (get_lflag('nosort') || $options =~ /\bFAST\b/);
	add_option("BOND", $options) 	  if get_lflag('bond');
	add_option("NOCONECT", $options)  if get_lflag('noconect');
	add_option("DUPCONECT", $options) if get_lflag('dupconect');
	add_option("MODEL", $options) 	  if get_lflag('pdb-model');
	
	make_atom_names_unique($mol) if get_flag('unique-atom-names', 'l');

	$reslen = 4 if get_flag('pdb-long-resname', 'l')  || $options && $options =~ /\bLONG_RESNAME\b/;
	
	# Generate bonds if specified in options
	molecule_check_and_fix_connectivity($mol, $options) if ($options =~ /\bBOND\b/);
	
	if (!defined $mol->{NUMATOMS}) {
		carp();
		die "Molecule NUMATOMS not defined\n";
	}
	
	# Make sure a single file does not contain duplicate compound names which 
	# creates problems when reading in Pymol
	my $n = $mol->{NAME} || 'mol';
	if (defined $fr->{NAMES}{$n}) {
	
		$n .= "_".($fr->{MOL_WRITE_COUNT}+1);
		$mol->{NAME} = $n;
	}
	++$fr->{NAMES}{$n};
	#print "Name: $mol->{NAME}\n";
	
	print $FH "HEADER    $mol->{NAME}\n" if $mol->{NAME};
	print $FH "COMPND    $mol->{NAME}\n" if $mol->{NAME};
	print $FH "TITLE     $mol->{TITLE}\n" if $mol->{TITLE};
	
	write_pdb_sdfdata("DATA", $mol->{SDF_DATA}, $FH);

	if ($options =~ /\bDUPCONECT\b/) {
		print $FH "REMARK    Multiple bonds are recorded as duplicate bond records\n";
		convert_aromatic_bonds_kekule($mol);
		$dupflag = 1;
	}
	
	if ($mol->{NUMATOMS} > 99999) {
		print $FH "REMARK    Molecule contains $mol->{NUMATOMS} atoms which is > 99,999. Extended format(s) are used for atom numbers\n";
		silico_msg('w', "Molecule contains $mol->{NUMATOMS} atoms which is > 99,999. Extended format(s) will be used for atom numbers\n");
		mol_warn($mol,'More than 99,999 atoms. Using PDB atom number format(s)');
	}
	
	write_pdb_cryst($mol, $FH);
	if ($options =~ /\bMODEL\b/) {
		print $FH "MODEL        ".($fr->{MOL_WRITE_COUNT}+1)."\n";	
	}

	# Print atom records
	my $i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		++$atomcount;
		my @f;
		
		# Default values
		$f[0] = "ATOM";			#Record type
		$f[3] = "";			#Alternate location indicator
		$f[5] = '';
		$f[6] = 1;			#Residue number
		$f[7] = "";			#Residue number extension
		$f[8] = 0;			#x
		$f[9] = 0;			#y
		$f[10] = 0;			#z
		$f[11] = 0;			#Occupancy
		$f[12] = 0;			#Temperature Factor
		$f[13] = "";			#SegID
		$f[15] = "";			#Formal charge
			
		$f[0] = $atom->{LABEL}	if defined $atom->{LABEL};	#Record type
		$f[1] = pdb_format_atom_number($atom->{NUM});	#Atom number
		$f[2] = pdb_format_atom_name($atom, $mol);		#Atom name
		$f[3] = $atom->{ALT}	if defined $atom->{ALT};	#Alternate location indicator
		$f[4] = pdb_format_subname($atom, $mol, $reslen);	#Residue name
		$f[5] = $atom->{CHAIN}  if defined $atom->{CHAIN};	#Chain
		$f[6] = $atom->{SUBID}	if defined $atom->{SUBID};      #Residue number
		$f[7] = $atom->{EXT}	if defined $atom->{EXT};	#Residue number extension
		$f[8] = $atom->{X}	if defined $atom->{X};		#x
		$f[9] = $atom->{Y}	if defined $atom->{Y};		#y
		$f[10] = $atom->{Z}	if defined $atom->{Z};		#z
		$f[11] = $atom->{OCC}	if defined $atom->{OCC};	#Occupancy
		$f[12] = $atom->{TEMP}  if defined $atom->{TEMP};	#Temperature Factor
		$f[13] = $atom->{SEGID} if defined $atom->{SEGID};	#SegID
		$f[14] = uc ($atom->{ELEMENT} || ''); 			#Element
		
		# PDB formal charge format is 1+, 2-, etc
		if (defined $atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE}) {
			$atom->{FORMAL_CHARGE} = sprintf "%d", $atom->{FORMAL_CHARGE};
			if ($atom->{FORMAL_CHARGE} < 0) {
				$f[15] = -$atom->{FORMAL_CHARGE}."-";
			} 
			if ($atom->{FORMAL_CHARGE} > 0) {
				$f[15] = $atom->{FORMAL_CHARGE}."+";
			} 
		} 
		
		if (length $f[3] > 1) {
			$f[3] = substr ($f[3],0,1);
			mol_warn($mol,'PDB alternate location indicator too long');
		}
		
		if (length $f[5] > 1) {
			$f[5] = substr ($f[5],0,1);
			mol_warn($mol,'PDB alternate chain indicator too long');
		}

		# Restart residue numbering if residue number is greater than 9999
		if ($f[6] > 9999) {
			$f[6] = pdb_format_residue_number($f[6]);
			mol_warn($mol,'PDB residue number too long. Formatting using hybrid format');
		}
		
		# The original (FORTRAN) format for pdb files is:
		# We are now allowing for 4 char residue names with flag
		# ( A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2 )
		
		printf $FH "%-6s%5s %-4s%1s%4s%1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n",@f;
	
		# Print TER records if atom->{TER} flag is set
		print $FH "TER\n" if $atom->{TER};
	}

	# Do CONECT records
	if ($options !~ /\bNOCONECT\b/) {
	
		my $i = 0;
		foreach my $atom (@{$mol->{ATOMS}}) {
		
			++$i;
				
			# Skip atoms with no bonds
			next if ($#{$atom->{CONNECT}} < 0);
			
			printf $FH "CONECT%5s", pdb_format_atom_number($atom->{NUM});
			
			my $bo;
			my $j = -1;
			my $c = 0;
			foreach my $connum (@{$atom->{CONNECT}}) {

				++$j;

				if (!defined $mol->{ATOMS}[$connum]) {
					silico_msg('e', "Atom number $i is connected to a nonexistent atom (number $connum)!\n",
							"Skipping.\n");
					next;
				}

				if ($dupflag) {
					$bo = $atom->{BORDERS}[$j];
					if ($bo > 3 || $bo == 0 ) {
						mol_warn("Bondorder $bo converted to single");
						$bo = 1;
					}
				} else {
					$bo = 1;
				}

				foreach (1..$bo) {
					++$c;
					if ($c > 4) {
						print  $FH "\n";
						printf $FH "CONECT%5s", pdb_format_atom_number($atom->{NUM});			
						$c = 1;
					}
					printf $FH "%5s",  pdb_format_atom_number($mol->{ATOMS}[$connum]{NUM});		
				}
			}
			print $FH "\n";
		}
	}

	if ($options =~ /\bMODEL\b/) {
		print $FH "ENDMDL\n";	
	} else {
		print $FH "END\n"; # To make Xplor happy!
	}
	
	mol_print_warnings($mol);
	#molecule_printout($mol);
	#die;
	
	return 1;
}



sub write_pdb_sdfdata {

	#<
	#? Print SDF_data to pdb file
	#. New version November 2024
	#. Note. There is currently no read routine for this format
	#; Requires: field_name, data_to_print, filehandle, data level
	#; Returns: nothing
	#>

	my $field = $_[0]; # Field name
	my $data = $_[1];  # Pointer to data
	my $FH = $_[2];    # Filehandle
	my $level = $_[3] || 0;

	return if !defined $data;
		
	my $s = "REMARK     ".('  ' x $level);
	print $FH $s."START SILICO DATA\n" if $level == 0;;
	
	my $type = ref($data) || "SCALAR";
	
	if ($type eq 'SCALAR') {
		
		my @ndata = formatdata($data, length($field)+3);
		
		print $FH $s."'$field' ".(shift(@ndata));
			
		foreach my $d(@ndata) {
			print $FH "\n".$s."  ".$d;
		}
		
		print $FH "\n";
		
	} elsif ($type eq 'ARRAY') {
	
		my $n = $#{$data}+1;
		
		print $FH $s."'$field' <ARRAY> of $n elements. Level $level\n";
		print $FH $s."{\n";
		
		my $i = 0;
		foreach my $d (@$data) {
			++$i;
			write_pdb_sdfdata("$field $i", $d, $FH, $level+1);
		}
		
		print $FH $s."}\n";
			
	} elsif ($type eq 'HASH') {
	
		my $n =  (keys %$data);
		
		print $FH $s."'$field' <HASH> of $n elements. Level $level\n";
		print $FH $s."{\n";
		
		my $i = 0;
		foreach my $key (keys %$data) {
		
			++$i;
			
			my $d = $data->{$key};
		
			write_pdb_sdfdata($key, $d, $FH, $level+1);
		}
		
		print $FH $s."}\n";
		
	} else {
		silico_msg('w', "Don't know how to handle SDF DATA field of type $type!\n");
	}
	
	print $FH $s."END SILICO DATA\n" if $level == 0;
}

sub formatdata {

	#<
	#? Format SDF_DATA for output to PDB file
	#. Long data is split up into lines.  All data is surrounded by single quotes.  A
	#  backslash at the end of the line indicates that the line is to be continued
	#; Requires: string, length of key
	#; Returns: array
	#>

	my $data = $_[0];
	my $lkey = $_[1] || 0;
	
	my $i = 0;
	my $maxlen = 59;
	my $len = $maxlen - $lkey;
	my @ndata;

	return "''" if !defined $data || $data eq '';
	
	chomp $data;
	
	my @f = split "\n", $data;
	
	foreach my $v (@f) {
	
		while (length($v) > $len) {
	
			push @ndata,"'".(substr($v, 0, $len, ''))."'\\";
			$len = $maxlen;
			++$i;	
		}
		push @ndata,"'".($v)."'";
	}
	
	return @ndata;
}

sub write_pdb_cryst {

	#<
	#? Write PDB crystal record
	#; Requires: molecule, filehandle
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $FH = $_[1];

	my $space_group;

	# Print out unit cell data
	return if !defined $mol->{CELL}[0];
	
	# In the printout, the last value (after $space_group) is the 
	# Z-value. This is the number of instances of the most numerous
	# chain in each unit cell (as distinct from each asymmetric unit).
	# Consequently, it will depend on the space group as well as on
	# the stoichiometry. If it is not set elsewhere, the value 1 will be used.

	if (defined $mol->{CELL}[0] &&
	    defined $mol->{CELL}[1] &&
	    defined $mol->{CELL}[2] &&
	    defined $mol->{CELL_ALPHA} &&
	    defined $mol->{CELL_BETA} &&
	    defined $mol->{CELL_GAMMA}) {
		
		$space_group = space_group_silico_to_pdb($mol->{SPACE_GROUP});
		if (!defined $space_group) {
			mol_warn($mol,"Setting space group to P 1");
			$space_group = 'P 1';
		}
		$mol->{CHAINS_PER_CELL} ||= 1;
		
		printf $FH "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
			"CRYST1",
			$mol->{CELL}[0],
			$mol->{CELL}[1],
			$mol->{CELL}[2],
			$mol->{CELL_ALPHA},
			$mol->{CELL_BETA},
			$mol->{CELL_GAMMA},
			$space_group,
			$mol->{CHAINS_PER_CELL};
			
	} elsif ( defined $mol->{CELL}[0] &&
		  defined $mol->{CELL}[1] &&
		  defined $mol->{CELL}[2] ) {

		mol_warn("Setting space group to P 1");

		$space_group = 'P 1';
		$mol->{CHAINS_PER_CELL} ||= 1;
		printf $FH "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
			"CRYST1",
			$mol->{CELL}[0],
			$mol->{CELL}[1],
			$mol->{CELL}[2],
			90, 90, 90,
			$space_group,
			$mol->{CHAINS_PER_CELL};
	} else { 
		
		 my $warn = "Bad unit cell data!\n".
			"\ta: ".($mol->{CELL}[0] || "UNDEF")." ".
			"b: ".($mol->{CELL}[1] || "UNDEF")." ".
			"c: ".($mol->{CELL}[2] || "UNDEF")." ".
			"alpha: ".($mol->{CELL_ALPHA} || "UNDEF")." ".
			"beta: ".($mol->{CELL_BETA} || "UNDEF")." ".
			"gamma: ".($mol->{CELL_GAMMA} || "UNDEF")."\n".
			"\tUnit cell data has been omitted";
			
			mol_warn($mol, $warn);
	}
}

sub space_group_silico_to_pdb {
	
	#<
	#? Turn a Silico space group code into a PDB space group string.
	#; Requires: string
	#; Returns: string
	#>
	
	my $sg = $_[0];
	my $file = $_[1] || "$ENV{SILICO_HOME}"."/data/pdb_space_groups.dat";
	
	my %hash = parse_pdb_space_group($file);
	
	return undef if !defined $sg;

	if (!defined $hash{$sg}) {
		silico_msg('e', "No PDB format available for space group code \"$sg\".\n",
				"Please check the file $file.\n");
		return undef;
	}
	
	return $hash{$sg};
}

sub space_group_pdb_to_silico {
	
	#<
	#? Turn a PDB space group string into a Silico space group code.
	#. This subroutine assumes that there is only one key (Silico space group code)
	#  for each value (PDB space group string).
	#; Requires: string
	#; Returns: string
	#>
	
	my $sg = $_[0];
	my $file = $_[1] || "$ENV{SILICO_HOME}"."/data/pdb_space_groups.dat";

		
	my %hash = parse_pdb_space_group($file);
	
	foreach my $key (keys %hash) {
		return "$key" if $hash{$key} eq $sg;
	}
	
	silico_msg('w', "Can not identify space group code for the space group \"$sg\".\n",
			"\tUsing space group P1.\n");

	return "0011";
}

sub parse_pdb_space_group {
	
	#<
	#? Parse the file pdb_space_group.dat
	#; Requires: File
	#; Returns: hash
	#>
	
	my $file = $_[0];
	
	my %hash;

	return %$Silico::pdb_space_group_hash if defined $Silico::pdb_space_group_hash;
	
	if (!open(FILE, "$file")) {
		file_read_error($file, 0);
		return %hash;
	}
	
	my $c = 0;
	while (<FILE>) {
		
		++$c;
		
		my @f;
		
		my $line = $_;
		chomp $line;
		
		next if $line =~ /^\s*#/;
		next if $line =~ /^\s*$/;
		
		@f = split /\t/, $line;
		
		if ($#f > 1) {
			silico_msg('w', "Error on line $c of file $file: too many elements!\n",
					"Line: \"$_\"\n",
					"Skipping.\n");
			next;
		} elsif ($#f < 1) {
			silico_msg('w', "Error on line $c of file $file: too few elements!\n",
					"Line: \"$_\"\n",
					"Skipping.\n");
			next;
		} else {
			$hash{"$f[0]"} = "$f[1]"
		}
	}

	%$Silico::pdb_space_group_hash = %hash;
	
	return %hash;
}

sub pdb_format_atom_name {

	#<
	#? Format atom name for pdb file
	#. Pads atom names with single character element names with a space
	# (Insight sometimes doesn't work otherwise)
	#; Requires: atom, $mol->{WARN}{FIELD_TRUNCATION} as a local variable
	#; Returns; string
	#>
	
	my $atom = $_[0];
	my $mol = $_[1];
	
	my $name = $atom->{NAME};
	
	# Default value
	$name = "X" if !defined $name;
	
	$name =~ s/ //g;
	
	if ($name !~ /^\d/ && length ($name) < 4  && length ($atom->{ELEMENT}) == 1) {
		$name = ' '.$name;
	}
	
	if (length $name > 4) {
		$name = substr ($name,0,4);
		mol_warn($mol, "Atom name was truncated");
	}

	return $name;
}

sub pdb_format_subname {

	#<
	#? Format substructure name to a length of four characters.  
	#. Length can be set to 3 or 4 characters
	#  Substructure name is truncated if it is too long
	#; Requires: atom, molecule, maxumum substructure name length
 	#; Returns; string
	#>
	
	my $atom = $_[0];
	my $mol = $_[1];
	my $maxlen = $_[2] || 3;
	
	my $subname = $atom->{SUBNAME};
	
	if ($subname =~ /^\s*$/) {
		$subname = 'UNK';
		mol_warn($mol, "Substructre renamed");
	}
	
	$subname =~ s/^\s*//;
	
	my $l = length ($subname);
	
	if ($l > $maxlen) {
		
		$subname = substr ($subname,0,$maxlen);
		mol_warn($mol, "Substructre renamed");
		$l = $maxlen;
	}

	if ($l == $maxlen) {
		# Do nothing
	} elsif ($l == $maxlen-1) {
		$subname = $subname." ";
	} elsif ($l == $maxlen-2) {
		$subname = $subname."  ";
	} elsif ($l == $maxlen-3) {
		$subname = $subname."   ";
	}
	
	return $subname;
}


sub pdb_format_residue_number {

	#<
	#? Format PDB residue numbers so that they change to base16 when greater than 9999 starting at A000
	#; Requires: decimal number
	#; Returns: string
	#>
	
	return $_[0] if $_[0] <= 9999;
	
	# Convert to base16 starting at A000
	return base16($_[0] - 10000 + 40960);	
}

sub pdb_format_atom_number {

	#<
	#? Format PDB atom numbers so that they change to base16 when greater than 9999 starting at A0000
	#. Then change to base32 when greater than FFFFF starting at G0000
	#. Base16 is to maximise compatiblity with OpenFF
	#; Requires: decimal number
	#; Returns: string
	#>

	return $_[0] if $_[0] <= 99999;
	
	# Use base16 starting at A0000 ( = 655360) if base10 has overflowed (FFFFF = 1048575)
	return base16($_[0] - 100000 + 655360) if $_[0] <= 1048575 + 100000 - 655360;
	
	# Use base36 starting at G0000 ( = 26873856) if base16 has overflowed
	return base36($_[0] + 26873856 - (1048576 + 100000 - 655360));
}

sub base16 {

        #<
        #? Convert a decimal number to base 36
        #. Source: https://www.perlmonks.org/?node_id=1197777
        #>
        
        my $val = $_[0];
	
	my $newstr = '';
        while ($val) {
                $newstr = substr('0123456789ABCDEF', $val % 16, 1) . $newstr;
                $val = int $val / 16;
        }
        return $newstr || '0';
}

sub base36 {

        #<
        #? Convert a decimal number to base 36
        #. Source: https://www.perlmonks.org/?node_id=1197777
        #>
        
        my $val = $_[0];
	
        my $newstr = '';
        while ($val) {
                #$newstr = substr('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ', $val % 36, 1) . $newstr;
		$newstr = substr('0123456789abcdefghijklmnopqrstuvwxyz', $val % 36, 1) . $newstr;
                $val = int $val / 36;
        }
        return $newstr || '0';
}


##################################################################
#
#	Misc PDB routines
#
##################################################################


sub read_pdb_connected {

	#<
	#? Read in datafile of amino acid connected atoms
	#. Requires; optional datafile (Default datafile is 'amino_acid_atoms.dat' which is found in the
	#  SILICO_DATA directory)
	#; Returns: nothing
	#; Sets $Silico::pdb_connected and $Silico::pdb_bondorder
	#>
	
	my $datafile = $_[0];
	
	$datafile ||= $Silico::data_dir."amino_acid_atoms.dat";
	
	open (DATAFILE, $datafile) || silico_msg('d', "Could not open file $datafile for reading!\n");
	while (<DATAFILE>) {

		# Skip comments and blank lines
		next if /^#/;
		next if !/\w/;

		my @f = split;

		my @array1 = ();
		my @array2 = ();
		for (my $i = 2; $i <= 5; ++$i) {
			
			if (!defined $f[$i]) {
				silico_msg('d', "Error in file $datafile!\n");
			}
			
			last if $f[$i] eq '-';
			
			push @array1,  $f[$i];
		}
		
		for (my $i = 6; $i <= 9; ++$i) {
		
			if (!defined $f[$i]) {
				silico_msg('d', "Error in file $datafile!\n");
			}
			
			last if $f[$i] eq '-';
			
			push @array2,  $f[$i];
		}

		# Make array of residue name, atom name
		my $hashname = "$f[0] $f[1]";
		@{$Silico::pdb_connected->{$hashname}} = @array1;
		@{$Silico::pdb_bondorder->{$hashname}} = @array2;
	}

	close DATAFILE;
}

sub pdb_molecule_sort_atoms {

	#<
	#? Routine to sort atoms in to a standard order
	#  Sorts by SEGID, CHAIN, residue number then by atom number.
	#  HETATMS are sorted to the end of the file
	#  Atoms with a ' ' chain specifier are sorted to the end of the file.
	#. Atoms are then renumbered.
	#; Requires: Molecule. (optional, flags to ignore chain and segid when sorting).
	#; Returns: Nothing
	#>

	my $mol = $_[0];
	my $ignore_chain = $_[1];
	my $ignore_segid = $_[2];
	
	my @atomorder = qw(N HN CA H C O CB HB HB1 HB2 HB3 CG CG1 CG2 OG1 SG HG HG1 HG2 HG3 HG HG1 HG11 HG12 HG13 HG2 HG21 HG22 HG23 HG3 CD CD1 CD2 ND2 OD1 OD2 SD HD HD1 HD2 HD3 HD11 HD12 HD21 HD22 CE CE1 CE2 CE3 NE OE1 OE2 HE HE1 HE2 HE3 CZ HZ HZ1 HZ2 HZ3  CH2 NH1 NH2 OH HH HH1 HH2 HH3);
	
	# Globals - so that they can be used inside pdb_sort
	$Silico::uc = !$ignore_chain;
	$Silico::us = !$ignore_segid;

	my $i = 0;
	foreach (@atomorder) {
		$Silico::ord->{$_} = $i;
		++$i;
	}

	# Save the old ordering
	$i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
		$atom->{OLDINDEX} = $i;
		$atom->{NAME} ||= 'X';
		$atom->{NAME} =~ s/ //g;
		++$i;
	}

	# Sort routine
	sub pdb_sort {

		# Sort by SEGID
		if ($Silico::us && (defined $a->{SEGID} || defined $b->{SEGID})) {
			return -1 if !defined $a->{SEGID};
			return 1 if !defined $b->{SEGID};
			return -1 if $a->{SEGID} lt $b->{SEGID};
			return 1 if $a->{SEGID} gt $b->{SEGID};
		}

		# Sort by chain
		if ($Silico::uc && (defined $a->{CHAIN} || defined $b->{CHAIN})) {

			# Sort so undefined chain comes after a defined chain
			return 1 if !defined $a->{CHAIN};
			return -1 if !defined $b->{CHAIN};

			# Sort so that chain ' ' comes after a defined chain
			return 1  if $a->{CHAIN} eq ' ' && $b->{CHAIN} ne ' ';
			return -1 if $a->{CHAIN} ne ' ' && $b->{CHAIN} eq ' ';

			return -1 if $a->{CHAIN} lt $b->{CHAIN};
			return 1 if $a->{CHAIN} gt $b->{CHAIN};
		}

		# Sort by substructure number
		if (defined $a->{SUBID} || defined $b->{SUBID}) {
			return -1 if !defined $a->{SUBID};
			return 1 if !defined $b->{SUBID};
			return -1 if $a->{SUBID} < $b->{SUBID};
			return 1 if $a->{SUBID} > $b->{SUBID};
		}

		# Sort by substructure number extension
		if (defined $a->{EXT} || defined $b->{EXT}) {
			return -1 if !defined $a->{EXT};
			return 1 if !defined $b->{EXT};
			return -1 if $a->{EXT} lt $b->{EXT};
			return 1 if $a->{EXT} gt $b->{EXT};
		}

		# Sort by atom name
		if (defined $Silico::ord->{$a->{NAME}} || defined $Silico::ord->{$b->{NAME}}) {
			return 1 if !defined $Silico::ord->{$a->{NAME}};
			return -1 if !defined $Silico::ord->{$b->{NAME}};
			return -1 if $Silico::ord->{$a->{NAME}} lt $Silico::ord->{$b->{NAME}};
			return 1 if $Silico::ord->{$a->{NAME}} gt $Silico::ord->{$b->{NAME}};
		}

		# If all else fails sort by atom number
		return 1  if !defined $a->{NUM};
		return -1 if !defined $b->{NUM};
		
		$a->{NUM} <=> $b->{NUM};
	}

	# Sort the atoms
	@{$mol->{ATOMS}} = sort pdb_sort (@{$mol->{ATOMS}});

	# Make translation table
	$i = 0;
	my @ctable;
	foreach my $atom (@{$mol->{ATOMS}}) {
		$ctable[$mol->{ATOMS}[$i]{OLDINDEX}] = $i;
		++$i;
	}

	# Translate connection table
	$i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
		$atom = $mol->{ATOMS}[$i];

		# Translate connection table
		for (my $j=0; $j<=$#{$atom->{CONNECT}}; ++$j) {

			#. Replace the old atom number with
			# the translated one
			$atom->{CONNECT}[$j] = $ctable[$atom->{CONNECT}[$j]];
		}
		++$i;
	}

	# Molecule now needs to be renumbered
	molecule_renumber($mol);
	
	# Rings are now bad
	$mol->{BAD_RINGS} = 1;
}

sub fix_pdb_chains {

	#<
	#? Renumber pdb chains based on changes in CHAIN, SEGID or a backwards step in SUBID
	#. If 'start' and the 'force flag' are not specified this routine will only renumber chains if
	#  required
	#. If 'start' is set and the 'force flag' is not set then the routine will renumber if
	#  'start' does not match the starting letter and if it is necessary
	#. If 'force flag' is set then the chains will be renumbered from 'A' by default.
	#. Note: This subroutine deletes the 'TERCOUNT' field can potentially cause problems if the pdb file has been modified
	#; Requires: molecule, starting chain letter (optional), flag to force renumbering
	#; Returns: last chain letter of current structure, flag to indicate that chains were renumbered
	#>

	my $mol = $_[0];
	my $start = $_[1];
	my $force = $_[2];

	my $renumber= 0; # Renumber all subsequent chain
	
	my $newchain;
	if (defined $start) {
	
		# Starting letter has been specfied
		
		$newchain = $start;
		$renumber = 1;
			
	} else {
	
		# No starting letter has been specfied
		
		if ($mol->{ATOMS}[0]{CHAIN} =~ /[A-Z]/i) {

			# The current starting letter is OK
			# Set the start to the chain of the first atom
			$newchain = $mol->{ATOMS}[0]{CHAIN};
			
		} else {
		
			# The current starting letter is bad
			# Renumber from A
			$newchain = 'A';
			$renumber = 1;
		}
	}
	
	# The force renumbering flag was set
	$renumber = 1 if $force;

	my $oldchain = $mol->{ATOMS}[0]{CHAIN} || '' ;
	my $oldsegid = $mol->{ATOMS}[0]{SEGID} || '' ;
	my $oldsubid = $mol->{ATOMS}[0]{SUBID} || 0;

	my $i = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		++$i;

		my $chain = $atom->{CHAIN} || '' ;
		my $segid = $atom->{SEGID} || '';
		my $subid = $atom->{SUBID} || 0;
		
		# Chain number remains the same
		if ($chain eq $oldchain) {
		
			# If segid, or subid changes but the chain remains the same
			# then force renumbering
			if ($subid < $oldsubid || $segid ne $oldsegid) {
				
				pdb_increment_chain($newchain);
				$renumber = 1;
			}

		} else {

			# Chain number changed but is a ' ' or something else
			# Or renumber is set
			if ($renumber || $chain !~ /[A-Z]/i) {

				increment_chain($newchain);
				$renumber = 1;
		
			} else {

				$newchain = $chain;
			}
		}
		
		$atom->{CHAIN} = $newchain;
		
		$oldchain = $chain;
		$oldsubid = $subid;
		$oldsegid = $segid;

		# Clean up TERCOUNT so that it does not make problems later
		delete $atom->{TERCOUNT};
	}

	return $newchain, $renumber;
}

sub pdb_increment_chain {

	#<
	#? Increment pdb chain
	#; Requires: Alphabetical character
	#; Returns: nothing
	#>
	
	# chain = $_[0]

	if ($_[0] eq ' ' || $_[0] eq '' || !defined $_[0] || $_[0] !~ /^[a-zA-Z]$/) {

		$_[0] = 'A';

	} else {

		++$_[0];
	
		# Cycle A-Z, a-z, etc
		$_[0] = 'a' if $_[0] eq 'AA';
		$_[0] = 'A' if $_[0] eq 'aa';
	}
}

sub printbox {

	#<
	#? Write out a pdb format box
	#. Writes eight atoms. One at each cornder of the box. Compatible with the box produced by DOCK
	#; Requires: maxx, maxy, maxz, minz, miny, minz, (optional) filename
	#  default 'box.pdb'
	#; Returns: undef if file open failed
	#>
	
	my $maxx = $_[0];
	my $maxy = $_[1];
	my $maxz = $_[2];
	my $minx = $_[3];
	my $miny = $_[4];
	my $minz = $_[5];
	my $filename = $_[6] || 'box.pdb';

	if (!open (OUTFILE, ">$filename")) {
		silico_msg('e', "Can not create or open file $filename for writing!\n");
		return undef;
	}

	print  OUTFILE "HEADER  box from $0\n";
	printf OUTFILE "REMARK    CENTER (X Y Z)  %8.3f %8.3f %8.3f\n", ($maxx-$minx)/2+$minx, ($maxy-$miny)/2+$miny, ($maxz-$minz)/2+$minz;
	printf OUTFILE "REMARK    DIMENSIONS (X Y Z) %8.3f %8.3f %8.3f\n", ($maxx-$minx), ($maxy-$miny), ($maxz-$minz);
	printf OUTFILE "ATOM      1  DUA BOX     1    %8.3f%8.3f%8.3f\n",$minx,$miny,$minz;
	printf OUTFILE "ATOM      2  DUB BOX     1    %8.3f%8.3f%8.3f\n",$maxx,$miny,$minz;
	printf OUTFILE "ATOM      3  DUC BOX     1    %8.3f%8.3f%8.3f\n",$maxx,$miny,$maxz;
	printf OUTFILE "ATOM      4  DUD BOX     1    %8.3f%8.3f%8.3f\n",$minx,$miny,$maxz;
	printf OUTFILE "ATOM      5  DUE BOX     1    %8.3f%8.3f%8.3f\n",$minx,$maxy,$minz;
	printf OUTFILE "ATOM      6  DUF BOX     1    %8.3f%8.3f%8.3f\n",$maxx,$maxy,$minz;
	printf OUTFILE "ATOM      7  DUG BOX     1    %8.3f%8.3f%8.3f\n",$maxx,$maxy,$maxz;
	printf OUTFILE "ATOM      8  DUH BOX     1    %8.3f%8.3f%8.3f\n",$minx,$maxy,$maxz;
	print  OUTFILE "CONECT    1    2    4    5\n";
	print  OUTFILE "CONECT    2    1    3    6\n";
	print  OUTFILE "CONECT    3    2    4    7\n";
	print  OUTFILE "CONECT    4    1    3    8\n";
	print  OUTFILE "CONECT    5    1    6    8\n";
	print  OUTFILE "CONECT    6    2    5    7\n";
	print  OUTFILE "CONECT    7    3    6    8\n";
	print  OUTFILE "CONECT    8    4    5    7\n";

	close OUTFILE;

	return 1;
}

sub molecule_pdb_rename_backbone {

	#<
	#? Rename pdb backbone atoms
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol = $_[0];

	my $aaflag;
	
	label_aa_backbone($mol);

	# Rename amino acid atoms
	foreach my $atom (@{$mol->{ATOMS}}) {

		$aaflag = 1 if  $atom->{FG}{AA_N};
		$atom->{NAME} = 'N'  if $atom->{FG}{AA_N};
		$atom->{NAME} = 'O'  if $atom->{FG}{AA_O};
		$atom->{NAME} = 'C'  if $atom->{FG}{AA_C};
		$atom->{NAME} = 'CA' if $atom->{FG}{AA_CA};
		$atom->{NAME} = 'CB' if $atom->{FG}{AA_CB};
		if($atom->{ELEMENT} eq 'H') {

			my $con;
			($con) = connected($atom, $mol);

			$atom->{NAME} = 'HA' if $con->{FG}{AA_CA};
			$atom->{NAME} = 'HB' if $con->{FG}{AA_CB};
			$atom->{NAME} = 'HN' if $con->{FG}{AA_N};
		}
	}
	$mol->{SYBTYPE} = 'BIOPOLYMER' if $aaflag;

}

sub molecule_pdb_rename_h {

	#<
	#? Rename protein hydrogens to conform to pdb nomenclature. Molecule is sorted if
	#  hydrogens are out of order in file
	#; Requires: molecule,  aminoacid datafile (default amino_acid_atoms.dat), 
	#  flag to delete unnamed hydrogens
	#; Returns: number of unassigned hydrogens
	#>

	my $mol = $_[0];
	my $datafile = $_[1] || silico_msg('d', "\$datafile is not defined\n");
	my $del = $_[2];
	
	my $missing_count;
	my $sort;
	
	silico_msg('c', "\n", "Renaming hydrogen atoms using file \"$datafile\"\n");
	
	my ($conatoms, $conhydrogens, $reshash) = molecule_pdb_rename_h_dataread($datafile);

	for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {
		
		my $atom = $mol->{ATOMS}[$i];

		# Skip hydrogens
		next if $atom->{ELEMENT} eq "H";
		
		# Strip out spaces from name
		my $aname = $atom->{NAME};
		$aname =~ s/ //g;
		my $sname = $atom->{SUBNAME};
		$sname =~ s/ //g;

		# Create hash of heavy atom name
		my $hashname = "$sname $aname";

		# Loop over connected
		my $con_hcount = -1;
		foreach my $connum (@{$atom->{CONNECT}}) {
		
			my $con = $mol->{ATOMS}[$connum];
			
			# Next if connected atom is not H
			next if ($con->{ELEMENT} ne 'H');
			++$con_hcount;

			# Find if hydrogens are a long way from their attached atoms
			# 20 spaces seems like a reasonable distance
			$sort = 1 if abs($atom->{NUM} - $con->{NUM}) > 20;

			# Skip hydrogens that are not connected within the same
			# residue (ie there is an error in the bonding)
			if ($con->{SUBID} != $atom->{SUBID}) {
				my $warn = 
					sprintf "Skipping hydrogen %s %s %s connected to %s %s %s\n",
					$con->{NUM},$con->{NAME},$con->{SUBID},
					$atom->{NUM},$atom->{NAME},$atom->{SUBID};
				$warn .= "		Distance ".distance($atom, $con)."\n";
				$warn =~ s/ +/ /g;
				mol_warn($mol, $warn);
				next;
			}

			if (defined ${$conhydrogens->{$hashname}}[$con_hcount]) {

				# PDB atom name is known: rename hydrogen atom 
				$con->{NAME} = ${$conhydrogens->{$hashname}}[$con_hcount];
				if ($con->{NAME} !~ /^\d/) {
					$con->{NAME} = ' '.$con->{NAME};
				}

			} elsif (defined $reshash->{$sname}) {

				# PDB residue name is known but atomname is not known: Change name to 'HX'
				my $warn = sprintf "Renaming H %s %s %s %s joined to %s %s %s %s",
					$con->{NUM},$con->{NAME},$con->{SUBNAME},$con->{SUBID},
					$atom->{NUM},$atom->{NAME},$atom->{SUBNAME},$atom->{SUBID};
				$warn .= "to HX. Num con H: ".($con_hcount+1)."\n";
				$warn =~ s/ +/ /g;
				mol_warn($mol, $warn);
				
				$con->{NAME} = ' HX';
				if ($Silico::debug) {
					silico_msg('g', "Converting atom element to F\n");
					$con->{NAME} = ' FX';
					$con->{ELEMENT} = "F";
				}
				
				++$missing_count;
			} else {

				# Residue name is not known
				silico_msg('n', "Ignoring hydrogen  $con->{NUM}, $con->{NAME}, $con->{SUBNAME} in unknown residue '$sname'\n");
			}
		}
	}

	# Delete any hydrogens we have called 'HX' if flag is set
	if ($del) {
		for (my $i=0; $i < $mol->{NUMATOMS}; ++$i) {

			my $atom = $mol->{ATOMS}[$i];
			if ($atom->{NAME} eq " HX") {
				
				my $warn = sprintf("Deleting atom  %s %s %s %s\n",
					$atom->{NUM},$atom->{NAME},$atom->{SUBNAME},$atom->{SUBID});
				$warn =~ s/ +/ /g;
				mol_warn($mol, $warn);
				molecule_delete_atom($mol, $i);
			}
		}
		molecule_pack ($mol);
	}

	if ($sort) {
		silico_msg('w', "Hydrogen atoms were out of order.  The atoms have been sorted\n");
		pdb_molecule_sort_atoms($mol);
	}

	return $missing_count;
}

sub molecule_pdb_rename_h_dataread {

	#<
	#? Read in data for molecule_pdb_rename_h
	#; Requires: datafile
	#; Returns: array, array, hash of residues
	#>
	
	my $datafile = $_[0];
	
	my $conatoms;
	my $conhydrogens;
	my $reshash;
	
	$datafile = $Silico::data_dir.$datafile;
	
	# Read in datafile of amino acid connectivities
	open (DATAFILE, $datafile) || silico_msg('d', "Could not open file $datafile for reading!\n");
	while (<DATAFILE>) {

		# Skip comments and blank lines
		next if /^#/;
		next if !/\w/;

		my @f = split;

		my @array1 = ();
		for (my $i=2; $i <= 5; ++$i) {
			last if $f[$i] eq '-';
			push @array1,  $f[$i];
		}

		my @array2 = ();
		for (my $i=2; $i <= 5; ++$i) {
			last if  $f[$i] eq '-';
			if ($f[$i] =~ /^H/ || $f[$i] =~ /^\dH/) {
				push @array2,  $f[$i];
			}
		}

		# Make a hash of residue names
		$reshash->{$f[0]} = 1;

		# Make array of residue name, atom name
		my $hashname = "$f[0] $f[1]";
		@{$conatoms->{$hashname}} = @array1;  # conatoms is not used in this version of the program
		@{$conhydrogens->{$hashname}} = @array2;
	}
	
	close DATAFILE;
	
	return ($conatoms, $conhydrogens, $reshash);
}

sub molecule_delete_alternate_atoms {

	#<
	#? Remove alternate atoms from a pdb file
	#; Requires: Molecule, aternate atom to RETAIN (taken from get flag or default value A)
	#; Returns: Molecule, modified and renumbered
	#>

	my $mol = $_[0];
	my $alt = $_[1] || get_flag("delete-alt-except", 'l') || 'A';

	my $a;
	my $delcount = 0;
	my $h;

	my $i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
		
		if ($atom->{ALT}) {
		
			$a = $atom->{ALT} || '';
			$a =~ s/ //;

			if ($a && ($a ne $alt)) {
				molecule_delete_atom($mol, $i);
				++$delcount;
				++$h->{$a};
			}
		}
		++$i;
	}

	molecule_pack ($mol);
	my @k = keys(%$h);
	if ($delcount) {
		silico_msg('n', "Molecule \"$mol->{NAME}\": Removed $delcount atoms in alternate position(s): @k\n");
	}
}

sub pdb_get_prot {

        #<
        #? Simple routine to extract pdb protein atoms
	#; Requires: molecule 
        #>

        my $mol = $_[0];

	my $newmol;

	$newmol = deep_copy($mol);
	$newmol->{ATOMS} = ();
	$newmol->{NUMATOMS} = 0;
	$newmol->{NUMBONDS} = 0;

	foreach my $atom (atoms($mol)) {

		next if $atom->{LABEL} !~ /^ATOM/;
		push @{$newmol->{ATOMS}}, $atom;
		++$newmol->{NUMATOMS};
	}

	return $newmol;
}



sub pdb_get_ligands {

        #<
        #? Extract ligands from pdb file
        #. Only the first example of each ligand  residue type is extracted
	#; Requires: molecule, min heavy atoms, max heavy atoms, min carbon atoms
        #>

        my $mol = $_[0];
        my $min  = $_[1] || get_flag("min-heavy", 'l') || 10;;  # Note minimum number of HEAVY ATOMS
        my $max  = $_[2] || get_flag("max-heavy", 'l') || 60;  # Note maximum number of HEAVY ATOMS
        my $minc = $_[3] || get_flag("min-carbons", 'l') || 5;  # Note maximum number of CARBON ATOMS

        my $ligands;
        my %sublist;

        # Druglike atoms
        my @atomlist = qw (C H N O S P F Cl Br I);

        my $l = molecule_split($mol);

        LIG: foreach my $lig (@$l) {

                # Remove non HETATM residues
                next if $lig->{ATOMS}[0]{LABEL} =~ /^ATOM/;
                my $h = mol_count_heavy_atoms($lig);

                # Remove ligands that are too large or too small
                next if $h < $min;
                next if $h > $max;
                $lig->{HEAVY} = $h;

                # Retain only first instance of each ligand;
                my $subname = $lig->{ATOMS}[0]{SUBNAME};
                if ($sublist{$subname}) {
                        silico_msg ('c', "Duplicate: Skipping duplicate ligand $subname\n");
                        next;
                }
                ++$sublist{$subname};

                # Add hydrogens
                my @g = mol_add_hydrogens($lig);
                print "$lig->{ATOMS}[0]{SUBNAME}. Added $g[0] hydrogens. Removed $g[1] hydrogens\n";

                my $f = molecule_formula($lig);
                $lig->{FORMULA} = formula_string($f);

                # Check that all atoms are in allowed list
                KEY: foreach my $key (keys %$f) {
                        foreach (@atomlist) {
                                next KEY if $key eq $_;
                        }
                        silico_msg ('c', "Nondruglike: Disallowed atom type $key\n");
                        next LIG;
                }

                if (!$f->{C} || $f->{C} < $minc) {
                        silico_msg ('c', "Nondruglike: Too few carbon atoms: $f\n");
                        next;
                }

                my $mw = molecule_mw($lig);
                $lig->{MW} = sprintf "%.1f", $mw;

                # Remove comments
                undef $lig->{SDF_DATA};

                # Name
                $lig->{NAME} = $subname;

                # May not be necessary
                pdb_molecule_sort_atoms($lig);
                molecule_renumber($lig);

                push @$ligands, $lig ;
        }

        return $ligands;
}

#
# PDBX format
#


sub write_pdbx {
	
	#<
	#? Write a pdb file
	#; Requires: Ensemble or molecule, filename, $options.
	#; Returns: Zero if file open failed otherwise returns 1.
	#; Options: BOND NOBOND QUIET
	#>

	my $molecules = ens($_[0]);
	my $outfile = $_[1];
	my $options = $_[2] || '';

	my $molcount = 0;

	print("warning: Subroutine is incomplete\n");

	$options = uc $options if $options;
	$outfile = get_ofilebase($molecules).".pdbx" if !$outfile;
	
	my $fr = open_molfile(undef, $outfile, undef, 'write', $options);
	if (!defined $fr) {
		silico_msg('e', "Could not open $outfile");
		return undef;
	}
	my $FH = $fr->{FH};
	
	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Writing pdbx file: $outfile");
		silico_msg('c', " using $options") if $options;
		silico_msg('c', "\n");
	}
	
	# Write out molecules
	foreach my $mol (@$molecules) {
	
		die if !defined $mol;
	
		++$molcount;
		if (!defined $FH) {
			silico_msg ('e',"Attempt to write to unopened filehandle. Perhaps you don't have write permission?\n");
			return undef;
		}
		
		write_pdbx_molecule($fr, $mol, $options);
	}
	
	close_molfile($fr);
	
	return 1;
}

sub write_pdbx_molecule {

	#<
	#? Write a single record from a pdb file.
	#; Options: BOND NOBOND 
	#; Requires: file record, molecule, options
	#. See read_pdb for general description
	#; Returns: 1 or undef if failed
	#>

	my $fr = $_[0]  || silico_msg('d', "File record not defined\n");
	my $mol = $_[1] || silico_msg('d', "Molecule not defined\n");
	my $options = $_[2] || '';
	
	my $FH = $fr->{FH};
	
	++$fr->{MOL_WRITE_COUNT};
	
	# Set standard output options based on flags
	make_atom_names_unique($mol) if get_flag('unique-atom-names', 'l');
	molecule_check_and_fix_connectivity($mol, $options) if ($options =~ /\bBOND\b/);
	
	if (!$mol->{PDBX_NAME}) {
		$mol->{PDBX_NAME} = $mol->{NAME} || "Mol_$fr->{MOL_WRITE_COUNT}";
		$mol->{PDBX_NAME} =~ s/\s/_/g;
		$mol->{PDBX_NAME} =~ s/^_*//g;
		$mol->{PDBX_NAME} =~ s/_*$//g;
	}
	
	print "name: $mol->{PDBX_NAME}\n";
	
	print $FH "# generated by Silico\n#\n";
	print $FH "data_$mol->{PDBX_NAME}\n#\n";
	print $FH "_entry.id $mol->{PDBX_NAME}\n#\n";
	print $FH "\n#\n";
	
	write_pdbx_cryst($mol, $FH);
	
	print $FH "#\n";
	print $FH "loop_\n";
	print $FH 
"#
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.pdbx_PDB_model_num\n";
	
	# Print atom records
	my $i = 0;
	foreach my $atom (@{$mol->{ATOMS}}) {
	
		my @f;
		++$i;
			
		$f[0] = $atom->{LABEL};		 #Record type - _atom_site.group_PDB
		$f[1] = $atom->{NUM};		 #Atom number - _atom_site.id
		$f[2] = $atom->{NAME};		 #Atom name - _atom_site.label_atom_id
		$f[3] = $atom->{ALT};		 #Alternate location indicator  - _atom_site.label_alt_id
		$f[4] = $atom->{SUBNAME};	 #Residue name
		$f[5] = $atom->{CHAIN};		 #Chain
		$f[6] = $atom->{SUBID};      	 #Residue number
		$f[7] = $atom->{X};		 #x _atom_site.Cartn_x
		$f[8] = $atom->{Y};		 #y _atom_site.Cartn_y
		$f[9] = $atom->{Z};		 #z _atom_site.Cartn_z
		$f[10] = $atom->{OCC};		 #Occupancy _atom_site.occupancy
		$f[11] = $atom->{TEMP};		 #Temperature Factor _atom_site.B_iso_or_equiv
		$f[12] = uc ($atom->{ELEMENT});  #Element - _atom_site.type_symbol
		
		
		# PDB formal charge format is 1+, 2-, etc _atom_site.pdbx_formal_charge
		if (defined $atom->{FORMAL_CHARGE} && $atom->{FORMAL_CHARGE}) {
			if ($atom->{FORMAL_CHARGE} < 0) {
				$f[13] = -$atom->{FORMAL_CHARGE}."-";
			} 
			if ($atom->{FORMAL_CHARGE} > 0) {
				$f[13] = $atom->{FORMAL_CHARGE}."+";
			} 
		} 
		
		$f[14] = $fr->{MOL_WRITE_COUNT}+1; #Model number
		
		# Remove leading and trailing spaces
		foreach (@f) {
			$_ ||= '';
			s/^\s*//;
			s/\s*$//;
		}
		
		$f[0] ||= "ATOM";
		$f[1] ||= $i;
		$f[2] ||= 'X$i';
		$f[3] ||= '.';
		$f[4] ||=  'UNK';
		$f[5] ||=  'A';
		$f[6] ||=  1;
		$f[7] ||=  0.0;
		$f[8] ||=  0.0;
		$f[9] ||=  0.0;
		$f[10] ||=  0.0;
		$f[11] ||=  0.0;
		$f[12] ||=  '.';
		$f[13] ||=  0;
		$f[14] ||=  '.';
		
		if ($atom->{NUM} == 1) {
			print $FH "#\n#Label     AtNum AtName Alt  SubNm Chn  SubID        X            Y            Z    Occ   Temp Elmt FlChg   Mdl\n#\n";
		}
		
		#           0    1    2    3    4   5   6  7      8      9      10    11    12   13    14   
		printf $FH "%6s  %8d  %5s  %2s  %5s  %2s  %d %12.3f %12.3f %12.3f %6.2f %6.2f   %2s    %2s    %2s\n",@f;
	}
	
	my $dupflag; # Need to fix this up.
	# Do CONECT records
	if (0 || $options !~ /\bNOCONECT\b/ && $mol->{NUMATOMS} <= 99999) {
	
		my $i = 0;
		foreach my $atom (@{$mol->{ATOMS}}) {
		
			++$i;
				
			# Skip atoms with no bonds
			next if ($#{$atom->{CONNECT}} < 0);
			
			printf $FH "CONECT%5s", $atom->{NUM};
			
			my $bo;
			my $j = -1;
			my $c = 0;
			foreach my $connum (@{$atom->{CONNECT}}) {

				++$j;

				if (!defined $mol->{ATOMS}[$connum]) {
					silico_msg('e', "Atom number $i is connected to a nonexistent atom (number $connum)!\n",
							"Skipping.\n");
					next;
				}

				if ($dupflag) {
					$bo = $atom->{BORDERS}[$j];
					if ($bo > 3 || $bo == 0 ) {
						mol_warn("Bondorder $bo converted to single");
						$bo = 1;
					}
				} else {
					$bo = 1;
				}

				foreach (1..$bo) {
					++$c;
					if ($c > 4) {
						print  $FH "\n";
						printf $FH "CONECT%5s", $atom->{NUM};			
						$c = 1;
					}
					printf $FH "%5s", $mol->{ATOMS}[$connum]{NUM};		
				}
			}
			print $FH "\n";
		}
	}

	if ($options =~ /\bMODEL\b/) {
		print $FH "ENDMDL\n";	
	} else {
		print $FH "END\n"; # To make Xplor happy!
	}
	
	mol_print_warnings($mol);
	#molecule_printout($mol);
	#die;
	
	return 1;
}

sub write_pdbx_cryst {

	#<
	#? Write PDBX crystal record
	#; Requires: molecule, filehandle
	#; Returns: nothing
	#>

	my $mol = $_[0];
	my $FH = $_[1];

	# Print out unit cell data
	return if !defined $mol->{CELL}[0];
	
	# In the printout, the last value (after $space_group) is the 
	# Z-value. This is the number of instances of the most numerous
	# chain in each unit cell (as distinct from each asymmetric unit).
	# Consequently, it will depend on the space group as well as on
	# the stoichiometry. If it is not set elsewhere, the value 1 will be used.
	
	if ( 	defined $mol->{CELL}[0] &&
		defined $mol->{CELL}[1] &&
		defined $mol->{CELL}[2] ) {

		if (	!defined $mol->{CELL_ALPHA} || 
			!defined $mol->{CELL_ALPHA} ||
			!defined $mol->{CELL_ALPHA}) {
			
			$mol->{CELL_ALPHA} = 90.0;
			$mol->{CELL_BETA}  = 90.0;
			$mol->{CELL_GAMMA} = 90.0;
			
			mol_warn("Setting space group to P 1");

			$mol->{SPACE_GROUP} = 'P 1';
		}
		
		$mol->{CHAINS_PER_CELL} ||= 1;
		
		my $space_group = space_group_silico_to_pdb($mol->{SPACE_GROUP});
		
		if (!defined $space_group) {
			mol_warn($mol,"Setting space group to P 1");
			$space_group = 'P 1';
		}
		
		print $FH, "#
			_cell.entry_id $mol->{PDBX_NAME}
			_cell.length_a $mol->{CELL}[0]
			_cell.length_b $mol->{CELL}[1]
			_cell.length_c $mol->{CELL}[2]
			_cell.length_alpha $mol->{CELL_ALPHA}
			_cell.length_beta $mol->{CELL_BETA}
			_cell.length_gamma $mol->{CELL_GAMMA}
			_symmetry.space_group_name_H-M '$space_group'\n";
	} else { 
		
		 my $warn = "Bad unit cell data!\n".
			"\ta: ".($mol->{CELL}[0] || "UNDEF")." ".
			"b: ".($mol->{CELL}[1] || "UNDEF")." ".
			"c: ".($mol->{CELL}[2] || "UNDEF")." ".
			"alpha: ".($mol->{CELL_ALPHA} || "UNDEF")." ".
			"beta: ".($mol->{CELL_BETA} || "UNDEF")." ".
			"gamma: ".($mol->{CELL_GAMMA} || "UNDEF")."\n".
			"\tUnit cell data has been omitted";
			
			mol_warn($mol, $warn);
	}
}

return 1;
