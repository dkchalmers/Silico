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
#! silico_io.pm
#? Silico general file handling routines and routines for getting
#  arguments and flags from the command line.
#. $Revision: 1.15-74-g1e475b1 $
#>

use strict;
package Silico;

sub open_file {

	#<
	#? General subroutine to open molecule files for reading.
	#. Note open_file and close_file are not designed to handle concurrently open files
	#. Example:
	#. if (!open_file(*INFILE, $infile)) {
	#.	print "Error (read_mol2): Can not open file $infile for reading!\n";
	#.	return undef;
	#. }
	#. Can handle compressed files on UNIX systems.
	#; Requires: filehandle, filename, 'single' flag to prevent opening files that are already open
	#; Returns: true if file opened ok or if the file was already open and the 
	#  'single' flag was set
	#>

	*FH = $_[0];
	my $infile = $_[1];
	my $single = $_[2]; # Keep track of open files

	use vars qw(
		$mae_prev_molecule 
		$silico_current_fh 
		$silico_current_filename
		$silico_infile_tokens 
		$silico_molcount 
		$silico_readcount 
		$silico_mol2_comments
		$silico_pdb_head
		$silico_pdb_sdata
	);
	
	if (!$infile) {
		silico_msg('e', "No input file was specified!\n");
		return undef;
	}

	use Config;
	
	# Keep track of open files and do
	# nothing if file is already open
	if ($single) {
		
		if ($silico_current_fh && $infile eq $silico_current_filename) {
			silico_msg('g', "Continuing read from file $infile\n");
			return 1;
		}
		if ($silico_current_fh && $infile ne $silico_current_filename) {
			silico_msg('w', "Opening file $infile although file $silico_current_filename was already open!\n");
		}
	}
	
	silico_msg('g', "Opening file $infile\n");

	$silico_current_fh = *FH;

	# Reset globals for the various read_* routines
	undef $mae_prev_molecule;
	undef $silico_infile_tokens;
	undef $silico_pdb_head;
	undef $silico_pdb_sdata;
	@$silico_mol2_comments = ();
	$silico_molcount = -1;
	$silico_readcount = 0;
	$silico_current_filename = $infile;

	if ($infile =~ /\.Z$/) {
		open (FH, "zcat $infile |") || return undef;
	} elsif ($infile =~ /\.gz$/ || $infile =~/.maegz$/) {
		open (FH, "gunzip -c $infile |") || return undef;
	} elsif ($infile =~ /\.bz2$/) {
		open (FH, "bunzip2 -c $infile |") || return undef;
	} else {
		open (FH, $infile) || return undef;
	}
}

sub close_file {

	#<
	#? Close the current input file opened using open_file
	#; Requires: nothing
	#; Returns: nothing
	#>
	
	silico_msg('g', "Closing file\n");

	use vars qw($silico_current_fh);

	if (!defined  $silico_current_fh) {
		silico_msg('e', "Call to close file with no open file\n");
		carp();
		return;
	}
	
	*FH = $silico_current_fh;
	close(FH);

	undef $silico_current_fh;
	undef $silico_current_filename;
}

sub read_mol_any {
	
	#<
	#? Read in any supported molecule file type.
	#. Uses the extension to determine file type.
	#. Reads following extensions:
	#, .arc (Mopac archive)
	#, .gro (Gromacs)
	#, .dat/.out (will determine the difference between Macromdel, Gaussian and Mopac files)
	#, .mae,.maegz Schroedinger maestro formats
	#, .mol2 (Sybyl)
	#, .mrk  (Merck format used with MMFF forcefield in CHARMm)
	#, .rtf (CHARMm topology definition format)
	#, .sdf .mol .rdf (MDL format)
	#, .pdb .ent .coor (text format from NAMD)
	#, .psf (CHARMm connectivity file) format
	#, .xyz (Tinker or other xyz) 
	#. Will also handle gzipped (.gz), bzip2'd (.bz2) and compressed (.Z) versions of these
	#  files
	#; Requires: file name. Optionally can pass start structure, max molecules, and option
	#  string to read routines although not all routines support 'start structure'
	#  and 'max molecules'
	#; Returns: undef if file open failed.
	#>

	my $name = $_[0];
	my $start = $_[1];
	my $max_molecules = $_[2];
	my $options = $_[3];
	
	my $mopac = 0;
	my $line;
	my $sname;
	my $success;

	return undef if (!defined $name);
	
	$sname = $name;
	$sname =~ s/\.Z$//;
	$sname =~ s/\.gz$//;
	$sname =~ s/\.bz2$//;

	# pdb
	if ($sname =~ /\.pdb$|\.ent$|.coor$/) {
		return (read_pdb ($name, $start, $max_molecules, $options));
	}
	# mol2
	if ($sname =~ /\.mol2$/) {
		return (read_mol2 ($name, $start, $max_molecules, $options));
	}
	# Merck
	if ($sname =~ /\.mrk$/) {
		require silico_merck;
		return (read_merck ($name, undef, undef, $options));
	}
	# Silico meta molecule file (mmm)
	if ($sname =~ /\.mmm$/) {
		require silico_mmm;
		return (read_mmm ($name, undef, undef, $options));
	}
	# sdf
	if ($sname =~ /\.sd$|\.sdf$|\.mol$/) {
		return (read_sdf ($name, $start, $max_molecules, $options));
	}
	# rdf
	if ($sname =~ /\.rdf$/) {
		return (read_rdf ($name, $start, $max_molecules, $options));
	}
	# mopac arc
	if ($sname =~ /\.arc$/) {
		require silico_mopac;
		return (read_mopac_cart ($name, undef, undef, $options));
	}
	# xyz
	if ($sname =~ /\.xyz_*[0-9]*$/) {
		require silico_xyz;
		return (read_xyz ($name, undef, undef, $options));
	}
	# Gromacs gro
	if ($sname =~ /\.gro$/) {
		require silico_gromacs;
		return (read_gromacs_gro ($name, undef, undef, $options));
	}
	# Gromacs itp
	if ($sname =~ /\.itp$/) {
		require silico_gromacs;
		return (read_gromacs_itp($name, undef, undef, $options));
	}
	# Gromacs top
	if ($sname =~ /.top$/) {
		require silico_gromacs;
		return (read_gromacs_top($name, undef, undef, $options));
	}
	# Gromacs trr, trj, xtc or cpt formats using trjconv wrapper script
	# Structure must be supplied with the -S or --struct flag
	if ($sname =~ /\.trr$/ || $sname =~ /\.trj$/  ||  $sname =~ /\.xtc$/ ||  $sname =~ /\.cpt$/ ) {
		require silico_gromacs;
		return (read_gromacs_trr ($name, undef, undef, $options));
	}
	if ($sname =~ /\.rtf$/) {
		require silico_charmm;
		return (read_charmm_rtf ($name, undef, undef, $options));
	}
	if ($sname =~ /\.psf$/) {
		require silico_charmm;
		return (read_charmm_psf ($name, undef, undef, $options));
	}
	if ($sname =~ /\.mae$/ || $name=~ /\.maegz$/) {
		require silico_mmod;
		return (read_maestro ($name, undef, $max_molecules, $options));
	}
	if ($sname =~ /\.cml$/ || $sname =~ /\.xml$/) {
		require silico_cml;
		return (read_cml ($name));
	}
	if ($sname =~ /\.smi$/) {
		require silico_smiles;
		return (read_smiles ($name));
	}
	if ($sname =~ /\.tab$/) {
		require silico_smiles;
		return (read_smiles ($name, undef, undef, "DTAB"));
	}
	if ($sname =~ /\.csv$/) {
		require silico_smiles;
		return (read_smiles ($name, undef, undef, "DCOMMA"));
	}
	# Macromodel/mopac/divcon dat/out
	# We must read the file because we can't tell by extension
	if ($sname =~ /\.dat$|\.out$/) {

		$success = open (IO, "$name");
		if (!$success) {
			silico_msg('e', "Can not open file $name for reading!\n");
			return undef;
		}

		# Read first 3 lines
		$line = '';
		for (1..3) {
			$line .= <IO>;
			chomp $line;
		}

		close IO;

		# Mopac dat file
		$mopac = 1 if ($line =~ /AM1|MNDO|PM3|am1|mndo|pm3/);
		# Mopac out file
		$mopac = 1 if ($line =~ /MOPAC/);
		$mopac = 1 if ($line =~ /SEILER/);
		
		if ($mopac) {
			silico_msg('n', "Assuming $name is a mopac .out file\n");
			require silico_mopac;
			return (read_mopac_cart ($name, undef, undef, $options));
		}

		# Divcon out file
		if ($line =~ /DIVCON/) {
			require silico_mopac;
			return (read_divcon_cart ($name, undef, undef, $options));
		}
		
		# Quick out file
		if ($line =~ /^  \*{10}/) {
			require silico_quick;
			return (read_quick_out ($name, undef, undef, $options));
		}

		silico_msg('n', "Assuming $name is a macromodel file\n");
		require silico_mmod;
		return (read_mmod ($name, $start, $max_molecules, $options));
	}
	
	# Read in the first line to decide what to do
	if (!open_file(*IO, $name)) {
		silico_msg('e', "Can not open file $name for reading!\n");
		return undef;
	}

	silico_msg('n', "Determining file type by reading first line\n");
	
	$line = <IO>;
	close IO;

	if (!defined $line) {
		silico_msg('e', "Error reading file \"$name\"!\n");
		return undef;
	}
	
	# Maybe PDB
	if ($line=~ /^HEADER/ || $line=~ /^ATOM/ || $line=~ /^HETATM/) {
		return (read_pdb ($name, undef, undef, $options));
	}
	
	# Maybe Gamess?
	if ($line=~ 'GAMESS') {
		require silico_gamess;
		return (read_gamess_cart ($name, undef, undef, $options));
	}
	
	# Maybe Gaussian?
	if ($line=~ 'Gaussian') {
		require silico_gaussian;
		return (read_gaussian_cart ($name, undef, undef, $options));
	}	
	
	# Maybe Sparky reslist?
	if ($line=~ 'Group   Atom  Nuc    Shift   SDev  Assignments') {
		require silico_nmr;
		return (read_sparky_reslist ($name, undef, undef, $options));
	}
	
	# Maybe Sparky peaklist?
	if ($line=~ /(\d)+\s+peaks/) {
		require silico_nmr;
		return (read_sparky_peaklist ($name, undef, undef, $options));
	}
	
	# Don't know
	silico_msg('e', "File $name is of an unknown type!\n",
			"Skipping.\n");
	return undef;
	
}


sub write_mol_any {
	
	#<
	#? Write any supported molecule file type.
	#; Requires: ensemble (or molecule), filebase (determined by subroutine 'obase' if not specified),
	#  format (determined by subroutine 'get_oformat' if not specified), (optional) options string
	#; Returns: undef if file open failed or other problem
	#>

	my $ensemble = ens($_[0]);
	my $filebase = $_[1];
	my $format = $_[2];
	my $options = $_[3];
	my $append = $_[4];

	my $outfile;
	
	if (!defined $ensemble->[0]) {
		silico_msg('e', "Molecules not defined. No file written.\n");
		use Carp;
		Carp::cluck if $Silico::debug;
		return undef;
	}

	# Modify filename and format based on any flags supplied
	($filebase, $format) = get_ofilename2($ensemble, $filebase, $format, $append);

	$outfile = "$filebase.$format";

	# cml
	if ($format eq 'cml') {
		require silico_cml;
		return (write_cml ($ensemble, $outfile, $options));
	}
	# pdb
	if ($format eq 'pdb' || $format eq 'ent' || $format eq 'coor' || $format eq 'pdbqt') {
		return (write_pdb ($ensemble, $outfile, $options));
	}
	# pdbx
	if ($format eq 'pdbx') {
		return (write_pdbx ($ensemble, $outfile, $options));
	}
	# mol2
	if ($format eq 'mol2') {
		return (write_mol2 ($ensemble, $outfile, $options));
	}
	# Merck
	if ($format eq 'mrk') {
		require silico_merck;
		return (write_merck ($ensemble, $outfile, $options));
	}
	# sdf
	if ($format eq 'sdf' || $format eq 'sdf_noparse' || $format eq 'mol') {
		return (write_sdf ($ensemble, $outfile, $options));
	}
	# smiles
	if ($format eq 'smi' ) {
		return (write_smiles ($ensemble, $outfile, $options));
	}
	if ($format eq 'tab' ) {
		$options ||= '';
		return (write_smiles ($ensemble, $outfile, "$options DTAB"));
	}
	if ($format eq 'csv' ) {
		$options ||= '';
		return (write_smiles ($ensemble, $outfile, "$options DCOMMA"));
	}
	# tinker xyz
	if ($format eq 'tinker') {
		require silico_tinker;
		$outfile = "$filebase.xyz";
		return (write_tinker_xyz ($ensemble, $outfile, $options));
	}
	if ($format eq 'mopac') {
		require silico_mopac;
		$outfile = "$filebase.dat";
		return (write_mopac_cart ($ensemble, $outfile, $options));
	}
	if ($format eq 'divcon') {
		require silico_mopac;
		return (write_divcon_cart ($ensemble, $outfile, $options));
	}
	if ($format eq 'mmod') {
		require silico_mmod;
		$outfile = "$filebase.dat";
		return (write_mmod ($ensemble, $outfile, $options));
	}
	if ($format eq 'mae') {
		require silico_mmod;
		return (write_mae ($ensemble, $outfile, $options));
	}
	if ($format eq 'maegz') {
		require silico_mmod;
		return (write_mae ($ensemble, $outfile, ($options||'')." GZ"));
	}
	if ($format eq 'gro') {
		require silico_gromacs;
		return (write_gromacs_gro ($ensemble, $outfile, $options));
	}
	if ($format eq 'itp') {
		require silico_gromacs;
		return (write_gromacs_itp ($ensemble, $outfile, $options));
	}
	if ($format eq 'top') {
		require silico_gromacs;
		return (write_gromacs_top ($ensemble, $outfile, $options));
	}
	if ($format eq 'rtf') {
		require silico_charmm;
		return (write_charmm_rtf ($ensemble, $outfile, $options));
	}
	if ($format eq 'psf') {
		require silico_charmm;
		return (write_charmm_psf ($ensemble, $outfile, $options));
	}
	if ($format eq 'xyz') {
		require silico_xyz;
		return (write_xyz ($ensemble, $outfile, $options));
	}
	if ($format eq 'zmat') {
		require silico_gaussian;
		return (write_gaussian_zmatrix ($ensemble, $outfile, $options));
	}

	# Don't know
	silico_msg('w', "Unknown file format '$format' requested.  Using mol2 format\n");
	require silico_mol2;
	return (write_mol2 ($ensemble, "$filebase.mol2", $options));
}

sub write_mol_any_check {

	#<
	#? Check that a specified output format is valid for write_mol_any
	#; Requires: format string
	#; Returns: 1 or 0
	#>
	
	my $output_format = $_[0];
	my @list = qw(ent gro mmod mol mol2 mopac mrk pdb psf sdf tinker);
	
	return 0 if !defined $output_format;
	
	foreach (@list) {
		return 1 if ($_ eq $output_format);
	}
	
	return 0;
}

sub check_grace_output_format {

	#<
	#? Checks an output format for compatibility with Grace
	#; Requires: Output format (text)
	#; Returns: Correctly formatted output format, extension
	#>
	
	my $format = $_[0];
	
	my $ext;
	my @flist;
	my $format_flag = 0;
	
	@flist = qw( ps eps mif svg pnm jpeg png metafile);
	
	$format = "jpeg" if $format eq "jpg";

	foreach (@flist) {$format_flag = 1 if $format eq $_}

	if ($format_flag == 0) {
		silico_msg('e', "Output file format $format is not supported by Grace!\n",
				"Please change to one of [".join(", ",@flist)."]\n");
		return undef;
	}
	
	if ($format eq 'metafile') {
		$format = 'Metafile';
		$ext = 'gmf';
	} else {
		if ($format eq 'jpeg') {
			$ext = 'jpg';
			$format = 'JPEG';
		} else {
			$ext = $format;
			if ($format eq 'ps') {
				$format = 'PostScript';
			} else {
				$format = uc($format);
			}
		}
	}
	
	return ($format, $ext);
}


##################################################################
#
# Routines for reading molecules one at a time
#
##################################################################

sub open_molfile {

	#<
        #? Open molecule file for reading one molecule at a time
        #. Uses the extension to determine file type.
        #; Requires: file record, file name, mode (read(default), write, append), format, option  string 
        #; Returns: file record or undef if open failed. 
        #; Options: GZ - gzip compression, BZ2 - bzip2 compression, QUIET - less output
        #>

	my $fr = $_[0];
        my $filename = $_[1];
	my $format = $_[2];
	my $mode = $_[3] || 'read';
        my $options = uc($_[4] || '');

	my $compress;
	my $fh;
	my $quiet;
	my $s = '';
	my $readwrite;
	
	use FileHandle;

	$format = lc $format if $format;
	$mode   = lc $mode if $mode;

  	$quiet = 1 if ($options =~ /\bQUIET\b/);
	
	# Close (and reopen) if filename or mode has changed
	# Note that you should check if $fr is defined before checking
	# hash elements, otherwise $fr WILL become defined

	if ((defined $fr->{FH} && defined $filename && defined $fr->{FILENAME} && $filename ne $fr->{FILENAME}) 
		|| (defined $fr->{FH} && defined $mode && defined $fr->{MODE} && $mode ne $fr->{MODE})) {
		
		silico_msg('g', "Closing file $fr->{FILENAME} and opening file $filename\n");
		close_molfile($fr);
		undef $fr;
	} 
	
	# Return existing filehandle without change
	if (defined $fr->{FH}) {
		return $fr;
	}

	$mode =~ s/\s//g if $mode;
	if ($mode eq 'read') {

		# Handle reading compressed files
		if ($filename =~ /\.Z$/) {
			$s = "zcat $filename |";
			$compress = 'compress';
		} elsif ($filename =~ /\.gz$/) {
			$s = "gunzip -c $filename |";
			$compress = 'gzip';
		} elsif ($filename =~/.maegz$/) {
			$s = "gunzip -c $filename -S .maegz |";
			$compress = 'maegz';
		} elsif ($filename =~ /\.bz2$/) {
			$s = "bunzip2 -c $filename |";
			$compress = 'bzip2';
		} else {
			$s = $filename;
		}

		$readwrite= 'reading';

	} elsif ($mode =~ 'write') {

		# Supplied flag
		my $c = get_lflag('compress') || '';

		if ($c eq 'none') {

			$filename =~ s/.maegz$/.mae/;
			$s = "> $filename";

		} elsif ($options =~ /\bGZ\b/ || $c eq 'gz' || $filename =~ /\.gz$/ || $filename =~ /\.maegz$/) {

			$filename =~ s/\.gz$//;
			if ($filename =~ /mae$/) {
				$filename =~ s/.mae$/.maegz/;
			} elsif ($filename =~ /maegz$/) { 
				$compress = 'maegz';
			} else {
				$filename .= ".gz";
				$compress = 'gzip';
			}
			$s = "| gzip -c > $filename";

		} elsif ($options =~ /\bBZ2\b/ || $c eq 'bz2' || $filename =~ /\.bz2$/) {

			$filename =~ s/\.bz2$//;
			$filename .= '.bz2';
			$s = "| bzip2 -c > $filename";
			$compress = 'bzip2';

		} else {
			$s = "> $filename";
		}
		$readwrite = 'writing';

	} elsif ($mode =~ 'append') {
		$s = ">> $filename";
		$readwrite = 'appending';
	} else {
		silico_msg 'd', "File mode '$mode' is not supported\n";
	}

	$fr->{FH} = $fh = new FileHandle;
	$fr->{FILENAME} = $filename;
	$_[1] = $filename; # Change the passed filename to match compression options
	$fr->{MOL_WRITE_COUNT} = -1;
	$fr->{MOL_READ_COUNT} = -1;
	$fr->{MODE} = $mode;
	$fr->{FORMAT} = $format || get_fileext($filename);
	$fr->{OPEN_COMMAND} = $s;
	$fr->{COMPRESS} = $compress;

	if (!$fh->open("$s")) {
		silico_msg('e', "Can not open file $filename for $readwrite\n");
		$fh->close();
		return undef;
	}

	return $fr;
}

sub mol_set_source_data {

	#< 
        #? Set up molecule source data such as SOURCE_FILE_TYPE, SOURCE_FILE_NAME, etc
        #; Requires: molecule, file_record, format (optional)
	#; Returns: nothing
	#>

	my $mol = $_[0] || croak();
	my $fr = $_[1] || croak();
	my $format = $_[2];

	$mol->{SOURCE_FILE_TYPE} = $format || $fr->{FORMAT};
	$mol->{SOURCE_FILE_NAME} = $fr->{FILENAME};
	$mol->{SOURCE_FILE_COMPRESS} = $fr->{COMPRESS};
	$mol->{NUM} = $fr->{MOL_READ_COUNT}+1;
}

sub read_mol_single {

	#< 
        #? Read a single molecule
        #. Uses the extension to determine file type.
        #; Requires: file record, filename (optional), options, start structure (optional), max molecules (optional)
	#; Returns: molecule or undef 
	#>

	#my $fr = $_[0];
	my $filename = $_[1];
	my $options = $_[2];
	my $start = $_[3] || get_flag('ss', 's');  
	my $max_molecules = $_[4] || get_flag('ms', 's'); 
	
	my $mol;
	my $fr;
	
	silico_msg ('d', "Neither file record nor filename provided\n") if !($_[0] || $filename);
		
	# Note that we must operate on the supplied variable '$_[0]'
	$_[0] = open_molfile($_[0], $filename, undef, 'read', $options);

	$fr = $_[0];
	
	return undef if $max_molecules && $fr->{MOL_READ_COUNT}+1 >= $max_molecules;
	$fr->{MOL_SKIP_COUNT} = 0 if $start;
	
	if (!defined $fr) {
		silico_msg('e', "Could not open file $filename\n");
		carp();
		$fr->{ERROR} = 1;
		return undef;
	}
	
	if (!defined $fr->{FORMAT}) {
		fr_printout($fr);
		silico_msg('d', "Format undefined\n");
	}

	while (1) {

		if ($fr->{FORMAT} eq 'mol2') {
			require silico_mol2;
			$mol =  read_mol2_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'pdb') {
			$mol =  read_pdb_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'pdbqt') {
			$mol =  read_pdb_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'gro') {
			require silico_gromacs;
			$mol =  read_gromacs_gro_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'sdf') {
			require silico_sdf;
			$mol =  read_sdf_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'mae' || $fr->{FORMAT} eq 'maegz') {
			require silico_mmod;
			$mol =  read_maestro_single($fr, $options); 
			
		} elsif ($fr->{FORMAT} eq 'smi') {
			require silico_smiles;
			$mol =  read_smiles_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'arc') {
			$mol =  read_mopac_arc_single($fr, $options);
	
		} elsif ($fr->{FORMAT} eq 'xyz') {
			require silico_xyz;
			$mol =  read_xyz_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'rtp') {
			require silico_gromacs;
			$mol =  read_rtp_single($fr, $options);

		} elsif ($fr->{FORMAT} eq 'in') {
			require silico_quick;
			$mol =  read_quick_in_single($fr, $options);
		
		} elsif ($fr->{FORMAT} eq 'out') {
			require silico_quick;
			$mol =  read_quick_out_single($fr, $options);
		
		
		} else {
			silico_msg('d', "Can not read file format '$fr->{FORMAT}' using single mode. Try using the --nosingle flag\n");
		}

		return undef if !defined $mol;

		if ($start && $fr->{MOL_SKIP_COUNT} < $start) {

			++$fr->{MOL_SKIP_COUNT};
			next;
		}

		#++$fr->{MOL_READ_COUNT}; # This is done in the called subroutines
		return $mol;
	}
}

sub write_mol_single {

	#< 
        #? Write a single molecule
        #; Requires: file_record, molecule, filename base (optional), file format (optional), filename append_string (optional), options
	#; Returns: 1 or undef 
	#>
	
	#my $fr = $_[0];
	my $mol = $_[1];
	my $filebase = $_[2];
	my $format = $_[3];
	my $append = $_[4];
	my $options = uc ($_[5] || '');
	
	my $ens;
	my $filename;
	my $fr;
	my @formats = qw (pdb sdf mol2 mae maegz gro xyz);
	my $mode = 'write';

	$mode = 'append' if $options =~/\bAPPEND\b/;
	
	$ens->[0] = $mol;

	# Set filebase and/or format if not provided
	if (!$filebase || !$format) {
	
		my ($filebase1, $format1) = get_ofilename2($ens, $filebase, $format, $append);
	
		$filebase = $filebase1 if (!$filebase || $append);
		$format = $format1 if !$format;
	}

	my $flag = 0;
	foreach (@formats) {
		next if $format ne $_;
		$flag = 1;
	}

	if ($flag == 0) {
		silico_msg('w', "Can not write format '$format'.  Using pdb\n");
		$format = 'pdb';
	}
	
	$filename = "$filebase.$format";
	
	# Open a new output file if none was supplied
	if (!defined $_[0]->{FH}) {
		# Note that we must operate on the supplied variable '$_[0]'
		$_[0] = open_molfile($fr, $filename, undef, $mode,  $options);
		silico_msg('d', "Could not open output file $filename\n") if !defined $_[0];
		#fr_printout($_[0]);
	}
	
	# $fr is defined here to ensure that we are operating on scalar from calling subroutine
	$fr = $_[0];
	
	if (!defined $fr) {
		carp();
		silico_msg('d', "File record undefined\n");
	}
	
	$fr->{FORMAT} = $format;
	# increment of MOL_WRITE_COUNT should be in write_single routines
	#++$fr->{MOL_WRITE_COUNT};

	my $val;
	if ($format eq 'pdb') {
		$val =  write_pdb_molecule($fr, $mol, $options);
	}
	if ($format eq 'mol2') {
		$val =  write_mol2_molecule($fr, $mol, $options);
	}
	if ($format eq 'gro') {
		require silico_gromacs;
		$val =  write_gromacs_gro_molecule($fr, $mol, $options);
	}
	if ($format eq 'sdf') {
		$val =  write_sdf_molecule($fr, $mol, $options);
	}
	if ($format eq 'mae') {
		$val =  write_mae_molecule($fr, $mol, $options);
	}
	if ($format eq 'maegz') {
		$val =  write_mae_molecule($fr, $mol, $options." GZ");
	}
	if ($format eq 'xyz') {
		require silico_xyz;
		$val =  write_xyz_molecule($fr, $mol, $options);
	}
	if ($format eq 'in') {
		require silico_quick;
		$val =  write_quick_in_single($fr, $mol, $options);
	}

	close_molfile ($fr) if $mode eq 'append';

	return $val;
}

sub close_molfile {

	#< 
        #? Close a molecule file that has been opened to read or write
        #; Requires: file record
	#>

	use FileHandle;
	my $fh = $_[0]->{FH};
	silico_msg('d', "Attempt to use undefined filehandle\n") if !$fh;
	$fh->close if $fh;
	undef $_[0];
}

sub silico_readline {

        #<
        #? Read a line from $fr->{FH}. Remove carriage returns and line feeds
	#. Lines can be stored in var $fr->{READLINE}
	#. $fr->{END} is set if the end of the file is reached.
        #; Requires: file record
	#; Returns: processed line
        #>

        my $fr = $_[0] || croak();

        my $line;
        my $FH = $fr->{FH};
	
	if (!defined $FH) {
		fr_printout($fr);
		carp();
		silico_msg('d', "Urk. Undefined filehandle\n");
	}

	# Stored lines if any
        if (defined $fr->{READLINE}[0]) {
                $line =  pop @{$fr->{READLINE}};
        } else {
                $line = <$FH>; 
        }
	
	# End of file
	if (!defined $line) {
		$fr->{END} = 1;
		#fr_printout($fr);
		return undef;
	}
                
	chomp $line;
        $line =~ s/\r$//; # Remove any trailing line feeds from pc text files
	
	#print "line '$line'\n";
	#fr_printout($fr);
	
	++$fr->{LINECOUNT};

        return $line;
}

sub silico_putline {

        #<
        #? Put a line back on 'the stack'
        #; Requires: line, file record
	#; Returns: nothing
        #>
	
	my $line = $_[0];
	my $fr = $_[1] || croak();
	
        push @{$fr->{READLINE}} ,$line;
	$fr->{LINECOUNT} ||= 0;
	--$fr->{LINECOUNT};
}

#
# Debuging
#

sub fr_printout {

	#< 
        #? Print out a silico file record for debugging purposes
        #; Requires: file record
	#>

	my $fr = $_[0];
	
	my $val;
	my @f;
	my @g;
	my @h;
	
	@f = caller(0);
	@g = caller(1);
	@h = caller(2);
	
	$f[1] =~ s/.*\///;
	
	print heading("File Record Printout: ".($fr || 'UNDEFINED')."\n");

	if ($g[1]) {
		$g[1] =~ s/.*\///;
		print "$g[3] : $f[1] line $f[2]\n";
	}
	if ($h[1]) {
		$h[1] =~ s/.*\///;
		print "$h[3] : $g[1] line $g[2]\n" if $h[3];
	}
	print "\n";
	
	return if !defined $fr;
	
	foreach (sort keys(%$fr)) {
		if (defined $_) {
			$val = $fr->{$_};
			$val = 'UNDEFINED' if (!defined $fr->{$_});
			printf "%-16s %-10s\n", $_, $val;
		}  else {
			print "UNDEFINED\n";
		}
	}
	print "\n";
}

##################################################################
#
#	AA sequence routines
#
##################################################################

sub read_seq_any {
  
	#<
	#? Read any sequence format
	#; Requires: filename, options
	#; Returns: sequences
	#>

	my $name = $_[0];
	my $options = $_[1];

	my $line;
	my $success;

	return undef if (!defined $name);

	$success = open (IO, "$name");
	if (!$success) {
		silico_msg('e', "Can not open file $name for reading!\n",
				"Skipping.\n");
		return undef;
	}

	$line = <IO>;
	close IO;

	# clustal
	if ($line =~ /^CLUSTAL/) {
		return (read_clustal ($name, $options));
	}

	# t-coffee score file
	if ($line =~ /^T-COFFEE,/) {
		return (read_tcoffee_score ($name, $options));
	}

	# Fasta
	if ($line =~ /^>/ || $name =~ /\.fasta/) {
		return (read_fasta($name, $options));
	}
	
	silico_msg('e', "Could not read sequence data from file $name.",
		"Available types are Clustal, T-Coffe score files, Fasta\n");

	return undef;
}

sub write_seq_any {
	
	#<
	#? Write any supported sequence file type.
	#; Requires: sequence ensemble, filebase, format
	#; Returns: undef if file open failed or other problem
	#>
	
	my $sequences = $_[0];
	my $filebase = $_[1];
	my $format = $_[2];
	my $options = $_[3] || '';
	
	my $outfile;
	
	if (!defined $sequences->[0]) {
		
		silico_msg('e', "Sequence ensemble not defined. No file written.\n");
		
		use Carp;
		Carp::cluck if $Silico::debug;
		
		return undef;
	}
	
	($filebase, $format) = get_ofilename2($sequences, $filebase, $format, undef, 'fasta');
	$outfile = "$filebase.$format";

	# clustal
	if ($format eq 'aln') {
		return (write_clustal ($sequences, $outfile, $options));
	}
	
	# fasta
	if ($format eq 'fasta') {
		return (write_fasta ($sequences, $outfile, $options));
	}
	
	# pir
	if ($format eq 'pir') {
		return (write_pir ($sequences, $outfile, $options));
	}
	
	# txt (space delimited text format suitable for excel
	if ($format eq 'txt') {
		$options .= "NUMRES=200 DELIMITER=' '";
		return (write_clustal ($sequences, $outfile, $options));
	}

	# No sensible format found
	silico_msg('w', "Unknown sequence output format '$format' requested. Using fasta format instead.\n");

	return (write_fasta ($sequences, "$filebase.fasta", $options));
}

sub write_seq_any_check {

	#<
	#? Check that a specified output format is valid for write_seq_any
	#; Requires: format string
	#; Returns: 1 or 0
	#>
	
	my $output_format = $_[0];
	my @list = qw(	aln fasta pir txt);
	
	return 0 if !defined $output_format;
	
	foreach (@list) {
		return 1 if ($_ eq $output_format);
	}
	return 0;
}

##################################################################
#
#	Input argument handling routines
#
##################################################################

sub get_std_output_options {

	#<
	#? Read a standard set of output options from the command line
	#. Deprecated subroutine;
	#; Requires: a list of options
	#; Returns: option string
	#>

	my $option_list = uc $_[0];
	
	my $options = '';

	$options .= 'NOSORT NOTYPE '	if ($option_list =~ 'FAST'		&& get_flags('fast', undef , 'Do not sort atoms'));
	$options .= 'NOSORT '		if ($option_list =~ 'NOSORT'		&& get_flags('nosort', undef , 'Do not sort atoms'));
	$options .= 'NOTYPE '		if ($option_list =~ 'NOTYPE'		&& get_flags('notype', undef, 'Do not type atoms'));
	$options .= ' BOND '		if ($option_list =~ 'BOND'		&& get_flags('bond', undef, 'Generate bonds'));
	$options .= 'NOBOND '		if ($option_list =~ 'NOBOND'		&& get_flags('nobond', undef, 'Do not generate any bonds'));
	$options .= 'NOCONECT '		if ($option_list =~ 'NOCONECT'		&& get_flags('noconect', undef, 'Do not write CONECT records'));

	return $options;
}

sub add_option {

	#<
	#? Add a text string to the list of command-line options
	#; Requires: Option to add (string), list of existing options (string)
	#; Returns: Modifies string provided to function
	#>
	
	my $option = uc ($_[0] ||'');
	my $options = uc ($_[1] ||'');
	
	# Get pointer to second argument
	my $addr = \$_[1];
	
	# Do nothing if the requested option is already being used
	return if ($options =~ /\b$option\b/);
	
	# Add a space on the end of the options string if the existing
	# string finishes in a non-whitespace character
	$options .= " " if ($options =~ /\S$/);
	
	# Add the requested option to the end of the options string
	$options .= "$option";
	
	# Change original variable
	${$addr} = $options
	
}

sub remove_option {

	#<
	#? Delete a text string from the list of command-line options
	#; Requires: Option to remove (string), list of existing options (string)
	#; Returns: Modifies string provided to function
	#>
	
	my $option = uc $_[0];
	my $options = uc $_[1];
	
	# Do nothing if the requested option is not in the options string
	return if ($options !~ /\b$option\b/);
	
	# Get pointer to second argument
	my $addr = \$_[1];
	
	
	# Remove the requested option from its positions in the list
	# (no matter how many times it appears)
	$options =~ s/\b$option\b//g;
	
	# Remove all spaces from the beginning
	$options =~ s/^\s+//;
	
	# Remove all spaces from the end
	$options =~ s/\s+$//;
	
	# Remove all extra spaces from the middle
	$options =~ s/\s{2,}/ /g;
	
	# Change original variable
	${$addr} = $options
}

sub setup_flags {

	#<
	#? Initiate flag by checking for menu or help options
	#>

	my $arg;
	my $flag;
	my $i;
	
	@Silico::Arguments = @ARGV;

	# Remember that we have done this
	$Silico::Flags_setup = 1;

	# Print out documentation and exit if help flag(s) are set
	my $eflag;
	foreach (@Silico::Arguments) {
		if ($_ eq '-help' || $_ eq '--help') {
			printdoc();
			$eflag = 1;
		}	
		if ($_ eq '--help-flags') {
			printdoc("$Silico::data_dir/stdflags.txt");
			$eflag = 1;
		} 
		if ($_ eq '--help-mopac') {
			printdoc("$Silico::data_dir/mopacflags.txt");
			$eflag = 1;
		}
		if ($_ eq '--help-asl') {
			printdoc("$Silico::data_dir/aslflags1.txt");
			printdoc("$Silico::data_dir/aslflags2.txt");
			$eflag = 1;
		}
	}
	exit if $eflag;
	
	# Find out if the 'menu' flag is set
	for ($i = 0; $i <= $#Silico::Arguments; ++$i) {

		$arg = $Silico::Arguments[$i];
		
		if ($arg eq '-m' || $arg eq '--menu') {
			
			$Silico::Menu = 1;
			
			splice (@Silico::Arguments, $i, 1);
			--$i;
		}
	}
}

sub get_default_flags_standard {

	#<
	#? Complete flag setup by setting silico-wide default flags
	#>
	
	print "## Silico standard flags ##\n" if $Silico::Menu;

	
		make_flag('o',		'output-format', 	"Output file format",	1);
		make_flag('O',		'output-file-name', 	"Output file name", 1);
		make_flag('compress', 	'output-file-compression',  	"Output file compression [gz, bz2, none]", 1);
		make_flag('nosort',	'nosort',		"Do not sort atoms");
		make_flag('notype',	'notype',		"Do not type atoms");
	my $fast = make_flag('fast',	'fast',			"Fast write. Equivalent to --nosort --notype");

		# Atoms
		make_flag(undef, 	'unique-atom-names', 	"Make molecule atom names unique");
	
		# Bonding
		make_flag('bond',	'bond',			"Generate bonds");
		make_flag('nobond',	'nobond',		"Do not generate bonds");
		make_flag('noconect',	'noconect',		"Do not generate CONECT records");
		
		# Structure reading
		make_flag('me',		'max-energy',		"Maximum relative energy for input molecule", 1, undef, undef, "decimal" );
		make_flag('ms',		'max-structures',	"Maximum number of structures to read", 1, undef, undef, "INT >=0");
		make_flag('ss',		'start-structure',	"Starting structure number", 1, undef, undef, "INT >=0");
		make_flag(undef,	'last-structure',	"Retain only last read structure");
		
		make_flag('stride', 	'stride',  		"Stride.  Read every nth structure", 1, undef, undef, "INT >=0");

		# Binsize for connecting atoms
		make_flag(undef, 	'binsize',		"Binsize for use in connecting atoms", 1, undef, undef, "DECIMAL > 0");
		
		# Gromacs specific flags
		make_flag(undef,	'coords',		"Structure input for reading gromacs trajectories (trr or cpt files)");
		
		# SDF files
		make_flag(undef, 	'sdf-2D',  		"Force SDF files to 2D format");
		make_flag(undef, 	'sdf-3D',  		"Force SDF files to 3D format");
		make_flag(undef, 	'sdf-squash-hcount',  	"Unset HCOUNT field in SDF files");
		
		# PDB files
		make_flag('longres', 	'pdb-long-resname',  	"Use 4 character residue names in pdb files");
		make_flag(undef, 	'pdb-model',  		"Write MODEL records between PDB structures (default END)");
		make_flag(undef,	'dupconect',		"Represent bond orders with duplicate bonds in pdb files");
		
	# Get non-molecule flags
	get_default_flags_minimal();
	
	# Die if mutually exclusive flags (-bond and -nobond) are used
	if (get_lflag('bond') && get_lflag('nobond')) {
	
		silico_msg('d', "-bond and -nobond options are mutually exclusive!\n",
				"Please review your command line options.\n");
	}
	
	if ($fast) {
		set_flag('nosort');
		set_flag('notype');
	}
	
	# Turn off Quiet and Very Quiet options if debugging is to be done
	if ($Silico::debug) {

		print "Entering debug mode.\n";

		if (get_lflag('quiet') || get_lflag('very-quiet')) {
			print "--quiet and --very-quiet options have been disabled.\n";
			delete_flag('quiet', 'l');
			delete_flag('very-quiet', 'l');
			$Silico::quiet =  0;
			$Silico::veryquiet = 0;
		}
		
		print "The output option \"QUIET\" will not be used.\n";
	}
}

sub get_default_flags_minimal {

	$Silico::quiet =	make_flag(undef, 	'quiet', 	"Suppress notes", 0, 0);
	$Silico::veryquiet =	make_flag(undef, 	'very-quiet', 	"Suppress notes and warnings", 0, 0);
	$Silico::TIMING	=  	make_flag('timing',	'timing',	"Timing mode", 0, 0);
	$Silico::debug	= 	make_flag('debug',	'debug',	"Debug mode", 0, 0);
	$Silico::warnings = 	make_flag('warnings',	'warnings',	"Print extended warning information", 0, 0);
	$Silico::memuse	= 	make_flag('memuse',	'memuse',	"Memory use mode", 0, 0);
	$Silico::temp_dir = 	make_flag(undef,	'tempdir',	"Temporary working directory", 1,);
	
	# Turn on Quiet option if Very Quiet is set
	if (get_lflag('very-quiet')) {
		set_flag('quiet');
		$Silico::quiet = 1;
	}

	#print "silico_warnings: $Silico::warnings\n";
	#print "silico_debug: $Silico::debug\n";
}

sub make_flag {

	#<
	#? Create a flag based on input options
	#; Use: make_flag (short_name, long_name, "description", expects_value, default_value, prompt, verify_string, 
	#  duplicate defnition allowed)
	#. 'short_name' is a short format switch starting with a single '-'
	#. 'long_name' is a long format switch starting with '--'
	#. 'description' is a text description of the flag
	#. 'expects_value' 0 - flag is a switch, 1 - expects scalar value, 2 - expects array (values separated by spaces).
	#. 'default_value' default value for flag (can be undef)
	#. 'option_to_prompt' (1 or 0) when set to 1 silico will prompt for data if it undefined
	#. 'verify_string' is used to check that the input data is of the correct type: integer,
	#  real number, string, etc.  See 'verify_input_data'
	#. 'duplicate_definititions_allowed'. This will cause make_flag to return an array of values for this flag
	#  This is required when we wish to set default values for a Silico standard flag
	#; Requires: Short name, Long name, Description, Expects value?, Default value, option to Prompt,
	#  verify string, duplicate definition allowed flag
	#. Note expects_value = 2 and duplicate_definitions_allowed returns an array
	#; Returns: value the flag is to take. 
	#>
	
	# Set the flag (hash pointer) up based on input values
	my $flag->{SHORT_NAME} = $_[0];
	$flag->{LONG_NAME} = $_[1];
	$flag->{DESCRIPTION} = $_[2];
	$flag->{HAS_VALUE} = $_[3];
	$flag->{DEFAULT} = $_[4];
	$flag->{PROMPT} = $_[5];
	my $verify_string = $_[6];
	$flag->{DUPLICATES_ALLOWED} = $_[7];

	# Run setup_flags if it has not already been done
	setup_flags() if (!defined $Silico::Flags_setup);
	
	# If the flag's short name is -m...
	if (defined $flag->{SHORT_NAME} && $flag->{SHORT_NAME} eq 'm') {
		silico_msg('d', "Attempt to make a flag with the reserved name -m!\n");
	}
	
	# If the flag's long name is --help...
	if (defined $flag->{LONG_NAME} && $flag->{LONG_NAME} eq 'help') {
		silico_msg('d', "Attempt to make a flag with the reserved name --help!\n");
	}

	# Check for duplicates
	foreach my $t (@$Silico::Flags) {

		next if $t->{DUPLICATES_ALLOWED};
	
		# "new" short name matches any existing flag's short name
		if (defined $flag->{SHORT_NAME} && defined $t->{SHORT_NAME} && $flag->{SHORT_NAME} eq $t->{SHORT_NAME}) {
			
			silico_msg('d',	"Attempt to make two flags with the same short name!\n",
					"Short name: -$flag->{SHORT_NAME}\n");
		}
		
		# "new" long name matches any existing flag's long name
		elsif (defined $flag->{LONG_NAME} && defined $t->{LONG_NAME} && $flag->{LONG_NAME} eq $t->{LONG_NAME}) {
			
			silico_msg('d',	"Attempt to make two flags with the same long name!\n",
					"Long name: --$flag->{LONG_NAME}\n");
		}
	}

	set_flag_from_arguments($flag);
	
	# If the current flag's value is not set, use the default, but only if
	# this flag is not to be prompted for.
	$flag->{VALUE} = $flag->{DEFAULT} if (!defined $flag->{VALUE} && !$flag->{PROMPT});
	
	#
	# There problem with the prompting, we are prompting for a value before executing the get_arguments subroutine
	# so we can not know if a value has been set.
	#

	# If the "prompt for every flag" (Menu) option is being used,
	# or this flag in particular must be prompted for...
	if ($Silico::Menu || ($flag->{PROMPT} && !defined $flag->{VALUE})) {

		# Print out a prompt.
		print "Enter value for flag: ";
		
		# Print out the flag's short name if defined (and a comma
		# if the long name is defined).
		if (defined $flag->{SHORT_NAME}) {
			print "-$flag->{SHORT_NAME}";
			print ", " if (defined $flag->{LONG_NAME});
		}
		
		# Print out the flag's long name if defined.
		if (defined $flag->{LONG_NAME}) {
			print "--$flag->{LONG_NAME}";
		}
		
		# Print out the flag's description if defined.
		print " $flag->{DESCRIPTION} " if (defined $flag->{DESCRIPTION});
		
		# If the flag expects to have a value...  
		if ($flag->{HAS_VALUE} && $flag->{HAS_VALUE} == 1) {
			
			# Print out what the value is currently (if e.g. it has
			# been set by the command line)
			if (defined $flag->{VALUE}) {
				print "[$flag->{VALUE}]";
			
			# If the value hasn't been set already, print out the
			# default value.
			} elsif (defined $flag->{DEFAULT}) {
				print "[$flag->{DEFAULT}]";
			
			# If neither Value nor Default has been set, print out
			# "not set".
			} else {
				print "[Not set]";
			}
			
		} elsif ($flag->{HAS_VALUE} && $flag->{HAS_VALUE} == 2) {
			
			# Print out what the value is currently (if e.g. it has
			# been set by the command line)
			if (defined $flag->{VALUE}) {
				my $string = @{$flag->{VALUE}};
				$string = substr(0, 20, $string) . "..." if length($string) > 23;
				print "[$string]";
			
			# If the value hasn't been set already, print out the
			# default value.
			} elsif (defined $flag->{DEFAULT}) {
				print "[$flag->{DEFAULT}]";
			
			# If neither Value nor Default has been set, print out
			# "not set".
			} else {
				print "[Not set]";
				@{$flag->{VALUE}}  = ();
			}
			
		# Otherwise, the flag is boolean:
		} else {
			
			# Print a "yes/no" prompt
			print "(y/n) ";
			
			# If the value or default value are true, print out
			# "yes".
			if ($flag->{VALUE} || $flag->{DEFAULT}) {
				print "[y]";
			
			# Otherwise, print out "no".
			} else {
				print "[n]";
			}
		}
		
		# Print a colon to finish the advisory statement.
		print ": ";
		
		# Get the value to set the flag to from standard input.
		my $string = <STDIN>;
		chomp $string;
		
		# If standard input returns an empty string:
		if ($string eq "") {
			
			# Set the value as the default value if there isn't already
			# a defined value.
			$flag->{VALUE} = $flag->{DEFAULT} if (!defined $flag->{VALUE});
			
		# Otherwise,
		} else {
			# If the current flag is expecting a scalar value...
			if ($flag->{HAS_VALUE}  == 1) {
				
				# Set the value of the flag based on the just-supplied
				# string
				$flag->{VALUE} = $string;
				
			# If the current flag is expecting an array
			} elsif ($flag->{HAS_VALUE}  == 2) {
				
				my $list;
				@$list = split " ", $string;
				$flag->{VALUE} = $list;
			
			# If the current flag is boolean...
			} else {
				
				while (1) {
					
					# Get the first character of the supplied string,
					# and make it lower case.
					$string = lc substr (0, 1, $string." ");
				
					# Exit the loop if the supplied string is a valid
					# value.
					last if $string eq "y";
					last if $string eq "n";
					last if $string eq "";
					
					# Prompt for the supply of a string.
					print "Enter value for flag: ";
		
					# Print out short name as part of the prompt (plus a
					# comma if the long name is used).
					if (defined $flag->{SHORT_NAME}) {
						print "-$flag->{SHORT_NAME}";
						print ", " if (defined $flag->{LONG_NAME});
					}
					
					# Print out long name as part of the prompt.
					if (defined $flag->{LONG_NAME}) {
						print "--$flag->{LONG_NAME}";
					}
		
					# Print out description as part of the prompt.
					print " $flag->{DESCRIPTION} " if (defined $flag->{DESCRIPTION});
					
					# Print out the available options.
					print "(y/n) ";
										
					# Print the default option (that will be used if
					# an empty string is supplied). True if value or
					# default value are true; false otherwise.
					if ($flag->{VALUE} || $flag->{DEFAULT}) {
						print "[y]";
					} else {
						print "[n]";
					}
					
					# Print a colon.
					print ": ";
		
					# Set a string based on the standard input.
					$string = <STDIN>;
					chomp $string;
					
					# Redo the loop.
				}
				
				# Set the flag's value to TRUE if the prompt returns yes.
				$flag->{VALUE} = 1 if $string eq "y";
				
				# Set the flag's value to FALSE if the prompt returns no.
				$flag->{VALUE} = undef if $string eq "n";
				
				# Set the flag's value to the default value if the empty string
				# is returned (and the flag isn't already set).
				$flag->{VALUE} = $flag->{DEFAULT}
					if ($string eq "" && !defined $flag->{VALUE});
			}
		}
	}
		
	# Check flag value is valid
	verify_input_data($flag->{VALUE}, $verify_string, "Flag: -$flag->{SHORT_NAME}. $flag->{DESCRIPTION}.", $flag->{DEFAULT})
		if (defined $verify_string && defined $flag->{VALUE});
	
	# Add the flag to the flag list
	push @$Silico::Flags, $flag;

	return $flag->{VALUE};
}


sub set_flag_from_arguments {

	#<
	#? Set flag values from arguments on the command line
	#. Flag and value are removed from the Silico::arguments array
	#; Requires: flag record
	#; Returns: nothing
	#>

	my $flag = $_[0];

	for (my $i = 0; $i <= $#Silico::Arguments; ++$i) {
	
		my $arg = $Silico::Arguments[$i];
		
		#  Long flags (--flag)
		if ($arg =~ /^\-\-/) {
		
			# If we have more than one equals sign in the flag, this is bad.
			if ($arg =~ /=.*=/) {
				silico_msg('d', "Long format flag supplied with more than one \"=\" sign! E.g. --Flag==20\n",
						"Please check your command line.\n");
			}
			
			# Get the "flag name" of the argument - everything between the opening "--"
			# and the equals sign.
			my $clipped = $arg;
			$clipped =~ s/^\-\-//;
			$clipped =~ s/=.*$//;
			
			# Skip if the current flag shouldn't have a long name.
			next if !defined $flag->{LONG_NAME};
			
			# Skip if the current flag's long name isn't the same
			# as this argument's long name.
			next if $flag->{LONG_NAME} ne $clipped;

			# Prevent duplicate values unless allowed
			if (!$flag->{DUPLICATES_ALLOWED} && defined $flag->{VALUE}) {
				silico_msg('d', "Attempted double definition of flag --$flag->{LONG_NAME}.\n");
			}
			
			# If the flag is expecting a scalar value...
			if ($flag->{HAS_VALUE}) {
			
				if ($arg !~ /=/) {
					silico_msg('d', "Long format flag does not contain '=' sign to set value.\n",
						"  Should be of the format '--flag=<val>'\n",
						"  Please check your command line.\n");
				}
				
				if ($flag->{HAS_VALUE} == 1) {

					# Set the flag's value to the original (i.e., unclipped) argument
					$flag->{VALUE} = $arg;

					silico_msg( 'd', "No value supplied for flag '".($flag->{SHORT_NAME} || $flag->{LONG_NAME})."'\n") if !defined ($flag->{VALUE});
				
					# Remove all up to and including the "=" sign
					$flag->{VALUE} =~ s/^.*=//;
					
				} elsif ($flag->{HAS_VALUE} == 2) {
				
					die ("Code not written here");
				}
			
			} else {

				# We have a Boolean flag. Set the value to 1 (i.e., true).
				$flag->{VALUE} = 1;
			}
			
			# Remove this entry from the Arguments array
			splice @Silico::Arguments, $i, 1;
			
			# Redo the loop from the same place
			--$i;
			
		# Short flags ("-")
		} elsif ($arg =~ /^\-/) {
		
			# Get the "flag name" of the argument - everything after the opening "-".
			my $clipped = $arg;
			$clipped =~ s/^\-//;

			# Skip if the current flag shouldn't have a short name.
			next if !defined $flag->{SHORT_NAME};
			
			# Skip if the argument doesn't match
			next if $flag->{SHORT_NAME} ne $clipped;
			
			# If the current flag is expecting a value...
			if ($flag->{HAS_VALUE}) {

				if (!$flag->{DUPLICATES_ALLOWED}) {

					# Prevent duplicates unless allowed
					if (defined $flag->{VALUE}) {
					
						silico_msg('d', "Attempted double definition of flag -$flag->{SHORT_NAME}.\n");
					
					} elsif ($flag->{HAS_VALUE} == 1) {
						
						# Remove this entry from the arguments array and set value
						splice @Silico::Arguments, $i, 1;
						$flag->{VALUE} = splice @Silico::Arguments, $i, 1;
						silico_msg( 'd', "No value supplied for flag '".($flag->{SHORT_NAME} || $flag->{LONG_NAME})."'\n") if !defined ($flag->{VALUE});

					} elsif ($flag->{HAS_VALUE} == 2) {
					
						splice @Silico::Arguments, $i, 1;
						my $val = splice @Silico::Arguments, $i, 1;
						@{$flag->{VALUE}} = split " ", $val;
						silico_msg( 'd', "No value supplied for flag '".($flag->{SHORT_NAME} || $flag->{LONG_NAME})."'\n") if !defined ($flag->{VALUE});
						#print "'@{$flag->{VALUE}}'\n";
					
					} else {

						silico_msg('e', "Error in flag_has_value: $flag->{HAS_VALUE}\n");
						carp();
						die;
					}

				} else {

					if ($flag->{HAS_VALUE} == 1) {

						# Remove entry from arguments array and set value
                                                splice @Silico::Arguments, $i, 1;
                                                push @{$flag->{VALUE}},  splice @Silico::Arguments, $i, 1;

					} elsif ($flag->{HAS_VALUE} == 2) {

						splice @Silico::Arguments, $i, 1;
						my $val = splice @Silico::Arguments, $i, 1;
						my @f = split " ", $val;
						@{$flag->{VALUE}} = (@{$flag->{VALUE}}, @f);

					} else {

						silico_msg('e', "Error in flag_has_value: $flag->{HAS_VALUE}\n");
                                                carp();
                                                die;
					}
				}
				
			# Otherwise, the current flag is boolean.
			} else {
				
				# Set the value of the flag to 1 (i.e., true)
				$flag->{VALUE} = 1;
				
				# Remove this entry from the arguments array.
				splice @Silico::Arguments, $i, 1;
			}
			
			# Decrement i in all cases, so the arguments loop starts again
			# from the same place (now that the array has been altered).
			--$i;
		}
	}
}

sub delete_flag {

	#<
	#? Delete a flag from the internal flags list.
	#. Searches by short name (most cases) or long name (if 'l' is supplied as
	#  a second argument).
	#; Requires: name of flag to delete (string), search type (optional string)
	#; Returns: nothing
	#>
	
	my $name = $_[0];
	my $type = lc $_[1];
		
	my $flag;
	my $i;
	
	if (defined $type && $type ne 'l' && $type ne 's' && $type ne '') {
		silico_msg('w', "Unrecognised search type \"$type\"!\n",
				"Assuming search is for short name.\n",
				"To delete a flag based on long name, pass 'l' as a second argument.\n");
	}

	for ($i = 0; $i <= $#{$Silico::Flags}; ++$i) {
	
		$flag = $Silico::Flags->[$i];
		
		if (($type eq 'l' && defined $flag->{LONG_NAME} && $flag->{LONG_NAME} eq $name)
			|| (defined $flag->{SHORT_NAME} && $flag->{SHORT_NAME} eq $name)) {
					
				splice @$Silico::Flags, $i, 1;
				--$i;
		}
	}
}

sub get_flag {

	#<
	#? Return the value held by a particular flag.
	#. Searches by short name (in most cases) or long name (if 'l' is supplied
	#  as a second argument).
	#; Requires: flag name (string), search type (optional string)
	#; Returns: flag value  (or pointer to array of values if duplicate definition is allowed);
	#>
	
	my $name = $_[0];
	my $type = $_[1] || silico_msg('d', "Flag type not provided to subroutine");

	$type = lc $type if $type;

	if (defined $type && $type ne 'l' && $type ne 's') {
		silico_msg('d', "Unrecognised search type \"$type\"!\n");
	}

	my $dup;
	@$dup = ();
	foreach my $f (@$Silico::Flags) {

		#printf "f_short: %s\n", ($f->{SHORT_NAME} || '-');
		#printf "f_long:  %s\n", ($f->{LONG_NAME} || '-');

		if ((defined $f->{SHORT_NAME} && $f->{SHORT_NAME} eq $name)
			|| ($type && $type eq 'l' && defined $f->{LONG_NAME} && $f->{LONG_NAME} eq $name)) {

			if ($f->{DUPLICATES_ALLOWED}) {

				die;
				push @$dup, $f->{VALUE};

			} else {

				return $f->{VALUE};
			}
		}
	}

	# Return array of values for flags where duplicate values are allowed
	return $dup if defined $dup->[0];;

	# It is possible that flags can be queried but have not been defined
	# in the main routine. 
	if (0) {
		silico_msg('w', "Attempt to use undefined flag $name\n");
	}
	
	return undef;
}

sub flag_is_boolean {

	#<
	#? Return true if a flag is a simple boolean switch
	#. Searches by short name (in most cases) or long name (if 'l' is supplied
	#  as a second argument).
	#; Requires: flag name (string), search type (optional string)
	#; Returns: true or false
	#>
	
	my $name = $_[0];
	my $type = $_[1];

	$type = lc $type if $type;
		
	if (defined $type && $type ne 'l' && $type ne 's' && $type ne '') {
		silico_msg('w', "Unrecognised search type \"$type\"!\n",
				"Assuming search is for short name.\n",
				"To provide a value based on long name, pass 'l' as a second argument.\n");
		carp();
	}

	foreach (@$Silico::Flags) {
		return $_->{HAS_VALUE} if ((defined $_->{SHORT_NAME} && $_->{SHORT_NAME} eq $name)
			|| ($type && $type eq 'l' && defined $_->{LONG_NAME} && $_->{LONG_NAME} eq $name));
	}
	
	return undef;
}

sub get_lflag {

	#<
	#? Get flag by long name
	#; Requires: flag name
	#. Calls get_flag
	#>
	
	return get_flag($_[0], 'l');
}

sub get_sflag {

	#<
	#? Get flag by short name
	#; Requires: flag name
	#. Calls get_flag
	#>
	
	return get_flag($_[0], 's');
}


sub set_flag {

	#<
	#? Set a flag to a value
	#. Searches by short name (default) or long name (if 'l' is supplied
	#  as a second argument).
	#; Requires: flag name (string), search type (optional string), value (default 1)
	#; Returns: nothing
	#>
	
	my $name = $_[0];
	my $type = $_[1] || 's';
	my $value = $_[2];
	
	$type = lc $type if $type;
	
	#This allows flags to be explicitly set to 0 or ''
	$value = 1 if !defined $value;
	
	if ($type ne 'l' && $type ne 's') {
		silico_msg('w', "Unrecognised search type \"$type\"!\n",
				"Assuming search is for short name.\n",
				"To provide a value based on long name, pass 'l' as a second argument.\n");
		$type = 's';
	}

	# Set flag if it exists already
	foreach (@$Silico::Flags) {
		if ((defined $_->{SHORT_NAME} && $_->{SHORT_NAME} eq $name)
			|| ($type eq 'l' && defined $_->{LONG_NAME} && $_->{LONG_NAME} eq $name)) {
			
			$_->{VALUE} = $value;
			return;
		}
	}
	
	# Make new flag if it does not exist.
	make_flag(undef, $name, undef, undef, $value) if $type eq 'l';
	make_flag($name, undef, undef, undef, $value) if $type eq 's';
}

sub set_sflag {

	#<
	#? Set flag by short name
	#; Requires: flag name, value
	#. Calls set_flag
	#>
	
	return set_flag($_[0], 's', $_[1]);

}

sub set_lflag {

	#<
	#? Set flag by long name
	#; Requires: flag name, value
	#. Calls set_flag
	#>

	return set_flag($_[0], 'l', $_[1]);
}

sub print_flags {

	#<
	#? Print out valuues of all flags
	#; Requires: optional filehadle (default STDOUT)
	#>

	my $FH = $_[0] || *STDOUT;
	
	return if  (get_lflag('very-quiet'));

	foreach my $f (@$Silico::Flags) {
	
		next if (!defined $f->{VALUE});
		
		print $FH "Using flag: ";

		print $FH "-$f->{SHORT_NAME}" if (defined $f->{SHORT_NAME});
		print $FH ", " if (defined $f->{SHORT_NAME} && defined $f->{LONG_NAME});
		print $FH "--$f->{LONG_NAME}" if (defined $f->{LONG_NAME});
		print $FH " ($f->{DESCRIPTION})" if (defined $f->{DESCRIPTION});
	
		if ($f->{HAS_VALUE}) {	
			
			my $ref = $f->{VALUE};
			if ($ref =~ 'ARRAY') {
				print $FH " with values:\n";
				foreach (@{$f->{VALUE}}) {
					print $FH "\t$_\n";		
				}
			} else {
				print $FH " with value '$f->{VALUE}'\n";
			}
		} else {
			print $FH "\n";
		}
	}
}

sub quiet {
	
	#<
	#? Check if we are being quiet (i.e., the --quiet option is used).
	#; Requires: nothing
	#; Returns: 1 (quiet), 0 (not quiet)
	#>
	
	return 1 if get_lflag('quiet');
	return 1 if get_lflag('very-quiet');
	return 0;

}

sub veryquiet {
	
	#<
	#? Check if we are being very quiet (i.e., the --very-quiet option is used).
	#; Requires: nothing
	#; Returns: 1 (very quiet), 0 (not very quiet)
	#>
	
	return 1 if (get_flag('very-quiet', 'l'));
	return 0;

}

sub get_arguments2 {

	#<
	#? Improved routine to get arguments from command line.
	#. UNIX:
	#. Accepts wildcards '*', '?' and [...].
	#  Skips any arguments starting with '-' which
	#  are used as flags. These can be found using 'get_flags'
	#. UNIX file read supports directories containing very large numbers
	#  of files where the standard UNIX Glob fails, if wildcard characters
	#  quoted (eg '*') on the command line.
	#. If a -help argument is found then the silico documentation is printed
	#  and the progarm exits
	#; Requires: arguments in @ARGV.
	#; Returns: argument list as an array.
	#>

	my $numargs = $_[0];
	my $use_default_flags = $_[1] || 'standard';
	
	my @arguments;
	my @files;
	my @listing;

	$use_default_flags = lc $use_default_flags;

	# Run setup_flags if it has not already been done
	setup_flags() if (!defined $Silico::Flags_setup);
	
	# Complete flag setup by checking for default flags
	if ($use_default_flags =~ /\bstandard\b/) {
		get_default_flags_standard();
	} else {
		get_default_flags_minimal();
	}

	# Check arguments
	my @arglist;
	foreach my $re (@Silico::Arguments) {
	
		if ($re =~ /^\-/) {
			silico_msg('d', "Attempt to use unrecognised flag '$re'.\n");
		}
		push @arglist, $re;
	}

	# Handle filenames
	foreach my $re (@arglist) {
	
		# If we have a wildcard then expand it
		if ($re =~ /\*/ || $re =~ /\?/) {
	
			my @f = split $Silico::FS, $re;
			undef @files;

			# Find out if we are starting from the root directory
			if ($f[0] eq '') {
				shift @f;
				$files[0] = '/';
			} else {
				$files[0] = '.';
			}
			
			foreach $re (@f) {

				my @newfiles;

				if ($re =~ /\*/ || $re  =~ /\?/ || $re  =~ /\[.*\]/) {

					$re =~ s/\./\\\./g; # Make a dot a dot (not a wildcard)
					$re =~ s/\*/.*/g; # Star wildcard
					$re =~ s/\?/.{1}/g; # Question mark wildcard
					print "Matching Perl regular expression: $re\n";
				
					# Contains wildcard
					foreach my $dir (@files) {

						# Get listing of subdirectories
						next if  !-d $dir;
						@listing = `ls $dir`;

						# Keep them if they match the current directory regex
						foreach (@listing) {
							next if !/^$re$/;
							chomp;
							push @newfiles,$dir.$Silico::FS.$_;
						}
					}
				} else {
					# Does not contain wildcard
					foreach my $dir (@files) {
						push @newfiles,$dir.$Silico::FS.$re;
					}
				}
				@files = @newfiles;
			}
			@arguments = (@arguments, @files);

		} else {
			# Otherwise just add file to arguments if no wildcard
			chomp $re;
			push  (@arguments, $re);
		}
	}

	if ($#arguments >= 0) {
		print "Matched ".($#arguments+1)." file";
		print "s" if ($#arguments != 0);
		print "\n";
	}

	# Change '//' to '/' and './' to ''
	my $arglist;
	foreach (@arguments) {

		s|//|/|g;
		s|^\./||g;
		
		push @$arglist, $_;
	}
	
	# We want to check for number of arguments.
	check_arguments($numargs, $arglist) if (defined $numargs);
	
	# Print the flags
	print_flags();

	return @arguments;
}

sub check_arguments {
	
	#<
	#? Check the argument list for how many arguments there are.
	#  If an error has been made in the instruction for number of arguments,
	#  the script will die. If the tests are not passed, the script will
	#  exit after running PRINTDOC. Otherwise, nothing is done: the program
	#  carries on.
	#; Requires: string specifying number of arguments; list of arguments
	#; Returns: nothing
	#>
	
	my $numargs = $_[0];
	my $arglist = $_[1];
	
	my @eq;
	my @gt;
	my @lt;
	my @sortgt;
	my @sortlt;
	
	my @f = split(/\s+/, $numargs);
	
	foreach my $f (@f) {
		
		# Remove finishing carriage returns and whitespace
		chomp $f;
		$f =~ s/ //g;
				
		# Put less-than, less-than-or-equal-to, greater-than,
		# greater-than-or-equal-to, and equal-to values in
		# appropriate lists.
		
		# Note:
		# Greater-than-or-equal-to and less-than-or-equal-to
		# values are combined into greater-than and less-than
		# respectively, according to >=a is equivalent to >(a-1)
		# and <=b is equivalent to <(b+1)
		
		if ($f =~ /^<=\d+$/) {
			$f =~ s/^<=//;
			push @lt, ($f+1);
		} elsif ($f =~ /^<\d+$/) {
			$f =~ s/^<//;
			push @lt, $f;
		} elsif ($f =~ /^>=\d+$/) {
			$f =~ s/^>=//;
			push @gt, ($f-1);
		} elsif ($f =~ /^>\d+$/) {
			$f =~ s/^>//;
			push @gt, $f;
		} elsif ($f =~ /^\d+$/) {
			push @eq, $f;
		} else {
			silico_msg('d', "Bad expected number of arguments!\n",
					"Expected number: $numargs\n",
					"Expected: one or more of \"num\", \">=num\", \">num\", \"<=num\", \"<num\".\n");
		}
	}
	
	# If there are different "equals" values, print an error
	foreach (@eq) {
		if ($_ != $eq[0]) {
			silico_msg('d', "Bad expected number of arguments!\n",
					"Expected number: $numargs\n",
					"Number of arguments can't be $eq[0] and $_ at the same time.\n");
		}
	}
	
	my $eq = $eq[0] if ($#eq >= 0);
	
	# Discard all but the highest greater-than value
	@sortgt = sort {$b <=> $a} @gt;
	my $gt = $sortgt[0];
	
	# Discard all but the lowest less-than value
	@sortlt = sort {$a <=> $b} @lt;
	my $lt = $sortlt[0];
	
	# Print errors if we've tried to do something impossible.
	
	# Case 1. $gt is at least $eq
	if (defined $gt && defined $eq && $gt >= $eq) {
		silico_msg('d', "Bad expected number of arguments!\n",
			"Expected number: $numargs\n",
			"Number of arguments can't be $eq and greater than $gt at the same time.\n");
	}
	
	# Case 2. $lt is at most $eq
	if (defined $lt && defined $eq && $lt <= $eq) {
		silico_msg('d', "Bad expected number of arguments!\n",
			"Expected number: $numargs\n",
			"Number of arguments can't be $eq and less than $lt at the same time.\n");
	}
	
	# Case 3. $gt is at least $lt
	if (defined $lt && defined $gt && $lt <= $gt) {
		silico_msg('d', "Bad expected number of arguments!\n",
			"Expected number: $numargs\n",
			"Number of arguments can't be greater than $gt and less than $lt at the same time.\n");
	}
	
	# Perform the tests
	if (defined $eq && $#{$arglist} != ($eq-1)
		|| defined $lt && $#{$arglist} >= ($lt-1)
		|| defined $gt && $#{$arglist} <= ($gt-1)) {

		# If the number of arguments is wrong (but greater than 0),
		# print out a message explaining that this is the problem
		silico_msg ('h', "Incorrect number of arguments specified. Exiting.\n");
	}
}

sub get_allflags {

	#<
	#? Return a list of all command line flags
	#; Requires: nothing
	#; Returns: flags in array
	#>
	
	my @list;
	my $arg;
	
	foreach $arg (@ARGV) {
		push (@list, $arg) if $arg =~ /^-/;
	}
	
	return @list;
}

##################################################################
#
#	Filename routines
#
##################################################################

sub get_filebase {
	
	#<
	#? Get file basename
	#. Removes only the final extension unless it is Z. .gz or .bz then the last
	#  two are removed
	#; Returns: string
	#>
	
	my $filename = $_[0];
	
	silico_msg('d', "Undefined input value\n") if !defined $filename;
	
	# Remove .gz or Z
	$filename =~ s/\.gz$//;
	$filename =~ s/\.Z$//;
	$filename =~ s/\.bz2$//;

	# Return unchanged filename if there are no dots in it
	return $_[0] if ($filename !~ /\./);

	my @f = split(/\./, $filename);
	return  join(".",@f[0..($#f-1)]);
}

sub get_filebase_short {
	
	#<
	#? Get file basename without any leading directories
	#. Removes the final extension unless it is .gz then the last
	#  two are removed
	#; Requires: filename
	#; Returns: string
	#>

	my $filename = $_[0];
	my @f;

	# Remove .gz or Z
	$filename =~ s/\.gz$//;
	$filename =~ s/\.Z$//;
	$filename =~ s/\.bz2$//;

	@f = split(/$Silico::FS/, $filename);
	$filename = $f[$#f];

	# Return unchanged filename if there are no dots in it
	return $filename if ($filename !~ /\./);

	@f = split(/\./, $filename);
	return  join(".",@f[0..($#f-1)]);
}

sub get_file_dir {
	
	#<
	#? Get file directory
	#; Requires: filename
	#; Returns: string
	#>

	my $filename = $_[0];
	
	# Return '' if there are no directories given in path
	return '' if ($filename !~ $Silico::FS);
	
	my @f = split(/$Silico::FS/, $filename);

	return  join(".",@f[0..($#f-1)]);
}

sub get_fileext {
	
	#<
	#? Get file extension name
	#. Note.  Strips off any .gz, .bz2 or .Z extension before returning file extension
	#; Requires: filename
	#; Returns: string
	#>

	my $filename = $_[0];
	
	if (!defined $filename) {
		#OK to return undefined value. Let the calling subroutine deal with it.
		#silico_msg('w', "Undefined input value!\n");
		#carp();
		return undef;
	}

	# Remove .gz , .bz2 or Z
	$filename =~ s/\.gz$//;
	$filename =~ s/\.bz2$//;
	$filename =~ s/\.Z$//;

	my @f = split(/\./, $filename);
	#print "f @f fn $filename\n";
	return $f[-1] if $f[-1] ne $filename;
	return "";	# Return nothing if there are no dots
}

sub get_file_compression {
	
	#<
	#? Get file compression method (.Z, .bz2 or .gz)
	#; Requires: filename
	#; Returns: string
	#>

	my $filename = $_[0];
	
	if (!defined $filename) {
		silico_msg('g', "Undefined input value!\n");
		return undef;
	}

	return 'gz'if $filename =~ /\.gz$/;
	return 'bz2' if $filename =~ /\.bz2$/;
	return 'Z' if $filename =~ /\.Z$/;
	return '';	
}


sub format_oname {

	#<
	#? Generate a numbered filename (various possible styles)
	#. Styles:
	#; molname: Molecule Name
	#; molname_i: Renames molecule name using 'insight safe' name
	#; named:    <base>_<name>.<ext>
	#; numbered: <base>_<num>.<ext>
	#; fnumbered: <base>_000<num>.<ext> (%05d)
	#; fnumbered_N: <base>_000<num>.<ext> (%0Nd) - ie number formated to N places
	#; Requires: style, molecule (requred for style 'name'), basename (required for numbered styles),
	#  number (required for numbered styles) or name (required for name style), extension (optional),
        #  compression (optional)
	#; Returns: filename (string)
	#>

	my $style = lc $_[0] || 'fnumbered';
	my $mol = $_[1];
	my $base = $_[2] || 'Mol';
	my $number_or_name = $_[3] || 0;
	my $ext = $_[4] || '';
	my $compression = $_[5] || '';
	
	my $numberz = 5;
	
	$ext = ".".$ext if $ext;
	$ext = ".".$compression if $compression;
	
	$base =~ s/ //g;
	$base =~ s/\///g;
	
	if (($style eq 'molname') && $mol) {
		return ($mol->{NAME} || 'Mol').$ext;
	}
	
	if (($style eq 'molname_i') && $mol) {
		return insight_rename(($mol->{NAME} || 'Mol'), $Silico::of_namehash).$ext;
	}
	
	if ($style eq 'named') {
		return ($base || 'Mol')."_".$number_or_name.$ext;
	}
	
	if ($style  =~ 'fnumbered') {
		$style =~ m/fnumbered_(.*)/;
		$numberz = $1 || 5;
		return $base.sprintf ("_%0$numberz"."d", $number_or_name).$ext;
	}
	
	if ($style eq 'numbered') {
		return $base.'_'.$number_or_name.$ext;
	}
	
	silico_msg ('w', "Output style '$style' does not exist or is not suitable.  Using 'fnumbered' format\n");
	return $base.sprintf ("_%0$numberz"."d", $number_or_name).$ext;
}

##################################################################
#
#	Text files
#
##################################################################

sub read_textfile {

	#<
	#? Read a file into a string
	#; Requires: filename
	#; Returns: text
	#>
	
	my $file = $_[0];
	
	open (INFILE, $file) || silico_msg('d', "Could not open file $file for reading\n");
	
	my $in = '';
	
	while (<INFILE>) {
		$in .= $_;
	}
	
	close INFILE;
	
	return $in;
}

sub write_textfile {

	#<
	#? Write a string to a file
	#; Requires: filename, text;
	#>
	
	my $file = $_[0];
	my $text = $_[1];
	
	open (OUTFILE, ">$file") || silico_msg('d', "Could not open file $file for writing\n");
	print OUTFILE $text;
	close OUTFILE;
}

##################################################################
#
#	Misc
#
##################################################################

sub slurp_file {

	#<
	#? Read a molecule file and split into separate molecules as
	#  TEXT files.
	#. Currently only splits pdb, mol2 and sdf files.  All others are
	#  returned unsplit
	#; Requires: filename
	#; Returns: an array containing each molecule as text
	#>

	my $name = $_[0];
	
	my $slurp;
	my $molcount = 0;
	my $modelcount = 0;

	my $sname = $name;
	$sname =~ s/\.Z//;
	$sname =~ s/\.gz//;
	$sname =~ s/\.bz2//;
	
	if (!open_file(*SLURP, $name)) {
		silico_msg('e', "Can not open file $name for reading!\n");
		return undef;
	}
	
	# pdb
	if ($sname =~ /\.pdb$|\.ent$|.coor$/) {
		while (<SLURP>) {
			# Increment moleucule count for each END or nnmr
			# MODEL statement
			if (/^MODEL/) {
				++$modelcount;
				++$molcount if $modelcount > 1;
			}
			$slurp->[$molcount] .= $_;
			++$molcount if /^END/;
		}
		
		return $slurp;
	}
	# mol2
	my $flag = 1;
	if ($sname =~ /\.mol2$/) {
		while (<SLURP>) {
			# Increment molecule count on <TRIPOS>MOLECULE of
			# comment. Set flag so that we do not start a new
			# molecule at each comment or <TRIPOS>MOLECULE line
			if (/^@<TRIPOS>MOLECULE/)  {

				++$molcount if !$flag;
				$flag = 0;
			}
			if (/^#/) {

				++$molcount if !$flag;
				$flag = 1;
			}

			$slurp->[$molcount] .= $_;
		}
	}
	# sdf
	if ($sname =~ /\.sdf$|\.mol$/) {
		while (<SLURP>) {
			# Increment model count if we find $$$$
			$slurp->[$molcount] .= $_;
			++$molcount if /\$\$\$\$/;
		}
	}

	# All others
	while (<SLURP>) {
		$slurp->[0] .= $_;
	}
	
	return $slurp;
}

sub overwrite_check {

	#<
	#? Check to see if an output file exists and prompt for response if it does
	#. Check is overridden by -force flag
	#; Requires: filename, flag to disallow option of changing filename
	#; Returns: new filename (if changed) or original filename.  The '-force' flag
	#  is set if 'Yes to all' is selected.
	#> 

	my $infile = $_[0];
	my $flag = $_[1];

	return $infile if get_flag('force', 's');
	
	my $file = $infile;
	
	while (1) {
	
		last if (!-e $file);
		last if $Silico::force;
	
		silico_msg('v', "File $file already exists! Do you want to continue anyway?\n");
		
		while (1) {
		
			my $answer;
			if ($flag) {
				$answer = prompt("[y]es (overwrite) / yes to [a]ll / [n]o (exit): ");
			} else {
				$answer = prompt("[y]es (overwrite) / yes to [a]ll / [n]o (exit) / [c]hange name: ");
			}
			$answer =~ s/ //g;
			$answer = substr($answer, 0, 1);
			$answer = lc($answer);
			
			if ($answer eq 'a') {
				$Silico::force = 1;
				return $file;
			} elsif ($answer eq 'y') {
				return $file;
			} elsif ($answer eq 'n') {
				silico_msg('v', "Exiting.\n");
				exit;
			} elsif (!$flag && $answer eq 'c') {
				last;
			}
		}
		
		while (1) {
			$file = prompt("Please enter alternative file name: ");
			$file = cleanline($file);
			last if (check_data_type($file, 'FILENAME') == 1);
			silico_msg('v', "File name $file is inappropriate.\n");
		}
	}
	
	return $file;
}


sub file_read_error {

	#<
	#? Print an error on failing to read a file
	#; Requires: filename, flag to write error message and die
	#; Returns: Nothing
	#> 

	carp() if $Silico::debug;
	
	if (!defined $_[0]) {
		silico_msg('d', "No filename supplied!\n");
	}
	
	if ($_[1]) {
		silico_msg('d', "Can not open file $_[0] for reading!\n");
	} else {
		silico_msg('e', "Can not open file $_[0] for reading!\n",
				"Skipping.\n");
	}
}

sub file_write_error {

	#<
	#? Print an error on failing to write a file
	#; Requires: filename, flag to write error message and die
	#; Returns: Nothing
	#> 
	
	carp() if $Silico::debug;

	if (!defined $_[0]) {
		silico_msg('d', "No filename supplied!\n");
	}
	
	if ($_[1]) {
		silico_msg('d', "Can not create or open file $_[0] for writing!\n");
	} else {
		silico_msg('e', "Can not create or open file $_[0] for writing!\n",
				"Skipping.\n");
	}
}

sub delimited_read_error {

	#<
	#? Print an error on failing to read delimited data from file
	#; Requires: filename, flag to write error message and die
	#; Returns: Nothing
	#> 
	
	carp() if $Silico::debug;

	if (!defined $_[0]) {
		silico_msg('d', "No filename supplied!\n");
	}

	
	if ($_[1]) {
		silico_msg('d', "Can not read delimited data from file $_[0]!\n");
	} else {
		silico_msg('e', "Can not read delimited data from file $_[0]!\n",
				"Skipping.\n");
	}
}

sub mol_read_error {
	
	#<
	#? Print an error on failing to read molecule data from a file
	#; Requires: filename, flag to write error message and die
	#; Returns: Nothing
	#> 

	carp() if $Silico::debug;

	if (!defined $_[0]) {
		silico_msg('d', "No filename supplied!\n");
	}

	if ($_[1]) {
		silico_msg('d', "Can not read molecule data from file $_[0]!\n");
	} else {
		silico_msg('e', "Can not read molecule data from file $_[0]!\n",
				"Skipping.\n");
	}
}

sub seq_read_error {
	
	#<
	#? Print an error on failing to read sequence data from a file
	#; Requires: filename, flag to write error message and die
	#; Returns: Nothing
	#> 

	carp() if $Silico::debug;

	if (!defined $_[0]) {
		silico_msg('d', "No filename supplied!\n");
	}
	
	if ($_[1]) {
		silico_msg('d', "Can not read sequence data from file $_[0]!\n");
	} else {
		silico_msg('e', "Can not read sequence data from file $_[0]!\n",
				"Skipping.\n");
	}
}


##################################################################
#
#	Output filename and format generation
#
##################################################################


sub set_default_oappend {

        #<
        #? Define a default oappend string
	#; Requires: string
        #>

        $Silico::output_append = $_[0];
}


sub get_ofilename {

	#<
	#? Get the output filebase and format of an ensemble
	#  checking flags and input format
	#. 'Base name' and 'output format' will be overridden by command line flags if they are supplied.
	#  The 'default format' will also be overridden by command line options.
	#. See also: get_ofilename2
	#; Requires: (ensemble, base name, format, append string, default_format) all optional
	#  splitflag to provide separate filename and extension, flag force use of supplied
    	#  names ignoring command line flags
	#; Returns: filename if 'splitflag' not set or filebase & format if it is set.
	#>

	my $ens = ens($_[0]);
	my $obase = $_[1];
	my $oformat = $_[2];
	my $append = $_[3];
	my $default_format = $_[4];
	my $splitflag = $_[5];
	
	my $ofilename;
	
	# Setting append to '' will prevent any append string being added
	# if a specific filebase has been provided
	$append = '' if defined $obase && !defined $append;
	$obase = get_ofilebase($ens, $obase, $append);
	$oformat = get_oformat($ens, $oformat, $default_format);
	$ofilename = $obase.".".$oformat;
	
	return $ofilename if !$splitflag;
	return $obase, $oformat;
}

sub get_ofilename2 {

	#<
	#? Split version of get_ofilename
	#. See get_ofilename
	#>
	
	$_[5] = 1;  # Set 'splitflag'
	return get_ofilename(@_);
}


sub get_ofilebase {

	#<
	#? Get output filebase. Takes command line flags into account
	#. Filebase is obtained from 
	#; 1. -O flag on command line
	#; 2. Name depending on SOURCE_FILE_NAME
	#; 3. 'outfile'
	#; Input: ensemble (optional), obase (optional), append_string (optional);
	#; Output: filebase
	#>
	
	my $ens = ens($_[0]);
	my $obase = $_[1];
	my $append = $_[2];

	# If output name is defined on command line then use this.  The append string is used
	# if one is explicitlty provided
	if (get_flag('output-file-name', 'l')) {

		$obase = get_filebase(get_flag('output-file-name', 'l'));
		# Add append string unless it is set to '' or ' '
		if ($append && !$append ne ' ') {
			$append =~ s/ //g;
			$obase .= '_'.$append;
		}
		return $obase;
	}

	# Otherwise generate a filename base depending on mol SOURCE_FILE_NAME
	$obase ||= get_filebase($ens->[0]{SOURCE_FILE_NAME}) if defined $ens->[0]{SOURCE_FILE_NAME};
	
	# Still could not generate a filename base.  Use 'outfile'
	if (!defined $obase) {
		silico_msg('w', "Unable to determine output filename base.  Using 'outfile'\n");
		$obase = 'outfile';
	}
		
	# Use either the explicitly supplied value for append or the default value
	if (!defined $append) {
		$append = $Silico::output_append if defined $Silico::output_append;
		$append = 'new' if !defined $append;
	}
		
	# Add append string unless it is set to '' or ' '
	if ($append ne '' && $append ne ' ') {

		$append =~ s/ //g;
		$obase .= '_'.$append;
	}
	
	return $obase;
}

sub get_oformat {

	#< 
	#? Get an appropriate output file format.  
	#. Selection priority:
	#; format ($_[1])
	#; -o <ext> flag
	#; -O <file.ext> flag
	#; format of source file
	#; provided default format ($_[2])
	#; mol2
	#.
	#; Input: ensemble or molecule (optional), default format (optional)
	#; Output: format string
	#>
	
	my $ens_mol = $_[0]; # Ensemble or molecule
	my $format =  $_[1];
	my $default = $_[2] || 'mol2'; # Default value
	
	my $mol;
	my $oformat;

	# Find out if we have an ensemble or a molecule
	if (ref $ens_mol eq 'ARRAY') {
		$mol = $ens_mol->[0];
	} elsif (ref $ens_mol eq 'HASH') {
		$mol = $ens_mol;
	}
	
	$oformat ||= get_sflag('o');
	$oformat ||= get_fileext(get_flag('output-file-name', 'l'));
	$oformat ||= $format;
	$oformat ||= get_fileext($mol->{SOURCE_FILE_NAME}) if $mol->{SOURCE_FILE_NAME};

	if (!defined $oformat) {
		silico_msg('w', "Unable to determine output format.  Using $default\n");
		$oformat = $default;
	}
	
	return $oformat;
}

sub get_tempname {

        #<
        #? Make a temporary filename base by concatenating program_name and a random number
        #  Uses $Silico::temp_dir (from flag --tempdir) or /usr/tmp (or /tmp) if this is not set
        #; Requires: temp directory
        #; Returns: filename base
        #>

	my $dir;
        my $hostname;
        my $progname;
     
	my $r;
	my @chars=('A'..'Z',0..9, 'a'..'z');

	if ($Silico::tempdir) {
		$dir = $Silico::tempdir;
	} elsif ($Silico::debug) {
		$dir = '.';
	} else {
        	$dir = '/usr/tmp' if -d '/usr/tmp';
        	$dir = '/tmp' if -d '/tmp';
	}

	$dir = $_[0] || $Silico::temp_dir || $dir;

        $progname = $0;
        $progname =~ s/.*\///g; # Delete all but last part of directory path
	
	foreach (1..6) {
		
		$r .= $chars[rand @chars];
	}

        return "$dir/$progname\_$r";
}



return 1;

