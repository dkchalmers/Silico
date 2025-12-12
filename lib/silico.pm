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
#! silico.pm
#? The parent silico module
#.
#. Silico: a perl molecular toolkit.
#.
#. Written by David Chalmers.  29.1.2000 and onward.
#.
#? Silico.pm is the Silico parent module.
#. Loads a library of of subroutines for handling molecular structures
#  and routines to read and write a number of standard molecular modelling
#  file formats.
#.
#. Silico.pm is usually loaded by the subroutine silico_setup which
#  defines the location of the Silico libraries and executables. 
#.
#. By default the following Silico modules are loaded:
#, silico_defnitions - info about atoms etc
#, silico_check - molecule checking routines
#, silico_control - program control and user interface routines
#, silico_data - internal data handling routines
#, silico_debug - debugging routines
#, silico_doc - automatic documentaton generation routines
#, silico_geom - molecular geometry routines
#, silico_hydrogens - Add and delete hydrogens from molecules
#, silico_io - general input/output routines
#, silico_mol2 - read/write Sybyl mol2 files
#, silico_molecules - general operations on molecules
#, silico_pdb - read/write Brookhaven pdb files
#, silico_prop - calculate molecular and atomic properties
#, silico_rings - identify rings in molecules
#, silico_sdf - read/write MDL sdf files
#, silico_statistics - statistical routines
#, silico_residue - residue routines
#, silico_asl - atom selection language
#
#. These are loaded by silico_io if required
#, silico_gromacs - Gromacs format files
#, silico_merck - Merck format files
#, silico_tinker - Tinker format files
#.
#. These have to be specified using a 'require' statement
#, silico_smiles - Read pseudo smiles and type atoms
#, silico_sequence - Extract protein sequences
#, silico_split - Split up molecules based on connectivities
#
#. The global variable $silico_version is also set.
#. $Revision: 1.27.2.2.2.6 $
#>

use strict;
package Silico;

################################
# GLOBAL VARIABLE DECLARATIONS #
################################

use vars qw(

	$Hostname
	$PROG
	$silico_version
);


$Silico::FS = "/"; # Directory symbol (used in some places)
$Silico::data_dir = $Silico::home_dir."data"."/"; # Data directory
$Silico::lib_dir2 = $Silico::lib_dir."local"."/"; # Local lib directory
push @INC, , substr($Silico::lib_dir2, 0, -1);

# Read the libraries

# Silico_definitions comes first because
# it contains more global variables
use silico_definitions;

# Then other packages in alphabetical order
use silico_control;
use silico_data;
use silico_debug;
use silico_doc;
use silico_check;
use silico_geom;
use silico_gromacs;
use silico_hydrogens;
use silico_hb;
use silico_io;
use silico_mol2;
use silico_molecules;
use silico_mmod;
use silico_mopac;
use silico_pdb;
use silico_prop;
use silico_residue;
use silico_rings;
use silico_sdf;
use silico_split;
use silico_statistics;
use silico_asl;

# Set the program name
$Silico::PROG = $0;
$Silico::PROG =~ s/.*\///;

# Find version and return
$silico_version = 1.15;

# Set the hostname
$Silico::Hostname = `hostname -s`;
chomp $Silico::Hostname;

return 1;
