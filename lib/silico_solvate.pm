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
#! silico_solvate.pm
#? Routines for solvation, density etc
#. $Revision: 1.6.2.2.2.7 $
#>

use strict;
package Silico;


sub calc_molecule_density {

	#<
	#! Calc_molecule_density
	#? Calculate the number of molecules required to achieve a certain density
	#; Requires: Solute molecule, solvent molecule, box size, required density
	#; Returns: Number of solvent molecules required
	
	my $mol = $_[0]; # Solute
	my $smol = $_[1] || croak(); # Solvent
	my $boxx = $_[2];
	my $boxy = $_[3];
	my $boxz = $_[4];
	my $density = $_[5] || 1;
	
	# Mass of solute
	my $formula;
	my $mw;
	if ($mol) {
		print "\n";
		$mw = sprintf "%0.2f", molecule_mw($mol);
		$formula = molecule_formula($mol);
		print "Solute molecular weight: $mw\n";
		print "Solute molecular formula: ";
		foreach my $el (sort keys %$formula) {
			print $el.$formula->{$el}." ";
		}
		print "\n";
	}

	# Calculate Mass of solvent
	my $solvformula = molecule_formula($smol);
	my $solvmw = sprintf "%0.2f", molecule_mw($smol);

	print "\n";
	print "Solvent molecular weight: $solvmw\n";
	print "Solvent molecular formula: ";
	foreach my $el (sort keys %$solvformula) {
		print $el.$solvformula->{$el}." ";
	}
	print "\n";
	
	my $vol = $boxx * $boxy * $boxz;
	my $tot_mass = $density * $vol * $Avogadro / 1e24;
	my $req_mass = $tot_mass - $mw;
	my $numsmols = int ($req_mass/$solvmw);
	
	if ($numsmols < 0) {
		my $n = $mol->{NAME} || 'Mol';
		print "\n";
		silico_msg('w', "A negative number of solvent molecules ($numsmols) is required for molecule '$n'!\n");
		return $numsmols;
	}
	
	$vol = sprintf "%0.2f", $vol;
	$density =  sprintf "%0.2f", $density;
	$req_mass =  sprintf "%0.2f", $req_mass;
	$tot_mass =  sprintf "%0.2f", $tot_mass;
	
	silico_msg('c',
		"Molecular weight: $solvmw\n",
		"Box Volume: $vol A^3\n",
		"Total mass required for density $density: $tot_mass amu\n",
		"Mass required: $req_mass amu\n",
		"Total molecules: $numsmols\n\n");
	
	return $numsmols;
}


sub calc_system_density {

	#<
	#? Calculate the density of a periodic system
	#; Requires: System molecule, box size (x, y, z), quiet flag
	#; Returns: Density
	
	my $mol = $_[0]; # System
	my $boxx = $_[1] || $mol->{CELL}[0];
	my $boxy = $_[2] || $mol->{CELL}[1];
	my $boxz = $_[3] || $mol->{CELL}[2];
	my $quiet = $_[4];

	my $mw = sprintf "%8.2f", molecule_mw($mol);
	my $formula = molecule_formula($mol);
	my $vol = $boxx * $boxy * $boxz;
	my $density = 1e24 * $mw/($vol * $Silico::Avogadro);
	
	if (!$quiet) {	
		silico_msg('c', "\nSystem molecular weight: $mw\n",
				"System molecular formula: ");
		foreach my $el (sort keys %$formula) {
			silico_msg('c', $el.$formula->{$el}." ");
			}
		silico_msg('c',"\n");
		my $d = sprintf "%0.3f", $density;
		silico_msg('c', "System density (assuming all atoms present): $d\n");		
	}
	
	return $density;
}

sub create_water_molecule {

	#<
	#? Create a template of a water molecule.
	#: Requires: Residue name (optional)
	#: Returns: molecule
	#>
	
	my $rname = $_[0];
	
	my ($res, $o, $h1, $h2);
	
	if (get_sflag('amber')) {

		$res = 'WAT';
		$o = 'O';
		$h1 = 'H1';
		$h2 = 'H2';

	} elsif (get_sflag('gromacs')) {

		$res = 'SOL';
		$o = 'OW';
		$h1 = 'HW1';
		$h2 = 'HW2';

	} else {

		# Default TIP3P
		$res = 'HOH';
		$o = 'OH2';
		$h1 = 'H1';
		$h2 = 'H2';
	}
	
	# Override standard names
	$res = $rname if $rname;
	
	# Create a new molecule and give it a name.
	my $mol = create_molecule('Water');
	
	# Add oxygen and two hydrogens. Each additional atom is added
	# to the neighbour hash by mol_add_atom.
	my $atom1 = mol_add_atom($mol, " $o", "O", 0, 0, 0, $res, 1);
	my $atom2 = mol_add_atom($mol, " $h1", "H", -0.667, -0.639, 0.220, $res, 1);
	my $atom3 = mol_add_atom($mol, " $h2", "H", 0.882, -0.446, 0.206, $res, 1);
	
	foreach my $atom ($atom1, $atom2, $atom3) {
		$atom->{CHAIN} = 'W';
		$atom->{SEGID} = 'WATR';
		$atom->{SUBNAME} = $res;
		$atom->{FG}{W} = 1;
	}
				
	if (!get_sflag('nobond')) {
		bond_create_atom($mol, $atom1, $atom2, 1);
		bond_create_atom($mol, $atom1, $atom3, 1);
	}
	
	return $mol;
}


sub molecule_solvate {

	#<
	#? Solvate a molecule
	#: Requires: molecule, solvent molecule, cell, number of molecules, bilayer flag, output format (pdb etc)
	#: Returns: molecule
	#>

	my $embed = $_[0];
	my $smol = $_[1];
	my $cell = $_[2];
	my $nummols = $_[3];
	my $bilayer = $_[4];
	my $oformat = $_[5];

	my $smols;
	$smols->[0] = $smol;
	my $num_components;
	$num_components->[0] = ($nummols);

	my ($mollist, $chain, $segid) = make_molecule_copies($smols, $num_components, 'S', 'SOL', get_sflag('random'), $oformat);

	my $molcount = $#{$mollist}+1;

	# Pack molecules
	silico_msg('c', heading("Packing $molcount solvent molecules\n"));

	my $model;
	$model = pack_molecules($mollist, $embed, $cell);

	# Set cell
	$model->{CELL} = $cell;
	$model->{CELL_ALPHA} = $model->{CELL_BETA} = $model->{CELL_GAMMA}= 90;

	mol_move_into_cell($model);

	return $model;
}



return 1;
