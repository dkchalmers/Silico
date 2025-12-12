#!/usr/bin/perl -w

#  COPYRIGHT NOTICE
#  
#  Silico - A Perl Molecular Toolkit
#  Copyright (C) 2008 David K. Chalmers and Benjamin P. Roberts,
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
#! silico_geom.pm
#? Silico routines to calculate or operate on molecular geometries
#  (distances, angles, etc)
#. $Revision: 1.15-66-g8b88391 $
#>

use strict;
package Silico;

##################
# LINEAR ALGEBRA #
##################

sub crossprod {
	
	#<
	#? Calculates the cross product of two vectors.
	#; Requires: two vectors
	#; Returns: cross product (an array {x, y, z})
	#>
	
	foreach (@_[0..5]) {
		silico_msg('d', "Error in input") if !defined $_;
	}

	my ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
	
	return ($y1 * $z2 - $y2 * $z1, $z1 * $x2 - $z2 * $x1, $x1 * $y2 - $x2 * $y1);
}

sub dotprod {

	#<
	#? Calculates the dot product of two vectors.
	#; Requires: two vectors
	#; Returns: dot product (a scalar)
	#>
	
	foreach (@_[0..5]) {
		silico_msg('d', "Error in input") if !defined $_;
	}
	
	return ($_[0] * $_[3] + $_[1] * $_[4] + $_[2] * $_[5]);
}

sub magnitude {

	#<
	#? Calculates the magnitude of a vector.
	#; Requires: vector as an array
	#; Returns: scalar
	#>
	
	foreach (@_[0..2]) {
		silico_msg('d', "Error in input") if !defined $_;
	}

	return  sqrt($_[0]**2 + $_[1]**2 + $_[2]**2);
}

sub unit_normal_vector {

	#<
	#? Return a unit vector normal to a given vector
	#. Vector will also be normal to X, Y or Z axis
	#; Requires: vector
	#; Returns: vector
	#>
	
	my ($x, $y, $z) = @_;
	
	my @v;
	
	silico_msg('d', "No normal to a zero vector\n") if $x == 0 && $y == 0 && $z == 0;
	
	if ($x >= $y && $x >= $z) {
		@v = (0,1,0);
	} elsif ($y >= $x && $y >= $z) {
		@v = (0,0,1);
	} else {
		@v = (1,0,0)
	}
	
	return unit_vector(crossprod(@v,  $x, $y, $z));
	
}

sub quaternion {
	
	#<
	#? Calculates a quaternion - a 'rotation matrix'
	#  when you want to rotate a point around a vector.
	#; Requires: Angle, rotation vector
	#; Returns: Rotation matrix as a 9-long (3x3) array
	#>
	
	my $angle = $_[0];
	my ($x1, $y1, $z1) = @_[1..3];
	
	my @mat;
	my ($q0, $q1, $q2, $q3);
	my $size;
	my ($xn, $yn, $zn);
	
	$size = sqrt($x1**2 + $y1**2 + $z1**2);
	
	# Die if the vector has no length
	silico_msg('d', "Zero rotation vector\n") if ($size == 0);
	
	# Make axis unit vector
	$xn = $x1/$size;
	$yn = $y1/$size;
	$zn = $z1/$size;
	
	# Calculate rotation matrix
	#
	# Given angle r in radians and unit vector u = ai + bj + ck or [a,b,c]', define:
	# q0 = cos(r/2),  q1 = sin(r/2) a,  q2 = sin(r/2) b,  q3 = sin(r/2) c
	# The rotation matrix is:
	# q0q0 + q1q1 - q2q2 - q3q3	2(q1q2 - q0q3)			2(q1q3  + q0q2)
	# 2(q1q2 + q0q3)		q0q0 - q1q1 + q2q2 - q3q3	2(q3q2 - q0q1)
	# 2(q3q1 - q0q2)		2(q3q2 + q0q1)			q0q0 - q1q1 - q2q2 + q3q3
	
	$q0 = cos($angle/2);
	$q1 = sin($angle/2) * $xn;
	$q2 = sin($angle/2) * $yn;
	$q3 = sin($angle/2) * $zn;
	
	@mat = (
		$q0*$q0 + $q1*$q1 - $q2*$q2 - $q3*$q3, 2*($q1*$q2 - $q0*$q3), 2*($q1*$q3 + $q0*$q2),
		2*($q1*$q2 + $q0*$q3), $q0*$q0 - $q1*$q1 + $q2*$q2 - $q3*$q3, 2*($q3*$q2 - $q0*$q1),
		2*($q3*$q1 - $q0*$q2), 2*($q3*$q2 + $q0*$q1), $q0*$q0 - $q1*$q1 - $q2*$q2 + $q3*$q3);
	
	return @mat;
}

sub scale_vector {
	
	#<
	#? Multiply a vector by a scalar
	#; Requires: scalar, vector
	#; Returns: vector
	#>
	
	my $s = shift;
	my @vec = @_;
	
	foreach (@vec) {
		
		$_ *= $s;
	}
	
	return @vec;
}

sub vector_add {
	
	#<
	#? Add two vectors (length 3)
	#; Requires: vectors, vector
	#; Returns: vector
	#>
	
	return ($_[0]+$_[3], $_[1]+$_[4], $_[2]+$_[5]);
}

sub unit_vector {
	
	#<
	#? Normalises a vector, providing a unit vector.
	#; Requires: array of three elements (x, y and z components)
	#; Returns: unit vector as an array
	#>
	
	my ($ax, $ay, $az) = @_;
	
	my $length;
	my ($ux, $uy, $uz);
	my @u;
	
	$length = magnitude ($ax, $ay, $az);
	
	silico_msg('d', "Error in input") if !defined $length;
	
	return (0, 0, 0) if ($length == 0);
	
	$ux = $ax / $length; 
	push @u, $ux;
	$uy = $ay / $length; 
	push @u, $uy;
	$uz = $az / $length; 
	push @u, $uz;
	
	return @u;
}

sub vector {
	
	#<
	#? Calculates a vector from two Cartesian points.
	#; Requires: coord1, coord2
	#; Returns: vector as an array
	#>
	
	my ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
	
	return ($x2-$x1, $y2-$y1, $z2-$z1);
}

sub atom_vector {

	#<
	#? Calculates a vector from two atoms (a1 - a2)
	#; Requires: coord1, coord2
	#; Returns: vector as an array
	#>
	
	return ($_[0]->{X} - $_[1]->{X}, $_[0]->{Y} - $_[1]->{Y}, $_[0]->{Z} - $_[1]->{Z});
}



################
# TRIGONOMETRY #
################



sub deg_to_rad {
	
	#<
	#? Convert degrees to radians
	#; Requires: scalar
	#; Returns: scalar
	#>
	
	$_[0]*$Pi/180;
}

sub rad_to_deg {
	
	#<
	#? Convert radians to degrees
	#; Requires: scalar
	#; Returns: scalar
	#>
	
	$_[0]*180/$Pi;
}

##################
# Periodic cells #
##################

sub mol_check_cell {

	#<
	#? Check if periodic cell exists and create a functional one if it does not.
	#; Requires: molecule, flag to print cell
	#; Returns: $cell
	#>
	
	my $mol = $_[0];
	my $print = $_[1];

	# Set molecule dimensions if no cell is defined or the cell is obviously rubbish
	if (!defined $mol->{CELL}  || $mol->{CELL}[0] < 2) {
		
		@{$mol->{CELL}} = molecule_dimensions($mol, 10);
		
		$mol->{CELL_ALPHA} = 90;
		$mol->{CELL_BETA} = 90;
		$mol->{CELL_GAMMA} = 90;
		$mol->{SPACE_GROUP} = "0011";
		$print = 1;
	}
	
	foreach my $val (qw(CELL_ALPHA CELL_BETA CELL_GAMMA)) {
	
		next if defined $mol->{$val};
		$mol->{$val} = 90;
		$mol->{SPACE_GROUP} = "0011";
		$print = 1;
	}
		
	pdb_cell_to_gromacs_bv($mol) if !$mol->{GROMACS_BOX_VECTORS} || !$mol->{BV_MATINV};
	
	silico_msg('c', "Setting system dimensions: abc: @{$mol->{CELL}} alpha: $mol->{CELL_ALPHA} beta: $mol->{CELL_BETA} gamma: $mol->{CELL_GAMMA} sg: $mol->{SPACE_GROUP}\n") if $print;
}


sub pdb_cell_to_gromacs_bv {

	#<
	#? Create gromacs box vectors
	#; Requires; molecule
	#>

        use POSIX qw(cos sin sqrt);

        my $mol = $_[0];
	
        # Extract cell lengths and angles
        my ($a, $b, $c) = @{ $mol->{CELL} }[ 0 .. 2 ];
        my $alpha = $mol->{CELL_ALPHA} || silico_msg('d', "Cell alpha not set\n");    # degrees
        my $beta  = $mol->{CELL_BETA}  || silico_msg('d', "Cell beta not set\n");
        my $gamma = $mol->{CELL_GAMMA} || silico_msg('d', "Cell gamma not set\n");

        # Convert angles to radians
        my $alpha_rad = $alpha * 3.141592653589793 / 180.0;
        my $beta_rad  = $beta * 3.141592653589793 / 180.0;
        my $gamma_rad = $gamma * 3.141592653589793 / 180.0;

        # Build GROMACS box vectors
        my $ax  = $a;
        my $ay  = 0;
        my $az  = 0;
        my $bx  = $b * cos($gamma_rad);
        my $by  = $b * sin($gamma_rad);
        my $bz  = 0;
        my $cx  = $c * cos($beta_rad);
        my $cy  = $c * (cos($alpha_rad) - cos($beta_rad) * cos($gamma_rad)) / sin($gamma_rad);
        my $cz2 = $c**2 - $cx**2 - $cy**2;
        my $cz  = $cz2 > 0 ? sqrt($cz2) : 0;

	foreach ($ax, $ay, $az, $bx, $by, $bz, $cx, $cy, $cz) {
	
		$_ = sprintf "%8.4f", $_;
	}

        # Store box vectors
	# Note box vectors are in nm
        # Note ideosyncratic ordering
	$mol->{GROMACS_BOX_VECTORS}= [ $ax/10, $by/10, $cz/10, $ay/10, $az/10, $bx/10, $bz/10, $cx/10, $cy/10 ];

        mol_bv_mat($mol);
	
	#silico_msg('c', "Creating Gromacs box vectors: @{$mol->{GROMACS_BOX_VECTORS}}\n\n");
}

sub mol_bv_mat {

	#<
	#? Return gromacs box vectors as a transformation matricies in Ang
	#; Requires: molecule 
	#; Returns: nothing
	#; Sets: mol->{BV_MAT},  mol->{BV_MATINV}
	#>

	my $mol = $_[0];
	
	my $bv = $mol->{GROMACS_BOX_VECTORS};
	
	croak("Box vectors not set") if !defined $bv;
	
	return if (defined $mol->{BV_MAT} && defined $mol->{BV_MATINV});
	
	# GROMACS box vectors (in nm converted to Ang)
	my @a = ($bv->[0] * 10, $bv->[3] * 10, $bv->[4] * 10);
	my @b = ($bv->[5] * 10, $bv->[1] * 10, $bv->[6] * 10);
	my @c = ($bv->[7] * 10, $bv->[8] * 10, $bv->[2] * 10);

	# Box matrix and its inverse
	my ($mat, $matinv);
	@$mat  = (@a, @b, @c);
	@$matinv = inverse3x3(@$mat);
	
	$mol->{BV_MAT} = $mat; 
	$mol->{BV_MATINV} = $matinv;
	
	return;
}

sub inverse3x3 {

	my ($a11, $a12, $a13, $a21, $a22, $a23, $a31, $a32, $a33) = @_;
	
	#print "x $a11, $a12, $a13, $a21, $a22, $a23, $a31, $a32, $a33\n";
	
	croak("Missing array element") if (!defined $a33);

	my $det = $a11 * ($a22 * $a33 - $a23 * $a32) - $a12 * ($a21 * $a33 - $a23 * $a31) + $a13 * ($a21 * $a32 - $a22 * $a31);

	croak("Matrix determinant is zero!") if abs($det) < 1e-12;

	my @inv = (
		($a22 * $a33 - $a23 * $a32) / $det,
		($a13 * $a32 - $a12 * $a33) / $det,
		($a12 * $a23 - $a13 * $a22) / $det,
		($a23 * $a31 - $a21 * $a33) / $det,
		($a11 * $a33 - $a13 * $a31) / $det,
		($a13 * $a21 - $a11 * $a23) / $det,
		($a21 * $a32 - $a22 * $a31) / $det,
		($a12 * $a31 - $a11 * $a32) / $det,
		($a11 * $a22 - $a12 * $a21) / $det
	);

	return @inv;
}

sub matvec {

	my ($m, $v) = @_;

	my @r;
	$r[0] = $m->[0] * $v->[0] + $m->[1] * $v->[1] + $m->[2] * $v->[2];
	$r[1] = $m->[3] * $v->[0] + $m->[4] * $v->[1] + $m->[5] * $v->[2];
	$r[2] = $m->[6] * $v->[0] + $m->[7] * $v->[1] + $m->[8] * $v->[2];

	return @r;
}


#############
# DISTANCES #
#############

sub distance {
	
	#<
	#? Find the distance between two atoms.
	#; Requires: atom1, atom2.
	#; Returns: distance.
	#>
	
	#my $atom1 = $_[0];
	#my $atom2 = $_[1];
	
	return sqrt ( ($_[0]->{X} - $_[1]->{X})**2 + ($_[0]->{Y} - $_[1]->{Y})**2 + ($_[0]->{Z} - $_[1]->{Z})**2);
}

sub distance_points {
	
	#<
	#? Find the distance between two sets of Cartesian points.
	#; Requires: x1, y1, z1, x2, y2, z2
	#;
	#>
	
	#my @point1 = @_[0..2];
	#my @point2 = @_[3..5];
	
	return sqrt(($_[0] - $_[3])**2 + ($_[1] - $_[4])**2 + ($_[2] - $_[5])**2);
}

sub distance_sq {

	#<
	#? Find the square of the distance between two atoms.
	#; Requires: atom1, atom2.
	#; Returns: square of distance.
	#>
	
	#my $atom1 = $_[0];
	#my $atom2 = $_[1];
	
	return ($_[0]->{X} - $_[1]->{X})**2 + ($_[0]->{Y} - $_[1]->{Y})**2 + ($_[0]->{Z} - $_[1]->{Z})**2;
}

sub distance_atoms_min {
	
	#<
	#? Find the shortest distance between one atom and an array of atoms
	#; Requires: atom, array of atoms
	#; Returns: shortest distance
	#>
	
	my $atom1 = shift;
	my @atoms2 = @_;

	my $d2;

	my $x = $atom1->{X};
	my $y = $atom1->{Y};
	my $z = $atom1->{Z};

	foreach my $atom2 (@atoms2) {

		my $v = distance_sq($atom1, $atom2);
		$d2 = $v if (!defined $d2 || $v < $d2);
	}

	return sqrt $d2;
}

sub distance_point_atoms_min {
	
	#<
	#? Find the shortest distance between a point and an array of atoms
	#; Requires: x, y, z, array of atoms
	#; Returns: shortest distance
	#>
	
	my $x = shift;
	my $y = shift;
	my $z = shift;
	my @atoms2 = @_;

	my $d2;

	foreach my $atom2 (@atoms2) {

		my $v =  ( ($x - $atom2->{X})**2 + ($y - $atom2->{Y})**2 + ($z - $atom2->{Z})**2);
		$d2 = $v if (!defined $d2 || $v < $d2);
	}

	return sqrt $d2;
}

sub general_distance_sq {
	
	#<
	#? Find the square of the distance between two atoms.
	#  Version that uses periodic informatin if available
	#; Requires: atom1, atom2, molecule
	#; Returns: distance_squared.
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $mol   = $_[2] || croak();

	if ($mol->{GROMACS_BOX_VECTORS}) {
	
		mol_bv_mat($mol) if (!$mol->{BV_MAT});
	
		my @r =  distance_periodic_xyz_mat($atom1->{X}, $atom1->{Y}, $atom1->{Z}, $atom2->{X}, $atom2->{Y}, $atom2->{Z}, $mol->{BV_MAT}, $mol->{BV_MATINV});
		
		return ($r[0]**2 + $r[1]**2 + $r[2]**2);
	
	} elsif ($mol->{CELL}) {
	
		return distance_periodic_xyz($atom1, $atom2, $mol->{CELL});
	
	} else {
	
		return distance_sq($atom1, $atom2);
	}
}

sub general_distance_xyz {
	
	#<
	#? Find the x, y, z components of vector between two atoms
	#  Version that uses periodic informatin if available
	#; Requires: atom1, atom2, molecule
	#; Returns: x, y, z
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $mol   = $_[2] || croak();

	if ($mol->{GROMACS_BOX_VECTORS}) {
	
		mol_bv_mat($mol) if (!$mol->{BV_MAT});
	
		return  distance_periodic_xyz_mat($atom1->{X}, $atom1->{Y}, $atom1->{Z}, $atom2->{X}, $atom2->{Y}, $atom2->{Z}, $mol->{BV_MAT}, $mol->{BV_MATINV});
		
	} elsif ($mol->{CELL}) {
	
		return distance_periodic_sq($atom1, $atom2, $mol->{CELL});
	
	} else {
	
		silico_msg('d', "Not implemented for non-periodic systems yet");
	}
}

sub distance_periodic_sq {
	
	#<
	#? Find the square of the distance between two atoms.
	#  For a system that has cubic, on-axis, periodic boundary conditions
	#; Requires: atom1, atom2, $cell->[cellx, celly, cellz]
	#; Returns: distance_squared.
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $cell  = $_[2] || croak();
	
	my $i = 0;
	my $rval = 0;
	
	foreach ('X', 'Y', 'Z') {
	
		my $c = $cell->[$i];
		my $val = abs($atom1->{$_} - $atom2->{$_});
		$val = $val - floor(($val+$c/2)/$c)*$c; 
		$rval += $val**2;
		++$i;
	}
	
	return $rval;
}

sub distance_periodic_sq_bv {
	
	#<
	#? Find the square of the distance between two atoms.
	#  For cuboid systems with gromacs box vectors
	#; Requires: atom1, atom2, $cell->[cellx, celly, cellz]
	#; Returns: distance_squared.
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $bv  = $_[2] || croak();
	
	my @r =  distance_periodic_points_xyz_bv($atom1->{X}, $atom1->{Y}, $atom1->{Z}, $atom2->{X}, $atom2->{Y}, $atom2->{Z}, $bv);
	
	return $r[0]**2 + $r[1]**2 + $r[2]**2;
}

sub distance_periodic_xyz_bv {
	
	#<
	#? Find the square of the distance between two atoms.
	#  For a system using Gromacs box vectors
	#; Requires: atom1, atom2, gromacs_box_vectors
	#; Returns: x, y, z
	#>
	
	my $atom1 = $_[0];
        my $atom2 = $_[1];
        my $bv   = $_[2] || croak();

	my @r =  distance_periodic_points_xyz_bv($atom1->{X}, $atom1->{Y}, $atom1->{Z}, $atom2->{X}, $atom2->{Y}, $atom2->{Z}, $bv);

	return @r;
}

sub distance_periodic_points_xyz_mat {

	#<
	#? Find the x y and components of the vector between two points
	#  For a system using Gromacs box vectors
	#; Requires: point1, point2, ptr to mat matrix, ptr to inverse matrix 
	#; Returns: x, y, z
	#>

	my (@r1, @r2, $mat, $matinv);
	($r1[0], $r1[1], $r1[2], $r2[0], $r2[1], $r2[2], $mat, $matinv) = @_;

	croak("Subroutine missing mat or matinv") if !($mat && $matinv);
	
	# Displacement in Cartesian
	my $dr;
	@$dr = ($r2[0] - $r1[0], $r2[1] - $r1[1], $r2[2] - $r1[2]);

	# Convert displacement to fractional coordinates
	my $s;
	@$s = matvec($matinv, $dr);

	# Apply minimum image convention
	for my $i (0 .. 2) {
		$s->[$i] -= sprintf("%.0f", $s->[$i]);    # equivalent to nearest integer
	}

	# Back to Cartesian
	my @dr_min = matvec($mat, $s);

	return @dr_min;
}



sub distance_periodic_xyz {
	
	#<
	#? Find the X, Y and Z components of distance between two atoms.
	#  For a system that has cubic, on-axis, periodic boundary conditions
	#. Note: values are always positive.
	#; Requires: atom1, atom2, $cell->[cellx, celly, cellz]
	#; Returns: periodic distances x, y, z
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $cell  = $_[2] || croak();
	
	my $i = 0;
	my @r;
	
	foreach ('X', 'Y', 'Z') {
	
		my $c = $cell->[$i];
		my $val = abs($atom1->{$_} - $atom2->{$_});
		$r[$i] = $val - floor(($val+$c/2)/$c)*$c; 
		++$i;
	}
	
	return @r;
}

sub distance_periodic_points_xyz {
	
	#<
	#? Find the X, Y and Z components of distance between two points
	#  For a system that has cubic, on-axis, periodic boundary conditions
	#. Note: values are always positive.
	#; Requires: 6 coordinates, cell
	#; Returns: periodic distances x, y, z
	#>
	
	my @p1 = @_[0..2];
	my @p2 = @_[3..5];
	my $cell = $_[6] || croak();
	
	my $i = 0;
	my @r;
	
	foreach (0..2) {
	
		my $c = $cell->[$i];
		my $d = abs($p1[$i] - $p2[$i]);
		$r[$i] = $d - floor(($d+$c/2)/$c)*$c; 
		++$i;
	}
	
	return @r;
}

sub distance_periodic_signed_xyz {
	
	#<
	#? Find the X, Y and Z components of distance between two atoms.
	#  For a system that has cubic, on-axis, periodic boundary conditions
	#. Note: values are signed.
	#; Requires: atom1, atom2, $cell->[cellx, celly, cellz]
	#; Returns: periodic distances x, y, z
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $cell  = $_[2] || croak();
	
	my @r;
	
	my $i = 0;
	foreach ('X', 'Y', 'Z') {
	
		my $c = $cell->[$i];
		my $val = $atom1->{$_} - $atom2->{$_};
		$r[$i] = $val - floor(($val+$c/2)/$c)*$c; 
		++$i;
	}
	
	return @r;
}

sub distance_periodic_xyz_old {
	
	#<
	#? Find the X, Y and Z components of distance between two atoms.
	#  For a system that has cubic, on-axis, periodic boundary conditions
	#; Requires: atom1, atom2, $cell->[cellx, celly, cellz]
	#; Returns: periodic distances x, y, z
	#>
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $cell  = $_[2] || croak("Cell not defined");
	
	my $cx2 = $cell->[0]/2;
	my $cy2 = $cell->[1]/2;
	my $cz2 = $cell->[2]/2;
		
	my $xd = abs($atom1->{X} - $atom2->{X});
	if ($cx2 > 0) {
		while ($xd > $cx2) {$xd -= $cell->[0]}
	} else {
		silico_msg('w', "Periodic cell is 0 or negative in X-dimension.\n",
				"Distance in infinite space will be used.\n");
	}
	
	my $yd = abs($atom1->{Y} - $atom2->{Y});
	if ($cy2 > 0) {
		while ($yd > $cy2) {$yd -= $cell->[1]}
	} else {
		silico_msg('w', "Periodic cell is 0 or negative in Y-dimension.\n",
				"Distance in infinite space will be used.\n");
	}
	
	my $zd = abs($atom1->{Z} - $atom2->{Z});
	if ($cz2 > 0) {
		while ($zd > $cz2) {$zd -= $cell->[2]}
	} else {
		silico_msg('w', "Periodic cell is 0 or negative in Z-dimension.\n",
				"Distance in infinite space will be used.\n");
	}
	
	return (abs($xd), abs($yd), abs($zd));
}


sub rms {
	
	#<
	#? Find the rms distance between two molecules (heavy atoms by default).
	#. Note. Not periodic
	#; Requires: Two molecules, flag to only calculate for heavy atoms
	#; Returns: rms or -ve numbers if something is wrong.
	#>
	
	my $mol1 = $_[0];
	my $mol2 = $_[1];
	my $heavy = $_[2] || 0;
	
	my $count = 0;
	my $sumsq = 0;
	my $warn;
	
	my $name1 = $mol1->{NAME} || 'Unnamed';
	my $name2 = $mol2->{NAME} || 'Unnamed';
	
	# Sanity
	if ($mol1->{NUMATOMS} < 1 || $mol2->{NUMATOMS} < 1) {
		silico_msg('e', "At least one molecule has no atoms!\n",
			"\tMolecule 1: '$name1' $mol1->{NUMATOMS} atoms; Molecule 2: '$name2' $mol2->{NUMATOMS} atoms\n");
		return -2;
	}
	
	# Check that both molecules have the same number of atoms
	if ($heavy) {
		if (mol_count_heavy_atoms($mol1) != mol_count_heavy_atoms($mol2)) {
			silico_msg('w', "Comparison molecules have different numbers of heavy atoms!\n",
				"\tMolecule 1: '$name1' $mol1->{NUMATOMS} atoms; Molecule 2: '$name2' $mol2->{NUMATOMS} atoms\n");
			return -1;
		}
	} else {
		if ($mol1->{NUMATOMS} != $mol2->{NUMATOMS}) {
			silico_msg('w', "Comparison molecules have different numbers of atoms!\n",
				"\tMolecule 1: '$name1' $mol1->{NUMATOMS} atoms; Molecule 2: '$name2' $mol2->{NUMATOMS} atoms\n");
			return -1;
		}
	}
	
	for (my $i=0; $i < $#{$mol1->{ATOMS}}; ++$i) {	
		next if ($heavy && $mol1->{ATOMS}[$i]{ELEMENT} eq 'H');
		if ($mol1->{ATOMS}[$i]{ELEMENT} ne $mol2->{ATOMS}[$i]{ELEMENT}) {
			print "$mol1->{ATOMS}[$i]{ELEMENT} ne $mol2->{ATOMS}[$i]{ELEMENT}\n";
			++$warn;
		}
		$sumsq += distance_sq ($mol1->{ATOMS}[$i], $mol2->{ATOMS}[$i]);
		++$count;
	}
	
	if ($warn) {
		silico_msg('w', "Atomic elements of $warn atoms did not match in RMS comparison!\n\tMolecules: $name1, $name2\n");
		molecule_printout($mol1);
		molecule_printout($mol2);
		return -3;
	}

	return	sqrt($sumsq/$count);
}

sub rms_atomlist {
	
	#<
	#? Find the rms distance between atomlists.
	#. Note: Not periodic
	#; Requires: Two atomlists, flag to return sum of squares rather than rms
	#; Returns: rms or -ve numbers if something is wrong.
	#>
	
	my $list1 = $_[0];
	my $list2 = $_[1];
	my $flag = $_[2];
	my $count = 0;
	
	my $sumsq = 0;
	
	if ($#{$list1} < 1 || $#{$list2} < 1) {
		silico_msg('e', "Lists have no atoms!\n");
		return -2;
	}
	
	# Check that both lists have the same number of atoms
	if ($#{$list1} != $#{$list2}) {
		silico_msg('e', "Comparison list have different numbers of atoms! ",$#{$list1}+1, $#{$list2}+1,"\n" );
		return -1;
	}
	
	for (my $i=0; $i < $#{$list1}; ++$i) {
		$sumsq += distance_sq ($list1->[$i], $list2->[$i]);
		++$count;
	}
	
	return  $sumsq if $flag;
	return	sqrt($sumsq/$count);
}

sub calc_binsize {

	#<
	#? Determine size of bins in cell given an approximate binsize
	#; Requires: molecule, approx binsize(default 2 ang)
	#; Returns:  binsize in each dimension (array), number of bins in each dimension (array)
	#; Sets: mol->{BINSIZE} and mol->{NUMBINS};
	#>
	
	my $mol = $_[0] || croak("Subroutine missing mol");
	my $approx_binsize = get_lflag('binsize') || $_[1] || 2;
	
	# Shortcut
	return ($mol->{BINSIZE}, $mol->{NUMBINS}) if (defined $mol->{BINSIZE} && defined $mol->{NUMBINS});

	use POSIX qw(floor);
	
	mol_check_cell($mol);
	
	my $cell = $mol->{CELL};
	my $binsize;
	my $numbins;
	my $max_numbins = 100000; # Code gets slow if numbins is too large. Overridden by binsize flag
	my $tot;
	
	#print "Cell: @$cell\n";

	# Set binsize and numbins in each dimension
	# binsize has a minimum value of 2 Ang
	while (1) {
	
		foreach my $i (0..2) {
	
			$binsize->[$i] = 0;
			my $extra = 0;
			while ($binsize->[$i] < 2) { 
		
				$numbins->[$i] = floor($cell->[$i]/($approx_binsize+$extra)) || 1;
				$binsize->[$i] = $cell->[$i]/$numbins->[$i];
				last if $numbins->[$i] == 1;
			
				#print "cell $cell->[$i] extra $extra n $numbins->[$i] b $binsize->[$i]\n";
			
				$extra += 0.5;
			}
		}
		
		$tot = $numbins->[0] * $numbins->[1] * $numbins->[2];
		
		# Finish if binsize is set externally
		last if get_lflag('binsize');
		
		# Finish if we dont have too many bins
		last if $tot <= $max_numbins;
		
		++$approx_binsize;
		#print "approx_binsize $approx_binsize\n";
	}
	
	printf "Binsize: %.2f %.2f %.2f  Numbins: @$numbins  Total bins: $tot\n", @$binsize;
	silico_msg('w', "This system is large, which might slow things down\n") if $binsize->[0] > 5;
	
	$mol->{BINSIZE} = $binsize;
	$mol->{NUMBINS} = $numbins;
	
	# Initialise bins to empty arrays if not defined
	if (!defined $mol->{BINS}) {
		my $bins;
		for my $i (0..$numbins->[0]-1) {
			for my $j (0..$numbins->[1]-1) {
				for my $k (0..$numbins->[2]-1) {
					my $a;
					@$a = ();
					$bins->[$i][$j][$k] = $a;
				}
			}
		}
		
		$mol->{BINS} = $bins;
	}

	return ($binsize, $numbins);
}

sub molecule_bin_atoms_old {

	#<
	#? Bin atoms 
	#. A faster bin atoms routine. This routine can optionally set atom->{BIN}
	#; Requires: molecule, atomlist (optional), binsize (optional), numbins array (optional), flag to add bin to atom (optional)
	#; Returns: array of binned atoms
	#>
	
	my $mol = $_[0] || croak("Subroutine missing mol");
	my $atomlist = $_[1];
	my $binsize = $_[2];
	my $numbins = $_[3];
	my $add_bin_to_atom = $_[4];
	
	use POSIX qw(floor);
	
	# Use atomlist if one has been provided
	my $atoms;
	if ($atomlist) {
		$atoms = $atomlist;
	} else {
		$atoms = $mol->{ATOMS};
	}
	
	if (!defined $binsize || !defined $numbins) {
		molecule_check_cell($mol);
		($binsize, $numbins) = calc_binsize($mol);
	}
	
	# Initialise bins to empty arrays if not defined
	my $bins = $mol->{BINS};
	if (!defined $bins) {
		for my $i (0..$numbins->[0]-1) {
			for my $j (0..$numbins->[1]-1) {
				for my $k (0..$numbins->[2]-1) {
					my $a;
					@$a = ();
					$bins->[$i][$j][$k] = $a;
				}
			}
		}
	}

	foreach my $atom (@$atoms) {

		my $bin = atom_calc_bin($atom, $mol);
		push @{$bins->[$bin->[0]][$bin->[1]][$bin->[2]]}, $atom;
		
		#silico_msg('e', "$atom->{X} $atom->{Y} $atom->{Z}\n") if ($bin->[0] < 0 || $bin->[1] < 0 || $bin->[2] < 0) ;
		#silico_msg('e', "$atom->{X} $atom->{Y} $atom->{Z}\n") if ($bin->[0] > $numbins-1 || $bin->[1] > $numbins-1 || $bin->[2] > $numbins-1) ;
		
		if ($add_bin_to_atom) {
			$atom->{BINX} = $bin->[0];
			$atom->{BINY} = $bin->[1];
			$atom->{BINZ} = $bin->[2];
		}
	}
	
	$mol->{BINS} = $bins;

	return $bins;
}

sub molecule_bin_atoms {

	#<
	#? Bin atoms 
	#. A faster bin atoms routine. This routine can optionally set atom->{BIN}
	#; Requires: molecule, flag to add bin to atom 
	#; Returns: array of binned atoms
	#>
	
	my $mol = $_[0] || croak("Subroutine missing mol");
	my $atomlist = $_[1];
	my $add_bin_to_atom = $_[2];
	
	use POSIX qw(floor);
	
	# Use atomlist if one has been provided
	my $atoms;
	if ($atomlist) {
		$atoms = $atomlist;
	} else {
		$atoms = $mol->{ATOMS};
	}
	
	mol_check_cell($mol);
		
	calc_binsize($mol);
	
	my $binsize = $mol->{BINSIZE};
	my $numbins = $mol->{NUMBINS};
	my $bins = $mol->{BINS};
	
	foreach my $atom (@$atoms) {

		my $bin = atom_calc_bin($atom, $mol);
		push @{$bins->[$bin->[0]][$bin->[1]][$bin->[2]]}, $atom;
		
		#silico_msg('e', "$atom->{X} $atom->{Y} $atom->{Z}\n") if ($bin->[0] < 0 || $bin->[1] < 0 || $bin->[2] < 0) ;
		#silico_msg('e', "$atom->{X} $atom->{Y} $atom->{Z}\n") if ($bin->[0] > $numbins-1 || $bin->[1] > $numbins-1 || $bin->[2] > $numbins-1) ;
		
		if ($add_bin_to_atom) {
			$atom->{BINX} = $bin->[0];
			$atom->{BINY} = $bin->[1];
			$atom->{BINZ} = $bin->[2];
		}
	}
	
	$mol->{BINS} = $bins;

	return $bins;
}

sub atom_calc_bin {

	#<
	#? Determine which bin an atom belongs to
	#. Triclinic version
	#; Requires: atom, molecule
	#; Returns:  bin (array);
	#>

    	my $atom    = $_[0];
	my $mol = $_[1];
	
	my $matinv = $mol->{BV_MATINV} || (molecule_printout($mol) && die);
	my $numbins = $mol->{NUMBINS} || croak("Subroutine issing numbins");

    	# Convert to fractional coordinates f = MATINV * r
    	my @f;
    	for my $i (0..2) {
        	$f[$i] = $matinv->[$i]*$atom->{X} + $matinv->[$i+3]*$atom->{Y} + $matinv->[$i+6]*$atom->{Z};
       	 	$f[$i] -= POSIX::floor($f[$i]);  # wrap into [0,1)
    	}

    	# Convert fractional coordinates to bin indices
    	my $bin;
    	for my $i (0..2) {
	
        	my $b = int($f[$i] * $numbins->[$i]);
       	 	$b = 0  if $b < 0;
        	$b = $numbins->[$i]-1 if $b >= $numbins->[$i];
        	$bin->[$i] = $b;
    	}

    	return $bin;
}


sub atom_calc_bin_ {

	#<
	#? Determine which bin an atom belongs to
	#; Requires: atom, binsize, number of bins in each dimension
	#; Returns:  bin (array);
	#>

	my $atom = $_[0];
	my $binsize = $_[1];
	my $numbins = $_[2];

	my $bin;
	my $i = -1;

	# Calculate the bin
	foreach ('X', 'Y', 'Z') {
		
		++$i;
		
		my $c = $atom->{$_}; 			# coordinate
		my $b = floor($c / $binsize->[$i]); 	# bin
		$bin->[$i] = $b % $numbins->[$i];	# use modulo to wrap bin into cell
	}    	

	return $bin;
}

sub find_closest_atoms_periodic {

        #<
        #? Find closest atoms in two molecules, taking into account
        #  periodic boundary conditions.
        #. Simple-minded algorithm for testing stuff
        #; Requires: molecule1, molecule2, maximum distances, flag to ignore hydrogens
        #; Returns: sets of atom pairs with distances
	#>

        my $mol1 = $_[0];               # Molecule 1
        my $mol2 = $_[1];               # Molecule 2
        my $maxdist = $_[2] || 10;      # Maximum distance
        my $ignh = $_[3];               # Ignore hydrogens

        my $self_comparison;
        # Determine if we are comparing a molecule to itself
        if ($mol1 == $mol2) {
                $self_comparison = 1;
        }

        my $cell = molecule_check_cell($mol2);

        my $short;
        @$short = ();
        my $flag = 1;
        my $count = 0;
        my $factor = 10;
        foreach my $atom1 (@{$mol1->{ATOMS}}) {

                next if $ignh && $atom1->{ELEMENT} eq 'H';

                ++$count;
                if ($count % ($factor) == 0) {
                        print "atom $count\n" if $count % $factor == 0;
                        $factor *= 10 if $count == $factor;
                }

                foreach my $atom2 (@{$mol2->{ATOMS}}) {

                        next if $self_comparison && ($atom1->{NUM} >= $atom2->{NUM});

                        next if $ignh && $atom2->{ELEMENT} eq 'H';

                        my $dist = sqrt(general_distance_sq($atom1, $atom2, $mol2));
			
			if ($dist <= $maxdist) {
				
				my $v;
				@$v = ($atom1, $atom2, $dist);
				push @$short, $v;
				@$short = sort {$b->[2] <=> $a->[2]} @$short;
							
				#foreach (@$short) {
				#	printf "%0.3f ", $_->[2];
				#}
				#print "\n";
			}
                }
        }
	
	@$short = sort {$a->[2] <=> $b->[2]} @$short;
        return $short;
}

sub find_closest_atoms_periodic_binned {

	#<
	#? Find closest atoms in two molecules, taking into account
	#  periodic boundary conditions.
	#. Faster algorithm 
	#; Requires: molecule1, molecule2, maximum distances, flag to ignore hydrogens
        #; Returns: sets of atom pairs with distances
	#>

	my $mol1 = $_[0]; 		# Molecule 1
	my $mol2 = $_[1]; 		# Molecule 2
	my $maxdist = $_[2] || 0.5;      # Maximum distance
	my $ignh = $_[3]; 		# Ignore hydrogens
	
	my $self_comparison;
	# Determine if we are comparing a molecule to itself
	if ($mol1 == $mol2) {
		$self_comparison = 1;
	}
	
	my $cell = molecule_check_cell($mol2);
	my ($binsize, $numbins) = calc_binsize($mol2);
	my $bins = molecule_bin_atoms($mol2, undef $binsize, $numbins);
	mol_bv_mat($mol2);
	
	my $factor = 1000;
	my $short;
	foreach my $atom1 (@{$mol1->{ATOMS}}) {
	
		next if $ignh && $atom1->{ELEMENT} eq 'H';

                my $x1 = $atom1->{X};
                my $y1 = $atom1->{Y};
                my $z1 = $atom1->{Z};

		my $bin = atom_calc_bin($atom1, $mol2);

                # Make lists of bins within +/- one bin of the current bin
                my $blist;
                foreach my $i (0..2) {
                        @{$blist->[$i]} = ($bin->[$i]-1..$bin->[$i]+1);
                        foreach my $n (@{$blist->[$i]}) {
                                $n = $n % $numbins->[$i]; # Modulo
                       }
                }

                foreach my $bx (@{$blist->[0]}) {
                        foreach my $by (@{$blist->[1]}) {
                                 foreach my $bz (@{$blist->[2]}) {

                                        foreach my $atom2 (@{$bins->[$bx][$by][$bz]}) {
			
						next if $self_comparison && ($atom1->{NUM} >= $atom2->{NUM});
			
						next if $ignh && $atom2->{ELEMENT} eq 'H';
			
						my $dist = sqrt(distance_periodic_sq($atom1, $atom2, $cell));
					
						if ($dist <= $maxdist) {
				
							my $v;
							@$v = ($atom1, $atom2, $dist);
							push @$short, $v;
							@$short = sort {$b->[2] <=> $a->[2]} @$short;
							
							#foreach (@$short) {
							#	printf "%0.3f ", $_->[2];
							#}
							#print "\n";
						}
					}
				}
			}
		}
	}
	
	@$short = sort {$a->[2] <=> $b->[2]} @$short;
	return $short;
}


sub clash_periodic {

	#<
	#? Find clashes between two molecules, taking into account
	#  periodic boundary conditions.
	#. Simple-minded algorithm for testing stuff
	#; Requires: molecule1, molecule2, maximum number of clashes, vdW scale factor for added atoms
	#; Returns: Number of clashes (up to maximum if this is set), pairs of violating atoms and distances

	my $mol1 = $_[0]; 		# Molecule 1
	my $mol2 = $_[1]; 		# Molecule 2
	my $maxclash = $_[2] || 1; 	# Maximum number of allowed clashes
	my $scale = $_[3] || 1; 	# vdW scaling factor for packed atoms
	
	my $clashcount = 0;
	my $violations;
	
	my $self_comparison;
	# Determine if we are comparing a molecule to itself
	if ($mol1 == $mol2) {
		$self_comparison = 1;
	}
	
	foreach my $atom1 (@{$mol1->{ATOMS}}) {
	
		my $d1 = $Silico::Vdw_Radii[$atom1->{ELEMENT_NUM}];
		
		if (!defined $d1) {
			atom_printout($atom1);
			silico_msg('d', "No vdW radius found for atom\n");
		}

		my $s = $atom1->{DONT_SCALE};	

		foreach my $atom2 (@{$mol2->{ATOMS}}) {
			
			next if $self_comparison && ($atom1->{NUM} >= $atom2->{NUM});
			
			# Dont scale atoms if scale factor < 1 and atoms are marked DONT_SCALE (i.e. are embedded molecule);
			my $s1;
			if (($s || $atom2->{DONT_SCALE}) && $scale < 0) {
				$s1 = 1;
			} else {
				$s1 = $scale;
			}
						
			my $cutoff = ($Silico::Vdw_Radii[$atom2->{ELEMENT_NUM}]+$d1)*$s1;

			my $distsq = general_distance_sq($atom1, $atom2, $mol2);
				
			next if ($distsq > $cutoff*$cutoff);
			++$clashcount;
			
			my $v;
			@$v = ($atom1, $atom2, $cutoff);
			push @$violations, $v;
			return ($clashcount, $violations) if defined $maxclash && ($clashcount >= $maxclash);
		}
	}
	return $clashcount, $violations;
}


sub clash_periodic_binned {

	#<
	#? Find clashes between two molecules, taking into account
	#  periodic boundary conditions
	#. A much faster routine than clash_periodic
	#; Requires: molecule1, molecule2, maximum number of clashes, vdW scale factor for added atoms
	#. mol2-> CELL, BINSIZE, NUMBINS and BINS must all be defined
	#; Returns: Number of clashes (up to maximum if this is set), pairs of violating atoms and distances
	#>

	my $mol1 = $_[0]; 		# Molecule 1
	my $mol2 = $_[1]; 		# Molecule 2 (Should be the larger of the two molecules);
	my $maxclash = $_[2]  || 1; 	# Maximum number of allowed clashes 
	my $scale = $_[3]  || 1; 	# vdW scaling factor for packed molecules
	my $dummy_rad = $_[4] || 0; 	# Radius of dummy atoms
	my $ignh = $_[5] || 0;		# Flag to ignore hydrogen atoms
		
	# Determine if we are comparing a molecule to itself
	my $self_comparison;
	if ($mol1 == $mol2) {
		$self_comparison = 1;
	}
	
	use POSIX qw(floor);
	
	my $cell = $mol2->{CELL} || croak("CELL not set");
	my $binsize = $mol2->{BINSIZE} || (molecule_printout($mol2) && croak("Binsize not set"));
	my $numbins = $mol2->{NUMBINS} || croak("Numbins not set");
	my $bins = $mol2->{BINS} || croak("Missing bins");
		
	my $cx = $cell->[0];
	my $cy = $cell->[1];
	my $cz = $cell->[2];
	my $cx2 = $cell->[0]/2;
	my $cy2 = $cell->[1]/2;
	my $cz2 = $cell->[2]/2;
	
	my $clashcount = 0;
	my $violations;

        foreach my $atom1 (@{$mol1->{ATOMS}}) {
	
		next if $ignh && $atom1->{ELEMENT} eq 'H';

		get_radius($atom1, $dummy_rad) if !defined $atom1->{RADIUS};
		
		my $d1 = $atom1->{RADIUS};
	
		if (!defined $d1) {
			atom_printout($atom1);
			silico_msg('d', "No vdW radius found for atom\n");
		}
		
		my $x1 = $atom1->{X};
                my $y1 = $atom1->{Y};
                my $z1 = $atom1->{Z};

		my $bin = atom_calc_bin($atom1, $mol2);
		
		# Make lists of bins within +/- one bin of the current bin
		my $blist;
		foreach my $i (0..2) {
			@{$blist->[$i]} = ($bin->[$i]-1..$bin->[$i]+1);
			foreach my $n (@{$blist->[$i]}) {
				$n = $n % $numbins->[$i]; # Modulo
			}
		}

		my $s = $atom1->{DONT_SCALE};
		
  		foreach my $bx (@{$blist->[0]}) {
		
                        foreach my $by (@{$blist->[1]}) {
			
                               	 foreach my $bz (@{$blist->[2]}) {
					
					foreach my $atom2 (@{$bins->[$bx][$by][$bz]}) {
					
						# Note
						# This could be sped up by not binning hydrogens
						# Needs modification to the mol2 binning routine
						next if $ignh & $atom2->{ELEMENT} eq 'H';
						
						next if $self_comparison && ($atom1->{NUM} >= $atom2->{NUM});
						
						# Dont scale atoms if scale factor < 1 and atoms are marked DONT_SCALE (i.e. are embedded molecule);
						my $s1;
						if (($s || $atom2->{DONT_SCALE}) && $scale < 0) {
							$s1 = 1;
						} else {
							$s1 = $scale;
						}
						
						get_radius($atom2, $dummy_rad) if !defined $atom2->{RADIUS};
						
						my $cutoff = ($atom2->{RADIUS}+$d1)*$s1;
						
						#print "$Silico::Vdw_Radii[$atom2->{ELEMENT_NUM}] $d1 $s1\n";
						
						my $xd = $x1 - $atom2->{X};
						my $xd2 = abs($xd - floor(($xd+$cx2)/$cx)*$cx);
						next if $xd2 > $cutoff;
						
						my $yd = $y1 - $atom2->{Y};
						my $yd2 = abs($yd - floor(($yd+$cy2)/$cy)*$cy);
						next if $yd2 > $cutoff;
					
						my $zd = $z1 - $atom2->{Z};
						my $zd2 = abs($zd - floor(($zd+$cz2)/$cz)*$cz);
						next if $zd2 > $cutoff;
						
						my $distsq = $xd2*$xd2+$yd2*$yd2+$zd2*$zd2;
						
						next if ($distsq > $cutoff*$cutoff);
						
						++$clashcount;
						
						#print "clashcount $clashcount maxclash $maxclash\n";
						
						my $v;
						@$v = ($atom1, $atom2, $cutoff);
						push @$violations, $v;
						return ($clashcount, $violations) if defined $maxclash && ($clashcount >= $maxclash);
                                        }
                                }
                        }
                }
	}
	
	return ($clashcount, $violations);
}


sub get_radius {

	my $atom = $_[0];
	my $dummy_rad = $_[1];

	$atom->{RADIUS} = $Silico::Vdw_Radii[$atom->{ELEMENT_NUM}];	
	if ($dummy_rad && $atom->{RADIUS} == 0) {
		$atom->{RADIUS} = $dummy_rad;
			
	}
}	
	


sub dont_scale_vdw {

	my $mol = $_[0];

	foreach my $atom (atoms($mol)) {
	
		$atom->{DONT_SCALE} = 1;
	}	
}

sub find_close_atoms_binned {

	#<
	#? Find atoms in mol2 which are close to mol1
	#; Requires: molecule, molecule, cutoff
	#; Returns: array of close atoms
	#>

	my $mol1 = $_[0]; # Molecule 1 or array of atoms
	my $mol2 = $_[1]; # Molecule 2. Should be the larger of the 2 molecules. Molecule_bin_atoms should have been run on this molecule
	my $cutoff = $_[2] || croak("Cutoff not set"); # Cutoff for 'close' atoms
	
	use POSIX qw(ceil);
	
	my $cutoffsq = $cutoff **2;
	my $cell = $mol2->{CELL} || croak("Cell not set");
	my $binsize = $mol2->{BINSIZE} || croak("binsize not set");
	my $numbins = $mol2->{NUMBINS} || croak("numbins not set");
	my $bins = $mol2->{BINS} || croak("bins missing");
	
	my $atoms1;
	if (ref $mol1 eq 'HASH') {
		@$atoms1 = $mol1->{ATOMS};
	} elsif (ref $mol1 eq 'ARRAY') {
		$atoms1 = $mol1;
	} else {
		croak("Incorrect thing passwd as mol1");
	}
	
	my $cx = $cell->[0];
	my $cy = $cell->[1];
	my $cz = $cell->[2];
	my $cx2 = $cell->[0]/2;
	my $cy2 = $cell->[1]/2;
	my $cz2 = $cell->[2]/2;
	
	# Bin maximum
	my $bxmax = $cell->[0]/$binsize->[0];
	my $bymax = $cell->[1]/$binsize->[1];
	my $bzmax = $cell->[2]/$binsize->[2];
	
	# Calculate how many bins to consider in each dimension
	my @bc;
	$bc[0] = ceil($cutoff/$binsize->[0]);
	$bc[1] = ceil($cutoff/$binsize->[1]);
	$bc[2] = ceil($cutoff/$binsize->[2]);
	
	my $hash;
	my $list;
        foreach my $atom1 (@$atoms1) {
	
		my $x1 = $atom1->{X};
                my $y1 = $atom1->{Y};
                my $z1 = $atom1->{Z};
		
		my $bin = atom_calc_bin($atom1, $mol2);
		
		my $blist;
		foreach my $i (0..2) {
			@{$blist->[$i]} = ($bin->[$i]-$bc[$i]..$bin->[$i]+$bc[$i]);
			foreach my $n (@{$blist->[$i]}) {
				$n = $n % $numbins->[$i]; # Modulo
			}
		}
		
  		foreach my $bx (@{$blist->[0]}) {
		
                        foreach my $by (@{$blist->[1]}) {
			
                               	 foreach my $bz (@{$blist->[2]}) {
					
					LOOP: foreach my $atom2 (@{$bins->[$bx][$by][$bz]}) {

            					my $xd = $x1 - $atom2->{X};
						my $xd2 = abs($xd - floor(($xd+$cx2)/$cx)*$cx);
						next if $xd2 > $cutoff;
						
						my $yd = $y1 - $atom2->{Y};
						my $yd2 = abs($yd - floor(($yd+$cy2)/$cy)*$cy);
						next if $yd2 > $cutoff;
					
						my $zd = $z1 - $atom2->{Z};
						my $zd2 = abs($zd - floor(($zd+$cz2)/$cz)*$cz);
						next if $zd2 > $cutoff;
						
						my $distsq = $xd2*$xd2+$yd2*$yd2+$zd2*$zd2;
						
						next if ($distsq > $cutoffsq);
						
						# Exclude atoms from mol1
						foreach my $atom3 (@$atoms1) {
							next LOOP if $atom3 == $atom2;
						}
						
						# Exclude duplicates
						push @$list, $atom2 if !$hash->{$atom2};
						++$hash->{$atom2};
                                        }
                                }
                        }
                }
	}
	return $list;
}




##########
# ANGLES #
##########

sub bndangle {
	
	#<
	#? Find the angle defined by three atoms
	#; Requires: atom1, atom2, atom3.
	#; Returns: angle in degrees.
	#>
	
	use POSIX;
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $atom3 = $_[2];
		
	#    set default value in case atoms are degenerate
	#    compute the value in degrees of the bond angle
	
	my $bndangle = 0;
	
	my $xab = $atom1->{X} - $atom2->{X};
	my $yab = $atom1->{Y} - $atom2->{Y};
	my $zab = $atom1->{Z} - $atom2->{Z};
	my $rab2 = $xab*$xab + $yab*$yab + $zab*$zab;
	my $xcb = $atom3->{X} - $atom2->{X};
	my $ycb = $atom3->{Y} - $atom2->{Y};
	my $zcb = $atom3->{Z} - $atom2->{Z};
	my $rcb2 = $xcb*$xcb + $ycb*$ycb + $zcb*$zcb;
	
	if ($rab2 != 0 && $rcb2 != 0) {
		my $dot = $xab*$xcb + $yab*$ycb + $zab*$zcb;
		my $cosine = $dot / sqrt($rab2*$rcb2);
		$cosine = 1 if ($cosine > 1);
		$cosine = -1 if ($cosine < -1);
		$bndangle = (180/$Silico::Pi) * POSIX::acos($cosine);
	}
	
	return $bndangle;
}

sub dihedral {
	
	#<
	#?   Find the value of the dihedral angle in the
	#    range from -180 to +180 degrees defined by four input atoms
	#;   Requires: atom1, atom2, atom3, atom4
	#;   Returns: dihedral in radians
	
	use POSIX;
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $atom3 = $_[2];
	my $atom4 = $_[3];
	
	# Compute the value in degrees of the dihedral angle
	my $xba = $atom2->{X} - $atom1->{X};
	my $yba = $atom2->{Y} - $atom1->{Y};
	my $zba = $atom2->{Z} - $atom1->{Z};
	my $xcb = $atom3->{X} - $atom2->{X};
	my $ycb = $atom3->{Y} - $atom2->{Y};
	my $zcb = $atom3->{Z} - $atom2->{Z};
	my $xdc = $atom4->{X} - $atom3->{X};
	my $ydc = $atom4->{Y} - $atom3->{Y};
	my $zdc = $atom4->{Z} - $atom3->{Z};
	my $xt = $yba*$zcb - $ycb*$zba;
	my $yt = $xcb*$zba - $xba*$zcb;
	my $zt = $xba*$ycb - $xcb*$yba;
	my $xu = $ycb*$zdc - $ydc*$zcb;
	my $yu = $xdc*$zcb - $xcb*$zdc;
	my $zu = $xcb*$ydc - $xdc*$ycb;
	my $rt2 = $xt*$xt + $yt*$yt + $zt*$zt;
	my $ru2 = $xu*$xu + $yu*$yu + $zu*$zu;
	my $rtru = sqrt($rt2 * $ru2);
	
	# Set default in case atoms are colinear or degenerate
	my $dihedral = 0;
	if ($rtru != 0) {
		my $cosine = ($xt*$xu + $yt*$yu + $zt*$zu) / $rtru;
		$cosine = 1 if ($cosine > 1);
		$cosine = -1 if ($cosine < -1);
		$dihedral = POSIX::acos($cosine);
		my $sign = $xba*$xu + $yba*$yu + $zba*$zu;
		$dihedral = -$dihedral if ($sign < 0);
	}
	
	return $dihedral;
}

sub two_bond_angle {
	
	#<
	#? Find the angle between two bonds (four atoms)
	#; Requires: atom1, atom2, atom3, atom4
	#; Returns: angle in radians.
	#>
	
	use POSIX;
	
	my $atom1 = $_[0];
	my $atom2 = $_[1];
	my $atom3 = $_[2];
	my $atom4 = $_[3];
	
	# set default value in case atoms are degenerate
	# compute the value in degrees of the bond angle
	
	my $angle = 0;
	
	my $xab = $atom1->{X} - $atom2->{X};
	my $yab = $atom1->{Y} - $atom2->{Y};
	my $zab = $atom1->{Z} - $atom2->{Z};
	my $rab2 = $xab*$xab + $yab*$yab + $zab*$zab;
	my $xcb = $atom3->{X} - $atom4->{X};
	my $ycb = $atom3->{Y} - $atom4->{Y};
	my $zcb = $atom3->{Z} - $atom4->{Z};
	my $rcb2 = $xcb*$xcb + $ycb*$ycb + $zcb*$zcb;
	
	my $dot;
	if ($rab2 != 0 && $rcb2 != 0) {
		$dot = $xab*$xcb + $yab*$ycb + $zab*$zcb;
		my$cosine = $dot / sqrt($rab2*$rcb2);
		$cosine = 1 if ($cosine > 1);
		$cosine = -1 if ($cosine < -1);
		$angle = POSIX::acos($cosine);
	}
	
	return $angle;
}

sub mol_rotate_euler{

	#<
	#? Rotates a molecule using Euler angles
	#. Uses the 'x convention'. See http://demonstrations.wolfram.com/EulerAngles/
	#; Requires: mol, phi, theta, psi (angles in radians)
	#; Returns: rotated molecule
	#>
	
	my ($mol, $phi, $theta, $psi) = @_;
	
	my $cphi   = cos($phi);
	my $ctheta = cos($theta);
	my $cpsi   = cos($psi);
	my $sphi   = sin($phi);
	my $stheta = sin($theta);
	my $spsi   = sin($psi);
	
	
	my $a = $cpsi*$cphi-$ctheta*$sphi*$spsi;
	my $b = $cpsi*$sphi+$ctheta*$cphi*$spsi;
	my $c = $spsi*$stheta;
	my $d = -$spsi*$cphi-$ctheta*$sphi*$cpsi;
	my $e = -$spsi*$sphi+$ctheta*$cphi*$cpsi;
	my $f = $cpsi*$stheta;
	my $g = $stheta*$sphi;
	my $h = -$stheta*$cphi;
	my $i = $ctheta;
	
	foreach my $atom (atoms($mol)) {
	
		my $x = $a * $atom->{X} + $b * $atom->{Y} + $c * $atom->{Z};
		my $y = $d * $atom->{X} + $e * $atom->{Y} + $f * $atom->{Z};
		my $z = $g * $atom->{X} + $h * $atom->{Y} + $i * $atom->{Z};
		
		$atom->{X} = $x;
		$atom->{Y} = $y;
		$atom->{Z} = $z;
	}
	
	return $mol;
}


###########
# VECTORS #
###########

sub bond_vector {
	
	#<
	#? Calculates a vector from two atoms.
	#; Requires: atom1 (origin), atom2 (vector end)
	#; Returns: vector as an array
	#>
	
	my $a1 = $_[0];
	my $a2 = $_[1];
	
	return ($a2->{X} - $a1->{X}, $a2->{Y} - $a1->{Y}, $a2->{Z} - $a1->{Z});
}

sub two_bond_normal {
	
	#<
	#? Return a normal vector to two bonds (four atoms)
	#; Requires: atom1, atom2, atom3, atom4
	#; Returns: vector (x, y, z)
	#>
	
	my ($a1, $a2, $a3, $a4) = @_;
	
	my @v1 = bond_vector($a1, $a2);
	my @v2 = bond_vector($a3, $a4);
	
	my @cross = crossprod(@v1, @v2);
	
	# if cross product is close to zero then change vector 2.
	# This will still produce a normal to vector 1 but in
	# an arbitrary direction
	
	if (magnitude(@cross) < 10e-8) {
		@cross = crossprod(@v1, $v2[0]+0.1, $v2[1]+0.1, $v2[2]-0.1);
	}
	
	# Normalise the vector
	@cross = unit_vector(@cross);
	
	return @cross;
}

sub calc_rotation_vector_onto_vector {

	#<
	#? Calculate rotation matrix and axis to rotate one vector onto another.
	#. See: http://immersivemath.com/forum/question/rotation-matrix-from-one-vector-to-another/
	#. Note: The implementation could be more efficient
	#; Requires: origin vector, target vector (3 coordinates)
	#; Returns: axis, rotation matrix
	#>

	my @veco = @_[0, 1, 2]; # Original vector
	my @vect = @_[3, 4, 5];  #Target vector
	
	use POSIX;
	
	# Centre of molecule
	my $axis;
	@$axis = my ($x, $y, $z) = unit_vector(crossprod(@veco, @vect)); # Rotation axis
	my $dp = dotprod(unit_vector(@veco), unit_vector(@vect));
	my $ang = POSIX::acos($dp); 
	my $cang = cos($ang); # Sin of angle
	my $sang = sin($ang); # Cos of angle
	my $d = 1-$cang;
	
	# Rotation matrix
	my $mat;
	@$mat =   ( $x**2*$d+$cang,      $x*$y*$d-$sang*$z,    $x*$z*$d+$sang*$y,
	            $x*$y*$d+$sang*$z,   $y**2*$d+$cang,       $y*$z*$d-$sang*$x,
		    $x*$z*$d-$sang*$y,   $y*$z*$d+$sang*$x,    $z**2*$d+$cang);
		   
	#print "c @c vec @vect axis @axis ang dp $dp $ang cang $cang sang $sang x $x y $y z $z d $d\n";
	#print "\nmat @mat[0, 1, 2]\n";
	#print   "mat @mat[3, 4, 5]\n";
	#print   "mat @mat[6, 7, 8]\n\n";
	
	return $mat, $axis;
}



sub angle_between_vectors {
	
	#<
	#? Angle between two vectors (radians)
	#; Requires: vector1, vector2
	#; Returns: angle
	#>
	
	use POSIX;
	
	my @vec1 = @_[0, 1, 2]; 
	my @vec2 = @_[3, 4, 5]; 
	
	my $dp = dotprod(unit_vector(@vec1), unit_vector(@vec2));
	
	return POSIX::acos($dp); 
}


###############
# TRANSLATION #
###############

sub molecule_translate {
	
	#<
	#? Translate a molecule by a vector (x,y,z).
	#; Requires: molecule, x, y, z.
	#>
	
	my $mol = $_[0];
	my @vec = @_[1..3];
	
	foreach my $atom (atoms($mol)) {
		
		$atom->{X} += $vec[0];
		$atom->{Y} += $vec[1];
		$atom->{Z} += $vec[2];
	}
}

sub molecule_trans_to_centre {
	
	#<
	#? Translate the centre of a molecule to 0,0,0 or to centre of cell
	#; Requires: molecule, specifier for origin or centre of cell ('cell'  or 'origin'), 
	#  specifier for whole molecule, largest substructure or hydrophobic centre ('all', 'largest' or 'hydrophobic')
	#; Returns: vector of translation
	#>
	
	my $mol = $_[0];
	my $mode  = $_[1] || 'cell';
	my $mode2 = $_[2] || 'all';
	my $cell = $_[3] || $mol->{CELL};
	
	my @dim; 
	my ($dx,$dy,$dz);
	
	if ($mode2 eq 'largest') {
		
		molecule_check_and_fix_connectivity($mol);
			
		my $newmols = molecule_split($mol);
		
		# Find biggest fragment
		my $l;
		foreach my $newmol (@$newmols) {
			$l = $newmol if (!defined $l->{ATOMS}[0] || $#{$newmol->{ATOMS}} > $#{$l->{ATOMS}});
		}
		
		@dim = molecule_centre($l);
		
	} elsif ($mode2 eq 'hydrophobic') {
	
		silico_msg('w', "Hydrophobic option not implemented. Using whole molecule\n");
		@dim = molecule_centre($mol);
	
	} elsif ($mode2 eq 'all') {
	
		@dim = molecule_centre($mol);
		
	} else {
	
	 	silico_msg('d', "Option 'all', 'largest' or 'hydrophobic' not set\n");
	}
	
	if ($mode eq 'origin') {
		
		$dx = -$dim[6];
		$dy = -$dim[7];
		$dz = -$dim[8];
	
	} elsif ($mode eq 'cell') {
	
	 	if (!defined $cell->[0]) {
			@$cell = (0,0,0);
			silico_msg('c', "Cell not defined. Translating to origin\n");
		}

		$dx = $cell->[0]/2 - $dim[6];
		$dy = $cell->[1]/2 - $dim[7];
		$dz = $cell->[2]/2 - $dim[8];
		
	} else {
		 silico_msg('d', "Option 'origin' or 'cell' not set\n");
	}
	
	#print "Translation vector ($dim[6], $dim[7], $dim[8])\n";
		
	foreach my $atom (atoms($mol)) {
		
		$atom->{X} += $dx;
		$atom->{Y} += $dy;
		$atom->{Z} += $dz;
	}
	
	return ($dx, $dy, $dz);
}



sub mol_move_into_cell {

	#<
	#? Translate a molecule into a unit cell using the centroid of the molecule as a reference
	#. Unit cell has coordinates starting at 0,0,0
	#; Requires: Molecule, cell, flag to move based on position the 'middle' atom of the molecule rather than centroid (much faster).
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $cell = $_[1] || $mol->{CELL} || croak("Missing cell");
	my $move_atom = $_[2];
	
	
	my @centre;
	if ($move_atom) {
	
		my $n = int($#{$mol->{ATOMS}}/2);
		my $a = $mol->{ATOMS}[$n];
		@centre = ($a->{X}, $a->{Y}, $a->{Z});
	} else {
	
		@centre = centroid($mol);
	}

	my $i = -1;
	my $mv;
	foreach my $v ('X', 'Y', 'Z') {

		++$i;
		my $g = 0;
		my $c = $cell->[$i]; # Cell dimension
	
		while ($centre[$i] < 0) {
			$centre[$i] += $c;
			++$g;
		}
		while ($centre[$i] >= $c) {
			$centre[$i] -= $c;
			--$g;
		}
		
		$mv->{$v} = $g*$c;
	}
	
	foreach my $atom (@{$mol->{ATOMS}}) {
		foreach my $v ('X', 'Y', 'Z') {
			$atom->{$v} += $mv->{$v};
		}
	}
}


sub mol_move_atoms_into_cell {

	#<
	#? Translate all atoms into the unit cell (ie into the range 0 -> cell_dim).
	
	#; Requires: Molecule, cell
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $cell = $_[1] || $mol->{CELL};
	
	foreach my $atom (atoms($mol)) {

		my $i = 0;
		foreach (qw(X Y Z)) {
		
			$atom->{$_} -= $cell->[$i] while ($atom->{$_} > $cell->[$i]);
			$atom->{$_} += $cell->[$i] while ($atom->{$_} <= 0);
			
			++$i;
		}
	}
}

sub aa_chirality {

	#<
	#? Find chirality of amino acid
	#; Requires: atoms N, CA, C, CB
	#; Returns: D, L or U if chiral volume is close to 0
	#>
	
	foreach (@_) {
		croak() if !defined $_;
	}
	
	my $chiral = 'U';
	
	my $cv = chiral_volume($_[0], $_[1], $_[2], $_[3]);
	
	$chiral = 'L' if ($cv > 0.1);
	$chiral = 'D' if ($cv < 0.1);
	
	#print "cv $cv c $chiral\n";
	
	return $chiral;	
}

sub chiral_volume {

	#<
	#? Calculate the chiral volume
	#. For atoms AB(C)D - B is the central atom
	#. v1 = BA
	#. v2 = BC
	#. v3 = BD
	#. vol = v1 . (v2 x v3)
	#; Requires: four atoms
	#; Returns: scalar
	#>
	
	foreach my $a (@_) {
		if (!defined $a || !defined $a->{X}) {
			atom_printout($a);
			print "----\n";
			croak();
		}
	}
	
	my @v1 = atom_vector($_[1], $_[0]);
	my @v2 = atom_vector($_[1], $_[2]);
	my @v3 = atom_vector($_[1], $_[3]);
	
	my @c = crossprod(@v2, @v3);
	return dotprod(@v1, @c);
}

############
# ROTATION #
############

sub molecule_rotate {
	
	#<
	#? Rotate a molecule using a 3x3 matrix .
	#; Requires: molecule, m1..m9.
	#>
	
	my $mol = $_[0];
	my @mat = @_[1..9];	# 3x3 rotation matrix
	
	silico_msg('d', "Error in matrix") if !defined $mat[8];
	
	foreach my $atom (atoms($mol)) {
		
		my $nx = $atom->{X}*$mat[0] + $atom->{Y}*$mat[1] + $atom->{Z}*$mat[2];
		my $ny = $atom->{X}*$mat[3] + $atom->{Y}*$mat[4] + $atom->{Z}*$mat[5];
		my $nz = $atom->{X}*$mat[6] + $atom->{Y}*$mat[7] + $atom->{Z}*$mat[8];

		$atom->{X} = $nx;
		$atom->{Y} = $ny;
		$atom->{Z} = $nz;
	}
}


sub molecule_rot_angle_x {
	
	#<
	#? Rotate a molecule around X axis by angle (degrees).
	#; Requires: molecule, angle (deg);
	#>
	
	my $mol = $_[0];
	my $angle_deg = $_[1];
	
	my $theta = deg_to_rad($angle_deg);
	my $ct = cos ($theta);
	my $st = sin ($theta);
	
	foreach my $atom (atoms($mol)) {
		
		#atom->{X} doesn't change
		
		my $oy = $atom->{Y}; # Original value of y
		my $oz = $atom->{Z}; # Original value of z
		
		$atom->{Y} = $oy*$ct + $oz*$st;
		$atom->{Z} = - $oy*$st + $oz*$ct;
	}
}

sub molecule_rot_angle_y {
	
	#<
	#? Rotate a molecule around Y axis by angle (degrees).
	#; Requires: molecule, angle (deg);
	#>
	
	my $mol = $_[0];
	my $angle_deg = $_[1];
	
	my $theta = deg_to_rad($angle_deg);
	my $ct = cos ($theta);
	my $st = sin ($theta);
	
	foreach my $atom (atoms($mol)) {
		
		#atom->{Y} doesn't change
		
		my $ox = $atom->{X}; # Original value of x
		my $oz = $atom->{Z}; # Original value of z
		
		$atom->{Z} = $oz*$ct + $ox*$st;
		$atom->{X} = - $oz*$st + $ox*$ct;
	}
}

sub molecule_rot_angle_z {
	
	#<
	#? Rotate a molecule around Z axis by angle (degrees).
	#; Requires: molecule, angle (deg);
	#>
	
	my $mol = $_[0];
	my $angle_deg = $_[1];
	
	my $theta = deg_to_rad($angle_deg);
	my $ct = cos ($theta);
	my $st = sin ($theta);
	
	foreach my $atom (atoms($mol)) {
		
		#atom->{Z} doesn't change
		
		my $ox = $atom->{X}; # Original value of x
		my $oy = $atom->{Y}; # Original value of y
		
		$atom->{X} = $ox*$ct + $oy*$st;
		$atom->{Y} = - $ox*$st + $oy*$ct;
	}
}

sub molecule_rotate_axis {
	
	#<
	#? Rotate atoms about an arbitrary axis using quaternions
	#; Require: molecule, angle (in radians),
	#  axis (x1, y1, z1), origin (x2, y2, z2) default 0,0,0
	#>
	
	my $mol = $_[0];
	my $angle = $_[1];
	my ($x1, $y1, $z1) = @_[2..4];
	my ($x2, $y2, $z2) = @_[5..7];
	
	my @mat;
	
	$x2 ||= 0;
	$y2 ||= 0;
	$z2 ||= 0;
	
	my $flag = 1 if $x2 || $y2 || $z2;
	
	# Compute the rotation matrix - a quaternion
	# This includes a normalisation of the vector [x,y,z] to a unit vector
	@mat = quaternion($angle, $x1, $y1, $z1);
	
	# Move molecule so that the centre of rotation is the origin
	molecule_translate($mol, -$x2, -$y2, -$z2) if $flag;
	
	# Rotate molecule
	molecule_rotate($mol, @mat);
	
	# Put molecule back
	molecule_translate($mol, $x2, $y2, $z2) if $flag;
}

################
# Random moves #
################


sub molecule_random_move_cartesian {

        #<
        #? Move and rotate a molecule (around a random axis) by a random amount
        #. Will use a periodic cell if supplied.
	#  Molecule is translated back into the cell 
        #; Require: molecule, translation vector (default cell size, or 10, 10, 10), 
	#  rotation vector (default random vector with up to 360 deg rotation), cell (optional)
        #>

        my $mol = $_[0];
        my $mvec = $_[1]; # Array: max X displacement, max Y displacement, max Z
        my $rvec = $_[2]; # Array: max X vector for rotation axis, max Y vector for rotation axis, max Z vector, max rotation angle (radians)
        my $cell = $_[3];

        if (!$mvec) {
                if ($cell) {
                        @$mvec = @$cell;
                } else {
			@$mvec = (10, 10, 10);
                }
        }
        if (!$rvec) {
                @$rvec = (1, 1, 1, 2*$Silico::Pi);
        }

        # Random rotation
        my $x1 = rand($rvec->[0]);
        my $y1 = rand($rvec->[1]);
        my $z1 = rand($rvec->[2]);
        my $angle = rand($rvec->[3]);

        # Random translation
        my $x2 = rand($mvec->[0]*2) - $mvec->[0];
        my $y2 = rand($mvec->[1]*2) - $mvec->[1];
        my $z2 = rand($mvec->[2]*2) - $mvec->[2];

        # Compute the rotation matrix - a quaternion
        # This includes a normalisation of the vector [x,y,z] to a unit vector
        my @mat = quaternion($angle, $x1, $y1, $z1);

        foreach my $atom (atoms($mol)) {

                my $nx = $atom->{X}*$mat[0] + $atom->{Y}*$mat[1] + $atom->{Z}*$mat[2];
                my $ny = $atom->{X}*$mat[3] + $atom->{Y}*$mat[4] + $atom->{Z}*$mat[5];
                my $nz = $atom->{X}*$mat[6] + $atom->{Y}*$mat[7] + $atom->{Z}*$mat[8];

                $atom->{X} = $nx + $x2;
                $atom->{Y} = $ny + $y2;
                $atom->{Z} = $nz + $z2;
        }

        mol_move_into_cell($mol, $cell, 1) if $cell;
}

sub select_rotate_bond {
	
	#<
	#? Select a random rotatable bond from a list of atoms
	#; Requires: molecule, atom list (uses heavy atoms if this is undefined)
	#; Returns: atom1, atom2; undef if failed to find a rotatable bond
	#>
	
	my $mol = $_[0];
	my $list1 = $_[1];
	
	# Use heavy atoms if list is not defined
	if (!defined $list1) {
	
		foreach my $atom (atoms($mol)) {
		
			next if $atom->{ELEMENT} eq 'H';
			push @$list1, $atom;
		}
	}
	
	# Find rings of up to size 20
	molecule_find_rings($mol, 20);
	
	my $count = 0;
	my ($atom1, $atom2);
	my @rlist;
	ATOMONE: while (1) {
		
		++$count;
		
		return undef if $count > 100;
		
		# Choose a random atom from list1 that has more than one connected atom
		while (1) {
			
			my $r = int rand ($#{$list1}+1);
			$atom1 = $list1->[$r];
			last if ($#{$atom1->{CONNECT}} > 0);
		}
		
		# Find all rotatable bonds connected to this atom
		foreach my $connum (@{$atom1->{CONNECT}}) {
			
			my $con = $mol->{ATOMS}[$connum];
			
			# Skip bonds that connect to a single atom
			next if ($#{$con->{CONNECT}} < 1);
			
			my $bo = find_bondorder($atom1->{NUM}-1, $connum, $mol);
			
			# Skip if the bond between these two atoms is double, triple, or aromatic
			next if ($bo == 2 || $bo == 3 || $bo == 4);
			
			# Skip if these two atoms share one or more rings
			next if shared_rings($mol, $atom1, $con) > 0;
					
			# Skip if one of these atoms is a secondary amide nitrogen and the other is an amide carbon
			next if ($atom1->{FG}{S_AMIDE_C} && $con->{FG}{S_AMIDE_N});
			next if ($atom1->{FG}{S_AMIDE_N} && $con->{FG}{S_AMIDE_C});
			
			push @rlist, $con;
		}
		
		# Repeat if there were no rotatable bonds
		redo ATOMONE if ($#rlist < 0);
		
		# Select random atom
		my $r = int (rand ($#rlist+1));
		
		$atom2 = $rlist[$r];
		last;
	}
	
	return $atom1, $atom2;
}

sub rotate_torsion_angle {
	
	#<
	#? Rotate a torsion angle
	#. Will do nothing if bond is in a ring
	#. The static atoms are overridden if an ANCHOR atom is set.
	#; Requires: molecule, 2 atoms (first one static), angle, anchor (optional), flag to rotate
	#  atoms already marked with PFLAG
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $atom2 = $_[1]; # Static
	my $atom3 = $_[2]; # Moves
	my $rot_angle = $_[3]; # Angle in radians
	my $anchor = $_[4]; # Anchor atom
	my $dont_partition = $_[5]; # Use existing PFLAGS to rotate atoms
	
	return if shared_rings($mol, $atom2, $atom3) > 0;
	
	my @axis = bond_vector($atom3, $atom2);
	
	# Mark atoms to be rotated with atom->{PFLAG} unless the dont_partition flag is set
	if (!$dont_partition) {
		
		my $val = mol_partition_on_bond($mol, $atom2, $atom3);
		return if $val < 0;
	}
	
	# Check to see if there is an anchor atom in the
	# atoms marked to be rotated.  If so we must
	# change the direction of rotation.
	my $swap;
	if ($anchor) {
		foreach my $atom (atoms($mol)) {
			
			if ($atom == $anchor && $atom->{PFLAG}) {
				$swap = 1;
				$rot_angle = -$rot_angle;
				last;
			}
		}
	}
	
	my @mat = quaternion($rot_angle, @axis);
	
	my $xo = $atom3->{X};
	my $yo = $atom3->{Y};
	my $zo = $atom3->{Z};
	
	foreach my $atom (atoms($mol)) {
		
		if ($swap) {
			next if $atom->{PFLAG};
		} else {
			next if !$atom->{PFLAG};
		}
		
		my $x = $atom->{X} - $xo;
		my $y = $atom->{Y} - $yo;
		my $z = $atom->{Z} - $zo;
		
		my $nx = $x*$mat[0] + $y*$mat[1] + $z*$mat[2];
		my $ny = $x*$mat[3] + $y*$mat[4] + $z*$mat[5];
		my $nz = $x*$mat[6] + $y*$mat[7] + $z*$mat[8];
		
		$atom->{X} = $nx + $xo;
		$atom->{Y} = $ny + $yo;
		$atom->{Z} = $nz + $zo;
	}
}

sub rotate_torsion_angle2 {
	
	#<
	#? Rotate a torsion angle
	#. Will do nothing if bond is in a ring
	#; Requires: molecule, vector atom1, vector_atom2, list of atoms to be rotated, angle in radians
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $atom2 = $_[1]; # Static
	my $atom3 = $_[2]; # Moves
	my $list = $_[3];
	my $rot_angle = $_[4]; # Angle in radians
	
	return if shared_rings($mol, $atom2, $atom3) > 0;
	
	my @axis = bond_vector($atom3, $atom2);
	my @mat = quaternion($rot_angle, @axis);
	
	my $xo = $atom3->{X};
	my $yo = $atom3->{Y};
	my $zo = $atom3->{Z};
	
	foreach my $atom (@$list) {
		
		my $x = $atom->{X} - $xo;
		my $y = $atom->{Y} - $yo;
		my $z = $atom->{Z} - $zo;
		
		my $nx = $x*$mat[0] + $y*$mat[1] + $z*$mat[2];
		my $ny = $x*$mat[3] + $y*$mat[4] + $z*$mat[5];
		my $nz = $x*$mat[6] + $y*$mat[7] + $z*$mat[8];
		
		$atom->{X} = $nx + $xo;
		$atom->{Y} = $ny + $yo;
		$atom->{Z} = $nz + $zo;
	}
}

sub set_torsion_angle {
	
	#<
	#? Set a torsion angle to a particular angle
	#. Will rotate whole molecule if bond is in a ring
	#. The static atoms are overridden if an ANCHOR atom is set.
	#; Requires: Molecule, 4 atoms (first two static), angle
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $atom1 = $_[1]; # Static
	my $atom2 = $_[2]; # Static
	my $atom3 = $_[3]; # Moves
	my $atom4 = $_[4]; # Moves
	my $angle = $_[5]; # Angle in radians
	my $anchor = $_[6];
	
	my $dihedral = dihedral($atom1, $atom2, $atom3, $atom4);
	
	# Return if the angle is already set.
	return if abs($dihedral - $angle) < 10e-4;
	
	my $rot_angle = $dihedral-$angle;
	
	rotate_torsion_angle($mol, $atom2, $atom3, $rot_angle, $anchor);
}

sub make_amides_trans {
	
	#<
	#? Sets the C-N-C(O)-C bond in all secondary amides to 180 degrees
	#; Requires: molecule, list of atoms to check (uses entire molecule if this is not defined)
	#; Returns: nothing
	#>
	
	my $mol = $_[0];
	my $atom_list = $_[1];
	
	use vars qw($Pi);
	
	if (!defined $atom_list) {
		@$atom_list = @{$mol->{ATOMS}};
	}
	
	label_amides($mol, $atom_list);
	
	foreach my $atom (@$atom_list) {
	
		my ($atom1, $atom2, $atom3, $atom4);
		
		next if !$atom->{FG}{S_AMIDE_C};
		
		$atom2 = $atom;
		
		foreach my $connum1 (@{$atom->{CONNECT}}) {
			
			my $con1 = $mol->{ATOMS}[$connum1];
			
			# Find secondary and tertiary amides
			if ($con1->{FG}{S_AMIDE_N}) {
				
				$atom3 = $con1;
				
				foreach my $connum2 (@{$con1->{CONNECT}}) {
					
					my $con2 = $mol->{ATOMS}[$connum2];
					
					next if $con2 == $atom2;
					next if $con2->{ELEMENT} ne 'C';
					
					$atom4 = $con2;
					last;
				}
			}
			
			if (!$con1->{FG}{S_AMIDE_O} && !$con1->{FG}{S_AMIDE_N}) {
				$atom1 = $con1;
			}
		}
		
		if (!defined $atom1 || !defined $atom2 || !defined $atom3 || !defined $atom4) {
			
			silico_msg('e', "AMIDE_C defined but AMIDE_N and/or AMIDE_O missing!\n",
					"Amide C: Atom number $atom->{NUM}\n");
			return;
		}
		
		silico_msg('c', "Making amide trans. Atoms: $atom1->{NUM} $atom2->{NUM} $atom3->{NUM} $atom4->{NUM}\n");
		set_torsion_angle($mol, $atom1, $atom2, $atom3, $atom4, $Pi);
	}
}

sub mol_partition_on_bond {
	
	#<
	#? Label atoms on one side of a bond and not the other
	#. Is expected to fail if the bond is within a ring
	#; Requires: molecule, bond atom 1 (static), bond atom 2 (to move)
	#; Returns: Time taken
	#>
	
	my $mol = $_[0];
	my $atom1 = $_[1]; # Static
	my $atom2 = $_[2]; # Moves
	
	my @list;
	my $starttime = (times)[0];
	
	foreach my $atom (atoms($mol)) {
		$atom->{PFLAG} = 0;
	}
	
	$atom1->{PFLAG} = 2;
	$atom2->{PFLAG} = 1;
	
	push @list, $atom2;
	
	while (1) {
		
		last if $#list < 0;
		
		my $atom = pop (@list);
		
		foreach my $connum (@{$atom->{CONNECT}}) {
			
			my $con = $mol->{ATOMS}[$connum];
			
			next if $con->{PFLAG};
			
			$con->{PFLAG} = 1;
			push @list, $con;
		}
	}
	
	my $timetaken = (times)[0] - $starttime;
	return $timetaken;
}

sub mol_partition_on_bond2 {
	
	#<
	#? Label atoms on one side of a bond (best non-recursive version)
	#. Will fail if start atom is in a ring
	#; Requires: molecule, bond atom 1 (static), bond atom 2 (to move),
	#  anchor atom
	#; Returns: List of stationary atoms, list of mobile atoms
	#>
	
	my $mol = $_[0];
	my $atom1 = $_[1]; #
	my $atom2 = $_[2]; #
	my $anchor = $_[3];
	
	#Delete any previous flags
	foreach my $atom (atoms($mol)) {
		$atom->{PFLAG} = 0;
	}
	
	# Establish whether atom1 and atom2 are in the same ring. We hope not.
	if (shared_rings($mol, $atom1, $atom2) > 0) {
		silico_msg('e', "Bond to be partitioned on is in a ring!\n",
				"Atom 1: $atom1->{NUM} (\"$atom1->{NAME}\")   Atom 2: $atom2->{NUM} (\"$atom2->{NAME}\")\n",
				"Nothing will be done.\n");
		return (undef, undef);
	}
	
	# Delete bond in question
	my $bo = find_bondorder2($atom1, $atom2, $mol);
	silico_msg('d', "Atoms are not bonded! Atoms: $atom1->{NUM}, $atom2->{NUM}\n") if $bo == -1;
	
	bond_delete_atom($mol, $atom1, $atom2);
	
	$anchor->{PFLAG} = 1;
	
	my @list;
	push @list, $anchor;
	
	while (1) {
		
		last if $#list < 0;
		
		my $atom = pop(@list);
		
		foreach my $connum (@{$atom->{CONNECT}}) {
			
			my $con = $mol->{ATOMS}[$connum];
			
			next if $con->{PFLAG};
			
			$con->{PFLAG} = 1;
			push @list, $con;
		}
	}
	
	my $list1;
	my $list2;
	foreach my $atom (atoms($mol)) {
		if ($atom->{PFLAG}) {
			push @$list1, $atom; # Anchored set
		} else {
			push @$list2, $atom; # Non-anchored set
		}
	}
	
	# Recreate bond
	bond_create_atom($mol, $atom1, $atom2, $bo);
	
	# Check to make sure both lists are actually occupied
	if ($#{$list1} < 0 || $#{$list2} < 0) {
		silico_msg('e', "Molecule was not effectively partitioned on this bond!\n",
				"Atom 1: $atom1->{NUM} (\"$atom1->{NAME}\")   Atom 2: $atom2->{NUM} (\"$atom2->{NAME}\")\n",
				"Nothing has been done.\n");
		return (undef, undef);
	}
	
	return $list1, $list2;
}

sub record_best_coords {
	
	#<
	#? Records the current coordinates of each atom as being the best ones.
	#; Requires: Molecule
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	
	# Set each atom's "best coordinates" to the current ones
	foreach my $atom (atoms($mol)) {
		$atom->{BEST_X} = $atom->{X};
		$atom->{BEST_Y} = $atom->{Y};
		$atom->{BEST_Z} = $atom->{Z};
	}
}

sub revert_to_best_coords {
	
	#<
	#? Sets atoms back to its best known position.
	#; Requires: Molecule, atom list (optional)
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $list = $_[1];
	
	#Set list to all atoms in molecule if list is not defined
	@$list = @{$mol->{ATOMS}} if !defined $list;
	
	# Set each atom's current coordinates to the recorded "best coordinates"
	foreach my $atom (@$list) {
		$atom->{X} = $atom->{BEST_X};
		$atom->{Y} = $atom->{BEST_Y};
		$atom->{Z} = $atom->{BEST_Z};
	}
}

###########
# SCALING #
###########

sub mol_rescale_bonds {
	
	#<
	#? Rescale all bond lengths in a molecule, based on either a scaling
	#  factor, or an optimum average C-C bond length.
	#; Requires: Molecule, one of either scaling factor or target length
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $sfactor = $_[1];
	my $length = $_[2];
	
	my @distances;
	my $factor;
	
	# Check to make sure one, and only one, of scaling factor or target length
	# have been passed
	unless (defined $sfactor xor defined $length) {
		silico_msg('d', "Wrong arguments supplied to subroutine!\n",
				"Either a scaling factor (argument 2) or a target length (argument 3)\n",
				"must be supplied.\n");
	}
		
	if (defined $length && !defined $sfactor) {
		
		# Calculate mean C-C bondlength
		foreach my $atom (atoms($mol)) {
			
			next if $atom->{ELEMENT} ne 'C';
			
			foreach my $connum (@{$atom->{CONNECT}}) {
				
				my $con = $mol->{ATOMS}[$connum];
				next if $con->{ELEMENT} ne 'C';
				
				push @distances, distance($con, $atom);
			}
		}
		
		if (!defined $distances[0]) {
			silico_msg('w', "Molecule \"$mol->{NAME}\" contains no carbon-carbon bonds!\n",
					"Skipping.\n");
			return;
		}
		
		# Calculate the average C-C bond length
		my $mean = sprintf "%1.3f", mean(@distances);
		
		# Print out said C-C bond length
		silico_msg('c', "Molecule: \"$mol->{NAME}\"   Mean C-C bondlength: $mean A\n");
		
		if ($mean > ($length - 0.05) && $mean < ($length + 0.05)) {
			
			silico_msg('c', "Average C-C bond length is within 0.05 Angstroms of desired length.\n",
					"Molecule will not be rescaled.\n");
			
			return;
		}
		
		$factor = $length/$mean;
			
	} else {
		$factor = $sfactor;
	}
		
	silico_msg('c', "Rescaling molecule by a factor of $factor.\n");
	
	foreach my $atom (atoms($mol)) {
		$atom->{X} *= $factor;
		$atom->{Y} *= $factor;
		$atom->{Z} *= $factor;
	}
}

#################
# SYSTEM CENTRE #
#################

sub centroid {
	
	#<
	#? Calculate the average atom position of a molecule or subset.
	#; Requires: Molecule, optional atom list
	#; Returns: array
	#>
	
	my $mol = $_[0];
	my $atomlist = $_[1];
	
	if (!defined $atomlist) {
		@$atomlist = @{$mol->{ATOMS}};
	} else {
		my $r = ref($atomlist);
		silico_msg('d', "Subroutine requires pointer to array, not '$r'")if $r ne "ARRAY";
	}

	my $count = 0;
	my $xsum = 0;
	my $ysum = 0;
	my $zsum = 0;
	
	foreach my $atom (@$atomlist) {
		
		++$count;
		
		$xsum += $atom->{X};
		$ysum += $atom->{Y};
		$zsum += $atom->{Z};
	}
	
	return undef if ($count == 0);
	return ($xsum/$count, $ysum/$count, $zsum/$count);
}

sub centroid_periodic {
	
	#<
	#? Calculate the average atom position of a molecule or subset considering periodic boundary conditions
	#. First atom of set is used as the reference atom
	#; Requires: Molecule, optional atom list, molecule cell
	#; Returns: array
	#>
	
	my $mol = $_[0];
	my $atomlist = $_[1];
	my $cell = $_[2] || $mol->{CELL} || silico_msg('d', "No cell inormation provided");
	
	if (!defined $atomlist) {
		@$atomlist = @{$mol->{ATOMS}};
	} else {
		my $r = ref($atomlist);
		silico_msg('d', "Subroutine requires pointer to array, not '$r'")if $r ne "ARRAY";
	}

	my $count = 0;
	my $xsum = 0;
	my $ysum = 0;
	my $zsum = 0;
	
	my $ref = $atomlist->[0];
	my @d = (0, 0, 0);
	
	foreach my $atom (@$atomlist) {
		
		++$count;
		
		my @r = distance_periodic_signed_xyz($atom, $ref, $cell);
		#my @r = ($atom->{X}-$ref->{X}, $atom->{Y}-$ref->{Y}, $atom->{Z}-$ref->{Z});
		
		foreach my $i (0..2) {
			$d[$i] += $r[$i];
		}
	}
	
	silico_msg('d', "No atoms passed to routine") if $count == 0;
	
	foreach my $i (0..2) {
		$d[$i] /= $count;
	}
	
	$d[0] += $ref->{X};
	$d[1] += $ref->{Y};
	$d[2] += $ref->{Z};
	
	return @d;
}

sub centroid_periodic_bv {
	
	#<
	#? Calculate the average atom position of a molecule or subset considering periodic boundary conditions
	#. First atom of set is used as the reference atom
	#; Requires: Molecule, optional atom list, box vector
	#; Returns: array
	#>
	
	my $mol = $_[0];
	my $atomlist = $_[1];
	my $bv = $_[2] || $mol->{GROMACS_BOX_VECTORS} ||  silico_msg('d', "No box vector provided");
	
	if (!defined $atomlist) {
		@$atomlist = @{$mol->{ATOMS}};
	} else {
		my $r = ref($atomlist);
		silico_msg('d', "Subroutine requires pointer to array, not '$r'")if $r ne "ARRAY";
	}

	my $count = 0;
	my $xsum = 0;
	my $ysum = 0;
	my $zsum = 0;
	
	my $ref = $atomlist->[0];
	my @d = (0, 0, 0);
	
	foreach my $atom (@$atomlist) {
		
		++$count;
		
		my @r = distance_periodic_xyz_bv($atom, $ref, $bv);

		foreach my $i (0..2) {
			$d[$i] += $r[$i];
		}
	}
	
	silico_msg('d', "No atoms passed to routine") if $count == 0;
	
	foreach my $i (0..2) {
		$d[$i] /= $count;
	}
	
	$d[0] += $ref->{X};
	$d[1] += $ref->{Y};
	$d[2] += $ref->{Z};

	print "d @d\n";
	
	return @d;
}

sub molecule_dimensions {
	
	#<
	#? Calculate the x, y, z dimensions of a molecule
	#; Requires: Molecule, margin to add (optional)
	#; Returns: Dimensions
	#>
	
	my $mol = $_[0];
	my $margin = $_[1] || 0;
	
	my $count = 0;
	my ($maxx, $maxy, $maxz);
	my ($minx, $miny, $minz);

	foreach my $atom (atoms($mol)) {
		
		$maxx = $atom->{X} if (!defined $maxx || $atom->{X} > $maxx);
		$maxy = $atom->{Y} if (!defined $maxy || $atom->{Y} > $maxy);
		$maxz = $atom->{Z} if (!defined $maxz || $atom->{Z} > $maxz);

		$minx = $atom->{X} if (!defined $minx || $atom->{X} < $minx);
		$miny = $atom->{Y} if (!defined $miny || $atom->{Y} < $miny);
		$minz = $atom->{Z} if (!defined $minz || $atom->{Z} < $minz);
		
		++$count;
	}
	
	if ($count == 0) {

		silico_msg('w', "Calculating centre of molecule with $count atoms\n");
		return (0,0,0);
	}
	
	return ($maxx - $minx+$margin, $maxy-$miny+$margin, $maxz-$minz+$margin);
}

sub molecule_centre {
	
	#<
	#? Calculate the x, y, z extents of a molecule and the centre of the box.
	#; Requires: Molecule
	#; Returns: Array ({maxx, maxy, maxz, minx, miny, minz, centx, centy, centz})
	#>
	
	my $mol = $_[0];
	
	my $count = 0;
	my ($maxx, $maxy, $maxz);
	my ($minx, $miny, $minz);
	my ($centx, $centy, $centz);
	
	foreach my $atom (atoms($mol)) {
		
		$maxx = $atom->{X} if (!defined $maxx || $atom->{X} > $maxx);
		$maxy = $atom->{Y} if (!defined $maxy || $atom->{Y} > $maxy);
		$maxz = $atom->{Z} if (!defined $maxz || $atom->{Z} > $maxz);

		$minx = $atom->{X} if (!defined $minx || $atom->{X} < $minx);
		$miny = $atom->{Y} if (!defined $miny || $atom->{Y} < $miny);
		$minz = $atom->{Z} if (!defined $minz || $atom->{Z} < $minz);
		
		++$count;
	}
	
	if ($count == 0) {

		silico_msg('w', "Calculating centre of molecule with $count atoms\n");
		return (0,0,0,0,0,0,0,0,0);
	}
	
	$centx = ($maxx + $minx) / 2;
	$centy = ($maxy + $miny) / 2;
	$centz = ($maxz + $minz) / 2;
	
	return ($maxx, $maxy, $maxz, $minx, $miny, $minz, $centx, $centy, $centz);
}

sub molecule_centre_of_mass {
	
	#<
	#? Calculate the centre of mass of a molecule or atom subset
	#; Requires: Molecule, (optional) atom list
	#; Returns: Centre of mass as an array ({x, y, z})
	#>
	
	my $mol = $_[0];
	my $atomlist = $_[1];
	
	my $c;
	my $i = 0;
	my @masswarn;
	my $total_mass = 0;
	my $xsum = 0;
	my $ysum = 0;
	my $zsum = 0;
	my (@xwarn, @ywarn, @zwarn);
	
	atomic_masses();
	
	# Die with an error if the @Atomic_Masses array is undefined
	# This indicates a serious problem in Silico
	if (!defined $Silico::Atomic_Masses[0]) {
		silico_msg('d', "Array \@Silico::Atomic_Masses is not properly defined!\n",
				"This is a serious error. You may wish to check Silico setup.\n");
	}
	
	# Set atomlist to be all atoms in molecule if atomlist is not defined
	@$atomlist = @{$mol->{ATOMS}} if (!defined $atomlist);
	
	ATOM: foreach my $atom (@$atomlist) {
		
		# Set the $mass variable to the mass of this element
		# If this is not available then set to 0
		my $mass = $Silico::Atomic_Masses[$atom->{ELEMENT_NUM}] || 0;
		
		# Add all atoms whose mass is 0 to an array for warning
		if ($mass == 0) {
			push @masswarn, "Atom: $atom->{NUM} ($atom->{NAME}), Element: $atom->{ELEMENT} ($atom->{ELEMENT_NUM})";
			next ATOM;
		}
		
		# Add all atoms with no X-coord to an array for warning
		if (!defined $atom->{X}) {
			push @xwarn, "Atom: $atom->{NUM} ($atom->{NAME}), Y-coord $atom->{Y}, Z-coord $atom->{Z}";
			next ATOM;
		}
		
		# Add all atoms with no Y-coord to an array for warning
		if (!defined $atom->{Y}) {
			push @ywarn, "Atom: $atom->{NUM} ($atom->{NAME}), X-coord $atom->{X}, Z-coord $atom->{Z}";
			next ATOM;
		}
		
		# Add all atoms with no Z-coord to an array for warning
		if (!defined $atom->{Z}) {
			push @zwarn, "Atom: $atom->{NUM} ($atom->{NAME}), X-coord $atom->{X}, Y-coord $atom->{Y}";
			next ATOM;
		}
		
		# Add a weighted component for this atom to variables
		# These describe the weighted average coordinate - the centre of mass
		$xsum += ($atom->{X} * $mass);
		$ysum += ($atom->{Y} * $mass);
		$zsum += ($atom->{Z} * $mass);
		
		# Add the mass to a variable for checking
		$total_mass += $mass;
		
		# Count the atoms for checking
		++$i;
	}
	
	# Print a warning showing all atoms in this molecule that are missing X-coords
	if ($xwarn[0]) {
		silico_msg('w', ($#xwarn + 1)." atoms are missing X-coordinates!\n",
				"Molecule: $mol->{NAME}");
		if ($xwarn[10]) {
			silico_msg('q', "The first 10 atoms missing X-coordinates are:\n");
			for ($c = 0; $c <= 9; ++$c) {
				silico_msg('q', "$xwarn[$c]\n");
			}
		} else {

			silico_msg('q', "Atoms missing X-coordinates are:\n");
			foreach (@xwarn) {
				silico_msg('q', "$_\n");
			}
		}
		silico_msg('q', "\n");
	}
	
	# Print a warning showing all atoms in this molecule that are missing Y-coords
	if ($ywarn[0]) {
		silico_msg('w', ($#ywarn + 1)." atoms are missing Y-coordinates!\n",
				"Molecule: $mol->{NAME}");
		if ($ywarn[10]) {
			silico_msg('q', "The first 10 atoms missing Y-coordinates are:\n");
			for ($c = 0; $c <= 9; ++$c) {
				silico_msg('q', "$ywarn[$c]\n");
			}
		} else {

			silico_msg('q', "Atoms missing Y-coordinates are:\n");
			foreach (@ywarn) {
				silico_msg('q', "$_\n");
			}
		}
		silico_msg('q', "\n");
	}
	
	# Print a warning showing all atoms in this molecule that are missing Z-coords
	if ($zwarn[0]) {
		silico_msg('w', ($#zwarn + 1)." atoms are missing Z-coordinates!\n",
				"Molecule: $mol->{NAME}");
		if ($zwarn[10]) {
			silico_msg('q', "The first 10 atoms missing Z-coordinates are:\n");
			for ($c = 0; $c <= 9; ++$c) {
				silico_msg('q', "$zwarn[$c]\n");
			}
		} else {

			silico_msg('q', "Atoms missing Z-coordinates are:\n");
			foreach (@zwarn) {
				silico_msg('q', "$_\n");
			}
		}
		silico_msg('q', "\n");
	}
	
	# Print a warning showing all atoms in this molecule that have no mass
	if ($masswarn[0]) {
		silico_msg('w', ($#masswarn + 1)." atoms are missing mass data!\n",
				"Molecule: \"$mol->{NAME}\"\n");
		if ($masswarn[10]) {
			silico_msg('q', "The first 10 atoms missing mass data are:\n");
			for ($c = 0; $c <= 9; ++$c) {
				silico_msg('q', "$masswarn[$c]\n");
			}
		} else {
			silico_msg('q', "Atoms missing mass data are:\n");
			foreach (@masswarn) {
				silico_msg('q', "$_\n");
			}
		}
		silico_msg('q', "\n");
	}
	
	# Print a warning if there are no atoms in this molecule
	if ($i == 0) {
		silico_msg('e', "No atoms are present in this molecule.\n",
				"Unable to calculate centre of mass.\n",
				"Skipping.\n");
		return undef;
	}
	
	# Print a warning if such atoms as are in the molecule are massless
	if ($total_mass == 0) {
		silico_msg('e', "This molecule has no molecular weight.\n",
				"Unable to calculate centre of mass.\n",
				"Skipping.\n");
		return undef;
	}
	
	# Divide the "total mass coordinate" by the total mass to get the weighted average
	my $xc = $xsum/$total_mass;
	my $yc = $ysum/$total_mass;
	my $zc = $zsum/$total_mass;
	
	return ($xc, $yc, $zc);
}

sub mol_orient_polarity {

	#<
	#? Orient molecule polarity dipole along axis
	#; Requres: molecule, align vector (optional)
	#; Returns: nothing
	#. Default aligns polar centre along +Y axis. -xy flag aligns polar end to +Z
	#. Nonpolar centre is translated to the origin
	#>
	

	my $mol = shift;
	
	my @target = @_;
	
	if (!(defined $target[0] && defined $target[1] && defined $target[2])) {
	
		if (get_sflag('xy')) {
			#print "Aligning along Z axis\n";
			@target = (0, 0, 1);
		} else {
			#print "Aligning along Y axis\n";
			@target = (0, 1, 0);
		}
	}

	my ($p, $n) = get_polar_nonpolar_centres($mol);
			
	# Translate nonpolar centre to origin
	molecule_translate($mol, -$n->[0], -$n->[1], -$n->[2]);
			
	# Rotate
	my ($rotmat, $axis) = calc_rotation_vector_onto_vector($p->[0] - $n->[0], $p->[1] - $n->[1], $p->[2] - $n->[2], @target);
	molecule_rotate($mol, @$rotmat);
	
	#print "Rotmat: @$rotmat\n";
}

sub get_polar_nonpolar_centres {

	#<
	#? Get molecule polar and nonpolar centres of mass
	#; Requres: molecule
	#; Returns: polar and nonpolar xyz coordinates
	#. Polar atoms are defined as N,O atoms. Everything else is nonpolar. 
	#  Hydrogens are ignored
	#>

	my $mol = $_[0];
	
	my $p_count;
	my $n_count;
	my $p;
	my $n;
	
	foreach my $atom (atoms($mol)) {
	
		my $el = $atom->{ELEMENT};
		next if $atom->{ELEMENT} eq 'H';
	
		if ($el eq 'N' || $el eq 'O') {
		
			my $i = 0;
			foreach (qw(X Y Z)) {
			
				$p->[$i] += $atom->{$_};
				++$i;
			}
			++$p_count;
		
		} else {
		
			my $i = 0;
			foreach (qw(X Y Z)) {
			
				$n->[$i] += $atom->{$_};
				++$i;
			}
			++$n_count;
		}
	}

	if ($p_count) {
	
		foreach my $i (0..2) {
			$p->[$i] /= $p_count;
		}
		
		$mol->{POLAR_CENTRE} = $p;
		
	} else {
	
		silico_msg('w', "No polar atoms in molecule\n");
		@{$mol->{POLAR_CENTRE}} = ();
	}
	
	if ($n_count) {
	
		foreach my $i (0..2) {
			$n->[$i] /= $n_count;
		}
		
		$mol->{NONPOLAR_CENTRE} = $n;
		
	} else {
	
		silico_msg('w', "No nonpolar atoms in molecule\n");
		@{$mol->{NONPOLAR_CENTRE}} = ();	
	}
	
	return $p, $n;
	
}


sub calc_distance_function {
	
	#<
	#? Calculates a repulsive potential energy term between two sets of atoms
	#. This is based on distance from others in the molecule.  This routine
	#  makes use of the cutoff distance in the global variable $cut and a minimum
	#  approach distance (clash). A fixed penalty is added for each atom clash.
	#; Requires: Molecule, atom list 1, atom list2
	#; Returns: Distance function, number of clashes
	#>
	
	my $mol = $_[0];
	my $list1 = $_[1];
	my $list2 = $_[2];
	my $clash = $_[3];
	my $cut = $_[4];
	
	my $cutoffsq = $cut ** 2; # $cut is a global variable set using the -cut flag
	my $clashsq = $clash ** 2; # $clash is a global variable set using the -clash flag
	my $clashcount = 0;
	my $clash_e = 0;
	my $distfunc = 0;
	my $ints = 0;
	my $starttime = (times)[0];
		
	# Define the function at the cutoff point (used to make the overall function go to zero at the cutoff)
	my $cutoff_func = 1/($cut+1)**2;
	
	foreach my $atom1 (@$list1) {
		
		my @connect;
		
		# Get coordinates for this atom
		my $x1 = $atom1->{X};
		my $y1 = $atom1->{Y};
		my $z1 = $atom1->{Z};
		
		foreach my $connum (@{$atom1->{CONNECT}}) {
			push @connect, $mol->{ATOMS}[$connum];
		}
		
		ATOMTWO: foreach my $atom2 (@$list2) {
			
			next ATOMTWO if $atom1 == $atom2;
			
			my $a = ($atom2->{X} - $x1) **2;
			next ATOMTWO if $a > $cutoffsq;
			my $b = ($atom2->{Y} - $y1) **2;
			next ATOMTWO if $b > $cutoffsq;
			my $c = ($atom2->{Z} - $z1) **2;
			next ATOMTWO if $c > $cutoffsq;
			
			my $distancesq = $a + $b + $c;
			
			next ATOMTWO if $distancesq > $cutoffsq;
			
			foreach my $con (@connect) {
				next ATOMTWO if $con == $atom2;
			}
			
			++$ints;
			
			if ($distancesq <= $clashsq) {
				$clash_e += 10;
				++$clashcount;
			}
			
			my $dist = sqrt ($distancesq);
			
			my $func = 1 / ($dist + 1) ** 2 - $cutoff_func;
						
			$distfunc += $func;
		}
	}
	
	my $time = sprintf "%.2f", (times)[0]-$starttime;
		
	silico_msg('g', "Function $distfunc interactions $ints time $time\n");
	
	return $distfunc + $clash_e, $clashcount;
}

sub calc_box_size {

	#<
	#? Determine an appropriate periodic cell size for a molecule
	#  NOTE: Current assumption is that alpha = beta = gamma = 90 deg!!!
	#. Uses specified values if supplied
	#. Will then use mol->{CELL} if defined
	#. Otherwise prompts for values
	#. $mol->{CELL} is set_box
	#: Requires: molecule, boxx, boxy, boxz, flag to translate to centre,
	#  additional margin around structure, flag to ignore dimensions contained
	#  in molecule file
	#: Returns: boxx, boxy, boxz
	#>

	my $mol = $_[0];
	my $cell;
	@$cell = ($_[1], $_[2], $_[3]);
	my $trans = $_[4];
	my $margin = $_[5] || 0;
	my $ignore = $_[6] || 0;
	
	silico_msg ('d', "Molecule contains no atoms\n") if !defined $mol->{ATOMS}[0];
	
	silico_msg ('c', heading ("Molecule dimensions\n"));
	my @dim = molecule_centre($mol);
	silico_msg('c', "\t   x         y         z\n",
		sprintf("max:\t%8.2f %8.2f %8.2f\nmin:\t%8.2f %8.2f %8.2f\ncentre:\t%8.2f %8.2f %8.2f\n\nedge:\t%8.2f %8.2f %8.2f\n\n",
		@dim, $dim[0]-$dim[3], $dim[1]-$dim[4], $dim[2]-$dim[5]));
	
	if ($trans) {
		silico_msg('c', "Moving centre of molecule ".($mol->{NAME} || 'Mol')." to (0,0,0)\n");
		molecule_trans_to_centre($mol, 'cell', 'all', $cell);
	}
	
	# Calculate dimensions
	#---------------------
	
	if (defined $cell->[0] && defined $cell->[1] && defined $cell->[2]) {
	
		silico_msg('c', "\nBox dimensions: $cell->[0] A, $cell->[1] A, $cell->[2] A\n\n");
		@{$mol->{CELL}} = @$cell;
		$mol->{CELL_ALPHA} = $mol->{CELL_BETA} = $mol->{CELL_GAMMA}= 90;
		return @$cell;		
	} 
	
	if (defined $mol->{CELL}[0]) {
	
	 	if ($ignore) {
	
			silico_msg('c', "\nIgnoring existing cell dimensions: @{$mol->{CELL}}\n\n");
		
		} else {
	
			silico_msg ('c',"File contains cell dimensions: @{$mol->{CELL}}. Using these.\n");
			$cell->[0] = $mol->{CELL}[0];
			$cell->[1] = $mol->{CELL}[1];
			$cell->[2] = $mol->{CELL}[2];
			
			silico_msg('c', "\nBox dimensions: $cell->[0] A, $cell->[1] A, $cell->[2] A\n\n");
			@{$mol->{CELL}} = @$cell;
			$mol->{CELL_ALPHA} = $mol->{CELL_BETA} = $mol->{CELL_GAMMA}= 90;
			return @$cell;
		}
	}
	
	silico_msg('c', "Using margin ($margin A)\n\n") if $margin;

	my $i = 0;
	foreach ('X', 'Y', 'Z') {
	
		# If we are adding a margin, then round the box size to
		# the nearest integer value.  Otherwise use the exact value
		my $size = sprintf "%.2f", $dim[$i]-$dim[$i+3];
		my $size_i;
		if ($margin) {
			$size_i = int($size+$margin+0.5);
		} else {
			$size_i = $size;
		}

		# Use any cell dimensions provided on command line
		# This is used if less than 3 dimensions are provided
		if (defined $cell->[$i]) {
			$size_i = $cell->[$i];
		}
		
		my $val;
		while (1) {
		
			$val = prompt("Box size $_ ($size_i A): ");
			chomp $val;
			$val =~ s/ //g;
			
			last if $val eq '';
			last if check_data_type($val, 'DECIMAL');
		}
			
		$cell->[$i] = $val || $size_i;
		++$i;
	}

	silico_msg('c', "\nBox dimensions: $cell->[0] A, $cell->[1] A, $cell->[2] A\n\n");

	$mol->{CELL} = $cell;
	$mol->{CELL_ALPHA} = $mol->{CELL_BETA} = $mol->{CELL_GAMMA}= 90;
	return @$cell;
}

sub makeboxmolecule {

	#<
	#? Make a box molecule using (on-axis cuboid) unit cell size
	#; Requires: molecule cell, flag to centre box on 0,0,0
	#; Returns: box molecule
	#>
	
	my $cell = $_[0];
	my $centre = $_[1];
	
	my @vals = (0,1);
	
	if (defined $centre) {
		@vals = (-0.5, 0.5);
	}
	
	my $mol = create_molecule('box');
	
	my $count = 0;
	foreach my $i (@vals) {
		foreach my $j (@vals) {
			foreach my $k (@vals) {
				++$count;
				mol_add_atom($mol, "Du$count", "Du", $i * $cell->[0], $j * $cell->[0], $k * $cell->[0], 'BOX', 1);
			}
		}
	}
	
	bond_create ($mol, 0, 1);
	bond_create ($mol, 0, 2);
	bond_create ($mol, 0, 4);
	bond_create ($mol, 1, 3);
	bond_create ($mol, 1, 5);
	bond_create ($mol, 2, 3);
	bond_create ($mol, 2, 6);
	bond_create ($mol, 3, 7);
	bond_create ($mol, 4, 5);
	bond_create ($mol, 4, 6);
	bond_create ($mol, 5, 7);
	bond_create ($mol, 6, 7);
	
	return $mol;
}

################################
# ORTHOGONAL LEAST-SQUARES FIT #
################################

# This section of code fits a plane to a set of atoms, based on minimising
# the sum of the squares of the distances of the atoms to the plane.
#
# Instructions provided by David Eberly, of Geometric Tools LLC.
# http://www.geometrictools.com

sub orthogonal_least_squares_plane_atoms {
	
	#<
	#? Fits a plane to a set of atoms, using the orthogonal least squares method.
	#; Requires: molecule, set of atoms
	#; Returns: One point in the plane, normal vector to the plane
	#>
	
	my $mol = $_[0];
	my $list = $_[1] || $mol->{ATOMS};
	
	# Test to see whether Math::MatrixReal is installed
	my $matflag = 0;
	foreach (@INC) {
		my $matpath = "$_/Math/MatrixReal.pm";
		my $mat = 0;
		++$mat if ( -f $matpath || -l $matpath );
		++$mat if ( -r $matpath );
		
		++$matflag if $mat == 2;
	}
	
	if (!$matflag) {
		silico_msg('d', "The Perl library Math::MatrixReal is not installed!\n",
				"This third-party library is necessary for running this script.\n",
				"For information on obtaining Math::MatrixReal, please visit\n",
				"http://www.cpan.org.\n");
	}
	
	# Use matrix mathematics
	require Math::MatrixReal;
	
	# Centroid, which can be shown to lie within the plane
	my ($cx, $cy, $cz) = centroid($mol, $list);
	
	my ($m11, $m12, $m13);
	my ($m21, $m22, $m23);
	my ($m31, $m32, $m33);
	foreach my $atom (@$list) {
		
		my $ax = $atom->{X} - $cx;
		my $ay = $atom->{Y} - $cy;
		my $az = $atom->{Z} - $cz;
		
		$m11 += $ax*$ax;
		$m12 += $ax*$ay;
		$m13 += $ax*$az;
		$m21 += $ay*$ax;
		$m22 += $ay*$ay;
		$m23 += $ay*$az;
		$m31 += $az*$ax;
		$m32 += $az*$ay;
		$m33 += $az*$az;
	}
	
	my $matrix = Math::MatrixReal->new_from_rows([[$m11,$m12,$m13], [$m21,$m22,$m23], [$m31,$m32,$m33]]);
	my ($values, $vectors) = $matrix->sym_diagonalize();
	
	my $ev1 = $values->element(1,1);
	my $ev2 = $values->element(2,1);
	my $ev3 = $values->element(3,1);
	
	my @vector1 = ($vectors->element(1,1),$vectors->element(2,1),$vectors->element(3,1));
	my @vector2 = ($vectors->element(1,2),$vectors->element(2,2),$vectors->element(3,2));
	my @vector3 = ($vectors->element(1,3),$vectors->element(2,3),$vectors->element(3,3));

	my @vector;
	if ($ev1 < $ev2 && $ev1 < $ev3) {
		@vector = unit_vector(@vector1);
	} elsif ($ev2 < $ev3) {
		@vector = unit_vector(@vector2);
	} else {
		@vector = unit_vector(@vector3);
	}
	
	return ($cx, $cy, $cz, @vector);
}

return 1;
