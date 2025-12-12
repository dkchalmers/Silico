
sub distance_periodic_bv {

	my $atom1 = $_[0];
        my $atom2 = $_[1];
        my $bv   = $_[2] || croak();

	# Example input:
	my @r1 = ($atom1->{X}, $atom1-{Y}, $atom1->{Z});
	my @r1 = ($atom2->{X}, $atom2-{Y}, $atom2->{Z});

	# GROMACS box vectors (in nm)
	# These can be read from a .gro file last line
	my @a = ($bv->[0], $bv->[3], $bv->[4]);
	my @b = ($bv->[5], $bv->[1], $bv->[6]);
	my @c = ($bv->[7], $bv->[8], $bv->[2]);

	# Displacement in Cartesian
	my @dr = ($r2[0] - $r1[0], $r2[1] - $r1[1], $r2[2] - $r1[2]);

	# Box matrix and its inverse
	my @M    = (@a, @b, @c);
	my @Minv = inverse3x3(\@M);

	# Convert displacement to fractional coordinates
	my @s = matvec(\@Minv, \@dr);

	# Apply minimum image convention
	for my $i (0 .. 2) {
		$s[$i] -= sprintf("%.0f", $s[$i]);    # equivalent to nearest integer
	}

	# Back to Cartesian
	my @dr_min = matvec(\@M, \@s);

	# Distance
	my $dist = sqrt($dr_min[0]**2 + $dr_min[1]**2 + $dr_min[2]**2);
}

sub inverse3x3 {

	my ($a11, $a12, $a13, $a21, $a22, $a23, $a31, $a32, $a33) = @_;

	my $det = $a11 * ($a22 * $a33 - $a23 * $a32) - $a12 * ($a21 * $a33 - $a23 * $a31) + $a13 * ($a21 * $a32 - $a22 * $a31);

	die "Matrix determinant is zero!" if abs($det) < 1e-12;

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

