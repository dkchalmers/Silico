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
#! silico_statistics.pm
#? Silico routines to calculate statistical quantities
#. $Revision: 1.9.2.3.2.1 $
#>

use strict;
package Silico;

sub maximum {
	
	#<
	#? Determines the maximum value within a set of data
	#; Requires: Array
	#; Returns: Scalar
	#>
	
	my $max;
	my $val;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A maximum value has not been determined.\n");
		return undef;
	}
	
	foreach $val (@_) {
		if (!defined $val) {
			silico_msg('w', "Undefined value in list!\n");
			next;
		}
		$max = $val if !defined $max;
		$max = $val if $val > $max;
	}
	
	return $max;
}

sub minimum {
	
	#<
	#? Determines the minimum value within a set of data
	#; Requires: Array
	#; Returns: Scalar
	#>
	
	my $min;
	my $val;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A minimum value has not been determined.\n");
		return undef;
	}
	
	foreach $val (@_) {
		if (!defined $val) {
			silico_msg('w', "Undefined value in list!\n");
			next;
		}
		$min = $val if !defined $min;
		$min = $val if $val < $min;
	}
	
	return $min;
}

sub range {
	
	#<
	#? Determines the range of a set of data.
	#; Requires: Array
	#; Returns: scalar
	#>
	
	my $min;
	my $max;
	my $range;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A range has not been calculated.\n");
		return undef;
	}
	
	$min = minimum(@_);
	$max = maximum(@_);
	
	if (defined $max && defined $min) {
		$range = $max-$min;
	} else {
		undef $range;
	}
	
	return $range;
}

sub maximum_minus_mean {
	
	#<
	#? Determine the difference between the mean of a set of data
	#  and the maximum of the same set.
	#; Requires: Array
	#; Returns: Scalar
	#>
	
	my $max;
	my $mean;
	my $val;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A difference has not been calculated.\n");
		return undef;
	}
	
	$max = maximum(@_);
	$mean = mean(@_);
	
	if (defined $max && defined $mean) {
		$val = $max - $mean;
	} else {
		undef $val;
	}
	
	return $val;
}

sub mean_minus_minimum {
	
	#<
	#? Determine the difference between the mean of a set of data
	#  and the minimum of the same set.
	#; Requires: Array
	#; Returns: Scalar
	#>
	
	my $min;
	my $mean;
	my $val;
	
	if (!defined $_[0]) {
		#silico_msg('e', "No input data was provided!\n",
		#		"A difference has not been calculated.\n");
		return undef;
	}
	
	$min = minimum(@_);
	$mean = mean(@_);
	
	if (defined $min && defined $mean) {
		$val = $mean - $min;
	} else {
		undef $val;
	}
	
	return $val;
}


sub mean {

	#<
	#? Calculates the mean of a set of data
	#; Requires: Array
	#; Returns: mean (scalar quantity)
	#>
	
	my $i = 0;
	my $mean;
	my $total;
	my $val;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A mean has not been calculated.\n");
		return undef;
	}
	
	foreach $val (@_) {
		if (!defined $val) {
			silico_msg('e', "Undefined value in list!\n");
			next;
		}
		++$i;
		$total += $val;
	}
	
	$mean = $total/$i;
	
	return $mean;
}

sub median {

	#<
	#? Calculates the median of a set of data
	#; Requires: Array
	#; Returns: median (scalar quantity)
	#>
	
	my $median;
	my $n;
	my @pair;
	my @sorted_input;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A median has not been calculated.\n");
		return undef;
	}

	$n = $#_ + 1;
	
	@sorted_input = sort {$a <=> $b} @_;
	
	# If n is odd...
	# This uses a bitwise AND to determine oddness. If $n is odd,
	# its last binary digit will be a 1. This, when bitwise ANDed to the
	# number 1, will return 1 (or, put another way, not 0).
	if (($n & 1) != 0) {
	
		$median = $sorted_input[(($n+1)/2)-1];
	
	} else {
	
		@pair = ($sorted_input[($n/2)-1], $sorted_input[($n/2)]);
		
		$median = mean(@pair);
	}
	
	return $median;
}

sub mode {

	
	my %counts;
	my $key;
	my $max_value;
	my @modes;
	my @sorted_values;
	my $val;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A mode has not been calculated.\n");
		return undef;
	}
	
	foreach $val (@_) {
		if (!defined $val) {
			silico_msg('w', "Undefined value in list!\n");
			next;
		}
		++$counts{$val};
	}
	
	@sorted_values = sort {$a <=> $b} values(%counts);
	
	$max_value = $sorted_values[$#sorted_values];
	
	return undef if ($max_value == 1);
	
	foreach $key (keys %counts) {
	
		push @modes, $key if ($counts{$key} == $max_value);
	
	}
	
	return @modes;

}

sub median_absolute_deviation {

	#<
	#? Calculates the median and median absolute deviation
	#. See: en.wikipedia.org/wiki/Median_absolute_deviation
	#; Requires: Array
	#; Returns: median (scalar quantity), median absolute deveiation (scalar)
	#>
	
	my @list = @_;
	
	my @deviations;
	my $median = median (@list);
	
	foreach (@list) {
	
		push @deviations, abs($median - $_);
	}
	
	return $median, median(@deviations);

}

sub standard_deviation_bias_corrected {

	#<
	#? Calculates the bias-corrected (n-1) standard deviation of a set of data
	#; Requires: Array
	#; Returns: standard deviation (scalar quantity)
	#>
	
	my $sd;
	my $variance;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A standard deviation has not been calculated.\n");
		return undef;
	}
	
	$variance = variance_bias_corrected(@_);
	
	$sd = sqrt($variance);
	
	return $sd;
}

sub standard_deviation_bias_uncorrected {

	#<
	#? Calculates the bias-uncorrected standard deviation of a set of data
	#; Requires: Array
	#; Returns: standard deviation (scalar quantity)
	#>
		
	my $sd;
	my $variance;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A standard deviation has not been calculated.\n");
		return undef;
	}
	
	$variance = variance_bias_uncorrected(@_);
	
	$sd = sqrt($variance);
	
	return $sd;
}

sub variance_bias_corrected {

	#<
	#? Calculates the bias-corrected (n-1) variance of a set of data
	#; Requires: Array
	#; Returns: variance (scalar quantity)
	#>
		
	my $i;
	my $mean;
	my $sqtotal;
	my $val;
	my $variance;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A variance has not been calculated.\n");
		return undef;
	}
	
	# The variance is zero if we only have one value.
	# This has been put in explicitly so we don't get "0 divided by 0"
	# errors.
	return 0 if !defined $_[1];
	
	$mean = mean(@_);
	
	foreach $val (@_) {
		
		if (!defined $val) {
			silico_msg('e', "Undefined value in list!\n");
			next;
		}
		
		++$i;
		
		$sqtotal += ($val - $mean)**2
	
	}
	
	$variance = $sqtotal / ($i-1);
	
	return $variance;
}

sub variance_bias_uncorrected {

	#<
	#? Calculates the bias-uncorrected (n) variance of a set of data
	#; Requires: Array
	#; Returns: variance (scalar quantity)
	#>
	
	my $i;
	my $mean;
	my $sqtotal;
	my $val;
	my $variance;
	
	if (!defined $_[0]) {
		#silico_msg('w', "No input data was provided!\n",
		#		"A variance has not been calculated.\n");
		return undef;
	}
	
	$mean = mean(@_);
	
	foreach $val (@_) {
	
		if (!defined $val) {
			silico_msg('e', "Undefined value in list!\n");
			next;
		}
		
		++$i;
		
		$sqtotal += ($val - $mean)**2
	
	}
	
	$variance = $sqtotal / $i;
	
	return $variance;
}

return 1;
