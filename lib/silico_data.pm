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
#! silico_data.pm
#? Various subroutines for handling text data and other data structures.
#. $Revision: 1.18.2.1.2.5 $
#>

use strict;
package Silico;

sub calc_human_readable_time {

	#<
	#? Calc_human_readable_time
	#. Takes a length of time (in seconds) and converts it into ddd:hh:mm:ss.ss
	#; Requires: Time in seconds
	#; Returns: Time as a text string
	#>

	my $time = $_[0];
	
	my $seconds = 0;
	my $minutes = 0; # Minutes rounded to nearest minute
	my $hours = 0;
	my $days = 0;
	my $weeks = 0;
	
	my $string;
	
	$minutes = int ($time / 60);
	$seconds = $time - ($minutes * 60);
	$seconds = sprintf ("%.2f", $seconds);
	$hours = int ($minutes / 60);
	$minutes = $minutes - ($hours * 60);
	$days = int ($hours / 24);
	$hours = $hours - ($days * 24);
	$weeks = int ($days / 7);
	$days = $days - ($weeks * 7);
	
	$string = "$seconds seconds" if ($minutes == 0 && $hours == 0 && $days == 0);
	
	if ($minutes > 0 || $hours > 0 || $days > 0 || $weeks > 0) {
		if ($seconds < 10) {
			$string = "$minutes"."m:0$seconds"."s";
		} else {
			$string = "$minutes"."m:$seconds"."s";
		}
		if ($minutes < 10) {
		$string = "0"."$string";
		}
	}

	if ($hours > 0 || $days > 0 || $weeks > 0) {
		$string = "$hours"."h:$string";
		if ($hours < 10) {
			$string = "0"."$string";
		}
	}
	
	if ($days > 0 || $weeks > 0) {
		$string = "$days"."d:$string";
	}
	
	if ($weeks > 0) {
		$string = "$weeks"."w:$string";
	}
	
	return $string;
}

sub ch {

	#<
	#? Check if a hashkey is defined; return 0 if not
	#; Requires: Hash, key
	#; Returns: Value of key, or 0 if undefined
	#>
	
	my $hash = $_[0];
	my $key = $_[1];
	
	
	if (defined $hash->{$key}) {return $hash->{$key}}
	
	return 0;
}

sub check_data_type {

	#<
	#? Checks an input string to see if it is a specified type
	#. Currently looks for: ALPHA (A-Z characters only), WORD (A-Z plus _),
	#  INT or INTEGER (integers), DECIMAL (real number in standard or mantissa format),
	#  FILENAME, BOOLEAN
	#. Also has numeric checks >, >=, <, <=
	#; Requires: String of data, string of type information
	#; Optional: Description, Suppression of standard warnings (will still warn if
	#  serious errors are encountered)
	#; Returns: -1 if a serious error; 0 if a data type problem; 1 if OK
	#>
	
	my $data = $_[0];
	my $type = $_[1];
	my $desc = $_[2]  || '';
	my $suppress = $_[3];
	
	$desc= "$desc " if $desc;
	
	my $ge;
	my $geflag;
	my $gt;
	my $gtflag;
	my $le;
	my $leflag;
	my $lt;
	my $ltflag;
	my $value;
	
	my $intflag = 0;
	
	my @f = split (" ", $type);
	
	$f[0] = uc $f[0];
	
	my $alpha = 1 if ($f[0] eq 'ALPHA');
	my $integer = 1 if ($f[0] eq 'INT' || $f[0] eq 'INTEGER');
	my $decimal = 1 if ($f[0] eq 'DECIMAL');
	my $word = 1 if ($f[0] eq 'WORD');
	my $filename = 1 if ($f[0] eq 'FILENAME');
	my $boolean = 1 if ($f[0] eq 'BOOL' || $f[0] eq 'BOOLEAN');
	
	my $text = 1 if ($word || $alpha || $filename);
	my $numeric = 1 if ($integer || $decimal);
	
	for (my $i = 1; $i <= $#f; ++$i) {
			
		my $f = $f[$i];
		
		if ($f eq '>' && !defined $gtflag) {
			$gtflag = 1;
			next;
		}
		
		if ($f eq '<' && !defined $ltflag) {
			$ltflag = 1;
			next;
		}

		if ($f eq '>=' && !defined $geflag) {
			$geflag = 1;
			next;
		}

		if ($f eq '<=' && !defined $leflag) {
			$leflag = 1;
			next;
		}
		
		if ($gtflag || $geflag || $ltflag || $leflag) {
			if ($f =~ /^-?\d+\.?\d*$/ || $f =~ /^-?\d*\.?\d+$/) {
				if ($gtflag) {
					$gt = $f;
					$gtflag = undef;
					next;
				} elsif ($geflag) {
					$ge = $f;
					$geflag = undef;
					next;
				} elsif ($ltflag) {
					$lt = $f;
					$ltflag = undef;
					next;
				} elsif ($leflag) {
					$le = $f;
					$leflag = undef;
					next;
				}
			} else {
				silico_msg('e', $desc."\"$f\", used in a numeric comparison, is not numeric!\n");
				return -1;
			}
		}
		
		if (	($f eq '>' && defined $gt) || ($f eq '>=' && defined $ge)
			|| ($f eq '<' && defined $lt) || $f eq '<=' && defined $le) {
			silico_msg('e', $desc."More than one numeric comparison of the same type was used!\n");
			return -1;
		}
	}
	
	if ($text && (defined $gt || defined $ge || defined $lt || defined $le)) {
	
		silico_msg('e', $desc."You've tried to use a numeric comparison when you want your data to be non-numeric!\n");
		return -1;
	}
	
	elsif ($text && $numeric) {
	
		silico_msg('e', $desc."You've tried to check for both numeric and text in one go!\n");
		return -1;
	}
	
	elsif ((defined $gt && defined $lt && $gt > $lt) || (defined $ge && defined $lt && $ge > $lt)
		|| (defined $gt && defined $le && $gt > $le) || (defined $ge && defined $le && $ge > $le)) {
			
		if ($gt >= $lt) {
			silico_msg('e', $desc."A number can't be larger than $gt and smaller than $lt at the same time!\n");
			return -1;
		} elsif ($ge >= $lt) {
			silico_msg('e', $desc."A number can't be at least $ge and smaller than $lt at the same time!\n");
			return -1;
		} elsif ($gt >= $le) {
			silico_msg('e', $desc."A number can't be larger than $gt and at most $le at the same time!\n");
			return -1;
		} elsif ($ge > $le) {
			silico_msg('e', $desc."A number can't be at least $ge and at most $le at the same time!\n");
			return -1;
		}
	}
	
	elsif (!defined $data) {return 0}
	
	elsif ($alpha) {
	
		if ($data !~ /^[A-Za-z]*$/) {
			silico_msg('w', "\"$data\" is not an alphabetical string.\n") unless $suppress;
			return 0;
		
		} else {
			return 1;
		}
	}
	
	elsif ($boolean) {
		
		return 1 if (standardise_boolean($data) == 1 || standardise_boolean($data) == 0);
		
		silico_msg('w', "\"$data\" is not a boolean value.\n") unless $suppress;
		return 0;
	}
	
	elsif ($numeric) {

		# Delete leading and trailing white space
		$data =~ s/^\s*//;
		$data =~ s/\s*$//;

		# Check that the data is actually an integer
		if ($integer) {
			if ($data =~ /^-?\d+$/) {
				$intflag = 1;
			} elsif ($data =~ /^-?\d*\.?\d+[Ee][+-]?\d+$/ || $data =~ /^-?\d+\.?\d*[Ee][+-]?\d+$/) {
				my @g = split /[Ee]/, $data;
				my $mantissa = $g[0];
				my $exponent = $g[1];
				$value = $mantissa * (10 ** $exponent);
				if ($value == 0) {
					$intflag = 1;
				} elsif ($value < 1) {
					$intflag = 0;
				} elsif (sprintf("%.12f", $value) == int($value)) {
					$intflag = 1;
				} else {
					$intflag = 0;
				}
			} else {
				$intflag = 0;
			}
						
			if (!$intflag) {
				silico_msg('w', $desc."Error: \"$data\" is not an integer.\n") unless $suppress;
				return 0;
			} else {
				return 1;
			}
		
		# Check that the data is actually a decimal
		} elsif ($decimal && $data !~ /^-?\d+\.?\d*$/ && $data !~ /^-?\d*\.?\d+$/
			&& $data !~ /^-?\d+\.?\d*[Ee][+-]?\d+$/
			&& $data !~ /^-?\d*\.?\d+[Ee][+-]?\d+$/) {
		
			silico_msg('w', "\"$data\" is not decimal.\n") unless $suppress;
			return 0;
		
		# Here follow a series of numeric checks
		} elsif (defined $gt && $data <= $gt) {
		
			silico_msg('w', $desc."Error: \"$data\" is less than or equal to $gt.\n") unless $suppress;
			return 0;
			
		} elsif (defined $ge && $data < $ge) {
		
			silico_msg('w', $desc."Error: \"$data\" is less than $ge.\n") unless $suppress;
			return 0;
			
		} elsif (defined $lt && $data >= $lt) {
		
			silico_msg('w', $desc."Error: \"$data\" is greater than or equal to $lt.\n") unless $suppress;
			return 0;
		
		} elsif (defined $le && $data > $le) {
		
			silico_msg('w', $desc."Error: \"$data\" is greater than $le.\n") unless $suppress;
			return 0;
			
		} else {
			return 1;
		}
	}
	
	elsif ($word) {
	
		if ($data !~ /^\w+$/) {
		
			silico_msg('w', $desc."Error: \"$data\" does not contain all word characters (A-Z, 0-9 or _).\n") unless $suppress;
			return 0;
			
		} else {
			return 1;
		}
	}
	
	elsif ($filename) {
	
		if ($data !~ /^[-\w.]+$/) {
		
			if (!$suppress) {
			
				silico_msg('w', $desc."\"$data\" does not contain only characters\n",
						"suitable to filenames (A-Z, a-z, 0-9, _, - or .).\n");
				return 0;
			}
		} else {
			return 1;
		}
	} else {
	
		silico_msg('e', $desc."Unsure how to perform the specified check!\n",
				"Data: \"$data\"\n",
				"Type to check: \"$type\"\n");
	
		return -1;
	}
}

sub en {
	
	#<
	#? Check that two numbers are equal.
	#  Equality means that either: neither number is defined, or they
	#  are equal to each other (as ascertained using NOT !=).
	#. Inequality is assumed if one is defined and the other is not,
	#  even if the defined number is 0.
	#; Requires: numbers
	#; Returns: 1 if they are, 0 otherwise
	#>
	
	my $number1 = $_[0];
	my $number2 = $_[1];
	
	if (defined $number1 xor defined $number2) {
		return 0;
	} elsif (defined $number1 && defined $number2 && $number1 != $number2) {
		return 0;
	} else {
		return 1;
	}
}

sub et {
	
	#<
	#? Check that two text strings are equivalent.
	#. Equivalence means that either: neither string is defined, or
	#  they are equal to each other (as ascertained using NOT ne).
	#. Non-equivalence is assumed if one is defined and the other is not,
	#  even if the defined string is the empty string.
	#; Requires: strings
	#; Returns: 1 if they are, 0 otherwise
	#>
	
	my $string1 = $_[0];
	my $string2 = $_[1];
	
	if (defined $string1 xor defined $string2) {
		return 0;
	} elsif (defined $string1 && defined $string2 && $string1 ne $string2) {
		return 0;
	} else {
		return 1;
	}
}

sub fnfp {
	
	#<
	#? Format Number For Printing.
	#; Take a numeric value, and print it. If the numeric value is
	#  not defined, return a string of the user's choice (default:
	#  'not defined').
	#. This can handle NaN strings, but unlike the subroutine ftfp,
	#  it does not put the strings in quotes, so whitespace is not
	#  obvious.
	#; Requires: string, optional 'not defined' string
	#; Returns: string
	#>
	
	my $input = $_[0];
	my $nd = $_[1] || 'not defined';
	
	my $output;
	
	$output = (defined $input ? "$input" : "$nd");
	
	return $output;
}

sub ftfp {
	
	#<
	#? Format Text For Printing.
	#; Take a text value, enclose it in double quotes, and print it.
	#  If the text value is not defined, return a string of the user's
	#  choice (default: 'not defined').
	#. Note that the 'not defined' string will not be printed in quotes.
	#  This allows the user to distinguish between a variable containing
	#  the text "not defined" and an undefined variable.
	#. This can handle numbers, but unlike the subroutine fnfp, it
	#  puts the data in quotes, which can appear ugly for numbers.
	#; Requires: string, optional 'not defined' string
	#; Returns: string
	#>
	
	my $input = $_[0];
	my $nd = $_[1] || 'not defined';
	
	my $output;
	
	$output = (defined $input ? "\"$input\"" : "$nd");
	
	return $output;	
}	

sub is_integer {

	#<
	#? Return 1 if a number is an integer
	#; Requires: Number, optional accuracy (defaults to 6 decimal places)
	#; Returns: 1 or 0
	#>
	
	my $val = $_[0];
	my $acc = $_[1] || 6;
	
	return 1 if ((int($val))*(10 ** $acc)) == (int($val*(10 ** $acc)));
	return 0
}

sub parse_comma_separated_numbers {
	
	#<
	#? Take a list of comma-separated numbers, which may include
	#  number ranges, and turn it into a list of individual
	#  numbers in array format.
	#; Requires: string
	#; Returns: array
	#>
	
	my $string = $_[0];
	
	my @largelist;
	my @list;
	my @split;
	my @vals;
	
	# Split on commas
	@split = split(/\s*,\s*/, $string);
	
	# Convert number ranges to lists of numbers
	foreach (@split) {
		@vals = ();
		@vals = number_range($_);
		push @largelist, @vals;
	}
	
	@list = remove_duplicates(@largelist);
	
	return @list;
}

sub remove_duplicates {
	
	#<
	#? Remove duplicate entries in a list.
	#; Requires: list
	#; Returns: pruned list
	#>
		
	my %exists;
	my @out;
	
	foreach (@_) {
		next if $exists{"$_"};
		
		push @out, $_;
		++$exists{"$_"};
	}
	
	return @out;
}

sub number_range {
	
	#<
	#? Turn a range of numbers (indicated by X-Y, where X and Y are integers)
	#  into a list of those numbers spanned by the range.
	#; Requires: integer range as text string
	#; Returns: List of integers
	#>
	
	my $val = $_[0];
	
	my $error = 0;
	my @f;
	my $high;
	my $i;
	my @list;
	my $low;
	my $remainder;
	my $rightmost;
	
	# First, remove any spaces in the string.
	$val =~ s/ //g;
	
	# If the value is merely a number, return it.
	if (check_data_type($val, 'DECIMAL', 1)) {
		@list = ($val);
		return @list;
	}
	
	# Now, we want to check: is the overall string of the correct kind?
	if ($val !~ /^-?\d+--?\d+$/) {
		silico_msg('e', "Argument is not a range of integers!\n",
				"Argument: $_[0]\n",
				"Argument should look like: [-]<number>-[-]<number>, where square brackets\n",
				"([]) indicate an optional character.\n",
				"Attempting to proceed with non-numeric value, \"NaN\".\n");
		@list = ("NaN");
		return @list;
	}
	
	# Temporarily replace dashes with other things.
	$val =~ s/^-/N/;
	$val =~ s/--/-N/;
	$val =~ s/-/:/;
	
	# Split on the colon
	@f = split(/:/, $val);
	
	if ($#f != 1) {
		silico_msg('e', "Argument contains too many numbers to be a range of integers!\n",
				"Argument: $_[0]\n",
				"Attempting to proceed with non-numeric value, \"NaN\".\n");
		@list = ("NaN");
		return @list;
	}
	
	$low = $f[0];
	$high = $f[1];
	
	# Replace N with dash again.
	$low =~ s/N/-/;
	$high =~ s/N/-/;
	
	# Sometimes, $high can be less than $low. Other times, definitely not.
	# Acceptable:
	#  - $high is a single digit, greater than the last digit of $low (e.g., 517-9)
	#  - $high is a pair of digits, greater than the last two digits of $low (e.g., 636-54)
	#  - and so on
	# All other cases are unacceptable, except for $high to be greater than $low.
	# Test this.
	
	if ($high < $low) {
		if ($high < 0) {
			$error = 1;
		} else {
			# Get that part of the low number which corresponds
			# to the high number.
			$rightmost = substr($low, -1*length($high));
			$remainder = $low - $rightmost;
			if ($high < $rightmost) {
				$error = 1;
			} else {
				$high += $remainder;
				silico_msg('w', "Assuming that the range has been abbreviated\n",
						"in the upper limit. Upper limit evaluated as: $high\n");
			}
		}
	}
	
	if ($error) {
		silico_msg('e', "Upper end of range is numerically less than lower end!\n",
				"Argument: $_[0]\n",
				"Attempting to proceed with non-numeric value, \"NaN\".\n");
		@list = ("NaN");
		return @list;
	}
	
	foreach ($i = $low; $i <= $high; ++$i) {
		push @list, $i;
	}
	
	return @list;
}		

sub cleanline {
	
	#<
	#? Take a line of text and remove leading and trailing spaces
	#; Requires: Text
	#; Returns: Text
	#>
	
	my $line = $_[0];
	
	chomp $line;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	
	return $line;
}

sub splitline {
	
	#<
	#? Take a line of text, clean it up and split it into an array
	#; Requires: Text
	#; Returns: Array
	#>
		
	my @f;
	my $line;
	
	$line = cleanline($_[0]);
	@f = split ' ', $line;
	
	return @f;
}

sub standardise_boolean {
	
	#<
	#? Turn a Silico boolean value into 1 or 0 depending on what it is.
	#; Requires: Text
	#; Returns: Number
	#>
		
	return 1 if $_[0] =~ /^true$/i;
	return 1 if $_[0] =~ /^t$/i;
	return 1 if $_[0] =~ /^yes$/i;
	return 1 if $_[0] =~ /^y$/i;
	return 1 if $_[0] eq '1';
	
	return 0 if $_[0] =~ /^false$/i;
	return 0 if $_[0] =~ /^f$/i;
	return 0 if $_[0] =~ /^no$/i;
	return 0 if $_[0] =~ /^n$/i;
	return 0 if $_[0] eq '0';
	
	silico_msg('w', "\"$_[0]\" is not a known boolean value!\n",
			"Returning the non-boolean value -1.\n");
	return -1;
}

sub verify_input_data {

	#<
	#? Verifies that input data is of a specified type
	#; Requires:
	#     Data
	#     Expected data type (text string suitable for check_data_type)
	#     Prompt message (eg. "Timestep")
	#     Optional default value after first time
	#; Returns: Modifies value of Data directly
	#>
	
	my $input = $_[0];
	my $type = $_[1];
	my $string = $_[2]; # Flag description
	my $default = $_[3];
	
	# Get memory address of input variable
	my $addr = \$_[0];
	
	my $dstatus;
	my $output;
	my $status;
	
	$dstatus = check_data_type($default, uc($type)) if (defined $default);
	if (defined $dstatus && $dstatus != 1) {
		if ($type =~ /^[aAeEiIoOuU]/) {
			silico_msg('d', "Default value (\"$default\") is not appropriate for $string!\n",
					ucfirst($string)." should be an ".lc($type).".\n");
		} else {
			silico_msg('d', "Default value (\"$default\") is not appropriate for $string!\n",
					ucfirst($string)." should be a ".lc($type).".\n");
		}
	}
	
	#$string = lc($string);
	
	while (1) {
		
		$status = check_data_type($input, uc($type), $string);
		
		if ($status == 1) {
			$output = $input;
			last;
		} elsif ($status == -1) {
			silico_msg('d',	"A serious problem was encountered!\n",
					"Check: ".ucfirst($string).", $type\n");
		}
		
		$input = prompt(ucfirst($string).": ");
		silico_msg('d', "Undefined input\n") if !defined $input;
		$input =~ s/ //g;
		
		if (check_data_type($input, uc($type), $string) == 1) {
			$output = $input;
			last;
		}

		silico_msg('v', "Inappropriate value for $string.\n");
		
		if ($type =~ /^[aAeEiIoOuU]/) {
			silico_msg('v', ucfirst($string)." should be an ".lc($type).".\n");
		} else {
			silico_msg('v', ucfirst($string)." should be a ".lc($type).".\n");
		}
		silico_msg('v', "Using default value: $default.\n") if (defined $default);
		silico_msg('v', "\n");
		
		if (defined $default) {
			$output = $default;
			last;
		}
	}
	
	# Change original variable
	${$addr} = $output;
}

return 1;

