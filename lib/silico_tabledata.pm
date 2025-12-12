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
#! silico_tabledata.pm
#? Handle data that is in tabular form
#. $Revision: 1.33.2.4.2.6 $
#>

use strict;
package Silico;

sub read_table_delimited {

	#<
	#? Read a character delimited text file
	#. Creates a data structure that looks like:
	#; $data->{VALUES}[<linennum>][<columnnum>] Values in array
	#; $data->{HEADERS}[<columnnum>] Array of all column headers
	#; $data->{KEYHASH}{} Hash of all column headers to alow translation to column numbers
	#; $data->{NUMLINES} Total number of lines
	#; $data->{NUMKEYS} Total number of columns
	#; Requires: filename, delimiter (strings 'space', 'tab', or 'comma' with the default 'tab'), has_header flag, options
	#; Returns: data structure or undef if file open failed
	#. Options:
	#; QUIET do not print messages
	#; HEADER first line is a header line
	#; CONVERT_EMPTY convert empty values to '' (not undef)
	#; CASE_INSENSITIVE convert all headers to upper case
	#>

	my $infile = $_[0];
	my $delimiter = $_[1] || "\t";
	my $has_header = $_[2];
	my $options = uc ($_[3] || '');
	
	$has_header = 1 if !defined $has_header;
	
	my %undefined;
	my $warn = 0;

	$delimiter = lc $delimiter;
	my $split = "\t" if ($delimiter eq 'tab'   || $delimiter eq "\t") ;
	$split = '\s+'   if ($delimiter eq 'space' || $delimiter eq " ");
	$split = ","     if ($delimiter eq 'comma' || $delimiter eq ",");
	
	eilico_msg('d', "Could not use delimiter '$delimiter'\n") if !defined $split;
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;
	
	if (!open_file(*INFILE, $infile)) {
		silico_msg('e', "Can not open file $infile for reading!\n");
		return undef;
	}
	
	if ($options !~ /\bQUIET\b/) {
		silico_msg('c', "Reading text file: $_[0] using '$delimiter' delimiter");
		silico_msg('c', " using $options") if $options;
		silico_msg('c',  "\n");
	}
	my $convert = 1 if ($options =~ /\bCONVERT_UNDEF\b/);
	$has_header = 0 if ($options =~ /\bNOHEADER\b/);
	
	my $ignore_hash = 1 if ($options =~ /\bIGNORE_HASH_LINES\b/); # Ignore lines starting with #
	my $ignore_at = 1 if ($options =~ /\bIGNORE_AT_LINES\b/);# Ignore lines starting with @
	
	my $data = new_table_delimited();	
	$data->{READ_DELIMITER} = $delimiter;
	$data->{SOURCE_FILE_NAME} = $infile;

	#molecule_printout($data);

	my $linecount = 0; # Numbers lines from 1
	my $numfields;
	while (my $line = <INFILE>) {
	
		my $f;
	
		chomp $line;
		
		$line =~ s/^\s*//g; # Remove leading whitespace
		
		next if $ignore_hash && $line =~ /^#/;
		next if $ignore_at && $line =~ /^@/;
		
		#print "split '$split' line: $line\n";
		
		if ($delimiter eq ',') {
			@$f = parse_csv_line($line);
		} else {
			@$f = split $delimiter, $line;	
		}

		#print "f ".(join "++", @$f)."\n";
		
		# Header line
		
		if ($linecount == 0) {

			$numfields = $#{$f};
			if ($Silico::debug) {
				print "Read ".($numfields+1)." fields from header\n";
				print "@$f\n";
			}

			if ($has_header) {
		
				$f->[0] =~ s/^#\s*//; # Remove leading hash if present
				$data->{HEADERS} = $f;
				# Remove leading and trailing whitespace
				foreach (@{$data->{HEADERS}}) {
					$line = uc $line if $options =~ /\bCASE_INSENSITIVE\b/;
					$line =~ s/^\s*//;
					$line =~ s/\s*$//;
				}
				++$linecount;
				next;

			} else {
				@{$data->{HEADERS}} = 1..($#{$f}+1);
				silico_msg('w', "File does not contain a header line.  Using column numbers for headers.\n");
				++$linecount;
			}
		}
		
		# Handle undefined values which are encoded as ''
		
		foreach my $j (0 .. $numfields) {
		
			if ($convert) {
				$line = '' if !defined $line;
			} else {
			
				if ($line eq '') {
					$line = undef;
					++$undefined{$j};
				}
			}
		}
		
		if ($#{$f} != $numfields) {
			silico_msg('w', "Line $linecount, Read $#{$f} fields but $numfields headers were read\n") if $warn == 0;
			++$warn;
		}
		
		@{$data->{VALUES}[$linecount-1]} = @$f;
		
		++$linecount;
	}
	
	#print "Read ".$linecount." lines of data\n" if $Silico::debug;
	
	silico_msg('w', "Above warning repeated $warn times\n\n") if $warn > 1;
	
	$data->{NUMLINES} = $linecount;
	$data->{NUMKEYS} = $#{$data->{HEADERS}} + 1;

	my $i = 0;
	my $keyhash;
	foreach my $key (@{$data->{HEADERS}}) {

		$keyhash->{$key} = $i;
		++$i;
	}

	$data->{KEYHASH} = $keyhash;

	close INFILE;
	
	foreach (keys %undefined) {
	
		silico_msg('w', "Read $undefined{$_} undefined values in column $_\n");
	}
	
	return $data;
}

sub parse_csv_line {

	# Source: https://www.oreilly.com/library/view/perl-cookbook/1565922433/ch01s16.html

    my $text = shift;      # record containing comma-separated values
    my @new  = ();
    push(@new, $+) while $text =~ m{
        # the first part groups the phrase inside the quotes.
        # see explanation of this pattern in MRE
        "([^\"\\]*(?:\\.[^\"\\]*)*)",?
           |  ([^,]+),?
           | ,
       }gx;
       push(@new, undef) if substr($text, -1,1) eq ',';
       return @new;      # list of values that were comma-separated
}

sub new_table_delimited {

	#<
	#? Create a new empty delimited table 
	#; Requires: nothing
	#; Returns: table
	#>

	my $table;

        $table->{VALUES}[0][0] = undef;
        @{$table->{HEADERS}} = ();
        $table->{NUMLINES} = 0;
        $table->{NUMKEYS}  = 0;
        $table->{KEYHASH} = {};
	$table->{READ_DELIMITER} = undef;
	$table->{SOURCE_FILE_NAME} = undef;

	return $table;
}

sub sdf_data_to_table_delimited {

	#<
	#? Copy ensemble sdf_data to table format
	#; Requires: Ensemble
	#; Returns: Table data structure
	#>

	my $mols = $_[0];

	my $columns;
	my $key;
	my $keys;
	my $mol;
	my $mol0;
	my $table = new_table_delimited;

	my $keynum = 0;
	my $row = 0;
	my $col;
	foreach $mol (@$mols) {

		foreach $key (sort keys %{$mol->{SDF_DATA}}) {

			# Add new keys if they are not already present
			if (!defined $keys->{$key}) {
				$columns->[$keynum] = $key;
				$keys->{$key} = $keynum;
				++$keynum;
			}

			$col = $keys->{$key};
			#print "key $key col $col\n";
			$table->{VALUES}[$row][$col] = $mol->{SDF_DATA}{$key};
		}
		++$row;
	}

	$table->{NUMLINES} = $row;
        $table->{NUMKEYS}  = $keynum;
	$table->{HEADERS} = $columns;
	$table->{KEYHASH} = $keys;
	$table->{READ_DELIMITER} = undef;
	$table->{SOURCE_FILE_NAME} = $mols->[0]{SOURCE_FILE_NAME};

	return $table;

}

sub table_delimited_inc_row {

	#<
	#? Return a row from a delimited data structure
	#; Requires: Data structure
	#; Returns: array or undef if finished
	#. Increments $data->{CURRENT_LINE} to keep track of row.  Reset $data->{CURRENT_LINE}
	#  to undef to start from the beginning again
	#>
	
	my $data = $_[0];
	
	if (!defined $data->{CURRENT_LINE}) {
		$data->{CURRENT_LINE} = 0;
	} else {
	
		++$data->{CURRENT_LINE};
	}
	
	return undef if $data->{CURRENT_LINE} > $#{$data->{VALUES}};
	
	return $data->{VALUES}[$data->{CURRENT_LINE}];
	
}

sub table_delimited_getrow {

	#<
	#? Return a row from a delimited data structure
	#; Requires: Data structure, array of column names
	#; Returns: array or undef if finished
	#. Increments $data->{CURRENT_LINE} to keep track of row.  Reset $data->{CURRENT_LINE}
	#  to undef to start from the beginning again
	#>
	
	my $data = shift;
	my @array;
	my $colnum;
	my $val;
		
	foreach (@_) {
	
		$colnum = $data->{KEYHASH}{$_};
		
		if (defined $colnum) {
			$val = $data->{VALUES}[$data->{CURRENT_LINE} || 0][$colnum];
			$val = '' if !defined $val;
			#print "1 colnum $colnum colname '$_' val '$val'\n";
			push @array, $val;
		} else {
		
			#print "2 colname '$_'\n";
			push @array, '';
		}
	}
	
	return @array;
}

sub table_delimited_getval {

	#<
	#? Return a value from the current line
	#; Requires: Data structure, column name
	#; Returns: value
	#. Uses $data->{CURRENT_LINE} to keep track of row
	#>
	
	my $data = $_[0];
	my $colname = $_[1];
	
	die if !defined $data || !defined $colname;
	
	if (!defined $data->{KEYHASH}{$colname}) {
	
		silico_msg ('w', "Column name not found: $colname\n");
		return undef;
	}
	
	return $data->{VALUES}[$data->{CURRENT_LINE} || 0][$data->{KEYHASH}{$colname}];
	
}


sub table_delimited_addcolumn {

	#<
	#? Add a new column to a 'delimited' data structure
	#; Requires: Data structure, new key name
	#; Returns: new column number (zero based)
	#>

	my $data = $_[0];
	my $newkey = $_[1];

	my $col = $#{$data->{HEADERS}} +1;

	# Check that key does not already exist
	if (defined $data->{KEYHASH}{$newkey}) {
	
		silico_msg('d', "Column already exists in data structure!\n");
	}

	$data->{HEADERS}[$col] = $newkey;
	$data->{KEYHASH}{$newkey} = $col;

	$data->{NUMKEYS} = $col+1;

	return $col;

}

sub write_table_delimited {

	#<
	#? Write a character delimited text file
	#; Requires: silico tabulated data structure, output filename,
	#  delimiter (default tab), options, list of keys to write out (default all).
	#  Missing values are written as a pair of single quotes ('')
	#; Returns: 1 if OK or undef if file open failed
	#. Options: None yet
	#; QUIET Do not print advisory messages
	#>

	my $data = shift;
	my $outfile = shift;
	my $delimiter = shift;
	my $options = shift || '';
	my @keys = @_;
	
	my $col;
	my @cols;
	my $key;
	my $line;
	my @pline;
	
	$delimiter ||= "\t";
	
	add_option("QUIET", $options) if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	silico_msg('c', "Writing $outfile\n") if ($options !~ /\bQUIET\b/);

	if (!open (OUT, ">$outfile")) {
		silico_msg('e', "Can not create or open file $outfile for writing\n");
		return undef;
	}
	
	# If no keys specified then write out all
	# otherwise get columns for each key
	if (!defined $keys[0]) {
		@keys = @{$data->{HEADERS}};
	}

	@cols = getcols($data, @keys);

	print heading("Headers\n");
	my $i = 1;
	foreach (@keys) {
		print "'$_' ";
		print "\n" if $i == 10;
		$i = 0 if $i == 10;
	}
	print "\n" if $i != 10;

	# Write header
	print OUT "#";
	print OUT join ($delimiter, @keys);
	print OUT "\n";

	foreach $line (@{$data->{VALUES}})  {

		@pline = ();

		#print "cols @cols\n";
		foreach $col (@cols) {
				
			if (!defined  $line->[$col]) {
				push @pline,"''";
			} else {
				push @pline, $line->[$col];
			}
		}
		
		if (!defined $pline[0]) {
			silico_msg('w', "No keys selected\n");
			return undef;
		}
		
		print OUT join ($delimiter, @pline);
		print OUT "\n";
	}

	close OUT;
	
	return 1;
}

sub smooth_table_delimited {

	#<
	#? Create running average of data in a delimited file over a window
	#. Window is the number of values on either side to average over.
	#  i.e. a window value of 1 will average over 3 values.
	#>

	my $data = $_[0];
	my $window = $_[1];
	
	my $line;
	my $numlines;
	my ($i, $j, $k);
	
	my $points = $window*2+1;

	silico_msg('c', "\nRunning average of data over window of $points points (2 * smooth + 1)\n\n");
	
	$numlines = $#{$data->{VALUES}};
	
	
	for ($i = $window; $i <= $numlines - $window; ++$i) {
	
		$line = $data->{VALUES}[$i];
		
		$j = 0;
		foreach (@$line) {
		
			my $c = 0;
			for ($k = $i-$window; $k <= $i+$window; ++$k) {
			
				next if !defined $data->{VALUES}[$k][$j];
				++$c;
				$data->{NEWLINES}[$i-$window][$j] += $data->{VALUES}[$k][$j];
			}
			
			$data->{NEWLINES}[$i-$window][$j] = $data->{NEWLINES}[$i-$window][$j]/$c;
			
			++$j;
		}
	}
			
	$data->{VALUES} = $data->{NEWLINES};

}

sub getlines {

	#<
	#? Return array of all lines in text data structure
	#; Requires: data structure
	#>
	
	print "This subroutine (getlines) is deprecated. Pleae use get_rows_all\n";
	
	return getrows_all(@_);
}

sub getcol {

	#<
	#? Find a single column number that corresponds to a particular key
	#. Case insensitive match
	#; Requires: data structure, list of keys
	#; Returns: array of column numbers
	#>

	my $data = $_[0];
	my $key = $_[1];
	
	my $i;
	my $col;
	my $dkey;
	
	$i = -1;
	foreach $dkey (@{$data->{HEADERS}}) {
		
		++$i;
		next if !(uc $key eq uc $dkey);
		$col = $i;
		last;
	}
	

	return $col
}

sub getcols {

	#<
	#? Find column numbers that correspond to a particular header
	#. Case insensitive match
	#; Requires: data structure, list of keys
	#; Returns: array of column numbers
	#>

	my $data = shift;
	my @keylist = @_;
	
	my $i;
	my @cols;
	my $dkey;
	my $key;

	foreach $key (@keylist) {
	
		$i = 0;
		foreach $dkey (@{$data->{HEADERS}}) {
		
			push @cols, $i if (uc $key eq uc $dkey);
			++$i;
		}
	}

	return @cols
}

sub getcols_not {

	#<
	#? Find column numbers except specified headers
	#. Inverse of getcol. Case insensitive match
	#; Requires: data structure, list of keys
	#; Returns: array of column numbers
	#>

	my $data = shift;
	my @keylist = @_;
	
	my @cols;
	my $dkey;
	my $i;
	my $key;

	$i = -1;
	DKEY: foreach $dkey (@{$data->{HEADERS}}) {
		
		++$i;
		foreach $key (@keylist) {
		
			next DKEY if (uc $key eq uc $dkey);
		}
	
		push @cols, $i;
	}

	return @cols;
}

sub getcols_all {

	#<
	#? Return column list of all keys sorted according to the order in @{$data->{HEADERS}}
	#; Requires: data structure
	#; Requires: pointer to array of column names
	#>
	
	my $data;
	
	my @cols;
	my $i;
	my $key;
	
	$i = 0;
	foreach $key (@{$data->{HEADERS}}) {
		push @cols, $i;
		++$i;
	}
	
	return @cols;
}

sub getrows_all {

	#<
	#? Get all rows as array of arrays. Values are in header order
	#; Requires: tabledata structure
	#; Returns: array of arrays
	#>	
	
	my $data = $_[0];

	return $data->{VALUES};
}

sub getheaders_all {

	#<
	#? Get all headers
	#; Requires: tabledata structure
	#; Returns: Pointer to arrray of headers
	#>
	
	my $data = $_[0];

	return $data->{HEADERS};
}



# Subroutines to handle binned data
# Used by program 'tabulate'

sub make_bins {

	#<
	#? Create a bin data structure
	#. Bin data structure contains the following fields:
	#, BINSIZE      - bin size
	#, RANGEMIN     - minimium bin value
	#, RANGEMAX     - maximum bin value
	#, UNDERFLOW    - count of values less than rangemin
	#, OVERFLOW     - count of values greater than rangemin
	#, COUNTS[i] - the number of values found for each key value for bin i
	#, NUMBINS      - number of bins
	#, NOT_NUMERIC  - count of non-numeric values found
	#. Requires: Binsize, range minimium, range maximum
	#. Returns: pointer to bin datastructure
	#>

	my $bins;

	$bins->{BINSIZE} = $_[0] || croak();	# Bin size
	$bins->{RANGEMIN} = $_[1];		# Bin minimum value
	$bins->{RANGEMAX} = $_[2];	  	# Bin maximum value
	
	croak() if !defined $bins->{RANGEMIN};
	croak() if !defined $bins->{RANGEMAX};

	$bins->{NUMBINS} = int (($bins->{RANGEMAX} - $bins->{RANGEMIN})/$bins->{BINSIZE} +0.5);

	$bins->{HEADERS} = {};

	return $bins;
}

sub bin_values {

	#<
	#? Bin values into 'bins' data structure
	#. Bin minimum and maximum values should be set in $bins->{$key}{MINVAL} and $bins->{$key}{MINVAL}
	#. Counts of overflows and underflows are kept in $bins->{$key}{OVERFLOW} and $bins->{$key}{UNDERLFLOW}
	#. Non-numeric data is not counted
	#. Requires: data value, key, pointer to bins data structure
	#. Returns: nothing
	#>

	my $val = $_[0];	# Data value to be 'binned'
	my $key = $_[1];	# Key to denote data set
	my $bins = $_[2];	# Pointer to bins data structure

	my $binval;

	$bins->{HEADERS}{$key} = 1;

	++$bins->{$key}{TOTAL_COUNTS};

	# Non-numeric or empty data
	if ($val !~ /^[\d.eE+-]*$/ || $val eq '') {

		++$bins->{$key}{NOT_NUMERIC};
		#print "nonumeric val $val key $key \n";
		return;
	}
	
	silico_msg('d', "Bin RANGEMIN not defined") if !defined $bins->{RANGEMIN};
	silico_msg('d', "Bin BINSIZE not defined") if !defined $bins->{BINSIZE};

	$binval = sprintf "%d", (($val - $bins->{RANGEMIN})/$bins->{BINSIZE})+0.5;

	$bins->{$key}{MINVAL} = $val if !defined $bins->{$key}{MINVAL};
	$bins->{$key}{MAXVAL} = $val if !defined $bins->{$key}{MAXVAL};

	$bins->{$key}{MAXVAL} = $val if $val > $bins->{$key}{MAXVAL};
	$bins->{$key}{MINVAL} = $val if $val < $bins->{$key}{MINVAL};

	if ($binval > $bins->{NUMBINS}) {
		++$bins->{$key}{OVERFLOW};
		return;
	}
	if ($binval < 0) {
		++$bins->{$key}{UNDERFLOW};
		return;
	}

	++$bins->{$key}{COUNTS}[$binval];

	return undef;
}


sub write_binfile {

	#<
	#? Write out a tab-delimited file of binned data
	#. Requires: output filename
	#. Returns: one or undef if file open failed.
	#>

	my $outfile = $_[0];
	my $bins = $_[1];
	
	my $i;
	my $keys;
	my $key;
	my $filebase;
	my $numbins;
	my $max = 1003;
	my $val;

	open (BINFILE, ">$outfile") or return undef;

	@$keys = sort keys %{$bins->{HEADERS}};
	silico_msg('c', "Writing binned data using keys: @$keys\n");

	# Print header
	print BINFILE "#".(join "\t", ("Count", @$keys))."\n";

	for ($i = 0; $i <= $bins->{NUMBINS}; ++$i) {

		print BINFILE  ($i*$bins->{BINSIZE}+$bins->{RANGEMIN})."\t";

		foreach $key (@$keys) {

			$val = $bins->{$key}{COUNTS}[$i] || "00";

			print BINFILE $val."\t";

		}

		print BINFILE "\n";
	}

	close BINFILE;

	return 1;
}



return 1;
