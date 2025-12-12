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
#! silico_gnuplot.pm
#? Silico wrapper routines to run gnuplot
#. $Revision: 1.1.2.2.2.1 $
#>

use strict;
package Silico;



##################################################################
#
#	Gnuplot routines
#
##################################################################

sub write_gnuplot_datafile {

	my $dat = $_[0];
	my $outfile = $_[1];
	
	my $headers = $dat->[0];
	my $values_array = $dat->[1];
	
	my $cw = 14;

  	open (DATFILE, ">$outfile") || file_write_error($outfile, 1);
	
	# Clean up headers
	foreach (@$headers) {
		$_ =~ s/^ *//;
		$_ =~ s/ *$//;
		$_ =~ s/ /--/g;
	}
			
	# Print headers
	foreach (@$headers) {
		printf DATFILE "%-".$cw."s  ", $_;
	}
	print DATFILE "\n";
				
	foreach my $values (@$values_array) {
		foreach (@$values) {
			printf DATFILE "%".$cw."f  ", $_;
		}
	
		print DATFILE "\n";
	}
			
	close DATFILE;		
}



sub run_gnuplot1 {

	#<
	#? Plots multiple data columns on separate graphs
	#; Requires: Gnuplot tab delimited data, optional title
	#; Returns: nothing
	#>
	
	my $datafile = $_[0]; # File with white space delimited data
	my $title = $_[1];
	my $first = $_[2] || get_flag('first-column', 'l'); # First column (numbered from 1)
	my $last  = $_[3] || get_flag('last-column', 'l'); # Last column (numbered from 1)
	
	my $aspectstring="portrait";
	my $colorstring="color";
	my $format = get_flag('plot-format', 'l') || 'postscript';
	my $gnuplot = get_flag('gnuplot-exe', 'l') || croak();
	my $keystring="unset key";
	my $line;
	my $time_unit = get_flag('time-unit', 'l') || croak();
	my $ratio=1.0;
	my $xinterval;

	my $out_ext= 'ps';
	$out_ext = $format if $format ne 'postscript';
	my $date = `date`;
	chomp $date;
	
	my $comfile = get_filebase($datafile).".inp";
	my $outfile = get_filebase($datafile).".".$out_ext;
	
	if (!defined $title) {
		$title  = "$outfile $date";
	}
	
	# Produce command file for gnuplot
	open (GPFILE, ">$comfile") || file_write_error($comfile, 1);

	print GPFILE "set size square\n";
	print GPFILE "set key\n";
	print GPFILE "set title '$title' noenhanced\n";
	print GPFILE "$keystring\n";

        if ($format eq 'postscript') {
		print GPFILE "set output \"$outfile\"\n";
                print GPFILE "set term postscript $aspectstring $colorstring \"Helvetica,12\"\n";

        } else {
                print GPFILE "set term $format\n";
        }
	
	if (!open (DATFILE, $datafile)) {
		file_read_error($datafile);
		return undef;
	}
	
	my $i = 0;
	my $start;
	my $finish;
	my @headers;
	while (my $line = <DATFILE>) {
		
		++$i;
		
		if ($i == 1) {
			@headers = split ' ', $line;
		} elsif ($i == 2) {
			my @a = split ' ', $line;
			$start = $a[0];
		}
		
		my @b = split ' ', $line;
		$finish = $b[0];
	}
		
	close DATFILE;
		
	# X values span
	my $xspan = abs($finish - $start);
	if ($xspan == 0) {
		silico_msg('w', "X range is zero\n");
		$xspan = 10;
	}
	# Get the common log of $xspan
	my $lxspan = log($xspan)/log(10);
	# Round down to an integer and subtract 1
	my $lxspan2 = int($lxspan) - 1;
	# Get the corresponding power of 10
	my $xspan2 = 10**($lxspan2);
	# Find how many times $xspan2 goes into $xspan
	my $xnumint = $xspan/$xspan2;
	
	# There is no way $xnumint will be more than 100.
	# Set various sizes of interval based on how many times
	# $xspan2 is needed to complete $xspan.
	if ($xnumint <= 6) {
		$xinterval = $xspan2;
	} elsif ($xnumint <= 12) {
		$xinterval = 2*$xspan2;
	} elsif ($xnumint <= 30) {
		$xinterval = 5*$xspan2;
	} elsif ($xnumint <= 60) {
		$xinterval = 10*$xspan2;
	} else {
		$xinterval = 20*$xspan2;
	}
	
	my $headcount = 0;
	foreach my $header (@headers) {

		$header =~ s/--/ /g; # Reinstate spaces

		++$headcount;
		next if ($first && $headcount < $first) || ($last && $headcount > $last);

		if($format ne 'postscript')  {
			$outfile = get_filebase($datafile)."_".(sprintf "%02d",$headcount).".".$out_ext;
			print GPFILE "set output \"$outfile\"\n";
		}
	
		print GPFILE "\n";
		print GPFILE "set xlabel \"Time ($time_unit)\"\n";
		print GPFILE "set xtics $xinterval\n";
		print GPFILE "set ylabel \'$header\'\n";
		print GPFILE "plot \"$datafile\" using 1:$headcount w lines\n";
	}
	close GPFILE;
	
	# Run gnuplot	
	silico_msg('c', "Running gnuplot\n");

	my $failure = system "$gnuplot $comfile";
	if ($failure) {
		silico_msg('e', "Execution of gnuplot failed!\n",
				"No graph created.\n");
	}

	if (!get_flag('noclean', 'l')) {
		silico_msg('c', "Removing files: $comfile $datafile\n");
		system "rm $comfile $datafile";
	}
}

sub run_gnuplot2 {

	#<
	#? Plots multiple data columns on a single graph
	#; Requires: Gnuplot tab delimited data, optional title
	#; Returns: nothing
	#>
		
	my $datafile = $_[0]; # File with white space delimited data
	my $title = $_[1];
	my $first = $_[2] || get_flag('first-column', 'l'); # First column (numbered from 1)
	my $last  = $_[3] || get_flag('last-column', 'l'); # Last column (numbered from 1)
	
	my $aspectstring="portrait";;
	my $colorstring="color";
	my $format = get_flag('plot-format', 'l') || 'postscript';
	my $gnuplot = get_flag('gnuplot-exe', 'l') || croak();
	my $keystring="set key";
	my $out_ext = 'ps';
	my $time_unit = get_flag('time-unit', 'l') || croak();
	
	my $comfile = get_filebase($datafile).".inp";
	my $outfile = get_filebase($datafile).".".$out_ext;
	
	# Produce command file for gnuplot
	open (GPFILE, ">$comfile") || file_write_error($comfile, 1);
 
 	if ($format eq 'postscript') {

		$outfile = get_filebase($datafile).".".$out_ext;
		print GPFILE "set output \"$outfile\"\n";
                print GPFILE "set term postscript $aspectstring $colorstring \"Helvetica,12\"\n";

        } else {
                print GPFILE "set term $format\n";
        }
	
	my $date = `date`;
	chomp $date;
	$title = $title || "$outfile $date";

	print GPFILE "set output \"$outfile\"\n";
	print GPFILE "set size square\n";
	print GPFILE "set key\n";
	print GPFILE "set title '$title' noenhanced\n";
	print GPFILE "$keystring\n";

	if ($format eq 'postscript') {
		print GPFILE "set term postscript $aspectstring $colorstring \"Helvetica,12\"\n";
	} else {
		print GPFILE "set term $format\n";
	}
	
	if (!open (DATFILE, $datafile)) {
		file_read_error($datafile);
		return undef;
	}
	
	my $line = <DATFILE>;
	my @f = split ' ',$line;
	close DATFILE;
	
	print GPFILE "set xlabel \"Time ($time_unit)\"\n";
	print GPFILE "plot ";
	
	my $headcount = 1;
	my $count=0;
	foreach my $header (@f[1..$#f]) {

		$header =~ s/--/ /g; # Reinstate spaces

		next if ($first && $headcount < $first) || ($last && $headcount > $last);
		++$headcount;

		if($format ne 'postscript')  {

			$outfile = get_filebase($datafile)."_".(sprintf "%02d",$headcount).".".$out_ext;
			print GPFILE "set output \"$outfile\"\n";
		}
	
		++$count;
		if ($headcount > 2) {
			print GPFILE ",";
		}
		
		print GPFILE "\"$datafile\" using 1 :$headcount title '$header' w lines";
		
	}
	print GPFILE "\n";
	
	close GPFILE;
	
	# Run gnuplot	
	silico_msg('c', "Running gnuplot\n");
	my $failure = system "$gnuplot $comfile";
	if ($failure) {
		silico_msg('e', "Execution of gnuplot failed!\n",
				"No graphs created.\n");
	}

	if (!get_flag('noclean', 'l')) {
		silico_msg('c', "Removing files: $comfile $datafile\n");
		system "rm $comfile $datafile";	
	}
}

sub gnuplot_plot_bins {

	#<
	#? Use gnuplot to plot histograms from bin datastructures
	#. Calls write_binfile
	#. Requires: output filename base,bins structure
	#. Returns: one or undef if file open failed.
	#>

	my $filebase = $_[0];
	my $bins = $_[1];
	
	my $binfile = $filebase.".bins";
	my $comfile = $filebase.".inp";
	my $outfile = $filebase.".ps";
	
	my $ratio=0.773;
	my $scale=0.8;
	my $xlabel;
	my $ylabel;
	
	my $keys;
	@$keys = sort keys %{$bins->{KEYS}};
	
	write_binfile($binfile, $bins);
	
	# Write Gnuplot input file
	open (GP, ">$comfile") || return undef;
	
	print GP "set term postscript portrait color\n";
	print GP "set output '$outfile'\n";
	
	print GP "set size ".($ratio/$scale).", $scale\n";
	
	my $i = -1;
	foreach my $key (@$keys) {
	
		++$i;
		
		my $step;
		my $xmax;
		my $xmin;
	
		# Skip plotting this one if no values were in range (ie maxval never set)
		if (!$bins->{$key}{MAXVAL}) {
			silico_msg('c', "No values for $key at bin size $bins->{BINSIZE}\n");
			next;
		}
		
		silico_msg('g', "key $key tc $bins->{$key}{TOTAL_COUNTS} rmax $bins->{RANGEMAX} rmin $bins->{RANGEMIN} maxval $bins->{$key}{MAXVAL} minval $bins->{$key}{MINVAL}\n");
				
		# Skip plotting if more than 50% of values were out of range
		if ((($bins->{$key}{UNDERFLOW} || 0) + ($bins->{$key}{OVERFLOW} ||0) + ($bins->{$key}{NOT_NUMERIC} || 0)) > $bins->{$key}{TOTAL_COUNTS}/2) {
			silico_msg('c', "More than 50% of values are outside range.  Skipping plot for $key at bin size $bins->{BINSIZE}\n");
			next;
		}
	
		$xlabel = $key || '';
		$ylabel = 'Count';
		print GP "set xlabel '$xlabel'\n";
		print GP "set ylabel '$ylabel'\n";
		
		$step = $bins->{BINSIZE} *5;

		$xmax = (int($bins->{RANGEMAX}));
		$xmin = (int($bins->{RANGEMIN}));
		
		while ($bins->{$key}{MAXVAL} < $xmax-$step && $xmax > $step) {
			#$xmax = $xmax/2;
			$xmax -= $step;
		}
		while ($bins->{$key}{MINVAL} > $xmin+$step && $xmin < $xmax-(2*$step)) {
			#$xmin = $xmin/2;
			$xmin += $step;
		}
		
		#$xmin = (int($bins->{RANGEMIN}));
		if ($xmax == $xmin) {
			$xmax +=5;
			$xmin -=5;
		}
		
		$xmin = 0 if ($xmin > 0 && $xmin < $xmax/2);
		
		print GP "plot [$xmin:$xmax] '$binfile' using 1:".($i+2)." w boxes\n";
	}
	
	close GP;
	
	my $failure = system ("gnuplot < $comfile");
	
	if ($failure) {
		silico_msg('e', "Execution of gnuplot failed!\n");
	}
	
	return 1;
}

return 1;
