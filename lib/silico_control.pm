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
#! silico_control.pm
#? Various subroutines for program control, user interface, etc.
#. $Revision: 1.15-75-gd957e88 $
#>

use strict;
package Silico;

sub get_hosttype {

	#<
	#? Find the host type
	#: Requires: nothing
	#: Returns: sgi, linux, linuxem64, Darwin (OS X)  or undef
	#>

	my $uname = `uname`;

	return 'sgi' if ($uname =~ /IRIX/);
	if ($uname =~ /Linux/) {

		$uname = `uname -m`;

		 if ($uname =~ /i686/) {
			return 'linux'
		 } elsif ($uname =~ /x86_64/) {
			return 'linuxem64';
		} else {
			return undef;
		}
	}
	return 'Darwin' if ($uname =~ /Darwin/);
	return 'undef';
}

sub print_finished_message {

	#<
	#? Prints out a message to indicate a script is finished
	#; Requires: Optional message to print (default finished)
	#, optional character to surround it with (default *), filehandle (optional)
	#; Returns: Nothing
	#>

	return if get_flag('very-quiet', 'l');
	
	my $message = $_[0] || "FINISHED";
	my $char = $_[1] || "*";
	my $FH = $_[2] || *STDOUT;

	my $space = ' ';
	
	my $length;
	
	$message =~ s/\s*$//; # Remove any trailing spaces
	
	$length = length($message);
	
	# Display a note that the script is finished
	print "\n";
	print "\n";
	print "               ".$char x ($length + 14)."\n";
	print "               ".$char x 2 .$space x ($length + 10).$char x 2 ."\n";
	print "               ".$char x 2 ."     ".$message."     ".$char x 2 ."\n";
	print "               ".$char x 2 .$space x ($length + 10).$char x 2 ."\n";
	print "               ".$char x ($length + 14)."\n";
	print "\n";
	print "\n";
}

sub print_timing_message {

	#<
	#? Prints out a message containing the total user time
	#  and total wall time used by the program
	#; Requires: Starting user time, filehandle (optional)
	#; Returns: Nothing
	#>
	
	return if get_flag('very-quiet', 'l');

	my $start_usertime = $_[0];
	my $FH = $_[1] || *STDOUT;
	
	my $usertime;
	my $walltime;
	
	$usertime = calc_human_readable_time(times - $start_usertime) if (defined $start_usertime);
	$walltime = calc_human_readable_time(time  - $^T);
	
	print $FH "Total user time (to nearest 1/100s) in $Silico::PROG: $usertime\n" if (defined $usertime);
	print $FH "Total wall time (to nearest second) in $Silico::PROG: $walltime\n" if (defined $walltime);
	print $FH "\n\n" if (defined $usertime || defined $walltime);
}

sub print_title {
	
	#<
	#? Print a program title banner.
	#; Requires: program name, version (optional), filehandle (optional)
	#>

	# Problem is that flags may not be set when print_title is called
	return if get_flag('very-quiet', 'l');

	my $name = ucfirst ($_[0]);
	my $program_version = $_[1] || '';
	my $FH = $_[2] || *STDOUT;
	
	my $hashes = 40;
	my $longest;
	my $string1;
	my $string1_length;
	my $string2;
	my $string2_length;
	
	$program_version =~ s/\s+\$$//;
	$program_version =~ s/.*: //;
	
	$string1 = "$name $program_version";
	chomp $string1;
	$string2 = "Silico version $silico_version";
	chomp $string2;
	$string1_length = length($string1);
	$string2_length = length($string2);
	
	# Use of conditional operator to determine which of strings 1 and 2 is longest
	# and then assign the length of it to $longest
	$longest = ($string1_length >= $string2_length) ? $string1_length : $string2_length;
	
	$hashes = $longest + 8 if ($longest > 32);

	print $FH "\n";
	print $FH "          ".("#" x $hashes)."\n";
	print $FH "\n";
	print $FH "              $string1\n";
	print $FH "\n";
	print $FH "              $string2\n" if $string2;
	print $FH "              Host: ".ucfirst($Silico::Hostname)."\n";
	print $FH "\n";
	print $FH "          ".("#" x $hashes)."\n";
	print $FH "\n";
}

sub silico_msg {
	
	#<
	#? Print out a note, warning, error or fatal-error message.
	#. Types are c, comment etc, g debug etc etc 
	#; Requires: Type (single-character string), strings
	#; Returns: Nothing
	#>

	# Note about the caller function: goes back up the stack. caller(0) is the current subroutine,
	# caller(1) is the routine from which this was called, caller(2) is the routine from which caller(1)
	# was called, and so on.
	
	# What is in Caller?
	# [0] = package name to which the routine belongs
	# [1] = the filename from where the routine has been called
	#       (_not the file where the routine is declared)
	# [2] = the line in [1] from where the routine has been called
	# [3] = the routine's name
	# [4] = (Boolean) whether the routine was called with arguments
	# [5] = what type of variable is waiting for the results - is it
	#       an array, a scalar or nothing?
	
	if (!defined $_[0]) {
		silico_msg('d', "No message type selected!\n",
			"Message type should be one of:\n",
			make_types_summary());
	} 

	my $c = lc($_[0]);
	
	# Skip notes and comments if we are using "--quiet" or "--very-quiet"
	return if (($c eq 'n' || $c eq 'c') and $Silico::quiet || $Silico::veryquiet);
	
	# Skip warnings and 'quiet comments' if we are using "--very-quiet"
	return if (($c eq 'w' || $c eq 'q') && $Silico::veryquiet);
	
	#
	# Note: Do not insert any code before this point.  This will slow the silico_msg routine down
	# and adversely impact program performance (because it is called many times).
	#
	
	# Skip debugging print statements if not in debug mode
	return if ($c eq 'g' && !$Silico::debug);
	return if ($c eq 'r'  && !$Silico::debug_rings);
	return if ($c eq 't'  && !$Silico::TIMING);
	
	if ($c eq 'c' ||	$c eq 'd' ||	$c eq 'dq' ||	$c eq 'e' ||	$c eq 'g' ||
		$c eq 'h' ||	$c eq 'n' ||	$c eq 'q' ||	$c eq 'r' ||	$c eq 't' ||
		$c eq 'v' ||	$c eq 'w' ) {

		# If we are not printing a new line at the end, make sure output is to be
		# flushed immediately. Otherwise, flush any remaining output and rebuffer.
		if ($_[$#_] !~ /\n\z/) {
			$| = 1;
		} else {
			print "\n" if ($| && ($c eq 'n' || $c eq 'w' || $c eq 'e' || $c eq 'd' || $c eq 'h'));
			$| = 0;
		}
		
		# Print prefixes
		my $printstring = '';
		$printstring = "Note" if ($c eq 'n');
		$printstring = "Warning" if ($c eq 'w');
		$printstring = "Error" if ($c eq 'e');
		$printstring = "\nFatal error" if ($c eq 'd' || $c eq 'dq' || $c eq 'h');

		#  Name of the subroutine from which silico_msg was invoked
		#  We don't want to use the subroutine if it is empty
		my $subroutine = (caller(1))[3] || '';
		$subroutine =~ s/Silico:://g;    # Get rid of the initial package name
		
		# Only print out the subroutine name if the subroutine name is
		# defined (i.e., not the main sub) and the message is something
		# other than a simple Note.
		if ($c eq 'n' || $c eq 'w' || $c eq 'e' || $c eq 'd' || $c eq 'dq' || $c eq 'h') {
			$printstring .= " ($subroutine)" if ($subroutine && $c ne 'n');
			$printstring .= ": " if defined $_[1];
		}

		# Print out the text of the message
		if ($#_ > 0) {
			for (my $i = 1; $i <= $#_; ++$i) {
				$printstring .= "$_[$i]";
			}
		}
		
		print $printstring;

		# Die
		if ($c eq 'd' || $c eq 'dq') {

			if ($c ne 'dq') {

				print "\n";
				use Carp qw(cluck confess);
				carp();

				# If we are debugging, print out the entire call stack
				if ($Silico::debug) {
					print "\n";
					confess("Backtrace");
				}
			}

			$| = 1;
			print "\n$Silico::PROG will now exit.";
			$| = 0;
			
			# Reset the text colour if we have changed it
			use Term::ANSIColor;
			print color('reset');
			die "\n";
		}
		
		if ($c eq 'h') {
			printdoc();
			die "\n";
		}
		return;
	}

	silico_msg('w', "Unrecognised message type '$_[0]'!\n",
		"Message type should be one of:\n",
		make_types_summary(),
		"Message will be printed as a comment (type \'c\').\n");

	silico_msg('c', @_[1..$#_]);
}

sub make_types_summary {

	#<
        #? Make a silico_msg type summary string
        #; Requires: Types hash
        #; Returns: string
        #>
	
	my %types;
	my $types_summary;
	
	# Types which print out a prefix
	$types{'n'} = "Note";
	$types{'w'} = "Warning";
	$types{'e'} = "Error";
	$types{'d'} = "Die";
	$types{'dq'} = "Die quietly";
	$types{'h'} = "Print documentation";

	# Types which don't print out a prefix
	$types{'g'}  = "Debugging comment";
	$types{'c'}  = "Comment";
	$types{'q'}  = "Quiet comment";
	$types{'v'}  = "Very quiet comment";
	$types{'r'}  = "Debug-rings comment";
	$types{'t'}  = "Timing comment";
		
	# Prepare a list (text string) containing the types
	foreach (sort keys %types) {
		$types_summary .= "'";
		$types_summary .= "$_";
		$types_summary .= "' ";
		$types_summary .= "(";
		$types_summary .= lcfirst $types{$_};
		$types_summary .= ")";
		$types_summary .= "\n";
	}

	return $types_summary;
}

sub prompt {
	
	#<
	#? Print a prompting message and read a reply from standard input.
	#; Requires: Text to prompt with
	#; Returns: Chomped reply from standard input
	#>
	
	my $i;
	my $response;
	
	if ($#_ > 0) {
		for ($i = 0; $i < $#_; ++$i) {
			print "$_[$i]";
		}
	}
	
	print "$_[$#_]";
	
	$response = <STDIN>;
	$response = '' if !defined $response;
	chomp $response;
	
	return $response;
}

sub sleep_or_stop {
	
	#<
	#? Check for sleep and stop files
	#. Waits while a sleep file is present in the working directory
	#. Returns 1 if stop file found
	#>
	
	my $sleepcount = 0;
	
	if (-e 'stop') {
		print "Found 'stop' file.\n";
		return 1;
	}
	while (-e 'sleep') {
		print "Sleeping\n" if $sleepcount == 0;
		sleep 10;
		++$sleepcount;
		$sleepcount = 0 if $sleepcount == 10;
	}
	
	return 0;
}

sub subname {

	#<
	#? Return the name of the subroutine in which this was called
	#; Requires: Nothing
	#; Returns: Name
	#>
	
	my $subname;
	
	$subname = (caller(1))[3];
	$subname =~ s/^.*:://g;
	
	return $subname;
}


sub mol_warn {

	#<
	#? Keep track of errors in $mol->{WARN}
	#; Requires: Molecule, warning string
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	my $string = $_[1];
}


sub mol_print_warnings {
	
	#<
	#? Print warnings from $mol->{WARN}
	#; Requires: Molecule
	#; Returns: Nothing
	#>
	
	my $mol = $_[0];
	
	print "\n";
	foreach my $warn (keys(%{$mol->{WARN}})) {
			
		my $v = $mol->{WARN}{$warn};
		$warn .= ". This warning was repeated $v times" if $v > 1;

		chomp $warn;
	
		silico_msg('w', $warn."\n");
	}
	print "\n";
	
	# Remove printed errors
	delete $mol->{WARN};

}

return 1;

