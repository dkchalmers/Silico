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
#! silico_doc.pm
#? Silico routines to print formatted documentation in text and html formats
#  which is generated from program comments that have been marked up using a simple
#  syntax
#.
#. Simplified markup syntax:
#, #< = Start marked up comments
#, #> = End marked up comments
#, #! = Heading/Title  (Heading 2)
#, #!1, #!2, #!3, #!4 Different heading levels
#, #?  Short description of program/subroutine used for indexing.
#, #. = New paragraph.
#, #, = Indent
#, #; = Hanging indent (long).
#, #: = Hanging indent (short, for numbering).
#, #- = List.
#, #c = Code block
#.
#.  Standard headings:
#, #F = Flags
#, #SF = Silico Flags
#, #GF = Gromacs Flags
#.
#. $Revision: 1.36.2.4.2.7 $
#>

use strict;
package Silico;

sub printdoc {
	
	#<
	#? Print marked up documentation for a file 
	#  The pager 'more' is used if the parent program has been called from 
	#  a shell
	#; Requires: filename (default executing program);
	#>
	
	my $file = $_[0] || $0; # Default is executing program

	# Change the output so that it goes to 'more' if stdout is a terminal
	if (parent_proc_is_shell()) {            
    		my $pager = 'more';  
    		open(STDOUT, "| $pager")    or silico_msg ('w', "can't open a pager: $!");
	}

	END { 
		close(STDOUT);
	}
	
	my $header = read_doc($file, 1);
	$header = format_doc_txt($header);

	print "$header";

	return;
}

sub parent_proc_is_shell {

	#<
	#? Find if the parent process is a shell
	#  Shells: bash, csh, ksh, sh, tcsh and zsh
	#; Requires: nothing
	#; Returns: true or false

	my @shells = qw "bash csh ksh sh tcsh zsh\n";
	my $ppid = getppid(); # Parent process ID
	my $out = `ps -p $ppid`; 
	my @f = split "\n", $out;

	foreach (@f) {

		chomp;
		my @f = split " ", $_;
		next if $ppid ne $f[0];
		my $prog = $f[-1];
		$prog =~ s/^-//;

		foreach (@shells) {
			next if $prog ne $_;
			return 1;
		}
	}
	return 0
}

sub read_doc {

	#<
	#? Read marked up documentation from a silico file
	#; Requires: Text file marked up using silico convention, flag to exit
	#  after encountering first end markup marker ('#>')
	#; Returns: Tagged silico text
	#. Marked up text is also added to a global hash %$FILEHASH
	#>

	use vars qw($FILEHASH);

	my $file = $_[0];
	my $headeronly = $_[1] || 1;

	open (INFILE, $file) || return undef;

	my $flag;
	my $lines;

	my $i = 0;
	while (<INFILE>) {
		
		next if !/^\s*#/;
		
		if (/^\s*#</) {
			$flag = 1;
			next;
		}
		
		if (/^\s*#>/) {
			last if $headeronly;
			$flag = 0;
			next;
		}

		$lines .= $_ if $flag;
	}

	close INFILE;
	
	$lines = tag_markup($lines);

	$FILEHASH->{$file} = $lines;
	return $lines;
}

sub read_doc_subroutine {

	#<
	#? Read marked up documentation from silico subroutines
	#. Marked up text is returned in a global hash %$SUBHASH
	#; Requires: Text file marked up using silico convention
	#; Returns: All lines read
	#>

	use vars qw($SUBHASH);

	my $file = $_[0];

	my $alllines;
	my $flag;
	my $lines;
	my $subflag;
	my $subname;
	
	open (INFILE, $file) || return undef;

	my $i = 0;
	while (<INFILE>) {

		# Find subroutines
		if (/\s*sub\s.*{/) {

			$subflag = 1;
			$lines = '';

			chomp;
			s/^\s*sub\s*//;
			s/\s.*//;

			$subname = sprintf "%-30s %s", ucfirst($_), ucfirst($file);

			next;
		}

		# Skip if we are not within a subroutine
		next if !$subflag;

		# Skip if we are not within a silico comment
		next if !/^\s*#/;

		# Set flag if we are at the start of a silico comment
		if (/^\s*#</) {

			$flag = 1;
			next;
		}
		
		# Unset flag if we are at the end of a silico comment
		# markup lines and save in a hash by subroutine name
		if (/^\s*#>/) {

			$lines = tag_markup($lines);
			$SUBHASH->{$subname} = $lines;
			$flag = 0;
			$subflag = 0;
			$alllines .= $lines;
			$lines = '';
			next;
		}

		# Skip if flag not set
		next if !$flag;

		$lines .= $_;
	}

	close INFILE;

	return $alllines;
}

sub format_doc_txt {

	#<
	#? Format a tagged string into ascii text
	#>

	my $string = $_[0];
	my $linelength = $_[1] || 80;	# Total line length. Default 80
	my $indent = $_[2] || 8;# Indent for second and subsequent lines

	my @f;
	my $line;
	my $para = '';
	my $ptype = '';
	
	return "" if !$string;
	
	# Split in to paragraphs (on <STD>, <HEAD.> or <INDENT>)
	@f = split /(?=<STD>)|(?=<STD_DESC>)|(?=<HEAD.>)|(?=<INDENT>)|(?=<HANG>)|(?=<HANG_SHORT>)|(?=<LIST>)|(?=<CODE>)/, $string;

	foreach (@f) {
	
		if (get_flag('clean', 'l')) {
			next if $_ =~ />\s*Created:/;
			next if $_ =~ />\s*Revision/;
		}
		
		# Clean out HTML links
		s/<LINK#.*?>\s*(\S+)/$1/g; # Remove links to anchors
		s/<LINK(.*?)>\s*(\S+)/$2 ($1)/g; # place http addresses in brackets
		s/<ANCH.*?>//g; # Remove anchor tags

		# This all needs to be cleaned up
		if (/<HEAD1>/) {
			s/<HEAD1> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = heading($line, '=');
			$line = "\n".$line if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<HEAD2>/) {
			s/<HEAD2> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = heading($line, '-');
			$line = "\n".$line if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<HEAD3>/) {
			s/<HEAD3> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = heading($line, '~');
			$line = "\n".$line if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<HEAD4>/) {
			s/<HEAD4> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = heading($line, ".");
			$line = "\n".$line if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<STD>/) {
			s/<STD> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = "\n".$line if $ptype ne '';
			$ptype = 'STANDARD';
		}

		if (/<STD_DESC>/) {
			s/<STD_DESC> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = "\n".$line if $ptype ne '';
			$ptype = 'STANDARD';
		}

		if (/<HANG>/) {
		
			s/<HANG> *//;
			$line = formatpara_txt($_, $linelength, 0, int($indent*2));
			$line = "\n".$line if ($ptype ne 'INDENT' && $ptype ne 'HANG' && $ptype ne 'HANG_SHORT');
			$ptype = 'HANG';
		}
		
		if (/<HANG_SHORT>/) {
		
			s/<HANG_SHORT> *//;
			$line = formatpara_txt($_, $linelength, 0, 4);
			$line = "\n".$line if ($ptype ne 'INDENT' && $ptype ne 'HANG' && $ptype ne 'HANG_SHORT');
			$ptype = 'HANG_SHORT';
		}

		if (/<INDENT>/) {
			s/<INDENT> *//;
			$line = formatpara_txt($_, $linelength, $indent, 0);
			$line = "\n".$line if ($ptype ne 'INDENT' && $ptype ne 'HANG' && $ptype ne 'HANG_SHORT');
			$ptype = 'INDENT';
		}

		if (/<LIST>/) {
			s/<LIST> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = "\n".$line if ($ptype ne 'LIST');
			$ptype = 'LIST';
		}

		if (/<CODE>/) {
			# Do nothing - use STD
                        s/<CODE> *//;
                        $line = formatpara_txt($_, $linelength, 0);
                        $line = "\n".$line if $ptype ne '';
			print "line $line\n";
                        $ptype = 'CODE';
                }

		# No formatting type
		if ( !defined $line) {
			$line =  "$_\n";
			$ptype = 'EMPTY';
		}

		$para .= $line;
	}
	return $para;
}

sub format_doc_html_header_footer {

	#<
	#; Format a tagged document as html and add header and footer files
	#>

	my $para = $_[0];
	my $header = $_[1];
	my $footer = $_[2];
	my $hplaceholders = $_[3];
	my $fplaceholders = $_[4];
	
	my $h = read_textfile($header);
	foreach (@$hplaceholders) {
		$h =~ s/PLACEHOLDER/$_/;
	}
	
	my $f = read_textfile($footer);
	foreach (@$fplaceholders) {
		$f =~ s/PLACEHOLDER/$_/;
	}
	
	$para =  format_string_html ($para);
	
	my $out = $h.$para.$f;

	return $out;
}


sub format_doc_html {

	#<
	#; Returns: Silico markup formatted as text
	#>

	my $para = $_[0];
	
	$para = format_string_html ($para);
	$para = wraphtml($para);

	return $para;
}

sub format_doc_markdown {

	#<
	#? Format a tagged string into markdown
	#; Requires: string, linelength
	#; Returns: Silico markup formatted as markdown
	#>

	my $string = $_[0];
	my $linelength = $_[1] || 80;

	my $para = '';
	my $ptype = '';
	
	return "" if !$string;
	
	# Split in to paragraphs (on <STD>, <HEAD.> or <INDENT>)
	my @f = split /(?=<STD>)|(?=<STD_DESC>)|(?=<HEAD.>)|(?=<INDENT>)|(?=<HANG>)|(?=<HANG_SHORT>)|(?=<LIST>)|(?=<CODE>)/, $string;

	foreach (@f) {
	
		my $line;
		
		if (get_flag('clean', 'l')) {
			next if $_ =~ />\s*Created:/;
			next if $_ =~ />\s*Revision/;
		}

		my @w = split " ", $_;
		my $tag = $w[0];
		$tag =~ s/[<>]//g;
		
		# Clean out HTML links and convert to markdown
		s/<LINK#(\d+)>\s*(\S+)/[$2](#anchor$1)/g; # Convert anchor links
		s/<LINK(.*?)>\s*(\S+)/[$2]($1)/g; # Convert http links
		s/<ANCH(\d+)?>/<a name="anchor$1"><\/a>/g; # Keep anchor tags as HTML (markdown doesn't support anchors well)
		s/<FILE(.*?)>\s*(\S+)/[$2]($1)/g; # Convert file links

		# Headings
		if (/<HEAD1>/) {
			s/<HEAD1> *//;
			$line = formatpara_txt($_, $linelength, 0);
			chomp $line;
			$line = "# " . $line;
			$line = "\n".$line."\n" if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<HEAD2>/) {
			s/<HEAD2> *//;
			$line = formatpara_txt($_, $linelength, 0);
			chomp $line;
			$line = "## " . $line;
			$line = "\n".$line."\n" if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<HEAD3>/) {
			s/<HEAD3> *//;
			$line = formatpara_txt($_, $linelength, 0);
			chomp $line;
			$line = "### " . $line;
			$line = "\n".$line."\n" if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<HEAD4>/) {
			s/<HEAD4> *//;
			$line = formatpara_txt($_, $linelength, 0);
			chomp $line;
			$line = "#### " . $line;
			$line = "\n".$line."\n" if $ptype ne '';
			$ptype = 'HEADER';
		}

		if (/<STD>/) {
			s/<STD> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = "\n".$line if $ptype ne '';
			$ptype = 'STANDARD';
		}

		if (/<STD_DESC>/) {
			s/<STD_DESC> *//;
			$line = formatpara_txt($_, $linelength, 0);
			$line = "\n".$line if $ptype ne '';
			$ptype = 'STANDARD';
		}

		if (/<HANG>/) {
			$line = format_string_html($_);
			$line = "\n".$line if ($ptype ne 'INDENT' && $ptype ne 'HANG' && $ptype ne 'HANG_SHORT');
			$ptype = 'HANG';
		}
		
		if (/<HANG_SHORT>/) {
			$line = format_string_html($_);
			$line = "\n".$line if ($ptype ne 'INDENT' && $ptype ne 'HANG' && $ptype ne 'HANG_SHORT');
			$ptype = 'HANG_SHORT';
		}

		if (/<INDENT>/) {
			$line = format_string_html($_);
			$ptype = 'INDENT';
		}

		if (/<LIST>/) {
			$line = format_string_html($_);
			$ptype = 'LIST';
		}

		if (/<CODE>/) {
			s/<CODE> *//;
			$line = formatpara_txt($_,  $linelength, 0);
			chomp $line;
			$line = "`$line`\n\n";
			$ptype = 'CODE';
		}

		# No formatting type
		if ( !defined $line) {
			$line =  "$_\n";
			$ptype = 'EMPTY';
		}

		if ($tag ne $ptype) {

			chomp $line;
			chomp $line;
			$line .= "\n\n"; 
		}

		$para .= $line;
	}
	return $para;
}


sub tag_markup {

	#<
	#? Convert silico markup to tags
	#. Tags are:
	#, #! -> < HEAD# > 
	#, #> -> < STD > 
	#, #? -> < STD_DESC > 
	#, #; -> < HANG > 
	#, #: -> < HANG_SHORT >
	#, #, -> < INDENT >
	#, #- -> < LIST >
	#, #c -> < CODE >
	#. Boilerplate flags #F and #SF cause the standard text to be inserted
	#. #INCLUDE <filename> causes the file from $SILICO_HOME/data to be inserted
	#. (Note that extra spaces have been added to the tags above to
	#  stop them being recognised by makedoc).
	#. The first part of a hanging indent can use underscores to represent spaces
	#>

	my $lines = $_[0];

	my $current;
	my $digit;
	my @f;
	my $flag;
	my @g;
	my $olines = '';
	
	return '' if !$lines;

	@f = split "\n", $lines;

	# Interpolate Flags and Standard Flags
	foreach (@f) {
	
		# Text '#F' is shorthand for "#. Flags"
		if (/^\s*#F\s*$/) {
			push @g, "#. Flags:";
		}
		
		# Text '#SF' is shorthand for the Silico Standard Flags
		elsif (/^\s*#SF\s*$/) {
		
			push @g, "#. Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, -bond, -nobond, -noconect,";
			push @g, "#  -ms, -me, -ss, --last-structure, -longres";
			push @g, "#  --timing, --quiet, --very-quiet and -debug.  These are available in most scripts.";
			push @g, "#  Some special options apply to SDF and PDB files.";
			push @g, "#. For more information rerun this program using '--help-flags'";
			push @g, "#  or see the silico documentation: file:/$ENV{SILICO_HOME}/doc.html'" if -e "$ENV{SILICO_HOME}/doc.html";
			#print  "$ENV{SILICO_HOME}/doc.html\n";
		}

		# Text '#INCLUDE <filename>' inserts markup text from file in $SILICO_HOME/data/
		elsif (/^\s*#INCLUDE\s+/) {
			
			my $filename = $_;
			$filename =~ s/#INCLUDE\s*//;
			$filename =~ s/\s*$//;
			
			my $success = open (INFILE, "$ENV{SILICO_HOME}/data/$filename"); 

			if (!$success) {

				 silico_msg('w', "Could not open INCLUDE file $ENV{SILICO_HOME}/data/$filename\n");
			} else {

				my $i;
				while (<INFILE>) {

					# Ignore all text after #>
					$i = 0 if /^#>/;

					# Ignore all text before #<
					if (/^#</) {
						$i = 1;
						next;
					}

					next if !$i;
					push @g, $_;
				}
			}
		}
		
		else {
			chomp $_;
			push @g, $_;
		}
	}
	
	foreach (@g) {
		
		# Replace < and > with &lt; and &gt; respectively
		# except following a '#'
		s/(?=[^#])</&lt;/g;
		s/(?=[^#])>/&gt;/g;
			
		# Header
		if ( /^\s*#!/) {
			$current = 'HEAD';

			if (/^\s*#!(\d)/) {

				$digit = $1;
				s/^\s*#!./<HEAD$digit>/;

			} else {
				
				# Default heading 2
				s/^\s*#!/<HEAD2>/;
			}
		}

		# Paragraph
		if ( /^\s*#\.\s*/) {
			$current = 'STANDARD';
			s/^\s*#\./<STD>/;
		}

		# Short description line
		if ( /^\s*#\?\s*/) {
			$current = 'STANDARD_DESC';
			s/^\s*#\?/<STD_DESC>/;
		}

		# Indent
		if (/^\s*#,\s*/) {
			s/^\s*#,/<INDENT>/;
			$current = 'INDENT';
		}

		# Hanging indent
		if (/^\s*#;\s*/) {
			s/^\s*#;/<HANG>/;
			$current = 'HANG';
		}
			
		# Short hanging indent
		if (/^\s*#:\s*/) {
			s/^\s*#:/<HANG_SHORT>/;
			$current = 'HANG_SHORT';
		}
		
		# List (line with no blank line following)
		if (/^\s*#-\s*/) {
			s/^\s*#-/<LIST>/;
			$current = 'LIST';
		}

		# Code block
		if (/^\s*#c\s*/) {
			s/^\s*#c/<CODE>/;
			$current = 'CODE';
		}
		# Clean up output #
		
		# Remove multiple spaces and tabs
		s/\s+/ /g;

		# Remove spaces after newlines

		# Remove leading hashes & spaces
		s/^\s*//;
		s/^#*//;
		s/^\s*//;

		#Remove trailing hashes and spaces
		s/#*\s*$//;

		# Remove 'underlines'
		#s/ *-{2,} *//;
		#s/ *={2,} *//;

		# Remove multiple spaces and tabs
		s/\s{2,}/ /g;
				
		# Remove dollar signs if a "revision" line
		s/\$Revision:/Revision:/;
		s/\$$//;
			
		# Add single space to end of line
		$_ .= ' ';

		$olines .= $_;
	}

	return $olines;
}
			
sub format_string_html {

	#<
	#? Create html from a tagged string
	#. Tags
	#; <LINK...>
	#; <ANCH...>
	#; <FILE>
	#; <HEADn>
	#; <STD>
	#; <STD_DESC>
	#; <HANG>
	#; <HANG_SHORT>
	#; <INDENT>
	#; <LIST>
	#; <CODE>
	#>

	my $string = $_[0];

	my @f;
	my $fstring = '';
	my @g;
	my $hangsize = 20;
	my $hangsize_short = 5;
	my $in;
	my $indentsize = 10;
	
	return "" if !$string;
	
	# Split in to paragraphs (on <STD>, <HEADER> or <INDENT>)
	@f = split /(?=<STD>)|(?=<STD_DESC>)|(?=<HEAD.>)|(?=<INDENT>)|(?=<HANG>)|(?=<HANG_SHORT>)|(?=<LIST>)|(?=<CODE>)|(?=<TABLE_ROW>)/, $string;

	foreach (@f) {
	
		chomp;
	
		if (get_flag('clean', 'l')) {
			next if $_ =~ />\s*Created:/;
			next if $_ =~ />\s*Revision/;
		}
		
		# HTML links
		s/<LINK(.*?)>\s*(\S+)/<A HREF="$1">$2<\/A>/g;
		s/<FILE(.*?)>\s*(\S+)/<A HREF="$1">$2<\/A>/g;
		s/<ANCH(.*?)>/<A NAME="$1"><\/A>/g;
	
		if (/<HEAD1>/) {
			s/<HEAD1> *//;
			$fstring .= "<H1> ";
			$fstring .= "$_ ";
			$fstring .= "</H1>\n";
			next;
		}
		
		if (/<HEAD2>/) {
			s/<HEAD2> *//;
			$fstring .= "<H2> ";
			$fstring .= "$_ ";
			$fstring .= "</H2>\n";
			next;
		}
	
		if (/<HEAD3>/) {
			s/<HEAD3> *//;
			$fstring .= "<H3> ";
			$fstring .= "$_ ";
			$fstring .= "</H3>\n";
			next;
		}
		
		if (/<HEAD4>/) {
			s/<HEAD4> *//;
			$fstring .= "<H4> ";
			$fstring .= "$_ ";
			$fstring .= "</H4>\n";
			next;
		}

		if (/<STD>/) {
			s/<STD> *//;
			$fstring .= "<p> ";
			$fstring .= "$_ ";
			$fstring .= "</p>\n";
			next;
		}

		if (/<STD_DESC>/) {
			s/<STD_DESC> *//;
			$fstring .= "<p> ";
			$fstring .= "$_ ";
			$fstring .= "</p>\n";
			next;
		}

		if (/<HANG>/ || /<HANG_SHORT>/) {
		
			my $h;
			$h = $hangsize if (/<HANG>/);
			$h = $hangsize_short if (/<HANG_SHORT>/);
			
			s/<HANG> *//;
			s/<HANG_SHORT> *//;
			
			# Approximate hanging indents for HTML
			# Get first word and remaining part of sentence
			@g = $_ =~ m/(\s*\S+)(.*)/;
				
			if (defined $g[0]) {
				my $i = int($h - 2*length($g[0]));
				$fstring .= $g[0];
				$fstring =~ s/_/ /;  # Remove first underscore
				if ($i > 0) {
					$fstring .= "&nbsp;" x $i;
				} 
				$fstring .= $g[1] if defined $g[1];
			} else {
				$fstring .= ' ';
			}
			$fstring .= "<br>\n";
			next;
		}
		
		if (/<INDENT>/) {
			s/<INDENT> *//;
			$fstring .= "&nbsp;" x $indentsize;
			$fstring .= "$_ ";
			$fstring .= "<br>\n";
			next;
		}
		if (/<LIST>/) {
			s/<LIST> *//;
			$fstring .= "$_ ";
			$fstring .= "<br>\n";
			next;
		}
		if (/<CODE>/) {
			s/<CODE> *//;
			$fstring .= "<span style=\"font-family: monospace;\">$_ </span> ";
			$fstring .= "<br>\n";
			next;
		}
		if (/<TABLE_ROW>/) {
			s/<TABLE_ROW> *//;
			s/\[/<td>/g;
			s/\]/<\/td>/g;
			$fstring .= "<tr> ";
			$fstring .= "$_ ";
			$fstring .= " </tr>\n";
			next;
		}
		
		# Treat as <STD>
		$fstring .= "<p> ";
		$fstring .= "$_ ";
		$fstring .= "</p>\n";
	}
	
	return $fstring;
}

sub wraphtml {
	
	#<
	#? Add header and footer to formatted html
	#>

	my $string = $_[0];
	
	my $fstring;
	
	$fstring .= '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">'."\n";
	$fstring .= "<HTML>\n";
	
	$fstring .= $string;
	
	$fstring .= "</HTML>\n";

	return $fstring;
}

sub formatpara_txt {

	#<
	#? Format a paragraph so that it fits onto a line with
	#  optional indent or hanging indent
	#>

	my $bigline = $_[0];
	my $linelength = $_[1] || 80;
	my $indent = $_[2] || 0;
	my $hang = $_[3] || 0;

	my $format;
	my $in = '';
	my $lines = '';
	my $linecount;
	my $len = 0;
	my $temp;
	my $word;
	my @words;
	
	# Revert &lt; and &gt; to < and >
	$bigline =~ s/&lt;/</g;
	$bigline =~ s/&gt;/>/g;
	
	# Split bigline into words
	@words  = split (/ /,$bigline);

	return '' if !defined $words[0];
	
	# Make indent
	if ($indent) {
		
		$in = ' ' x $indent;
		$lines .= $in;
		$len = $indent;
	}
	
	if ($hang) {
	
		$in = ' ' x $hang;
		# Add a space if word is longer than indent
		if (length($words[0]) >= length($in)) {
			$words[0] = $words[0]." ";
		}
		$format =  "%-".($hang)."s";
		$temp .= sprintf $format, $words[0];
		$temp =~ s/__/&-&/g; # Replace double underscores with placeholder
		$temp =~ s/_/ /g;    # Replace single underscores with spaces
		$temp =~ s/&-&/_/g;  # Replace placeholders with underscores
		$lines .= $temp;
		shift @words;
		$len = length($lines);
	}
	
	$linecount = -1;
	foreach $word (@words) {
		
		++$linecount;

		# Put space back
		$word .= ' ';

		# Calculate new length of line
		$len += length($word);

		if ($len >= $linelength) {

			if ($linecount > 0) {
				$lines .= "\n".$in.$word;
				$len = length $in.$word;
			} else {
				$lines .= "\n".$word;
				$len = length $word;
			}
			next;
		}
		$lines .= $word;
	}

	# Add final carriage return
	$lines .= "\n";
	return $lines;
}


sub heading {

	#<
	#? Format a heading
	#; Requires: String, optional underline character (defaults to -)
	#; Returns: Formatted string
	#>

	my $string = $_[0];
	my $under = $_[1] || '-';
	
	my @f;
	my $l = 0;;
	my $out;
	
	$string = '' if !defined $string;
	
	@f = split "\n", $string;

	$out =  "\n";
	chomp $string;
	$out .=  "$string\n";
	
	foreach (@f) {
	
		s/ *$//; # Remove any trailing spaces
		$l = length($_) if length($_) > $l;
	}
		
	$out .=  $under x ($l);
	$out .= "\n";

	return $out;
}



return 1;
