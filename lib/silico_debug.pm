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
#! silico_debug.pm
#? Silico general debugging routines.
#. $Revision: 1.35.2.3.2.4 $
#>

use strict;
package Silico;

sub ensemble_printout {

	#<
	#? Print out all the key value pairs for a molecule.
	#. Only the first few molecules are printed out
	#  unless the optional second argument is set to 'all'.
	#  This is passed on to molecule_printout.
	#; Requires: ensemble, 'all' (optional), filehandle (optional)
	#
	#>

	my $molecules = $_[0];
	my $flag = uc ($_[1] || '');
	my $FH = $_[2] || *STDOUT;
	my $total_molecules = 5; # print $FH out this many molecules
	
	(print $FH "Ensemble_printout: Molecule is undefined\n") && carp() && return if !defined $molecules;

	if (defined $flag && $flag eq "ALL") {
		$total_molecules = 1000;
	}
	
	print $FH "\nEnsemble Printout\n";
	print $FH "=================\n";
	print $FH "Location $molecules\n";
	print $FH "Number of molecules ".($#{$molecules}+1)."\n";

	my $molcount = -1;
	foreach my $mol (@$molecules) {
		++$molcount;

		if ($molcount > $total_molecules) {
			print $FH "** Stopping.  Use option 'all' to print $FH all molecules **\n";
			last;
		}

		print $FH "** Molecule $molcount\n";
		molecule_printout ($mol, $flag, $FH);
	}

	print $FH "=================\n";
}
 
sub molecule_printout {

	#<
	#? Print out several atoms from a molecule record.
	#. Prints out all molecule records and those for the first 5 (MacOS)
	#  or 10 atoms.  Takes optional second argument 'last' prints last 
	#  10 atoms.  'all' prints all atoms.
	#; Requires: molecule, 'all' (optional), filehandle (optional).
	#>

	my $mol = $_[0];
	my $flag = $_[1];
	my $FH = $_[2] || *STDOUT;
	my $total_atoms = 9;
	
	$flag ||= '';
	$flag = uc $flag;
	
	carp();

	# print $FH out fewer atoms on a Mac because mine is slow!
	$total_atoms = 5 if ($^O eq "MacOS");

	if (defined $flag && $flag eq "ALL") {
		$total_atoms = 10000;  
		$total_atoms = $#{$mol->{ATOMS}} + 1 if (defined $mol->{ATOMS}[0]);
	}

	print $FH "\n";
	print $FH "Molecule Printout\n";
	print $FH "-----------------\n";
	if (!defined $mol) {
		print $FH "** Molecule not defined\n";
		print $FH "\n";
		return;
	}
	print $FH "** Molecule location $mol\n";
	print $FH "\n";

	if (ref $mol eq 'ARRAY') {

		print "** Molecule is incorrectly an array\n";
		print "@$mol\n";
		croak();
	}

	foreach my $key (sort keys(%$mol)) {

		# Key name
		printf $FH "\t%-15s\t",$key;

		# Key value
		if (defined $mol->{$key}) {
			print $FH "\t'$mol->{$key}'";
		} else {
			print $FH  "\tNOT DEFINED";
		}
		
		# Hashes
		if (ref ($mol->{$key}) eq 'HASH') {
		
			print $FH "\n";
			foreach (sort(keys %{$mol->{$key}})) {
			
				print $FH "\t              < '$_'";
				if (defined $mol->{$key}{$_}) {
					print $FH " '$mol->{$key}{$_}' >\n";
				} else {
					print $FH " NOT DEFINED >\n";
				}
			}
			next;
		}
		
		# Arrays
		if (ref ($mol->{$key}) eq 'ARRAY') {
			
			printf $FH " (%d)", $#{$mol->{$key}};
			print $FH " <";
			for (my $i=0; $i >= 0 && $i < 10 && $i <=$#{$mol->{$key}}; ++$i) {
				if (!defined $mol->{$key}[$i]) {
					print $FH "\tNOT DEFINED";
					next;
				}
				print $FH " '$mol->{$key}[$i]'";
			}
			print $FH " >";
			
			
		}
		print $FH "\n";
	}
	
	print $FH "Atom Records\n";
	print $FH "------------\n";

	my $i = -1;
	foreach my $atom (@{$mol->{ATOMS}}) {
		
		++$i;
		last if ($flag ne 'LAST' && $i > $total_atoms);
		next if ($flag eq 'LAST' && $i < ($#{$mol->{ATOMS}} - 9));
		print $FH "Atom ".($i+1)."\n";
		atom_printout ($atom, $FH);
	}
	
	if (!$flag && $i == 10) {
		print "Truncating atom printout after $i atoms\n";
	}

}

sub atom_printout {

	#<
	#? Print out a formated atom record for debugging purposes
	#. Reauires: atom, filehandle (optional)
	#>

	my $atom = $_[0];
	my $FH = $_[1] || *STDOUT;

	if (!defined $atom) {

		print $FH "ATOM NOT DEFINED\n";
		return;
	}
		
	print $FH "Location: $atom\n";
	
	foreach my $key (sort keys(%$atom)) {

		printf $FH "\t%-15s", $key;

		if (!defined $atom->{$key}) {
			print $FH "\tUNDEFINED\n";
			next;
		}

		if (ref ($atom->{$key}) eq '') {
			print $FH "\t'$atom->{$key}'\n";
			next;
		}

		# Arrays
		if (ref ($atom->{$key}) eq 'ARRAY') {
			print $FH "\t[";
			my $i = 0;
			foreach my $val (@{$atom->{$key}}) {
				if (++$i > 10) {
					print $FH " ...";
					last;
				}
				defined $val ? print $FH " '$val'" : print $FH " UNDEF";
			}
			print $FH " ]\n";
			next;
		}

		# Hashes
		if (ref ($atom->{$key}) eq 'HASH') {
			print $FH "\t{";
			my $i = 0;
			foreach my $val (keys %{$atom->{$key}}) {
				if (++$i > 10) {
					print $FH " ...";
					last;
				}

				print $FH "  ";
				if (defined $val) {
					print $FH "'$val'";
				} else {
					print $FH "NOT DEFINED";
				}

				print $FH " '";
				if (defined $atom->{$key}{$val}) {
					print $FH $atom->{$key}{$val};
				} else {
					print $FH "NOT DEFINED";
				}
				print $FH "'";
			}
			print $FH " }\n";
			next;
		}

		print $FH "SOMETHING ELSE\n";
	}
}

sub ensemble_print_connect {

	#<
	#? Print out connection table
	#. Has not been used for some time.  Check before using.
	#; Requires: ensemble, output filename (optional).
	#>

	my $ensemble = $_[0];
	
	my $molcount =0;
 
	if ($_[1]) {
		open (OUT, ">$_[1]");
	} else {
		open (OUT, ">&STDOUT");
	}
	
	foreach my $mol (@$ensemble) {
	
		print OUT "Molecule $molcount\n";

		for (my $i=0; $i <$mol->{NUMATOMS}; ++$i) {

			my $atom = $mol->{ATOMS}[$i];
			
			print OUT "$i  connect\t";
			foreach (@{$atom->{CONNECT}}) {
				print OUT "$_\t";
			}
			print OUT "\n";

			print OUT "$i  order\t";
			foreach (@{$atom->{BORDERS}}) {
				print OUT "$_\t";
			}
			print OUT "\n";

			print OUT "$i  submol\t$atom->{SUBMOL}[$i]\n" if (defined $mol->{SUBMOL}[$i]);
		}
		++$molcount;
	}
}

#
#  Memory usage
#

sub memuse {

    	#<
	#? Print out memory usage found using the SZ record of the 'ps' command
	#; Requires: string to identify calling location.
	#>

	return if ! $Silico::memuse;
	
	my $tag = $_[0] || '';
	my  $s = `ps -lp $$`;
	my @f = split "\n", $s;
	my $i = 0;
	my @g = split " ", $f[0];
	foreach (@g) {
		last if m/SZ/;
		++$i;
	}
	@g = split " ", $f[1];
	my $v = sprintf "%3.2f", $g[$i]/1024;
	print "Memuse $tag: $v MB\n";
}

return 1;
