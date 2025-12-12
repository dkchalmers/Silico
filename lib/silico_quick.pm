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
#! silico_quick.pm
#? Silico Quick routines.
#. $Revision: 1.1.2.11 $
#>

use strict;

package Silico;

##################################################################
#
#	Quick read routines
#
##################################################################

sub read_quick_out {

	#<
	#? Read in a Quick out file
	#; Requires: filename, (optional) unused, unused, options string
	#; Returns: ensemble or undef if read failed.
	#>

	my $infile  = $_[0];
	my $options = $_[3] || '';

	add_option("QUIET", $options)    if (quiet());
	remove_option("QUIET", $options) if $Silico::debug;

	my $fr;
	$fr = open_molfile($fr, $infile, undef, undef, $options);
	if (!defined $fr) {
		silico_msg('e', "Can not open file $infile for reading\n");
		return undef;
	}
	my $FH = $fr->{FH};

	silico_msg('c', "Reading Quick Cartesian file: $infile\n") if $options !~ /\bQUIET\b/;

	my $mol = read_quick_out_single($fr, $options);

	my $mols;
	$mols->[0] = $mol;
	return $mols;
}

sub read_quick_out_single {

	#
	#? Read a single structure from a Quick .out file and return ensemble
	#; Requires: File record, options
	#; Returns: molecule. Undef if read failed.
	#

	my $fr      = $_[0];
	my $options = $_[1];

	my $mol;
	my $head;
	my $comments;
	my $keywords;
	my %element_count;

	my $FH = $fr->{FH};
	
	my $en = '';
	my $opt;
	my $scf_converged = 0;
	my $scf_failure_count = 0;

	# Find start
	my $line;
	
	%$mol = ();
	
	while ($line = <$FH>) {
	
		if ($line =~ /^  KEYWORD=/) {
			$keywords = $line;
			$keywords =~ s/  KEYWORD=\s*//;
			$keywords =~ s/\s*$//;
			chomp $keywords;
			
			$opt = 1 if $keywords =~ /OPTIMIZE/;
		}
	
		if ($line =~ /^ TOTAL ENERGY/) {
		
			chomp $line;
			my @f = split ' ', $line;
			$en = $f[3];
		}
		
		
		if (!$opt && $line =~ /-- INPUT GEOMETRY --/ ) {
			read_quick_coords($mol, $FH);
		}
		
		if ($line =~ /^ @ Begin Energy Calculation/) {
		
			$scf_converged = 0;
		}
		
		if ($line =~ /REACH CONVERGENCE AFTER/) {
		
			$scf_converged = 1;
		}
		
		if ($line =~ /NO CONVERGENCE/) {
		
			++$scf_failure_count;
		}
		
		#
		# For optimization jobs, file will read 'OPTIMIZED GEOMETRY INFORMATION' 
		# or GEOMETRY INFORMATION (NOT OPTIMIZED)
		#
		last if ($opt && $line =~ /GEOMETRY INFORMATION/) || (!$opt && $line =~ /@ End Energy calculation/);
	}
	
	if (!defined $line) {
		silico_msg('w', "Premature end to file\n");
		return undef;
	}
	
	if (!$opt) {
	
		$mol->{SDF_DATA}{OPT} = 'No';
		
	} elsif ($line =~ /NOT OPTIMIZED/) {
		$mol->{SDF_DATA}{OPT} = 'Incomplete';
		silico_msg('c', "Optimization incomplete\n");
	} else {
		$mol->{SDF_DATA}{OPT} = 'Complete';
	}
	
	if ($scf_converged) {
		$mol->{SDF_DATA}{SCF_CONVERGED} = 'Yes';
	} else {
		$mol->{SDF_DATA}{SCF_CONVERGED} = 'No';
		silico_msg('c', "SCF convergence failed");
	}
	
	$mol->{SDF_DATA}{SCF_FAILURE_COUNT} = $scf_failure_count if $scf_failure_count;
	$mol->{SDF_DATA}{QUICK_KEYWORDS} = $keywords || '';
	
	if ($opt) {
	
		# Skip 3 lines
		<$FH>;
		<$FH>;
		<$FH>;

		# Read in coordinates
		read_quick_coords($mol, $FH);

		# Skip 2 more lines
		<$FH>;
		<$FH>;
		read_quick_forces($mol, $FH);
	}
	
	# Get energy from last printed total energy to get around bug in QUICK with energies < -1000
	$mol->{SDF_DATA}{ENERGY} = $en;
	$mol->{SDF_DATA}{ENERGY_UNITS} = "Hartrees";
	
	quick_finalise($mol, $fr, $options);

	return $mol;
}

sub read_quick_coords {

	my $mol = $_[0];
	my $FH = $_[1];

	my $i = 0;
	while (my $line = <$FH>) {
		
		if (!defined $line) {
			silico_msg('w', "Premature end to file\n");
			return undef;
		}
			
		last if $line !~ /\S/; # End of molecule on blank line

		chomp $line;
			
		my @l = split ' ', $line;

		my $atom;

		$atom->{NUM}     = $i + 1;    #Atom number
		$atom->{ELEMENT} = $l[0];     #Atom element
		$atom->{SUBNAME} = "MOL";     #Residue name
		$atom->{SUBID}   = 1;         #Residue number
		$atom->{X}       = $l[1];     #x
		$atom->{Y}       = $l[2];     #y
		$atom->{Z}       = $l[3];     #z
			
		#print "$atom->{X} $atom->{Y} $atom->{Z}\n";

		push @{ $mol->{ATOMS} }, $atom;
		++$mol->{NUMATOMS};
		++$i;
	}
}

sub read_quick_forces {

	my $mol = $_[0];
	my $FH = $_[1];

	my $i = 0;
	while (my $line = <$FH>) {

		if (!defined $line) {
			silico_msg('w', "Premature end to file\n");
			return undef;
		}
		last if $line !~ /\S/;
		
		chomp $line;
		my @l = split ' ', $line;

		my $atom = $mol->{ATOMS}[$i];

		$atom->{FX} = $l[1];    #force-x
		$atom->{FY} = $l[2];    #force-y
		$atom->{FZ} = $l[3];    #force-z

		++$i;
	}
}


sub read_quick_in_single {

	#<
	#? Read a quick .in file
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $fr      = $_[0];
	my $options = $_[1];

	my $mol;
	my $comments;
	my $keywords;

	my $FH = $fr->{FH};

	# Read in z-matrix type Cartesian coordinates
	$keywords = <$FH>;
	$comments = <$FH>;

	# Read in keywords and next two comment lines
	if (!(defined $keywords && defined $comments)) {
		silico_msg('e', "Error reading file $fr->{FILENAME}!\n");
		return undef;
	}

	chomp $keywords;

	$mol->{SDF_DATA}{QUICK_KEYWORDS} = $keywords;

	# Read in body of input file
	read_quick_coords($mol, $FH);
	quick_finalise($mol, $fr, $options);

	return $mol;
}

sub quick_finalise {

	#<
	#? Fix up atomic elements from quick file read and perform checks
	#; Requires: molecule
	#; Returns: nothing
	#>

	my $mol     = $_[0];
	my $fr      = $_[1];
	my $options = $_[2] || '';

	#fr_printout($fr);
	
	if ($#{ $mol->{ATOMS} } < 0) {
	
		silico_msg('e', "The file $fr->{FILENAME} does not contain any atom records!".
			" Perhaps Quick failed to run.\n", "Skipping.\n");
		return undef;
	}

	# Check to see if we read in a z-matrix
	if (       $mol->{ATOMS}[0]{X} == 0
		&& $mol->{ATOMS}[0]{Y} == 0
		&& $mol->{ATOMS}[0]{Z} == 0
		&& $mol->{ATOMS}[1]{Y} == 0
		&& $mol->{ATOMS}[1]{Z} == 0
		&& $mol->{ATOMS}[2]{Z} == 0)
	{

		silico_msg('w', "You appear to be reading in a quick z-matrix file as Cartesian coordinates!\n");
		molecule_printout($mol);
	}

	# Check and fix atomic elements and set atom names
	my $i = 0;
	my %element_count;
	foreach my $atom (atoms($mol)) {

		# Checks and calculated fields
		$atom->{NAME} = $atom->{ELEMENT} . (++$element_count{ $atom->{ELEMENT} });    #Atom name

		# Set atom element number (a hash would be better)
		my $aname = ucfirst lc $atom->{ELEMENT};
		for (my $j = 1 ; $j <= $#Atomic_Elements ; ++$j) {
			if ($aname =~ /^$Atomic_Elements[$j]/) {
				$atom->{ELEMENT_NUM} = $j;
				last;
			}
		}

		if (!defined $atom->{ELEMENT_NUM}) {
			$atom->{ELEMENT}     = "Du";
			$atom->{ELEMENT_NUM} = 0;
		}
		++$i;
	}

	$mol->{NUMATOMS} = $i;
	$mol->{NUMBONDS} = 0;
	
	molecule_generate_bonds($mol) if ($options !~ /\bNOBOND\b/);
	
	$mol->{SOURCE_FILE_NAME} = $fr->{FILENAME};
	$mol->{SDF_DATA}{SOURCE} = $fr->{FILENAME};
}

##################################################################
#
#	Quick write routines
#
##################################################################

sub write_quick_in {

	#<
	#? Write a Quick input format file
	#. Currently prints out only the first molecule
	#  of the ensemble
	#; Requires: Ensemble or molecule, filename.
	#; Returns: undef if file open fails otherwise returns 1.
	#>

	my $molecules = ens($_[0]);
	my $outfile   = $_[1];
	my $keywords  = $_[2] || $molecules->[0]{SDF_DATA}{QUICK_KEYWORDS} || silico_msg('d', "Keywords not provided\n");

	my $atom;
	my $i;
	my $mol;
	my $name;
	my $success;

	$outfile = get_ofilebase($molecules) . ".in" if !$outfile;

	$mol  = $molecules->[0];
	$name = $mol->{NAME} || "Mol";
	chomp($name);

	my $fr;
	$fr = open_molfile($fr, $outfile, undef, 'write');
	if (!defined $fr) {
		silico_msg('e', "Can not open file $outfile for writing\n");
		return undef;
	}
	my $FH = $fr->{FH};

	#silico_msg('c', "Writing Quick Cartesian file: $outfile\n");

	print $FH "$keywords\n";
	print $FH "\n";

	foreach my $atom (atoms($mol)) {

		printf $FH "%-8s %10.6f %10.6f %10.6f\n", $atom->{ELEMENT}, $atom->{X}, $atom->{Y}, $atom->{Z};
	}

	close($FH);

	return 1;
}


##################################################################
#
#	Run Quick
#
##################################################################


sub setup_quick_flags {

	#<
	#? Setup Quick flags
	#; Requires: Nothing
	#; Sets: quick flags
	#; Returns: quick keywords
	#>
	
	my $method  = make_flag('method',  'method',        "Method [HF, DFT]", 1, "DFT");
	my $functional = make_flag('f',    'functional',    "Method [B3LYP, PBE0, etc.]", 1, "B3LYP");
	my $basis   = make_flag('b',       'basis',         "Basis set [6-311G, 6-311G*, cc-pvDZ, etc]", 1, "6-311G");
	my $charge  = make_flag('c',       'formal-charge', "Formal charge", 1, 0);
	              make_flag('sc',      'calc-formal-charge', "Calculate formal charge");
	my $scf     = make_flag('scf',     'scf',           "SCF option [100]", 1, 100);
	my $densrms = make_flag('densrms', 'densrms',       "DENSRMS option [1.0d-6]", 1, "1.0d-6");
	my $cutoff  = make_flag('cutoff',  'cutoff',        "CUTOFF option [1.0d-6]",  1, "1.0d-10");
	my $noopt   = make_flag('noopt',   'no-optimze',    "Do not use OPTIMIZE keyword");
	my $steps   = make_flag('steps',   'optimize-steps', "OPTIMIZE steps (0 uses Quick default value)", 1, 50);
	              make_flag('dry',     'dry-run',        "Do not actually submit jobs");
	
	# Run flags
	my ($nproc, undef, undef) = split ' ',  `cat /proc/cpuinfo | grep processor | wc`;
	
	my $cuda =    make_flag('cuda',     'cuda',          "Use CUDA version");
	my $mpi =     make_flag('mpi',      'mpi',           "Use MPI version");
	my $ncpu =    make_flag('n',        'ncpu',          "Number of MPI cpus", 1, int($nproc/2));
	my $noclean = make_flag('noclean',  'noclean',       "Do not delete output files");
	
	my $quick_exe = `which quick`;
	chomp $quick_exe;
	$quick_exe  =  make_flag('x', 'quick-exe', "Location of Quick executable [quick]", 1, $quick_exe);
	
	$quick_exe .= ".cuda" if $cuda;
	$quick_exe .= ".MPI" if $mpi;

	if (!-x $quick_exe) {
		silico_msg('d', "$quick_exe is not a valid executable file!\n");
	}
	
	if (!$mpi) {
		$Silico::quick_command = $quick_exe;
	} else {
		$Silico::quick_command = "mpirun -n $ncpu $quick_exe";
	}

	my $keywords = $method;
	
	$keywords .= " $functional" if $method ne 'HF';
	$keywords .= " BASIS=$basis" if $basis;
	$keywords .= " DENSRMS=$densrms" if $densrms;
	$keywords .= " CUTOFF=$cutoff"  if $cutoff;
	$keywords .= " SCF=$scf"  if $scf;
	$keywords .= " OPTIMIZE" if !$noopt;
	$keywords .= "=$steps" if !$noopt &&$ steps;
	
	$Silico::quick_keywords = $keywords;

	return $keywords;
}


sub run_quick_batch {

	#<
	#? Run quick minimisation or quick ESP charge calculation on molecule.
	#  Jobs are submitted as a batch
	#. Keywords should be in $mol->{SDF_DATA}{KEYWORDS}
	#; Requires: Molecule
	#; Returns: molecule or undef if failed
	#>

	my $mols = ens($_[0]);
	
	my $outlist;
	my $timeout = 2;
	
	#
	# Write input files
	#
	
	my $filebase = get_ofilebase($mols->[0]);
	my $fileext = get_oformat($mols->[0]);
	my $keywords = get_lflag("");
	
	my $i = 0;
	my $rmols;
	
	foreach my $mol (@$mols) {
	
		++$i;
		$mol->{SDF_DATA}{Run_Quick_Num} = $i;
		$mol->{SDF_DATA}{FORMULA} = formula_string(molecule_formula($mol, undef, undef, 1));
	
		print heading( "$mol->{SDF_DATA}{Run_Quick_Num}. $mol->{NAME}\n");
		
		my $n = sprintf "%04d", $mol->{SDF_DATA}{Run_Quick_Num};
		my $infile = $mol->{QUICK_INFILE} = quotestuff("$filebase\_$n\_quick.in");
		my $ofile = $mol->{QUICK_OUTFILE} = quotestuff("$filebase\_$n\_quick.out");
		
		# Calculate formal charge
		my $charge = get_sflag('c');
		if (get_lflag('calc-formal-charge')) {
			silico_msg('c', "Calculating charge\n");
			$charge = molecule_formal_charge($mol);
			
		}
		
		silico_msg('c', "Setting total charge to $charge\n");
		
		# Set keywords 
		my $k2 = "$Silico::quick_keywords CHARGE=".$charge;
		print "Keywords: $k2\n";
		$mol->{SDF_DATA}{QUICK_KEYWORDS} = $k2;
		
		if (-e $ofile && file_modified_since($ofile) <= 1) {
			print "File $ofile exists and looks like job is still running. Skipping\n";
			next;
		}
		
		if (-e $ofile && quick_job_is_complete($ofile) == 1) {
			print "Output file $ofile is complete. Skipping\n";
			next;
		}
		
		if (-e $ofile && quick_job_is_complete($ofile) == -1) {
			print "Output file $ofile is present but job failed. Skipping\n";
			next;
		}
		
		push @$rmols, $mol;
	}
	
	#
	# Run quick jobs
	#
	
	print "\n --- Jobs ---\n";
	$i = 0;
	foreach my $mol (@$rmols) {
	
		++$i;
		
		printf "%4d ", $mol->{SDF_DATA}{Run_Quick_Num};
		
		print "\n" if $i%20 == 0;
	}
	print "\n" if $i%20 != 0;
	print "\n";
	
	foreach my $mol (@$rmols) {
		
		submit_quick_job($mol) if !get_lflag('dry-run');
	}
	
	# Wait for outputs
	
	my $found;
	my $failed;
	my $died;
	my $count = 0;
	my $completed;
	while (1) {

		my $flag = 1;
		foreach my $mol (@$mols) {
		
			my $outfile = $mol->{QUICK_OUTFILE};
			my $ret;

			next if $found->{$outfile};
			next if $failed->{$outfile};
			next if $died->{$outfile};
			
			if (quick_job_is_complete($outfile) == 1) {
				$found->{$outfile} = 1;
				print "Found: $outfile\n";
				push @$completed, $mol;
				$count = 0;
				next;
			}
			
			if (quick_job_is_complete($outfile) == -1) {
				$failed->{$outfile} = 1;
				print "Job failed: $outfile\n";
				$count = 0;
				next;
			}
			
			# Check if it is a long time since log was updated
			if (-e $outfile && file_modified_since($outfile) > $timeout) {
				print "Job died: $outfile\n";
				$died->{$outfile} = 1;
				next;
			}
			
			$flag = 0;
		}

		last if $flag;
		
		print "Waiting\n" if ++$count%30 == 0;
		sleep 1;
	}
	
	return $completed;
}


sub submit_quick_job {

	#<
	#? Submit Quick job in background
	#. Keywords should be in $mol->{SDF_DATA}{KEYWORDS}
	#; Requires: Molecule,  filebase, location of quick_exe, flag to not clean up quick files after running
	#; Returns: molecule or undef if failed
	#>

	my $mol       = $_[0] || die();
	my $qinfile   = $_[1] || $mol->{QUICK_INFILE} || die();
	my $qoutfile  = $_[2] || $mol->{QUICK_OUTFILE} || die();
	my $keywords  = $_[3] || $mol->{SDF_DATA}{QUICK_KEYWORDS};
	my $quick_exe = $_[4] || get_lflag('quick-exe');
	
	system("rm -f $qoutfile");

	my $ens;
	$ens->[0] = $mol;
	if (!write_quick_in($ens, $qinfile)) {
		silico_msg('e', "Could not write $qinfile!\n");
		return undef;
	}

	my $run = "$Silico::quick_command $qinfile";

	run_quick_host($run, $Silico::hostlist, 1);
	
	return 1;
}

sub quotestuff {

	# Remove metacharacters that cause problems

	foreach ('\(', '\)') {
		$_[0] =~ s/$_//;
	}
	
	$_[0] =~ s/\*/x/;
	
	
	return $_[0];

}

sub process_quick_out {

	my $mols = $_[0];
	
	foreach my $mol (@$mols) {
	
		my $qinfile = $mol->{QUICK_INFILE};
		my $qoutfile = $mol->{QUICK_OUTFILE};
		my $keywords = $mol->{SDF_DATA}{QUICK_KEYWORDS};
	
		if (!-e $qoutfile) {
			silico_msg('w', "Quick failed to produce output file: $qoutfile.\n\n");
			return;
		}

		my $newens = read_quick_out($qoutfile, undef, undef, "NOBOND")
		  || silico_msg('d', "Failed to read molecule data from file $qoutfile!\n");
	  
		return undef if !$newens;

		my $newmol = $newens->[0];
	
		# Copy sdf_data from quick output molecule to original
		my @list = qw(ENERGY ENERGY_UNITS OPT SOURCE SCF_CONVERGED SCF_FAILURE_COUNT);

		foreach (@list) {
			$mol->{SDF_DATA}{$_} = $newmol->{SDF_DATA}{$_} if defined $newmol->{SDF_DATA}{$_};
		}
	
		my $string = get_quick_sdf_data_name();
		$mol->{SDF_DATA}{$string} = $mol->{SDF_DATA}{ENERGY};
	
		# Copy coordinates and forces if OPTIMIZE was specified
		if($keywords =~ /\bOPTIMIZE\b/) {
	
			my $i = 0;
			foreach my $atom (atoms($mol)) {
	
				foreach (qw(X Y Z)) {
			
					$atom->{$_} = $newmol->{ATOMS}[$i]{$_};
				$atom->{"F$_"} = $newmol->{ATOMS}[$i]{"F$_"};
				}
				++$i;
			}
		}
	}


	return $mols;
}

sub quick_job_is_complete {

	my $file = $_[0];
	
	return 0 if !-e $file;

	open (INFILE, $file) || silico_msg('d', "Could not open input file $file\n");
	
	while (<INFILE>) {
	
		#print;
		return 1 if $_ =~ /Normal Termination. Task Finished/;
		return -1 if $_ =~ /Error Termination. Task Failed/;
		return -1 if $_ =~ /Residue conversion error/;
		return -1 if $_ =~ /DL-FIND ERROR:/;
	}
	
	close (FILE);
	
	return 0;
}

sub run_quick_host {

	my $cmd = $_[0];
	my $hostlist = $_[1];
	my $background = $_[2];
	
	my $hosttype = find_host();
	my $ret;

	$ret = run_quick_local($cmd, $hostlist, $background) if $hosttype eq 'local';
	$ret = run_quick_slurm_cluster($cmd)	if $hosttype eq 'slurm_cluster';

	return $ret;
}

sub run_quick_local {

	my $cmd = $_[0];
	my $hostlist = $_[1];
	my $background = $_[2];
	
	my $host;
	if (get_lflag('cuda')) {
		($host, undef) = find_available_localhost_gpu($hostlist);
	} else {
		($host, undef) = find_available_localhost_cpu($hostlist);
	}
	
	print "Run: $cmd\n";
	
	$cmd = quotestuff($cmd);

	run_ssh($host, $cmd, $background, 1);
	
	return "$host";
}

sub run_quick_slurm_cluster {
	die();
}

sub get_quick_sdf_data_name {

	my $basis = get_lflag('basis');
	my $method = get_lflag('functional') || get_lflag('$method');
	my $string = uc("QUICK_ENERGY_$method\_$basis");
	
	return $string;
}



return 1;
