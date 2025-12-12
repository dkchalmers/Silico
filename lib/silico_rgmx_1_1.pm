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
#! silico_rgmx.pm
#? Silico routines to run Gromacs
#. $Revision: 1.1.2.13 $
#>

use strict;
package Silico;

sub rgmx_setup {

	#<
        #? Set up flags for Gromacs jobs
        #; Requires: Nothing
        #; Returns: Nothing. Sets flags
        #>
	
	#Gromacs options
	make_flag('g_exe',      'gmx-exe', "Gromacs executable", 1, 'gmx');
	make_flag('g_nt',       'gmx-num-threads', "Gromacs -nt flag", 1);
	make_flag('g_nogpu',    'gmx-no-gpu',  "Do not use a GPU");
	make_flag('g_maxwarn',  'gmx-maxwarn', "Gromacs grompp maxwarn flag", 1);
	make_flag('g_pin',      'gmx-pin',     "Gromacs pin threads flag");
	make_flag('g_pinstride','gmx-pinstride',"Gromacs pinstride for use with gmx-pin", 1, 1);
	make_flag('g_pinoffset','gmx-pin-offset-array', "Space delimited array of pinoffset values", 1, "0 24");
	make_flag('g_notune',   'gmx-notune-pme',"Do not tune pme");
	make_flag('g_f',        'gmx-flags',   "Additional gromacs flags", 1);
	make_flag('g_v',        'gmx-verbose', "Print out gromacs timings with gmx -v flag");
	
	# Gromacs local job options
	make_flag('h',          'hostlist', "Local machine host list. Space separated string (use quotes)", 1);
	make_flag('maxj',       'max-jobs', "Local machine maximum number of GPU jobs(gmx, desmond, etc) allowed per GPU", 1, 2, undef, "INTEGER > 0");

	# SLURM job options
	make_flag('s_qos',      'slurm-qos',       "SLURM quality of service flag (qos)", 1, 'rtq');
	make_flag('s_time',     'slurm-run-time',  "SLURM max job run time (hrs:mins:secs)",  1);
	make_flag('s_nnode',    'slurm-num-node',  "SLURM number of nodes to request", 1);
	make_flag('s_ncpu',     'slurm-num-cpu',   "SLURM number of cpus to request", 1);
	make_flag('s_ngpu',     'slurm-num-gpu',   "SLURM number of GPUs to request",  1);
	make_flag('s_account',  'slurm-account',   "SLURM account", 1);
	make_flag('s_part',     'slurm-partition', "SLURM partition", 1);
	make_flag('s_option',   'slurm-option', "Additional SLURM option", 1);
}


sub run_gromacs {

	#<
        #? Run a gromacs job
        #; Requires: gromacs exec command, output file base name, list of hosts (optional), flag to run job in background
        #; Returns: job identifier
        #>

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];
	my $hostlist = $_[2];
	my $background = $_[3];
	
	my $hosttype = find_host_queue();
	my $ret;

	if ($hosttype eq 'local') {
	
		$ret = run_gmx_local($gmx_cmd, $outbase, $hostlist, $background);
		
	} elsif ($hosttype eq 'slurm_cluster') {
	
		$ret = run_gmx_slurm_single($gmx_cmd, $outbase);
		
	} else  {
		silico_msg('d', "Unknown hosttype: $hosttype\n");
	}

	return $ret;
}

sub run_grompp_gromacs {

	#<
        #? Run a grompp & gromacs job
        #; Requires: grommpp exe cmd, gromacs exec command, postprocess cmd, output file base name, 
	#  list of hosts (optional), flag to run job in background
        #; Returns: job identifier
        #>

	my $grompp_cmd = $_[0];
	my $gmx_cmd = $_[1];
	my $postprocess_cmd = $_[2];
	my $outbase = $_[3];
	my $hostlist = $_[4];
	my $background = $_[5];
	
	my $hosttype = find_host_queue();
	my $ret;
	
	if ($hosttype eq 'local') {

		$ret = run_grompp_gmx_local($grompp_cmd, $gmx_cmd, $postprocess_cmd, $outbase, $hostlist, $background) ;
	
	} elsif ($hosttype eq 'slurm_cluster') {
		
		$ret = run_grompp_gmx_slurm_single($grompp_cmd,$gmx_cmd, $postprocess_cmd, $outbase);
	
	} else {
		silico_msg('d', "Unknown hosttype: $hosttype\n");
	}

	return $ret;
}

sub get_gmx_cmd {

	#<
        #? Formulate Gromacs execution command with flags
        #; Requires: initial gromacs exec command, number of threads, gpuid, output filebase 
        #; Returns: command string
        #>

	my $cmd = $_[0] || get_lflag('gmx-exe') || 'gmx';
	my $nt = $_[1]  || get_lflag('gmx-num-threads') || '';
	my $gpuid = $_[2];
	my $outbase = $_[3] || '';

	my $pin = get_lflag('gmx-pin');
	my $xflags = get_lflag('gmx-flags') || '';
	
	$cmd .= " -nt $nt" if $nt;
	$cmd .= " -gpu_id $gpuid" if defined $gpuid;
	$cmd .= " $xflags";
	
	if ($pin) {
	
		my $pinstride = get_lflag('gmx-pinstride') || 1;
	
		if (!defined $Silico::pinoffset_array) {
			
			my $pinoffset = get_lflag('gmx-pinoffset-array') || '';
			my @arr = split " ", $pinoffset;
			print "arr: @arr";
			
			if ($#arr == -1) {
				print "Guessing pinoffsets to be (0, $nt)\n";
				@$Silico::pinoffset_array  = (0, $nt);
			} else {
				print "Setting pinnoffset array to: @arr\n";
				@$Silico::pinoffset_array  = @arr;
			}
		}
	
		$cmd .= " -pin on -pinstride $pinstride -pinoffset $Silico::pinoffset_array->[0] ";
		
		push @$Silico::pinoffset_array, shift @$Silico::pinoffset_array;
	}
	
	$outbase .= '_' if $outbase;
	my $stderr_string = " >& $outbase"."mdrun_stderr.out ";
	
	$cmd .= $stderr_string;
	
	$cmd =~ s/  / /g;
	$cmd =~ s/^ //g;
	
	#print "GMX command: $cmd\n";

	return $cmd;
}

sub run_gmx_local {

	#<
        #? Run a gromacs job on local workstation cluster
        #; Requires: gromacs exec command, output file base name
        #; Returns: job identifier
        #>
	
	#carp();

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];
	my $hostlist = $_[2];
	my $background = $_[3];
	
	my $host;
	my $jobid;
	
	($host, my $gpuid, my $nt) = find_available_localhost_gpu($hostlist);

	$gmx_cmd = get_gmx_cmd($gmx_cmd, $nt, $gpuid, $outbase);
		
	$jobid = "$host\_$gpuid";

	#print "Run: $gmx_cmd\n";

	run_ssh($host, $gmx_cmd, $background, 1);
	
	return $jobid;
}

sub run_grompp_gmx_local {

	#<
        #? Run a grompp then gromacs job on local workstation cluster
        #; Requires: grompp command, gromacs comand, postprocess command, output file base name, hostlist, flag run job in background
        #; Returns: job identifier
        #>
	
	#carp();

	my $grompp_cmd = $_[0];
	my $gmx_cmd = $_[1];
	my $postprocess_cmd = $_[2] || '';
	my $outbase = $_[3];
	my $hostlist = $_[4];
	my $background = $_[5];
	
	my $host;

	my $fnt = get_sflag('g_nt');
	
	($host, my $gpuid, my $nt) = find_available_localhost_gpu($hostlist);

	my $jobid = "$host\_$gpuid";
	my $stopfile = "stop.$jobid";

	$nt = $fnt if $fnt;

	$gmx_cmd = get_gmx_cmd($gmx_cmd, $nt, $gpuid, $outbase);
	
	# Get filebase
	$grompp_cmd =~ m/(\S+).tpr/;
	my $filebase = $1;
	silco_msg('d', "Could not get filebase") if !$filebase;
	
	my $cmdfile = "$filebase\_run.in";
	
	open (FILE, ">$cmdfile") || silico_msg('d', "Could not write file $cmdfile");
	
	print FILE "#\n#Host: $host GPUId: $gpuid\n#\n";
	print FILE "$grompp_cmd\n";
	# Two tests and sleep to allow for slow write
	print FILE "test ! -f $filebase.tpr && sleep 5 && test ! -f $filebase.tpr &&touch $stopfile && exit\n";
	print FILE "$gmx_cmd\n";
	print FILE "$postprocess_cmd\n" if $postprocess_cmd;
	
	close (FILE);
	
	my $cmd;
	$cmd = "sh $cmdfile";

	run_ssh($host, $cmd, $background, 1);
	
	return $jobid;
}

sub run_gmx_slurm_single {

	#carp();

	#<
        #? Run a gromacs job on a Slurm cluster, one slurm job per gmx job)
        #; Requires: gromacs exectue command, output file base name, list of hosts (optional), flag to run job in background
        #; Returns: job identifier
        #>
	
	#carp();

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];
	
	my $nt =     get_sflag('g_nt');
	my $time =   get_lflag('slurm-run-time') || "7:0:0"; # 7 days
	my $nnode =  get_sflag('slurm-nnode') || 1;
	my $s_ncpu = get_sflag('slurm-ncpu') || 1;
	my $ngpu =   get_sflag('slurm-ngpu') || 1;
	my $mem  =   get_sflag('slurm-mem') || 10000;
	
	$gmx_cmd = get_gmx_cmd($gmx_cmd, $nt, undef, $outbase);
	
	my $s_nt = get_sflag('s_ncpu');
	
	my $ntasks = $s_nt || $nt || 24;

	my $partition = find_slurm_partition();

   	my $memstring = '';
        $memstring = "#SBATCH --mem=$mem" if defined $mem;
		
	my $slurmstring = 
	
"#!/bin/bash
		
#SBATCH --time=$time
#SBATCH --partition=$partition		
#SBATCH --nodes=$nnode
#SBATCH --ntasks=$ntasks
#SBATCH --gres=gpu:$ngpu
$memstring

$gmx_cmd
rm $outbase.jobid\n";


	print "* Submitting slurm job *\n";
	#print "$slurmstring\n";
	
	my $slurmfile = "$outbase.slurm";
	open(OUT, ">$slurmfile") || silico_msg('d', "Could not open $slurmfile for writing\n");
	print OUT $slurmstring;
	close OUT;
	
	system "sbatch $slurmfile | sed 's/.*job //' > $outbase.jobid ";
	
	# To be done
	return "UNK";
}

sub run_grompp_gmx_slurm_single {

	#<
        #? Run a gromp && gromacs job on a slurm cluster, one slurm job per gmx job)
        #; Requires: gromacs exectue command, output file base name, list of hosts (optional), flag to run job in background
        #; Returns: job identifier
        #>

	my $grompp_cmd = $_[0];
	my $gmx_cmd = $_[1];
	my $postprocess_cmd = $_[2] || '';
	my $outbase = $_[3];
	
	my $nt =     get_sflag('g_nt');
	my $time =   get_lflag('slurm-run-time') || "7:0:0"; # 7 days
	my $nnode =  get_sflag('slurm-nnode') || 1;
	my $s_ncpu = get_sflag('slurm-ncpu') || 1;
	my $ngpu =   get_sflag('slurm-ngpu') || 1;
	my $mem  =   get_sflag('slurm-mem') || 10000;
	my $option = get_sflag('s_option');
	
	# Get filebase
	$grompp_cmd =~ m/(\S+).tpr/;
	my $filebase = $1;
	silco_msg('d', "Could not get filebase") if !$filebase;
	
	$gmx_cmd = get_gmx_cmd($gmx_cmd, $nt, undef, $outbase);
	my $s_nt = get_sflag('s_ncpu');
	
	my $ntasks = $s_nt || $nt || 24;
	
	my $ostring = '';
	$ostring = "#SBATCH $option\n" if $option;

	my $partition = find_slurm_partition();
	
	my $slurmstring = 

"#!/bin/bash

#SBATCH --time=$time
#SBATCH --partition=$partition		
#SBATCH --nodes=$nnode
#SBATCH --ntasks=$ntasks
#SBATCH --gres=gpu:$ngpu
#SBATCH --mem=$mem

$ostring

$grompp_cmd
test ! -f $filebase.tpr && { touch stop; exit; }
$gmx_cmd
$postprocess_cmd
rm $outbase.jobid\n";

	print "* Submitting slurm job *\n";
	#print "$slurmstring\n";
	
	my $slurmfile = "$outbase.slurm";
	open(OUT, ">$slurmfile") || silico_msg("d", "Could not open $slurmfile for writing\n");
	print OUT $slurmstring;
	close OUT;
	
	system "sbatch $slurmfile | sed 's/.*job //' > $outbase.jobid ";
	my $v = `cat $outbase.jobid`;
	chomp $v;

	return $v;
}

sub run_gmx_slurm_multi {

	#<
        #? Run a gromacs job on a slurm cluster. All jobs in one batch submission)
        #; Requires: gromacs exec command, output file base name
        #; Returns: job identifier
        #>

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];

	my $qos = get_sflag('s_qos');
	my $nt = get_lflag('slurm-number-threads');
	my $time = get_lflag('slurm-run-time') || 12;
	my $nnode = get_sflag('s_node')|| 1;
	my $ngpu = get_sflag('ngpu') || 1;
	my $account = get_sflag('s_account') || 'mo99';
	my $gmx_module =  get_sflag('s_module') || silico_msg('d', "Gromacs module not set");
	my $option = get_sflag('s_option');
	
	$gmx_cmd =~ s/gmx/gmx_mpi/;
	$gmx_cmd = get_gmx_cmd($gmx_cmd, $nt, undef, $outbase);
	
	my $ostring = '';
	$ostring = "#SBATCH $option\n" if $option;

	my $partition = find_slurm_partition();

	my $slurmstring = 

"#!/bin/bash	

#SBATCH --time=$time    
#SBATCH --account=$account          
#SBATCH --partition=$partition
#SBATCH --nodes=$nnode
#SBATCH --gres=gpu:V100:$ngpu 
#SBATCH --ntasks=$nt
#SBATCH --qos=$qos 

# Environment
env

$option

module load gromacs/$gmx_module
mpiexec --bind-to none $gmx_cmd
rm $outbase.jobid\n";

	print "* Submitting slurm job *\n\n";
	print "$slurmstring\n";

	my $slurmfile = "$outbase.slurm";
	open(OUT, ">$slurmfile") || silico_msg('d', "Could not open $slurmfile for writing\n");
	print OUT $slurmstring;
	close OUT;
	system "sbatch $slurmfile | sed 's/.*job //' > $outbase.jobid ";
	my $v = `cat $outbase.jobid`;
	chomp $v;
	
	return $v;
}

sub check_log {

	#<
	#? Check that log file is complete
	#; Requires: output filebase
	#; Returns: undef (job is running), 0 (start new job) or 1 (restart job)
	#>

	my $outbase = $_[0];
	
	my $restart = 0;

	# Check status of log or slurm file
	my $log = "$outbase.log";
	
	my $hosttype = find_host_queue();
	
	if (-e $log) {

		print "\n";

		if (mdrun_log_is_complete($log)) {
		
			silico_msg('c', "Log file '$log' is already complete. Skipping\n");
			return undef;
		}

		if ($hosttype eq 'slurm_cluster') {
		
			if (-e "$outbase.jobid" && $hosttype eq 'slurm_cluster') {
			
				my $jid = `cat $outbase.jobid`;
				chomp $jid;
				my @f = split "\n", `squeue`;
				
				foreach (@f) {
			
					next if /JOBID/; # Header
					my $n;
					($n) = $_ =~ /^\s*(\d+)/;
					
					next if $n != $jid;
					
					silico_msg('c', "File $outbase.jobid exists with jobid '$jid', which matches a job in slurm queue. Skipping\n");
					return undef;
				}
			}

		} else {

			if (file_modified_since($log) < 10 && !get_sflag('force')) {

				silico_msg('c', "Log file '$log' is recently modified. Job may be still running. Skipping\n", 
					"	Use -force option or delete the log file to override\n");
				return undef;
			}
		}

		silico_msg('c', "Log file '$log' incomplete. Restarting job\n");

		# Set restart to complete job
		$restart = 1;
	} else {
	
		if (-e "$outbase.gro") {
		
			silico_msg('n', "Output file $outbase.gro present although no log file was present. Deleting it\n");
			system "rm $outbase.gro";
		}
	
	}
	
	return $restart;
}

sub edit_mdpfile {

	#<
	#? Edit template mdp file
	#; Requires: template file, output file
	#; Returns: nothing
	#>

	my $mdp_template  = $_[0];
	my $mdp_outfile   = $_[1];
	my $lambda = $_[2];
	my $lambda_values = $_[3];
	my $clambda_values = $_[4];
	my $vlambda_values = $_[5];
	my $rlambda_values = $_[6];
	
	#croak();
	
	my $min = get_lflag('minimise');
	my $alg = get_lflag('minimisation-algorithm') || 'steep';
	
	my $err = "*** lambda values not set! ****";

	open(IN,  $mdp_template) || silico_msg('d', "Could not open $mdp_template for reading\n");
	open(OUT, ">$mdp_outfile")   || silico_msg('d', "Could not open $mdp_outfile for writing\n");
	
	my ($t1, $t2, $t3, $t4);
	
	$t1 = $t2 = $t3 = $t4 = $err;

	$t1 = "@$lambda_values"   if defined $lambda_values;
	$t2 = "@$clambda_values"  if defined $clambda_values;
	$t3 = "@$vlambda_values"  if defined $vlambda_values;
	$t4 = "@$rlambda_values"  if defined $rlambda_values;
	
	
	while (<IN>) {
	
		if($min && /^\s*integrator/) {
			print OUT "integrator  =  $alg\n";
			next;
		}
		s/\$LAMBDA\$/$lambda/;
		s/\$LAMBDA_VALUES\$/$t1/;
		s/\$COUL_LAMBDA_VALUES\$/$t2/;
		s/\$VDW_LAMBDA_VALUES\$/$t3/;
		s/\$RESTRAINT_LAMBDA_VALUES\$/$t4/;
		
		print OUT $_;
	}

	close(OUT);
	close(IN);
}

sub print_grompp_log_extract_and_exit {

	my $outbase = $_[0];

	sleep 1;
	print heading("Extract of grompp log file");
	open (ERR, "$outbase\_grompp_stderr.out") || silico_msg('d', "Grompp failed. Could not open error output file $outbase\_grompp_stderr.out\n");
	my $i = 0;
	while (<ERR>) {
		++$i if /^GROMACS/;
		next if $i < 2;
		print;
	}
	silico_msg('d', "Grompp failed");
	exit;
}

sub mdrun_log_is_complete {

        #<
        #? Check if gmx  mdrun log file is complete
        #; Requires: logfile
        #; Returns: 1 if complete or 0
        #>

        my $log = $_[0];
	
	if (!-e $log) {
	
		silico_msg('w', "File $log not found\n");
		return 0;
	}

        if (!open(LOG, "<$log")) {
	
		carp();
		print "File open error: $!\n";
		silico_msg('w', "Could not open log file '$log' for reading");
		return 0;
	}

        my $terminated = 0; # Run has been terminated
        my $finished = 0;   # Has finished flag

        while (<LOG>) {

                $terminated = 0 if /^Restarting from checkpoint/;
                if (/^Received the.*stopping within.*steps/) {
			$terminated = 1;
			$finished = 0;
		}
                $finished = 1 if /^Finished mdrun/;
        }

        print "Logfile: terminated $terminated.\n" if $terminated;

        close(LOG);

        return 1 if !$terminated && $finished;

        return 0;
}

sub is_intel {

	my $string = `cat /proc/cpuinfo |grep GenuineIntel`;

	return 1 if $string;
}

sub gently_randomise_velocities {

	my $mol = $_[0];
	my $fraction = $_[1] || die;
	
	my $num = int($fraction * $mol->{NUMATOMS});
	
	print "Exchanging velocities of $num atoms\n";
	
	for my $i (0..$num) {
	
		my $atom1 = $mol->{ATOMS}[int(rand($mol->{NUMATOMS}))];
		next if !defined $atom1;
		
		my $el = $atom1->{ELEMENT};
		
		redo if !($el eq 'C' || $el eq 'O' || $el eq "N" || $el eq 'H');
	
		while (1) {
		
			my $atom2 = $mol->{ATOMS}[int(rand($mol->{NUMATOMS}))];
			
			next if $atom2->{ELEMENT} ne $el;
			
			foreach ('VX', 'VY', 'VZ') {
			
				my $v = $atom1->{$_};
				$atom1->{$_} = $atom2->{$_};
				$atom2->{$_} = $v;
			}
			last;
		}
	}
	return $mol;
}

sub stop_file_exists {

	#<
	#? Check if a 'stop' or 'stop.<jobid>' file exists
	#. File is removed if found
	#; Requires: optional die flag
	#; Returns: true or false
	#>

	my $jobid = $_[0];
	my $nodie = $_[1];
	
	my $sfile = 'stop';
	if ($jobid) {
	
	 	$sfile = "stop.$jobid";
	}

	if (-e $sfile) {
		sleep 10;
		unlink($sfile);
		if (!$nodie) {
			silico_msg('dq', "Found '$sfile' file. Grompp may have failed. Removing file and finishing\n");
		}
		silico_msg('c', "Found $sfile' file. Grompp may have failed.\n");
		return 1;
	}
	if (-e 'stop') {
		sleep 10;
		unlink('stop');
		if (!$nodie) {
			silico_msg('dq', "Found 'stop' file. Removing file and finishing\n");
		}
		silico_msg('c', "Found 'stop' file but ignoring it.\n");
		return 1;
	}
	
	return 0;
}




return 1;
