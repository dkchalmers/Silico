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
#. $Revision: 1.1.2.33 $
#>

use strict;
package Silico;


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

	$ret = run_gmx_local($gmx_cmd, $outbase, $hostlist, $background) if $hosttype eq 'local' ;
	$ret = run_gmx_slurm_single($gmx_cmd, $outbase)	if $hosttype eq 'slurm_cluster';

	silico_msg("Error in hosttype") if !defined $ret;

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

	$ret = run_grompp_gmx_local($grompp_cmd, $gmx_cmd, $postprocess_cmd, $outbase, $hostlist, $background) if $hosttype eq 'local' ;
	$ret = run_grompp_gmx_slurm_single($grompp_cmd,$gmx_cmd, $postprocess_cmd, $outbase) if $hosttype eq 'slurm_cluster';

	silico_msg("Error in hosttype") if !defined $ret;

	return $ret;
}

sub run_gmx_local {

	#<
        #? Run a gromacs job on local workstation cluster
        #; Requires: gromacs exec command, output file base name
        #; Returns: job identifier
        #>

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];
	my $hostlist = $_[2];
	my $background = $_[3];
	
	my $host;
	my $jobid;
	
	if (get_sflag("nogpu")) {
	
		($host, my $nt) = find_available_localhost_cpu($hostlist, 10);
		
		$gmx_cmd .= " -nt $nt -nb cpu -bonded cpu -pme cpu -update cpu ";
		$gmx_cmd =~ s/  / /g;

		$jobid = "$host\_nogpu";
		
	} else {
	
		($host, my $gpuid, my $nt) = find_available_localhost_gpu($hostlist);

		$gmx_cmd .= " -nt $nt ";
		$gmx_cmd .= " -gpu_id $gpuid ";
		$gmx_cmd =~ s/  / /g;
	
		if (get_sflag('pin')) {
	
			# This does not take into account that we might 
			# be overloading the gpu and there might
			# actuall be free cores
			my $offset = $nt * $gpuid;
	
			$gmx_cmd .= "-pin on -pinoffset $offset ";
		}
		
		$jobid = "$host\_$gpuid";
	} 
	
	my $stderr_string;
	if (get_sflag('print')) {
	
		$gmx_cmd .= "-v ";
		$stderr_string = '';
	} else {
		$stderr_string = ">& $outbase\_mdrun_stderr.out ";
	}
	
	$gmx_cmd .= $stderr_string;
	
	print "Run: $gmx_cmd\n";

	run_ssh($host, $gmx_cmd, $background, 1);
	
	return $jobid;
}

sub run_grompp_gmx_local {

	#<
        #? Run a grompp then gromacs job on local workstation cluster
        #; Requires: grompp command, gromacs comand, postprocess command, output file base name, hostlist, flag run job in background
        #; Returns: job identifier
        #>

	my $grompp_cmd = $_[0];
	my $gmx_cmd = $_[1];
	my $postprocess_cmd = $_[2] || '';
	my $outbase = $_[3];
	my $hostlist = $_[4];
	my $background = $_[5];
	
	my $host;
	my $jobid;

	my $fnt = get_sflag('nt');
	
	if (get_sflag("nogpu")) {
	
		($host, my $nt) = find_available_localhost_cpu($hostlist, 10);

		$nt = $fnt if $fnt;
		
		$gmx_cmd .= " -nt $nt -nb cpu -bonded cpu -pme cpu -update cpu ";
		$gmx_cmd =~ s/  / /g;

		$jobid = "$host\_nogpu";
		
	} else {
	
		($host, my $gpuid, my $nt) = find_available_localhost_gpu($hostlist);

		$nt = $fnt if $fnt;

		$gmx_cmd .= " -nt $nt ";
		$gmx_cmd .= " -gpu_id $gpuid ";
		$gmx_cmd =~ s/  / /g;
	
		if (get_sflag('pin')) {
	
			# This does not take into account that we might 
			# be overloading the gpu and there might
			# actuall be free cores
			my $offset = $nt * $gpuid;
	
			$gmx_cmd .= "-pin on -pinoffset $offset ";
		}
		
		$jobid = "$host\_$gpuid";
	} 
	
	my $stderr_string;
	if (get_sflag('print')) {
	
		$gmx_cmd .= "-v ";
		$stderr_string = '';
	} else {
		$stderr_string = ">& $outbase\_mdrun_stderr.out ";
	}
	
	$gmx_cmd .= $stderr_string;
	
	# Get filebase
	$grompp_cmd =~ m/(\S+).tpr/;
	my $filebase = $1;
	silco_msg('d', "Could not get filebase") if !$filebase;
	
	my $cmdfile = "$filebase\_run.in";
	
	open (FILE, ">$cmdfile") || silico_msg('d', "Could not write file $cmdfile");
	
	print FILE "$grompp_cmd\n";
	print FILE "test ! -f $filebase.tpr && touch stop && exit\n";
	print FILE "$gmx_cmd\n";
	print FILE "$postprocess_cmd\n" if $postprocess_cmd;
	
	close (FILE);
	
	my $cmd;
	$cmd = "sh $cmdfile";

	run_ssh($host, $cmd, $background, 1);
	
	return $jobid;
}

sub run_gmx_slurm_single {

	#<
        #? Run a gromacs job on a slurm cluster, one slurm job per gmx job)
        #; Requires: gromacs exectue command, output file base name, list of hosts (optional), flag to run job in background
        #; Returns: job identifier
        #>

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];
	
	my $nt    = get_lflag('slurm-number-threads');
	my $time  = get_lflag('slurm-run-time') || "7:0:0"; # 7 days
	my $nnode = get_sflag('s_nnode') || 1;
	my $ngpu  = get_sflag('s_ngpu') || 1;
	my $mem   = get_slfag('s_mem') || 10000;
	
	if (get_sflag("no_gpu")) {
		$gmx_cmd .= "-nb cpu";
	} 
	
	my $slurmstring;
	
	if (!get_sflag('multidir')) {
	
		#
		# NOT multidir
		#
		
		my $ntasks = 1;
		my $tpn = 1;
		my $nt = get_sflag('nt') || 12;
		my $partition = find_slurm_partition();
		
		$slurmstring = "#!/bin/bash
		
#SBATCH --nodes=$nnode
#SBATCH --partition=$partition
#SBATCH --time=$time
#SBATCH --ntasks-per-node=$tpn 
#SBATCH --cpus-per-task=$nt
#SBATCH --gres=gpu:$ngpu
#SBATCH --mem=$mem

$gmx_cmd -nt $nt
rm $outbase.jobid\n";

	} else {
	
		#
		# Multidir
		#
	
		#my $ntasks = $#dirs+1;
		# This needs to be fixed
		#
		my $ntasks = $nt;
		#
		#
		my $tpn = $ntasks/$nnode;
		my $gpustring = '';
		
		if (!$nt) {
		
			{
				use integer;
				$nt = 48/$tpn;
				$nt ||= 1;
			}
		}
		
		if ($tpn != int($tpn)) {
			silico_msg('d', "Number of tasks-per-node (ntasks: $ntasks nnode: $nnode tpn: $tpn) is not an integer.\n");
		}
		
		# Set gputasks if they don't distribute evenly onto GPUs
		if ($tpn%$ngpu != 0) {
		
			my @s;
		
			my $g = 0;
			foreach (1..($tpn*2)) {
				push @s, $g;
				++$g;
				$g = 0 if $g == $ngpu;
			}
			
			@s = sort(@s);
			$gpustring = join '', @s;
		}
		
		if ($gpustring) {
			
			$gpustring = "-nb gpu -pme gpu -gputasks ".$gpustring;
		}
		
		print "nt: $nt time: $time nnode: $nnode ngpu: $ngpu ntasks $ntasks tasks-per-node: $tpn gpustring: $gpustring\n";
	
		$slurmstring = 
"#!/bin/bash
		
#SBATCH --nodes=$nnode
#SBATCH --partition=CLUSTER
#SBATCH --time=$time
#SBATCH --ntasks-per-node=$tpn
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:$ngpu
#SBATCH --mem=$mem

$gmx_cmd $gpustring
rm $outbase.jobid\n";
	
	}

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
        #? Run a gromp && gromacs job on slurm cluster, one slurm job per gmx job)
        #; Requires: gromacs exectue command, output file base name, list of hosts (optional), flag to run job in background
        #; Returns: job identifier
        #>

	my $grompp_cmd = $_[0];
	my $gmx_cmd = $_[1];
	my $postprocess_cmd = $_[2] || '';
	my $outbase = $_[3];
	
	my $nt = get_sflag('nt') || 24;
	my $time = get_lflag('slurm-run-time') || "7:0:0"; # 7 days
	my $nnode = get_sflag('slurm-nnode') || 1;
	my $ngpu = get_sflag('slurm-ngpu') || 1;
	my $mem = get_sflag('slurm-mem');
	
	# Get filebase
	$grompp_cmd =~ m/(\S+).tpr/;
	my $filebase = $1;
	silco_msg('d', "Could not get filebase") if !$filebase;
	
	if (get_sflag("no_gpu")) {
		$gmx_cmd .= "-nb cpu";
	} 
	if (get_sflag('multidir')) {
		silco_msg('d', "Multidir not supported");
	}
	
	my $ntasks = 1;
	my $tpn = 1;
	
	my $memstring = '';
	$memstring = "#SBATCH --mem=$mem" if defined $mem;
		
	my $slurmstring = "#!/bin/bash
		
#SBATCH --nodes=$nnode
#SBATCH --partition=CLUSTER
#SBATCH --time=$time
#SBATCH --ntasks-per-node=$tpn 
#SBATCH --cpus-per-task=$nt
#SBATCH --gres=gpu:$ngpu
$memstring

$grompp_cmd
test ! -f $filebase.tpr && { touch stop; exit; }
$gmx_cmd -nt $nt
$postprocess_cmd
rm $outbase.jobid\n";

	print "* Submitting slurm job *\n";
	#print "$slurmstring\n";
	
	my $slurmfile = "$outbase.slurm";
	open(OUT, ">$slurmfile") || silico_msg("d", "Could not open $slurmfile for writing\n");
	print OUT $slurmstring;
	close OUT;
	
	system "sbatch $slurmfile | sed 's/.*job //' > $outbase.jobid ";
	
	# To be done
	return "UNK";
}

sub run_gmx_slurm_multi {

	#<
        #? Run a gromacs job on Slurm cluster. All jobs in one batch submission)
        #; Requires: gromacs exec command, output file base name
        #; Returns: job identifier
        #>

	my $gmx_cmd = $_[0];
	my $outbase = $_[1];

	my $qos = get_sflag('s_qos');
	my $nt = get_lflag('slurm-number-threads') || 12;
	my $time = get_lflag('slurm-run-time') || 12;
	my $nnode = get_sflag('s_node')|| 1;
	my $ngpu = get_sflag('ngpu') || 1;
	my $account = get_sflag('s_account') || 'mo99';
	my $partition = get_sflag('s_part') || 'm3g';
	my $gmx_module =  get_sflag('s_module') || silico_msg('d', "Gromacs module not set");
	
	$gmx_cmd =~ s/gmx/gmx_mpi/;

	my $slurmstring = 
	
"#!/bin/bash

#SBATCH --account=$account          
#SBATCH --time=$time       
#SBATCH --partition=$partition
#SBATCH --nodes=$nnode
#SBATCH --gres=gpu:V100:$ngpu 
#SBATCH --ntasks=$nt
#SBATCH --qos=$qos 

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
	
	# To be done
	return "UNK";
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
