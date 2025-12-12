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
#! silico_schedule.pm
#? Silico routines to schedule jobs on local workstations
#. $Revision: 1.1.2.31 $
#>

use strict;

package Silico;

#
#
#

make_flag('run-local', undef , "Ignore slurm and run as a local job");


sub proc_hostlist {

	my $h = $_[0] || '';

	my $hostlist;
	@$hostlist = ();
	my $myname = `hostname -s`;
	chomp $myname;

	# Spilt into @$hostlist
	foreach (split(' ', $h)) {

		# Change 'localhost' to current machine name
		# to avoid issues with ssh
		s/localhost/$myname/;

		push @$hostlist, $_;
	}

	# set @$hostlist if not set
	if (!$hostlist->[0]) {
		$hostlist->[0] = $myname;
	}

	# Randomise array order so that we dont get to instances
	# of the program following each other in the same host order.
	# To reduce the likelihood of overlading hosts due to race 
	# condition
	$hostlist = fisher_yates_shuffle($hostlist);
	
	silico_msg('c', "\nHostlist: @$hostlist\n");

	$Silico::hostlist = $hostlist;

	return $hostlist;
}

sub fisher_yates_shuffle {

	# From the O'Reilly cookbook
	# fisher_yates_shuffle( \@array ) : generate a random permutation
	# of @array in place

	my $array = $_[0];

	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}

	return $array;
}

sub get_hostlist {

	if (!$Silico::hostlist) {

		$Silico::hostlist = `hostname -s`;
		chomp $Silico::hostlist;
	}
	return $Silico::hostlist;
}

sub find_hostname {

	#<
	#? Find information about host from hostname
	#; Returns host name
	#>

	my $name = `hostname -s` || "unknown";
	chomp $name;

	return $name
}

sub find_host_queue {

	#<
	#? Find information job submission queue
	#; Returns slurm_cluster or 'local
	#>

	return $Silco::host_queue if defined $Silco::host_queue;
	$Silco::host_queue = 'local';

	if (get_sflag('run-local')) {
		print "Flag -run-local is set. Running locally\n\n";

	} elsif (host_is_slurm_cluster()) {
		$Silco::host_queue = 'slurm_cluster';
	}

	return $Silco::host_queue;
}

sub host_is_slurm_cluster {

	#<
	#? Is host a slurm cluster by checking if squeue command exists
	#; Returns true or false
	#>

	return $Silico::is_slurm_cluster if defined $Silico::is_slurm_cluster;

	my $val = `which squeue 2> /dev/null`;

	if ($val) {

		print "\nDetected slurm. Submitting jobs via squeue. Use the flag --run-local to override this and run jobs locally\n\n";
		$Silico::is_slurm_cluster = 1;
		return 1;
	}

	$Silico::is_slurm_cluster = 0;
	return 0;
}

sub find_slurm_partition {

	return $Silico::slurm_partition if defined $Silico::slurm_partition;

	my $part;

	if (get_sflag('s_part')) {

		$part = get_sflag('s_part');

	} else {

		my $val;
		$val = `scontrol show partition | grep PartitionName`;
		$val =~ s/PartitionName=//g;

		print "Available partitions. Please select one using the -s_part flag\n";
		print "$val\n";
		silico_msg('d');
	}

	print "Using Slurm partition '$part'\n\n";

	$Silico::slurm_partition = $part;
	return $part;
}

sub find_available_localhost_cpu {

	#<
	#? Find cpus available on local workstations
	#; Requires: hostlist (array), sleep time (optional)
	#; Flags: max-jobs, number-threads
	#; Returns: (host, numthreads)
	#>

	my $hostlist = $_[0] || $Silico::hostlist;
	my $waittime = $_[1] || 5;

	# Jobs to avoid
	my @proclist = qw (gmx desmond gdesmond quick python);

	silico_msg('d', "Hostlist not specified") if !defined $hostlist;

	my $maxjob = get_lflag('max-jobs')       || 1;
	my $nt     = get_lflag('number-threads') || 4;

	# We will rotate hostlist
	my @hostlistcopy = @$hostlist;

	my $hostname = find_hostname;

	my $count = 0;
	while (1) {

		++$count;
		foreach my $host (@hostlistcopy) {

			my $nproc = find_nproc($host);

			my $cmd = "ps -e -o pcpu,cmd";

			my $out = run_ssh($host, $cmd);
			my @procs = split "\n", $out;

			my $cputot = 0;
			my $nprocs;
			foreach my $p (@procs) {

				foreach my $pl (@proclist) {

					if ($p =~ /$pl/) {

						#print "pl $pl p $p\n";
						my @f    = split " ", $p;
						my $pcpu = shift @f;
						$cputot += $pcpu;
						++$nprocs;
					}
				}
			}

			print "Host: $host. CPUtot = $cputot%. nprocs = $nprocs\n" if $count % 4 == 0;

			my $np = $nproc * 100;
			if ($cputot <= $np) {

				print "Host: $host. CPUs available ($cputot, $np).\n";
				return ($host, $cputot);
			}

			# Rotate hostlist
			push(@$hostlist, (shift @$hostlist));

			#print "hostlist @$hostlist\n";
		}

		print "Waiting for CPUs\n" if $count % 4 == 0;
		sleep $waittime;
	}
}

sub find_available_localhost_gpu {

	#<
	#? Find gpus available on local sysgtems
	#; Requires: hostlist (array), maximum number of jobs to run on each gpu is provided as a flag
	#; Uses: List @procnames which is a list of process names to avoid when scheduling
	#; Flags: max-jobs
	#; Returns: host, gpuid, nprocgpu (or nothing)
	#>

	my $hostlist = $_[0] || $Silico::hostlist;

	silico_msg('d', "Hostlist not specified") if !defined $hostlist;

	my $maxjob = get_lflag('max-jobs') || 1;

	if (!defined $hostlist->[0]) {
		$hostlist->[0] = 'localhost';
	}

	# We will rotate hostlist
	my @hostlistcopy = @$hostlist;

	my $hostname = find_hostname;

	while (1) {

		my $count;
		foreach my $host (@hostlistcopy) {

			# Wait for job to start if this is the last host we submitted a job on
			sleep 2 if ($Silico::lasthost && $Silico::lasthost eq $host);

			my $ngpus = find_ngpu($host);

			# Rotate hostlist if no gpus
			if (!$ngpus) {
				push(@$hostlist, (shift @$hostlist));
				print "Host $host: no GPUs found\n";
				next;
			}

			my $nproc    = find_nproc($host);
			my $nprocgpu = int($nproc / $ngpus);

			silico_msg('d', "Number of requred CPUs calculated to be $nprocgpu") if $nprocgpu <= 0;

			my $jobcounts;
			foreach my $gpuid (0 .. $ngpus - 1) {

				my ($jobcount, $avoid, $s) = countjobs($host, $gpuid);

				my $rec;
				$rec->{GPUID} = $gpuid;
				$rec->{COUNT} = $jobcount;
				$rec->{AVOID} = $avoid;
				$rec->{PROCS} = $s || '';

				push @$jobcounts, $rec;

				if ($jobcount == 0) {
					last;
				}
			}

			@$jobcounts = sort { ($a->{AVOID} * 100 || $a->{COUNT}) <=> ($b->{AVOID} * 100 || $b->{COUNT}) } @$jobcounts;
			
			my $use = $jobcounts->[0];

			print "Host $host ($nproc, $ngpus) ";
			
			#print heading("Jobcounts");
			#ensemble_printout($jobcounts);
			
			#print heading("Jobcounts");
			#molecule_printout($use);

			if ( !$use->{AVOID} && $use->{COUNT} < $maxjob) {

				# Rotate hostlist if we put a job on the last gpu
				push(@$hostlist, (shift @$hostlist)) if $use->{GPUID} == $ngpus - 1;

				$Silico::lasthost = $host;

				print "GPU $use->{GPUID}. Jobs $use->{COUNT}";
				print " ($use->{PROCS})" if $use->{PROCS};
				print ". Available\n";

				return ($host, $use->{GPUID}, $nprocgpu);
			}

			print "GPU $use->{GPUID}. $use->{COUNT} Jobs ($use->{PROCS}). None\n";

			# Rotate hostlist
			push(@$hostlist, (shift @$hostlist));
		}

		++$count;
		my $sleep = 2;
		print "Waiting for GPUs\n\n" if ($count * $sleep) % 15 == 0;
		sleep $sleep;
	}
}

sub countjobs {

	my $host  = $_[0];
	my $gpuid = $_[1];

	my $out = run_nvidia_smi($host, "-i $gpuid");

	# List of GPU intensive processes we want to completely avoid standing on
	my @procnames1 = qw(desmond spdyn quick.cuda quick.cuda.MPI python);

	# List of processes we can share with
	my @procnames2 = qw(gmx mdrun mdrun_mpi);

	my @lines = split "\n", $out;

	my $count    = 0;
	my $jobcount = 0;
	my $avoid    = 0;
	my $s        = '';

	foreach my $l (@lines) {

		++$count if $l =~ /====/;
		next     if $count < 2;
		next     if $l =~ /====/ || $l =~ /----/;

		foreach my $procname (@procnames1) {

			next if $l !~ /$procname/;

			$avoid = 1;
			++$jobcount;

			if (!$s) {
				$s = $procname;
			} else {
				$s = "$s, $procname";
			}
		}

		foreach my $procname (@procnames2) {

			next if $l !~ /$procname/;

			++$jobcount;

			if (!$s) {
				$s = $procname;
			} else {
				$s = "$s, $procname";
			}
		}
	}

	#print "jobcount: $jobcount s: $s\n";

	return $jobcount, $avoid, $s;
}

sub find_nproc {

	#<
	#? Find number of cpus on host using ssh
	#; Requires: hostname
	#; Returns: Count of CPUs
	#>
	
	my $host = $_[0];

	return $Silico::nprocs->{$host} if defined $Silico::nprocs->{$host};

	my $nproc = run_ssh($host, "nproc");
	chomp $nproc;
	silico_msg('d', "'nproc' failed") if $nproc <= 0;

	# Divide by 2 to allow for virtual cores
	$nproc = int($nproc/2);

	$Silico::nprocs->{$host} = $nproc;

	return $nproc;
}

sub find_ngpu {

	#<
	#? Find number of gpus on host
	#; Requires: hostname
	#; Returns: Count of GPUs
	#>

	my $host = $_[0];

	return $Silico::ngpus->{$host} if defined $Silico::ngpus->{$host};
	
	# Check that nvidia-smi command exists
	my $g = run_ssh($host, 'which nvidia-smi', undef, undef, 1); # Supress error message
	if ($g eq '') {
	
		print "Command nvidia-smi not present on host $host\n";
		$Silico::ngpus->{$host} = 0;
		return 0;
	}

	$g = run_nvidia_smi($host, '-L');
	return undef if !defined $g;

	my @gpus  = split "\n", $g;
	my $ngpus = $#gpus + 1;

	$Silico::ngpus->{$host} = $ngpus;

	return $ngpus;
}

sub run_nvidia_smi {

	my $host  = $_[0];
	my $flags = $_[1];

	chomp $host;
	croak() if !$host;

	#print "host '$host'\n";

	my $g = run_ssh($host, "nvidia-smi $flags");

	if (!defined $g) {
		print "Host $host: No output from nvidia-smi. Is nvidia-smi installed?\n";
		return '';
	}

	if ($g !~ /GPU/) {
		print "Host $host: Nvidia-smi did not find any GPUs\n";
		return '';
	}

	return $g;
}

sub run_ssh {

	#<
	#? Run a command by ssh
	#. Command must be a single commad (i.e. no semicolons)
	#; Requires: host, commmand, flag to put job in background, flag to print output from command, flag to supress error message
	#; Returns: output from command (if background flag not set). Undef if command failed
	#>

	my $host       = $_[0];
	my $command    = $_[1];
	my $background = $_[2];
	my $print      = $_[3];
	my $nowarn     = $_[4];

	croak() if !defined $host;

	if (!$Silico::thishost) {

		$Silico::thishost = find_hostname;
	}

	if ($host ne 'localhost' && $host ne $Silico::thishost) {

		my $dir = `pwd`;
		chomp $dir;

		$command = qq{"cd $dir; $command"};
		$command = "ssh $host $command";
	}

	# Execute command using system and return
	if ($background) {

		$command = "$command &";

		print "Run (run_ssh): $command\n" if $print;
		my $failure = system $command;
		if ($failure) {
			silico_msg('w', "Command $command failed\n");
			return undef;
		}
		return '';
	}
	
	#print "command $command\n";

	# Run command using backticks and return output
	my $output = qx{$command};

	if ($? != 0) {

		silico_msg('w', "Command '$command' failed with exit code $?\n") if !$nowarn;
		return '';
	}

	return $output;
}

sub file_modified_since {

	#<
	#? Check when file last modified (minutes)
	#; Requires: filename
	#; Returns: time since modified in minutes
	#>

	my $file = $_[0];

	if (open(FILE, $file)) { 

		my $file_modified_timestamp = (stat(FILE))[9];

		# Time since modified in minutes
		my $m = (time() - $file_modified_timestamp) / 60;

		return $m;
	} 
	
	silico_msg('w', "Could not open file '$file' to determine modification time\n");
	
	return 0;
}


return 1;
