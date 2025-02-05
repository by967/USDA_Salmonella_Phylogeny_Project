#!/usr/bin/perl
# perl Download_Raw_Reads.pl
# run this from your directory with a TXT file (SRR.txt) containing a list of raw reads you want to download 
# ly276@cornell.edu
# 08/21/2022

use strict;
use warnings;
use Parallel::ForkManager;  # you may need to install this module: `cpanm Parallel::ForkManager`

my @txt = glob "*SRR*"; # searches for all files in the current directory that contain "SRR" in their names

foreach my $filename (@txt) {
    open INPUT, "<", $filename;

    my $pm = Parallel::ForkManager->new(20);  # number of parallel processes

    while (defined(my $line = <INPUT>)) {
        chomp($line);

        my $pid = $pm->start and next;  # forks child process
        system("fastq-dump --gzip --split-files $line");
        $pm->finish;  # ends child process
    }

    close INPUT;
    $pm->wait_all_children;  # wait for all child processes to finish
}


