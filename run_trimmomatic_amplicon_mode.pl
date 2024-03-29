#!/usr/bin/perl

#Runs trimmomatic-0.33 in batch on 16s rRNA amplicon data
#Written by Tom de Man

use warnings;
use strict;
use File::Basename;

my $path;

my @files=@ARGV;

my %paired_files;
foreach my $file (@files){
	my ($file_name,$dir)=fileparse($file);
	$path = $dir;
	if($file_name =~ /(.+)_R([1|2])_/){
		$paired_files{$1}[$2-1]=$file;
	#attempt different naming scheme
	}elsif($file_name =~ /(.+)_([1|2])/){
		$paired_files{$1}[$2-1]=$file;
	}else{
		warn "Input file does not contain '_R1_' or '_R2_' in name: $file";
	}
}

foreach my $name (sort keys %paired_files){
	unless(defined($paired_files{$name}[0]) && defined($paired_files{$name}[1])){
		warn "Couldn't find matching paired end files for file starting with: $name";
		next;
	}
	print "trimming your data....\n";
	print "----------------------\n";
	print "$paired_files{$name}[0]"." <--> "."$paired_files{$name}[1]"."\n";

	my $cmd="trimmomatic PE -trimlog out.trimlog -phred33 -threads 2 $paired_files{$name}[0] $paired_files{$name}[1] $path$name"."_R1_paired_trimmed.fastq $path$name"."_R1_single_trimmed.fastq $path$name"."_R2_paired_trimmed.fastq $path$name"."_R2_single_trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:225";
	print $cmd,"\n";
	die if system($cmd);
}
