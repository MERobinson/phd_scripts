#!/usr/bin/perl
# Written: Mark Robinson Jan 2015
# Last updated: Mark Robinson Jan 2015
# USAGE:- perl read_coverage.pl
#
# Quick script to check read depth
# 
# Input: one or more .sgr files
# Output: to STDOUT (terminal)
################################################################################
use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
################################################################################
# User specified variables:
# $indir = input directory containing sgrs
# $percent = percent of chromosome to skip at ends (0.1 = 10%)
# $bin = binning of sgr files
# $dicty = whether sgr's are dicty genome or not (1 = true)
################################################################################
my $indir = "/Volumes/HARWOOD-WD3/Scripts/QC/sgr_read_coverage/indir";
my $bin = 10;
my $dicty = 1;
################################################################################

opendir(DIR, $indir) || die "Unable to access file at $indir: $!\n";
my @files = readdir(DIR);

foreach my $infile (@files){
    if (($infile !~ /^\.+/) && ($infile =~ /.*\.sgr/)){
    	open(IN, "$indir/$infile") || die "Unable to open $infile: $!";
    	print "Processing $infile \n";
    	my (%pos_count,%val_count);
    	LABEL: while (<IN>) {
    		chomp;
    		my @line = split('\t',$_);
    		my $chrn = $line[0];
    		my $pos = $line[1];
    		if ($dicty) {
    			if ($chrn eq "chr2" && $pos >= 3016080 && $pos <= 3768650) { ## Dicty chr2 duplication removal
    				next LABEL;
    			}
    		}
    		$pos_count{$chrn}++;
    		$val_count{$chrn}+=$line[2];
    	}

    	# calc coverage and output to STDOUT
    	for my $chrom (sort keys %val_count) {
    		my $coverage = nearest(0.001,$val_count{$chrom}/($pos_count{$chrom}*$bin));
    		print "Coverage for chromosome $chrom : $coverage\n";
    	}
    	close(IN);
    }
}
close(DIR);


