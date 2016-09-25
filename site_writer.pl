#!/usr/bin/perl
# Written: Nick Kent, 23rd Aug 2011
# Last updated: Mark Robinson, 2015
# USAGE:- perl Site_writer.pl [-options] <sgr_files> <site_files>
#
# Altered version of Nick's SiteWriter_Full.pl script to allow runnning on the command line and
# different normalisation/output options.
#
# Input: -comma seperated list of nucleosome maps (.sgr format)
#        -sites file (.txt format) containing: chrn\tfeature_name\tposition\torientation\n
#
# Output:
#
############################  Modules  ########################################

use strict;
use warnings;
use Getopt::Long;
use Cwd qw();
use File::Basename;
use Math::Round qw(:all);
use Data::Dumper;

############################  Settings  ########################################

my $bin = 5; # binning
my $outdir = Cwd::cwd();
my $window = 1000; # bp around each site to plot
my $scale = 1; # optional scaling of output values
my ($help);

my $help_string = "\nUsage: $0 [-options] <sgr_files> <site_files>\nOptional settings:
        --bin|-b = binning desired (default = 5bp)
        --window|-w = window to sum values within (default = 1000bp)
        --scale|-s = scaling of output data (default = 1.0)
        --out|-o = output directory (default = cwd)\n";

GetOptions ("bin=i" => \$bin,
            "out=s" => \$outdir,
            "window=i" => \$window,
            "help" => \$help,
            "scale=f" => \$scale)
or die ("\nError in command line arguments\n$help_string\n");

die $help_string if ($help);

my $sgr_list = shift @ARGV or die "Insufficient arguements: $help_string";
my $site_list = shift @ARGV or die "Insufficient arguements: $help_string";
my @sgr_files = split(',',$sgr_list);
my @site_files = split(',',$site_list);

################################################################################

print "\nStart:\t",`date`."\n";

my (%nuc_map,@conditions);

# store nuc maps
for my $sgr_file (@sgr_files) {
    
    my ($filename,$dir,$suffix) = fileparse($sgr_file,".sgr");
    die "$sgr_file is not an .sgr file \n" unless ($suffix);
    
    print "Reading in nucleosome map from: $filename\n";
    
    push(@conditions,$filename);
    
    open( my $in, '<', $sgr_file ) || die "Unable to open $sgr_file: $!\n";
    
    while(<$in>) {
        
        chomp;
        my ($chrn,$pos,$val) = split('\t');
        $nuc_map{$chrn}{$pos}{$filename} = $val;
    }
}

# loop through each site file in turn
for my $site_file (@site_files) {
    
    my (%cfd,%feature_count);
    
    my ($filename,$dir,$suffix) = fileparse($site_file,".txt");
    die "$site_file is not an .sgr file \n" unless ($suffix);
    
    my @strings = ("$filename\_sitewriter.txt","$filename\_individual_sites.txt");
    my %handles = filehandles(\@strings,\@conditions);
    
    # header for heatmapping file
    my @header = ( map {$_*$bin} (-$window/$bin)..($window/$bin) );
    for my $cond (@conditions) {
        my $fh = $handles{$cond}[1];
        print ($fh "gene_id".join("bp\t",@header)."\n");
    }
    
    print "Reading in sites from: $filename\n";
    
    open( my $in, '<', $site_file ) || die "Unable to open $site_file: $!\n";
    
    while(<$in>) {
        
        chomp;
        my ($chrn,$feature,$pos,$strand) = split('\t');
        $pos=nearest($bin,$pos);
        
        if (exists $nuc_map{$chrn}{$pos}) {
            
            my @range = (map {$_ * $bin} (($pos-$window)/$bin)..(($pos+$window)/$bin) );
            
            COND: for my $condition (sort keys %{$nuc_map{$chrn}{$pos}}) {
                
                my (@feature_sum);
                
                my $fh = $handles{$condition}[1];
                
                while (my ($index,$temp_pos) = each (@range)) {
                    
                    print "Feature $feature at $chrn $pos has missing values, skipped\n", next COND unless (exists $nuc_map{$chrn}{$temp_pos}{$condition});
                    $feature_sum[$index] += ($nuc_map{$chrn}{$temp_pos}{$condition} * $scale);
                    
                }
                
                $feature_count{$condition}++;
                
                @feature_sum = reverse @feature_sum if ($strand eq "R" || $strand eq "-");
                
                print ($fh "$feature\t".join("\t",@feature_sum)."\n");
                
                $cfd{$condition}[$_] += $feature_sum[$_] for (0..$#feature_sum);
            }
        }
    }
    
    for my $cond (sort keys %cfd) {
        
        my $fh = $handles{$cond}[0];
        
        print ($fh "Bin\tVal\n");
        
        while ( my ($index,$bin) = each (@header)) {
            
            my $val = $cfd{$cond}[$index]/$feature_count{$cond};
            
            print ($fh "$bin\t$val\n");
        }
    }
}

sub filehandles { # sub to create arrays of annon filehandles
    
    my ($strings,$conditions) = ($_[0],$_[1]);
    
    my %handles;
    
    for my $condition (@$conditions) {
        
        while (my ($index,$string) = each @$strings ) {
            
            my $fh = IO::File->new(">$outdir/$condition\_$string");
            $handles{$condition}[$index] = $fh;
            
        }
    }
    
    return %handles;
}

print "\nEnd:\t",`date`."\n";
exit;