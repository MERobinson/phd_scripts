#!/usr/bin/perl

# Written: Mark Robinson, April 2015
# Last updated: MER Jan 2016
# USAGE: PODS_2.0.pl -options [-c canonical_peaks] [-g gff_file] <list,of,SAM,files>
#
# Input:-Comma seperated list of SAM files from MNase-seq experiment
#       -[optional] Canonical peaks file (narrowPeak format)
#       -[optional] GFF file for gene coordinates
#       -If accession numbers used in alignment in stead of chr names change sub at end to match desired organism/chr numbers
#
# Output:-Info file (including read counts/coverage, peak numbers and average param values)
#        -Fragment size distribution (.txt format)
#	 -Canonical nucleosome positions (narrowPeak format)
#        -Individual + average nucleosome maps (.bed format)
#	 -Peak parameter files for each sample (.txt format)
#
# See guide for full details, briefly: nucleosomes are scored with a dynamic read extension from the dyad position,
# normalised to total read number and average maps generated for each condition. Consensus peaks are called via
# identification of peak shape in average nuc map, peak parameters are then quantified for each consensus peak from
# individual replicate maps. Parameters measured = position, occupancy, distribution and size.

################################ Load Modules ################################

use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
use POSIX;
use Getopt::Long;
use Cwd qw();
use File::Basename;
use IO::File;
use constant PI => 4 * atan2(1, 1);

################################ Checking Settings ################################

# default option settings
my $bin = 5; # binning to use for sgr file
my $min_size = 120; # minimum read length to map
my $max_size = 180; # maximum read length to map
my $window = 50; # window around peaks to plot distribution
my ($up_limit,$start_limit,$term_limit,$down_limit) = (1,3,3,1); # number of defined nuc category (upstream, at TSS, at TTS, downstream)
my $sigma = 3;
my $times_of_sigma = 3;
my $extend = 0.33; # read extension around dyad (if <= 1 read is multiplied by this factor, of > 1 read is extended by this many bp)
my $PE; # if data is paired end sequenced (default = SE)
my $chr_regex = 'chr[0-9A-za-z]+'; # regular expression string for searching for chromosome names (default = 'chr[a-zA-Z1-9]+')
my $outdir = Cwd::cwd();
my ($help,$print_transformed,$canonical_file,$gff_file);

# getting user specified options ### NEEDS UPDATING
GetOptions ("bin=i" => \$bin,
            "min=i" => \$min_size,
            "max=i"   => \$max_size,
            "canonical=s" => \$canonical_file,
            "gff=s" => \$gff_file,
            "extension=f"  => \$extend,
            "window=i" => \$window,
            "up_limit=i" => \$up_limit,
            "start_limit=i" => \$start_limit,
            "term_limit=i" => \$term_limit,
            "down_limit=i" => \$down_limit,
            "PE" => \$PE,
            "regex=s" => \$chr_regex,
            "help" => \$help,
            "transformed" => \$print_transformed,
            "out=s" => \$outdir)
or die("\nError in command line arguments\n");

my $help_string = "\nUsage: $0 -options [-c canonical_peaks] <list,of,SAM,files>\nOptional settings:
        --bin|-b = Binning desired (default = 5 bp)
        --min|-mi = Minimum fragment length to include (after rounding to closest bin) (default = 120 bp)
        --max|-ma = Maximum fragment length to include (after rounding to closest bin) (default = 180 bp)
        --canonical|-c = canonical nucleosome positions (narrowPeak format file)
        --gff|-g = gff file (default = none) 
        --window|-w = window around each peak to record score distribution (default = 50 bp)
        --up_limit|-u = number of defined nucleosome categories upstream of each gene to detect (default = 1)
        --start_limit|-s = number of defined nucleosome categories downstream of each TSS to detect (default = 3)
        --term_limit|-te = number of defined nucleosome upstream of each TTS to detect (default = 3)
        --down_limit|-d = number of defined nucleosome downstream of each gene to detect (default = 1)
        --extension|-e = Read extension length (if > 1 will be treated as set length, otherwise treated as fraction of read length) (default = 0.33)
        --PE|-p = Flag paired end reads (default = off)
        --regex|-r = Regular expression string for determining chromosomes included (default = chr[a-zA-Z0-9]+)
        --transformed|-tr = Print transformed average map (default=off)
        --out|-o = output directory (default = cwd)
        
        NB: unless a canonical peaks file is provided (with flag -c) one will be created from input files";

die $help_string if ($help);

die "\nIncorrect number of arguments supplied\n$help_string\n" if (@ARGV != 1);
my $sam_files = $ARGV[0] or die "\nArguments not supplied\n$help_string\n";

my $date = strftime "%d:%m:%Y", localtime;
my $date_time = localtime();

$| = 1; # wasn't burffering properly (not entirely sure why!)

################################ Global Variables ################################

# all global variables except those in settings
my ($rep_count,%peak_count,%nuc_map,%input_files,$info_out,@norm_array,%chrn_lengths,%peaks_hash,%gene_hash,%gene_ref);

################################ Check Input Files ################################

{
    print "\nStart: $date_time\n";

    # check number of conditions and reps, and filetype
    $rep_count = my @temp = split(',',$sam_files);

    # parse input files using sub (creates HoA with format: %input_files{file_no.(int)}[condition(string),file path(string)]
    my @suffixes = (".sam");
    %input_files = inputfiles($sam_files,\@suffixes);
    
    # open output info file
    my $infofile = "Run_info_".$0."_$date.txt";
    open($info_out, '>', "$outdir/$infofile") || die "Unable to open output file: $!\n";

    # print command to info file
    my $commandline = join " ", $0, @ARGV;
    print ($info_out "Command: $commandline\n");

}

############################## Read in Canonical Peaks #############################################

if ($canonical_file) {
    
    print "Reading in canonical peak positions\n";

    open (my $canonical_peaks,"<",$canonical_file) || die "Unable to open $canonical_file: $!\n";

    while(<$canonical_peaks>) {

        chomp;

        my ($chrn,$start,$end,$name,$score,$strand,$signal_val,$pval,$qval,$dyad_pos) = split("\t");

        next unless( $chrn =~ /$chr_regex/ );  

        # %peaks_hash{chr}{pos}[leading_edge,trailing_edge,[categories],[genes]]
        $dyad_pos = $start+$dyad_pos; # dyad pos is stored as offset from start, set back to actual coord
        $peaks_hash{$chrn}{$dyad_pos} = [$start,$end];
    }
}

############################## Read in GFF file #############################################

if ($gff_file) {
    
    print "Reading in gene coordinates from: $gff_file\n";

    my @suffixes = (".gff", ".gff3");
    my ($filename,$dir,$suffix) = fileparse($gff_file,".gff");
    die "$gff_file is not an accepted filetype file \n" unless ($suffix);
    
    open( my $in, '<', $gff_file ) || die "Unable to open $gff_file: $!\n";
    
    while(<$in>) {
        
        chomp;
        my ($chrn,$source,$feature,$start,$stop,$score,$strand,$frame,$attribute) = split('\t');
        $chrn = dicty_chr_names($chrn) if ($chrn =~ /^DDB/);
        
        next unless (defined $feature && $feature eq "gene");
        
        # calc length and round coord to nearest bins
        my $length = $stop-$start;
        ($start, $stop) = (nearest($bin,$start),nearest($bin,$stop));
        
        # split up attribute into gene name and id
        my $gene_id = $1 if ($attribute =~ /ID=([A-Za-z_0-9]+)\;/);
        my $gene_name = $1 if ($attribute =~ /Name=([A-Za-z_0-9]+)\;/);
        
        $gene_hash{$gene_id} = [$chrn,$start,$stop,$length,$strand,$gene_name];
        
        # create look-up hash for the genes mapping to each bin
        push(@{$gene_ref{$chrn}{$_}},$gene_id) for (map {$_*$bin} $start/$bin..$stop/$bin);
        
    }
}

############################## Generating Nucleosome Maps ##########################################

{
    
    print "Reading SAM files:\n";
    
    # use sub to get frag. dist. filehandles for individual rep files
    my @strings = "fragment_size_distribution.txt";
    my @indices = 1..$rep_count;
    my @handles = filehandles(\@strings,\%input_files,\@indices); 

    for my $file_no (sort keys %input_files) {
        
        my $condition = $input_files{$file_no}[0];
        
        print "\t$condition\n";
        
        print ($info_out "\n$condition:\n");
            
        open (my $file_in, '<', "$input_files{$file_no}[1]") || die "Unable to open $input_files{$file_no}[1]: $!";
        
        my ($unaligned_count,%size_hash,$read_count,$aligned_count,$mapped_count,$genome_size,$mapped_bp) = (0);

        LINE: while(<$file_in>) {
                
            chomp;

            # Initialise chromosome sizes from header
            if ($_ =~ /^@/) {
                
                my ($tag,$chr,$length) = split('\t');
                
                next LINE unless ($tag eq '@SQ');
                
                if ( $chr =~ /($chr_regex)/ ) {
                    
                    my $chrn = $1;
                    
                    $chrn = accession_chr($chrn) if ($chrn =~ /^NC/);
                    
                    next LINE unless ($chrn);
                    
                    if ($length =~ /^LN:(\d+)/) {
                        
                        my $len = ceil($1/$bin);
                        $chrn_lengths{$chrn} = $len;
                        $genome_size += $len;
                    
                    }
                }
                next LINE;
            }
            
            # check it recognised chromosomes
            die "Chromosome names not recognised - check chr regex\n" if (!defined $genome_size);
            
            # count all reads/bins depending on if sgr or sam
            $read_count++;
                
            my @line = split('\t');
            
            # Check for chromosome alignment and skip if unaligned
            my $chrn;
            if ($line[2] =~ /($chr_regex)/) { $chrn = $1; }
            $chrn = accession_chr($chrn) if ($chrn && $chrn =~ /^NC/);
            
            $unaligned_count++, next LINE if ($line[2] =~ /\*/);
                
            next LINE if (!defined $chrn);
                
            # set fragment length depending on whether PE or SE
            my $frag_size = ($PE) ? $line[8] : length($line[9]);
            my $adj_frag_size = nearest($bin,$frag_size);
                
            # increment count of size and mapped reads
            $size_hash{abs($adj_frag_size)}++;
            $aligned_count++;
                
            # Check fragment length is within size window for frag sizes window
            if ($adj_frag_size >= $min_size && $adj_frag_size <= $max_size) {
                
                # Determine dyad position
                my $dyad = ($line[3] + ($frag_size * 0.5));
                
                # Increment count at all bins within $extend bp of central dyad (rounded) (if within nuc size range)
                my $read_extension = ($extend > 1) ? $extend : ($frag_size * $extend);
                my $start = nearest( $bin, ($dyad - ($read_extension/2)) );
                my $end = nearest( $bin, ($dyad + ($read_extension/2)) );
                $nuc_map{$chrn}{$_}[$file_no][0]++ for ( map { $_ * $bin} ( ($start/$bin)..($end/$bin) ) );
                
                # Add size of current fragment size to sum of sizes for all bins covered by this read
                $nuc_map{$chrn}{$_}[$file_no][1]+= $frag_size for ( map {$_ * $bin} ( ($start/$bin)..($end/$bin) ) );
                
                # Increment count of reads within window
                $mapped_count++;
                $mapped_bp += ($end-$start)/$bin;
            }
        }
        
        # print some info (incl size distribution) to info outfile
        $genome_size *= $bin;
        my $coverage = $mapped_bp/$genome_size;
        print ($info_out "\tMapped genome size: $genome_size\n");
        print ($info_out "\tTotal reads input: $read_count\n\tReads mapping to the genome: $aligned_count\n");
        print ($info_out "\tUnaligned reads: $unaligned_count\n\tReads within specified size range: $mapped_count\n\tAverage per bin coverage: $coverage\n\n");
        
        # print frag size distribution
        my $fh = $handles[$file_no][0];
        print ($fh "Size(bp)\tRaw_count\tFrequency\n"); # header
        
        for my $size (sort {$a <=> $b} keys %size_hash) {
                
            my $count = $size_hash{$size};
            my $freq = nearest(0.00001,$count/$aligned_count);
            print ($fh "$size\t$count\t$freq\n");
            
        }
        
        # normalisation factor
        my $norm = $genome_size/$mapped_bp;
        $norm_array[$file_no] = $norm;

    }
}

######################## Normalise & Print Maps and Sum for Average ##################################

{
    
    my @strings = ("$min_size\-$max_size\_e$extend\_b$bin\_nuc_map.sgr");
    my @indices = 1..$rep_count;
    my @handles = filehandles(\@strings,\%input_files,\@indices);

    print "\nNormalising and printing nucleosome maps\n";
    
    for my $chrn (sort keys %nuc_map) {

        for my $pos ( map {$_ * $bin} (0..$chrn_lengths{$chrn}) ) {
            
            for my $file_no ( sort {$a <=> $b} keys %input_files ) {

                for my $val_no (reverse (0..1)) { # 0 = nucleosome occupancy map, 1 = nucleosome size map (need to do size first before changing read numbers)
                
                    my $fh = $handles[$file_no][0];
                    my $score = 0;

                    if ( exists $nuc_map{$chrn}{$pos}[$file_no][0] ) {
                        
                        # set normalisation factor depending on whether were dealing with size or occupancy
                        my $norm = ($val_no == 0) ? $norm_array[$file_no] : 1/$nuc_map{$chrn}{$pos}[$file_no][0];
                        
                        # normalise current value
                        $score = nearest(0.01,$nuc_map{$chrn}{$pos}[$file_no][$val_no]*$norm);
                    
                    } 

                    # replace un-normalised values in maps
                    $nuc_map{$chrn}{$pos}[$file_no][$val_no] = $score;

                    if ($val_no == 0) {
                        # print individual files
                        print ($fh "$chrn\t$pos\t$score\n");
                        # sum for average map
                        $nuc_map{$chrn}{$pos}[0][0] += $score unless($canonical_file);
                    }
                }
            }
        }
    }
    
}

# Add average conditions to list of conditions
$input_files{0} = ["canonical_nuc_map"] unless($canonical_file);

################################## Print Canonical Map ###########################################

# Print average maps per condition
unless($canonical_file) {
    
    print "Generating consensus nucleosome map\n";

    open ( my $outfile, ">", "$outdir/Averaged_canonical_nuc_map.sgr")
        || die "Unable to open output file: $!\n";

    for my $chrn (sort keys %nuc_map) {

        for my $pos ( sort {$a<=>$b} keys %{$nuc_map{$chrn}} ) {
            
            my $score = nearest(0.01,$nuc_map{$chrn}{$pos}[0][0]/$rep_count); # normalise
            $nuc_map{$chrn}{$pos}[0][0] = $score; # replace values with normalised values
            print ($outfile "$chrn\t$pos\t$score\n");

        }
    }
}

############################ Transform Data ##################################

## Perform first derivative of Gaussian and Laplacian of Gaussian convolutions
{
    my (@FDoG,@LoG,@current);

    # calculate weights
    for my $k ( -$sigma*$times_of_sigma..$sigma*$times_of_sigma ) {
        
        my $weight = exp( -( ($k**2)/(2*$sigma**2) ) ); # Gaussian
        push( @FDoG, ( -($k/$sigma**2) ) * $weight ); # First derivative of Gaussian
        push( @LoG, ((($k**2)/($sigma**4)) - (1/($sigma**2))) * $weight ); # Laplacian of Gaussian
        
    }

    my $factor = ( 1/sqrt(2*PI*$sigma*$sigma) );

    print "\n";

    # perform convolution on each bin
    for my $chrn (sort keys %nuc_map) {
        
        push (@current,"$chrn");
        print "Performing convolutions: @current\r";
        
        for my $x (  map {$_ * $bin} (0..$chrn_lengths{$chrn}) ) {
            
            my $x_index = $x/$bin;
            
            for my $file_no (sort keys %input_files) {
                
                my $fx = $nuc_map{$chrn}{$x}[$file_no][0];
                my ($dg,$lg);
                
                if ( $x_index >= $sigma*$times_of_sigma && $x_index <= $chrn_lengths{$chrn}-($sigma*$times_of_sigma) ) {
                    
                    my $weight_index = 0;
                    
                    for my $k ( map { $_ * $bin } ($x_index-$sigma*$times_of_sigma..$x_index+$sigma*$times_of_sigma) ) {
                        
                        my $val = $nuc_map{$chrn}{$k}[$file_no][0] || 0;
                        $dg += $val * $FDoG[$weight_index];
                        $lg += $val * $LoG[$weight_index];
                        $weight_index++;
                        
                    }
                    
                    $dg *= $factor;
                    $lg *= $factor;
                    
                    # store transformed data
                    my $size = $nuc_map{$chrn}{$x}[$file_no][1];
                    $nuc_map{$chrn}{$x}[$file_no] = [$fx,$size,$dg,$lg];
                    
                } else {
                    
                    # NB: this is a very lazy way of dealing with edges - just throws data away,
                    #     isn't a problem for standard nuc maps given unreliable data at chr ends anyway
                    #     but might want to change if doing something with other size classes
                    $nuc_map{$chrn}{$x}[$file_no] = [0,0,0,0];
                }
            }
        }
    }
    print "\n";
}

## Print transformed maps if $print_transformed flagged
if ($print_transformed) {
    
    print "\nPrinting transformed maps\n";
    
    for my $file_no (sort keys %input_files) {
        
        my $condition = $input_files{$file_no}[1];
        
        my $FDoG_outfile = $condition."_FDoG.sgr";
        open (my $FDoG_out, '>', "$outdir/$FDoG_outfile") || die "Unable to open output file $FDoG_outfile\n";
        
        my $LoG_outfile = $condition."_LoG.sgr";
        open (my $LoG_out, '>', "$outdir/$LoG_outfile") || die "Unable to open output file $LoG_outfile\n";
        
        for my $chrn (sort keys %nuc_map) {
            
            for my $pos ( sort {$a <=> $b} keys %{$nuc_map{$chrn}} ) {
                
                print ($FDoG_out "$chrn\t$pos\t$nuc_map{$chrn}{$pos}[$file_no][2]\n");
                print ($LoG_out "$chrn\t$pos\t$nuc_map{$chrn}{$pos}[$file_no][3]\n");
                
            }
        }
    }
}

############################ Identify Consensus Peaks ##################################

# Loop through transformed WT av. data and identify peaks using inflection points

unless($canonical_file) {
    
    print "\nDetecting consensus nucleosome positions\n";
    
    open ( my $sgrfile, ">", "$outdir/Canonical_nucleosome_positions.sgr") 
        || die "Unable to open output file: $!\n";

    open ( my $narrowPeak, ">", "$outdir/Canonical_nucleosome_positions.narrowPeak") 
        || die "Unable to open output file: $!\n";

    my ($peaks_count);

    for my $chrn (sort keys %nuc_map) {
        
        my ($start,$maxima,$end,$minima,$read_sum,$width);
        
        BIN: for my $pos ( sort {$a <=> $b} keys %{$nuc_map{$chrn}} ) {
                
            next unless (exists $nuc_map{$chrn}{$pos-(3*$bin)}[0][0] && exists $nuc_map{$chrn}{$pos+(3*$bin)}[0][0]); # skip chromosome ends
            
            unless (defined $nuc_map{$chrn}{$pos}[0][2] && defined $nuc_map{$chrn}{$pos-$bin}[0][2]) {
                print "Chrn: $chrn, Pos: $pos, Values: ".$nuc_map{$chrn}{$pos}[0][2].", ".$nuc_map{$chrn}{$pos-$bin}[0][2]."\n";
            }

            # identify min-extremum point (neg iflection in FDoG) = start/end of potential peak area
            if ( $nuc_map{$chrn}{$pos}[0][2] <= 0 && $nuc_map{$chrn}{$pos-$bin}[0][2] > 0 ) {
                
                # if every other peak point also identified, store as a peak
                if ( defined $start && defined $end && defined $maxima && defined $minima ) {
                    
                    my $av_occ = $read_sum/$width; # average coverage per bin for peak region
                    
                    if ( $av_occ > 2 && $av_occ < 50 ) { # eliminate any peaks that are <2x or >50x average coverage
                        
                        $peaks_count++;

                        # %peaks_hash{chr}{pos}[leading_edge,trailing_edge,[categories],[genes]]
                        $peaks_hash{$chrn}{$maxima} = [$start,$end,[],[]];
                        
                        # print peaks in .sgr format for visualisation
                        for my $coord ( map {$_ * $bin} (($start/$bin)..($end/$bin)) ) {
                            
                            my $val = ($coord == $maxima) ? $nuc_map{$chrn}{$maxima}[0][0] : 1;
                            print ($sgrfile "$chrn\t$coord\t$val\n");

                        }

                        # print narrowpeak file
                        my $point_source = $maxima - $start;
                        print ($narrowPeak "$chrn\t$start\t$end\tCN$peaks_count\t$av_occ\t.\t$read_sum\t-1\t-1\t$point_source\n");
                    }
                }
                
                # reset other points and continue
                ($minima,$read_sum,$width) = (1,0,0);
                $start = undef, $end = undef, $maxima = undef, next;
                
            }
            
            # identify leading peak edge (neg inflection in LoG) (only take first possible site in current candidate region)
            if ( $nuc_map{$chrn}{$pos}[0][3] <= 0 && $nuc_map{$chrn}{$pos-$bin}[0][3] > 0 && defined $minima  && !defined $start ) {
                
                $start = $pos, next BIN;
            }
            
            # keep track of total occupancy
            $read_sum += $nuc_map{$chrn}{$pos}[0][0], $width++ if (defined $start);
            
            # find max extremum point (positive inflection in FDoG)
            if ( $nuc_map{$chrn}{$pos}[0][2] >= 0 && $nuc_map{$chrn}{$pos-$bin}[0][2] < 0 && defined $start && defined $minima ) {
                
                $maxima = $pos, next BIN;
            }
            
            # identify trailing edge (pos inflection in LoG) (reset for each possible site in current candidate region)
            $end = $pos, next BIN if ( $nuc_map{$chrn}{$pos}[0][3] >= 0 && $nuc_map{$chrn}{$pos-$bin}[0][3] < 0 && defined $start);

        }
    }
}

############################## Determine Nucleosome Categories ######################################

if ($gff_file) {

    my (%gene_nuc_count);

    print "Determining nucleosome categories\n";

    my ($prev_chrn,$prev_nuc,$prev_gene);

    for my $chrn (keys %peaks_hash) {
        
        NUC: for my $nuc (sort {$a <=> $b} keys %{$peaks_hash{$chrn}}) {
            
            # skip and delete first nucleosome on each chromosome
            if (!defined $prev_chrn || $chrn ne $prev_chrn) {
                ($prev_nuc,$prev_chrn) = ($nuc,$chrn);
                delete($peaks_hash{$chrn}{$nuc}), next NUC;
            }
            
            # check if nucleosome is genic
            if (exists $gene_ref{$chrn}{$nuc}) {
                
                # loop through assocaited genes (in case of overlapping)
                for my $gene_id (@{$gene_ref{$chrn}{$nuc}}) {
                    
                    # push/unshift nucleosome position onto list of genic nucleosomes for current gene
                    my $strand = $gene_hash{$gene_id}[4];
                    push(@{$gene_hash{$gene_id}[7]}, $nuc) if ($strand eq "+");
                    unshift(@{$gene_hash{$gene_id}[7]}, $nuc) if ($strand eq "-");
                    
                    # track current nucleosome, gene and chrn for next iteration
                    ($prev_chrn,$prev_nuc,$prev_gene) = ($chrn,$nuc,$gene_id);
                }
                
                next NUC;
            
            } else {
                
                my $mapped_flag;
                
                # If within extension limit of the last gene seen
                if ( $prev_gene && ( ($nuc - $gene_hash{$prev_gene}[2]) < 1000 ) && ($chrn eq $gene_hash{$prev_gene}[0]) ) {
                    
                    # push/unshift nucleosome position onto list of upstream/downstream nucleosomes for current gene
                    my $strand = $gene_hash{$prev_gene}[4];
                    push(@{$gene_hash{$prev_gene}[8]}, $nuc) if ($strand eq "+");
                    unshift(@{$gene_hash{$prev_gene}[6]}, $nuc) if ($strand eq "-");
                    
                    $mapped_flag++;
                }
                
                # Then check if within extension limits of next gene
                WALK: for (my $index = $nuc; $index <= ($nuc+1000); $index += $bin) {
                    
                    if (exists $gene_ref{$chrn}{$index}) {
                        
                        my $gene_id = $gene_ref{$chrn}{$index}[0];
                        
                        my $strand = $gene_hash{$gene_id}[4];
                        push(@{$gene_hash{$gene_id}[6]}, $nuc) if ($strand eq "+");
                        unshift(@{$gene_hash{$gene_id}[8]}, $nuc) if ($strand eq "-");
                        
                        $mapped_flag++, last WALK;
                    }
                }
            }
            ($prev_nuc,$prev_chrn) = ($nuc,$chrn);
        }
    }

    # loop through each gene and determine nuc categories
    for my $gene_id (sort keys %gene_hash) {
        
        my ($chrn,$start,$stop,$length,$strand,$name,@mapped_nuc) = @{$gene_hash{$gene_id}};
        
        ### Upstream ###
        if (defined $mapped_nuc[0]) {
            
            while ( my ($index, $nucleosome) = each(@{$mapped_nuc[0]}) ) {
                
                # determine nucleosome category based on desired number of positions and number of upstream nuc
                my $nuc_category = ( $index >= (@{$mapped_nuc[0]}-$up_limit) ) ? ($up_limit-($#{$mapped_nuc[0]}-$index)) : 0;
                
                # %peaks_hash{chr}{pos}[leading_edge,trailing_edge,[categories],[genes]]
                push( @{$peaks_hash{$chrn}{$nucleosome}[2]}, $nuc_category );
                push( @{$peaks_hash{$chrn}{$nucleosome}[3]}, $gene_id );
            }
        }
        
        ### Genic ###
        if (defined $mapped_nuc[1]) {
            
            # to deal with genes with less than specified number of nucleosomes, need to change limits to be the smaller out of 1/2 nuc or the limit
            my $nuc_count = @{$mapped_nuc[1]};
            my $upper_limit = ( ceil($nuc_count/2)>$start_limit ) ? $start_limit : ceil($nuc_count/2);
            my $lower_limit = ( floor($nuc_count/2)>$term_limit ) ? $term_limit : floor($nuc_count/2);
            
            while ( my ($index, $nucleosome) = each (@{$mapped_nuc[1]}) ) {
                
                # set nucleosome category
                my $nuc_category;
                $nuc_category = ($index >= $upper_limit) ? ($up_limit+$start_limit+1) : ($up_limit+1+$index);
                $nuc_category = ((($nuc_count-1) - $lower_limit) < $index) ? ($up_limit+$start_limit+1+($index-(($nuc_count-1)-$lower_limit))) : $nuc_category;
                
                # %peaks_hash{chr}{pos}[leading_edge,trailing_edge,[categories],[genes]]
                push( @{$peaks_hash{$chrn}{$nucleosome}[2]}, $nuc_category );
                push( @{$peaks_hash{$chrn}{$nucleosome}[3]}, $gene_id );
            }
        } 
            
        
        ### Downstream ###
        if (defined $mapped_nuc[2]) {
            
            while ( my ($index, $nucleosome) = each (@{$mapped_nuc[2]}) ) {
                
                # set nucleosome category
                my $nuc_category = ($index >= $down_limit) ? ($up_limit+$start_limit+1+$term_limit+$down_limit+1) : ($up_limit+$start_limit+1+$term_limit+$index+1);
                
                # %peaks_hash{chr}{pos}[leading_edge,trailing_edge,[categories],[genes]]
                push( @{$peaks_hash{$chrn}{$nucleosome}[2]}, $nuc_category );
                push( @{$peaks_hash{$chrn}{$nucleosome}[3]}, $gene_id );
            }
        } 
    }
}

############################ Identify Replicate Peaks and Measure Peak Parameters ##################################

{
    
    print "Measuring peak parameters from replicates\n";
    
    my @conditions = 1..$rep_count;
    my @strings = ("peak_positions_$date.sgr","peak_parameters_$date.txt");
    my @handles = filehandles(\@strings,\%input_files,\@conditions);
    
    COND: for my $file_no (@conditions) {
        
        # print param file header
        my $param_fh = $handles[$file_no][1];
        my $sgr_fh = $handles[$file_no][0];
        my @header = map {$_ * $bin} ( -($window/$bin)..($window/$bin) );
        print ($param_fh "chr\tstart\tend\tdyad\tcanonical\tpos\tocc\tdist\tsize\twidth\tgene_id\tnuc_category\t".join("\t",@header)."\n");
        
        for my $chrn (sort keys %peaks_hash) {  

            # loop through reference peak positions for this condition
            CANONICAL: for my $canonical_pos (sort {$a <=> $b} keys %{$peaks_hash{$chrn}}) {
                
                next unless(defined $nuc_map{$chrn}{$canonical_pos}[$file_no][0]);

                # %peaks_hash array layout: [leading_edge,trailing_edge,]
                my $canonical_start = $peaks_hash{$chrn}{$canonical_pos}[0];
                my $canonical_end = $peaks_hash{$chrn}{$canonical_pos}[1];
                
                # loop through this peak region and find local max extremum points in FDoG
                my @maximas;
                for my $current_pos ( map {$_ * $bin} ( ($canonical_start/$bin)..($canonical_end/$bin) ) ) {
                    
                    next CANONICAL unless (exists $nuc_map{$chrn}{$current_pos-$bin}[$file_no][2] && exists $nuc_map{$chrn}{$current_pos}[$file_no][2]);
                    push @maximas, $current_pos if ($nuc_map{$chrn}{$current_pos}[$file_no][2] >= 0 && $nuc_map{$chrn}{$current_pos-$bin}[$file_no][2] < 0);
                }
                
                # Find closest maxima to av dyad
                my $pos_diff = 200; # maximum cut-off for pos change
                my $maxima;
                for my $temp_maxima (@maximas) {
                    if (abs($temp_maxima-$canonical_pos) < $pos_diff) {
                        ($maxima,$pos_diff) = ($temp_maxima,($temp_maxima-$canonical_pos))
                    }
                }

                # initialise param variables
                my ($pos,$occ,@cumulative,$grad,$size,$width,$leading,$trailing);

                # If failed to find peak in rep, use canonical co-ord to measure occ and size (skip pos)
                if (abs($pos_diff) > 199) {
                    ($pos,$grad,$maxima) = ('','',$canonical_pos);
                } else {
                    LEADING: for my $temp_pos ( map {$_ * $bin} (($canonical_start-50)/$bin..$maxima/$bin) ) {
                        $leading = $temp_pos if ( $nuc_map{$chrn}{$temp_pos}[$file_no][3] <= 0 && $nuc_map{$chrn}{$temp_pos-$bin}[$file_no][3] > 0 );
                    }
                    TRAILING: for my $temp_pos ( map {$_ * $bin} ($maxima/$bin..($canonical_end+50)/$bin) ) {
                        $trailing = $temp_pos if ( $nuc_map{$chrn}{$temp_pos}[$file_no][3] >= 0 && $nuc_map{$chrn}{$temp_pos-$bin}[$file_no][3] < 0 );
                    }
                    $pos = $pos_diff;
                }

                $leading = $canonical_start unless ($leading);
                $trailing = $canonical_end unless ($trailing);

                # Get score distribution around peak
                my @score_dist;
                my @temp_array = map {$_ * $bin} ( (($maxima-$window)/$bin)..(($maxima+$window)/$bin) );
                while ( my ($index, $temp_pos) = each (@temp_array) ) {
                    $score_dist[$index] = $nuc_map{$chrn}{$temp_pos}[$file_no][0] // '';
                } 

                # Measure parameters
                @temp_array = map {$_ * $bin} ( ($leading/$bin)..($trailing/$bin) );
                while ( my ($index, $temp_pos) = each ( @temp_array ) ) {
                    
                    unless (defined $nuc_map{$chrn}{$temp_pos}[$file_no][0]) {
                        print "Chr: $chrn, Pos: $temp_pos, File: $file_no\n"
                    }

                    $occ += $nuc_map{$chrn}{$temp_pos}[$file_no][0];
                    $cumulative[$index] = $occ; 
                    $width += $bin;
                    $size += $nuc_map{$chrn}{$temp_pos}[$file_no][1];
                    
                }                

                # Average/normalise/etc.
                $size /= ($width/$bin);
                
                # calculate gradient
                if (defined $pos && $cumulative[-1] > 0) {
                    
                    $cumulative[$_] /= $cumulative[-1] for (0..$#cumulative);
                    my $summit_index = ($maxima-$leading)/$bin;
                
                    # gradient calc works better over 5 bins but if peak too narrow, also try 3
                    if ( $summit_index >= 3 && ( (($width/$bin)-$summit_index) >= 3 ) ) {

                        $grad = nearest(0.001, ( (($cumulative[$summit_index+2]-$cumulative[$summit_index-2])/5) *100) );
                        
                    } elsif ( $summit_index >= 2 && ( (($width/$bin)-$summit_index) >= 2 ) ) { 
                
                        $grad = nearest(0.001, ( (($cumulative[$summit_index+1]-$cumulative[$summit_index-1])/3) *100) );

                    } 
                }
                $grad = $grad // '';

                # Print parameters for individual rep file
                # %peaks_hash{chr}{pos}[leading_edge,trailing_edge,[categories],[genes]]
                for my $index (0..$#{$peaks_hash{$chrn}{$canonical_pos}[2]}) {

                    my $nuc_category = $peaks_hash{$chrn}{$canonical_pos}[2][$index] // '';
                    my $gene_id = $peaks_hash{$chrn}{$canonical_pos}[3][$index] // '';
                    # if position change and gene defined, check orientation of gene and adjust change
                    if ($gene_id && $pos) {
                        #$gene_hash{$gene_id} = [$chrn,$start,$stop,$length,$strand,$gene_name]
                        $pos = -$pos if ($gene_hash{$gene_id}[4] eq "-");
                    }
                    print ($param_fh "$chrn\t$leading\t$trailing\t$maxima\t$canonical_pos\t$pos\t$occ\t$grad\t$size\t$width\t$gene_id\t$nuc_category\t".join("\t",@score_dist)."\n");
                }
                
                # print peaks in .sgr format for visualisation
                for my $coord ( map {$_ * $bin} (($leading/$bin)..($trailing/$bin)) ) {
                    
                    my $val = ($coord == $maxima) ? $nuc_map{$chrn}{$maxima}[$file_no][0] : 1;
                    print ($sgr_fh "$chrn\t$coord\t$val\n");

                }

                # Count peaks
                $peak_count{$file_no}++;
            }
        }
    }
}

print ($info_out "\nReplicate peaks identified matching canonical:\n");
for my $file_no (sort keys %peak_count) {
    next if ($file_no < 2);
    my $condition = $input_files{$file_no}[1];
    print ($info_out "\t$condition = $peak_count{$file_no}\n");
}

print "DONE!\n\nEnd: $date_time\n";
exit;

############################# Subroutines #############################

# sub to parse files and return HoA with basename and file path
sub inputfiles { 
    
    my ($file_list,$suffixes,$sample_count) = ($_[0],$_[1]);

    my @files = split(',',$file_list);

    for my $file (@files) {
            
        my ($condition,$dir,$suffix) = fileparse($file,@$suffixes);
        die "Input file $condition isn't an accepted filetype: $!\n" unless $suffix;
        $sample_count++;
        $input_files{$sample_count} = [$condition,$file];

    }

    return %input_files;
}

# sub to create arrays of filehandles (expects a hash of files like from inputfiles function)
# third argument to specify which files in file hash to include
sub filehandles { 
    
    my ($strings,$file_hash,$indices,@handles) = ($_[0],$_[1],$_[2]);
    my %conditions = map { $_ => 1 } @$indices;

    FILE: for my $file_no (sort {$a <=> $b} keys %$file_hash) {
        
        next FILE unless (exists($conditions{$file_no}));

        my $condition = $$file_hash{$file_no}[0];
        
        while (my ($index,$string) = each @$strings ) {
            
            my $fh = IO::File->new(">$outdir/$condition\_$string");
            $handles[$file_no][$index] = $fh;
            
        }
    }
    
    return @handles;
}

# sub to switch chr accession no. to names 
sub accesion_chr {
    
    my $id = $_[0];
    my $chrn;
    
    # yeast
    $chrn = "chr1" if ($id =~ /^NC_001133/);
    $chrn = "chr2" if ($id =~ /^NC_001134/);
    $chrn = "chr3" if ($id =~ /^NC_001135/);
    $chrn = "chr4" if ($id =~ /^NC_001136/);
    $chrn = "chr5" if ($id =~ /^NC_001137/);
    $chrn = "chr6" if ($id =~ /^NC_001138/);
    $chrn = "chr7" if ($id =~ /^NC_001139/);
    $chrn = "chr8" if ($id =~ /^NC_001140/);
    $chrn = "chr9" if ($id =~ /^NC_001141/);
    $chrn = "chr10" if ($id =~ /^NC_001142/);
    $chrn = "chr11" if ($id =~ /^NC_001143/);
    $chrn = "chr12" if ($id =~ /^NC_001144/);
    $chrn = "chr13" if ($id =~ /^NC_001145/);
    $chrn = "chr14" if ($id =~ /^NC_001146/);
    $chrn = "chr15" if ($id =~ /^NC_001147/);
    $chrn = "chr16" if ($id =~ /^NC_001148/);
    $chrn = "chrM" if ($id =~ /^NC_001224/);
    
    return $chrn;
}

sub dicty_chr_names {
    
    my $chrn = $_[0];
    
    $chrn = "chr1" if ($chrn=~ /^DDB0232428$/);
    $chrn = "chr2" if ($chrn=~ /^DDB0232429$/);
    $chrn = "chr3" if ($chrn=~ /^DDB0232430$/);
    $chrn = "chr4" if ($chrn=~ /^DDB0232431$/);
    $chrn = "chr5" if ($chrn=~ /^DDB0232432$/);
    $chrn = "chr6" if ($chrn=~ /^DDB0232433$/);
    $chrn = "chr3F" if ($chrn=~ /^DDB0215151$/);
    $chrn = "chr2F" if ($chrn=~ /^DDB0215018$/);
    $chrn = "chrBF" if ($chrn=~ /^DDB0220052$/);
    $chrn = "chrR" if ($chrn=~ /^DDB0237465$/);
    $chrn = "chrM" if ($chrn=~ /^DDB0169550$/);
    
    return $chrn;
}
