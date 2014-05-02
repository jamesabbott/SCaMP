#!/usr/bin/perl

############################################################
#
# $HeadURL$
# $Author$
# $Revision$
# $Date$
#
############################################################

# Script to align paired metagenomic reads against the human genome, and remove the reads which map
# successfully. If either one read or the other of a pair maps, we will discard the read-pair, only
# retaining those for which both reads are unmapped

#$ -j y
#$ -cwd
#$ -m a
#$ -pe smp 8
#$ -R y
#$ -l h_rt=08:00:00,h_vmem=2G
#$ -t 1-73

use warnings;
use strict;

use File::Copy;

{
    if ( !$ENV{'SGE_TASK_ID'} ) {
        die "This script must be run as an array job...";
    }

    my $task = $ENV{'SGE_TASK_ID'} - 1;
    print "Exec host is " . $ENV{'HOSTNAME'} . "\n";

    # directory containing qc_trimmed output from stage 1
    my $trimmed_data = "/data/florinash/trimmed";
    # directory to write resulting bam files to
    my $out_dir      = "/data/florinash/human_filtered_reads/";
    # indexed reference database to align against
    my $ref_db       = "/reference_data/Homo_sapiens/v38_p0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";

    opendir DIR, $trimmed_data or die "Error opening $trimmed_data:$!";
    my @samples = grep !/\.\.?\z|stdout/, readdir DIR;
    close DIR;

    die "Incorrect number of samples ... task id = $task;
        sample number = $#samples \n " if ( $#samples < $task );

    my $sample = $samples[$task];
    chomp $sample;
    print "Sample = $sample\n";

    # get the filenames of the paired fastq files
    opendir SAMPLE, "$trimmed_data/$sample/" or die "Could not open $trimmed_data/$sample: $!";
    my $pair1 = ( grep ( /_1.fq.gz/, readdir SAMPLE ) )[0];
    rewinddir SAMPLE;
    my $pair2 = ( grep ( /_2.fq.gz/, readdir SAMPLE ) )[0];
    close SAMPLE;

    # copy fastq files for sample to local storage...
    mkdir '/local_scratch/florinash/'
      or die "Could not create /local_scratch/florinash: $!"
      unless ( -d '/local_scratch/florinash' );

    my $scratch_dir = "/local_scratch/florinash/$sample";
    if ( !-d $scratch_dir ) {
        mkdir $scratch_dir or die "mkdir failed: $!";
    }

    copy( "$trimmed_data/$sample/$pair1", "$scratch_dir" ) or die "$pair1 Sample copy failed: $!";
    copy( "$trimmed_data/$sample/$pair2", "$scratch_dir" ) or die "$pair2 Sample copy failed: $!";

    my $cmd =
"/usr/biosoft/cluster/bwa/current/bwa mem -t 8 -M $ref_db $scratch_dir/$pair1 $scratch_dir/$pair2 | /usr/biosoft/cluster/samtools/current/samtools view -Su - > $scratch_dir/$sample.human.bam";
    system($cmd) == 0 or die "Error executing bwa: $!";

    # SAM flag 0x0004 indicates the read is unmapped, while 0x0008 means the mate is unmapped,
    # therefore we want those with 0x0012 (12), where both apply...
    $cmd = "/usr/biosoft/cluster/samtools/current/samtools view -b -f 12 -F 256 $scratch_dir/$sample.human.bam > $scratch_dir/$sample.unaligned.bam";
    system($cmd) == 0 or die "Error executing bwa: $!";

    # Modify sample id's to remove barcodes, since we won't require these anymore...
    my $orig_sample = $sample;
    $sample=~s/_[ACTG]+$//;

    $cmd =
     "/usr/bin/java  -Xmx4G -jar /usr/biosoft/cluster/picard/current/SamToFastq.jar " . 
 	"INPUT=$scratch_dir/$orig_sample.unaligned.bam FASTQ=$scratch_dir/$sample". "_1.fq SECOND_END_FASTQ=$scratch_dir/$sample" . "_2.fq";
     system($cmd) == 0 or die "Error executing command: $!";

    if ( !-d "$out_dir/$sample" ) {
        mkdir "$out_dir/$sample" or die $!;
    }

    foreach my $fastq ($sample . ".hsfilt_1.fq", $sample . ".hsfilt_2.fq") {
	my $cmd = "gzip $scratch_dir/$fastq";
	system($cmd)==0 or die "Error executing command: $cmd";
	copy( "$scratch_dir/$fastq.gz", "$out_dir/$sample/$fastq.gz" )
	  or die "Error copying $fastq.gz: $!";
    }
    `rm -rf $scratch_dir`;
}

