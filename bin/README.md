************************************************************

  $HeadURL$ 
  $Author$ 
  $Revision$
  $Date $ 

************************************************************

This is a set of scripts for analysis of metagenomic samples. 

  * count_reads.pl - counts no. of read pairs in .fq.gz files
  * qc_trim.pl - quality trims reads and runs fastqc on trimmed data
  * human_filter.pl - Removes reads mapping to the human genome from the data


Overall Workflow for carrying out metagenomic analysis
======================================================

1. Quality trimming/fastqc analysis 
-----------------------------------

The bin/qc_trim.pl script needs to be modified to set the '$orig_dir' variable
to point to a directory containing demultiplexed, paired reads. The outputs
will be written to the location specified by '$out_dir'. It may be necessary
to modify the adapter sequences in %barcodes according to sample preparation
methodology.

It will also be necessary to adjust the task range specified at the top
of the script (#$ -t 1-73) to reflect the number of samples in the analysis.

Submit the script using 'qsub bin/qc_trim.pl'.

$out_dir will contain one directory per-sample after the run, containing
trimmed, gzipped fastq files, and fastqc outputs.

2. Filter out human reads
------------------------- 

The bin/human_filter.pl script will carry out an alignment using bwa mem vs an appropriately
indexed human genome database. Read-pairs which are unmapped are then extracted, and resplit into
fastq files. This stage also removes the barcode tags from the sample names since they are no
longer required.

The following may require editing:

    # directory containing qc_trimmed output from stage 1
    my $trimmed_data = "/data/florinash/trimmed"; 
    # directory to write resulting fastq files 
    my $out_dir      = "/data/florinash/human_alignments_test/";
    # indexed reference database to align against
    my $ref_db       = "/reference_data/Homo_sapiens/v38_p0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";

Again, the '#$ -t' line will need altering to reflect sample numbers. Submit
the script to the cluster using:

qsub bin/human_filter.pl


