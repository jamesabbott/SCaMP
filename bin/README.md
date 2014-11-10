<!---
************************************************************

  $HeadURL$ 
  $Author$ 
  $Revision$
  $Date $ 

************************************************************
-->
# Florinash Metagenome Pipeline

This is a set of scripts for analysis of metagenomic samples. 

## Reference Databases

Each database used by for filtering by these scripts requires the following:

1. A fasta formatted sequence of sequences, named 'database.fa'
2. A fasta index,generatd with 'samtools faidx' ('database.fa.fai')
2. BWA indices appropriate for searching with BWA MEM
3. A tab-delimited file ('database.tax.dat') relating taxonomy data to each sequence:
   The format should be as follows:

Accession   Scientific name  Strain  NCBI TaxId
i.e.

`chr1    Homo sapiens            9606`

This example has no strain defined, so this field is simply left blank

###Building Databases

A series of standard databases can be downloaded and formatted using the
`build_ref_db` script. Databases will be written to subdirectories within
/data/florinash/reference_data, to a date-stamped or versioned directory.

`build_ref_db` will download the database using the defiinition from the DATA
block at the bottom of the script, then will uncompress it, convert to fasta
format while generating the required taxonomy .dat file, then index for bwa
alignment. 

e.g. To download and format the 'parasites' database, run the command:

`build_ref_db --db parasites --download`

###Database Configurations

`build_ref_db` contains pre-build configurations in YAML format at the bottom
of the script. Three different was of defining databases are currently
available. 

All database configurations require a `format` parameter, which should be
either `embl` or `genbank`.

#### ftp

Allows downloading a database from a directory on an ftp server

Configuration parameters:

  * `ftphost`: URI or ftp server
  * `ftpdir`: directory on ftp server where database is found
  * `filepattern`: regular expression defining list of files to be downloaded

#### entrez

Allows a database to build using an entrez query

Configuration parameters:

  * `entrez_query`: query to carry out to select records to include in database

#### entrez_accession

Allows a database to be built from a list of accessions

Configuration parameters:

  * `accessions`: list of accessions to include in database


## Overall Workflow for carrying out metagenomic analysis

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

2. Read Filtering
-----------------

Quality-trimmed sequence reads next require filtering to remove reads relating
to the host organism, and other categories of organism which may be of interest
i.e. parasites. This can be accomplished with the `align_reads` script, which
is designed to be run on the cluster.

`align_reads` can function in two modes.

  a) Filtering: to filter out reads which align to the specified reference
database, returning a pair of gzipped fastq files containing only read pairs
which do not align to the reference, along with a report
(`sample.db.filtered_reads.txt`) listing the IDs of the reads which aligned
with the reference and were removed. A second report
(`sample.db.mapping_report.txt`) provides a taxonomic breakdown of the excluded
reads. This reports the NCBI taxonomy id, species name and number of read-pairs
aligned to the organism.

  b) Alignment: Returns the resulting aligned reads in bam format.

`align_reads` needs to be provided with the paths to the input and output
directories to use, which should be accessible on the cluster as well as on
codon. Each input directory should contain a series of per-sample directories,
each containing the sequence reads in a pair of gzip compressed fastq files,
named `*_1.fq.gz` and `*_2.fq.gz`. The output directory is populated in a
similar fashion when the job runs, consequently the outputs from one filtering
stage can be used directly as inputs to the next stage, allowing multiple
filtering steps to be easily carried out. 

A list of available databases can be found in the %dbs hash defined on line 89
of `align_reads`. At the time of writing, these were 'human', 'viral',
'plants', and 'parasites'.

`align_reads` is intended to be submitted directly to SGE to carry out the
filtering operations. The number of tasks making up the job is defined by the
'-t' directive at the head of the script. 

A typical series of filtering jobs will look like:

* `qsub bin/align_reads --in_dir /data/florinash/trimmed/ --out_dir /data/florinash/filtered/human --db human --filter`
* `qsub bin/align_reads --in_dir /data/florinash/filtered/human --out_dir /data/florinash/filtered/parasites --db parasites --filter`
* `qsub bin/align_reads --in_dir /data/florinash/parasites/ --out_dir /data/florinash/filtered/viral --db viral --filter`
* `qsub bin/align_reads --in_dir /data/florinash/plants/ --out_dir /data/florinash/filtered/plants--db plants --filter`


3. Taxonomic Classification
---------------------------

TODO: Write this!

4. Assembly
-----------

Assembly needs to be carried out on ax3 due to the high memory demands of
IDBA-UD. This is currently setup in ~jamesa, so will need separate installation
for other users, and the paths in the script modifying according (specifically,
the idba_dir definition on line 50). The number of jobs to run will need
modifying in the `#PBS -J` directive on line 15.

### 1st-Pass Assembly

Fastq files should be copied to `${SCRATCH}/fastq`, and an output directory
created in `${SCRATCH}/idba`. The jobs can then be submitted with `qsub
assemble.ax3`. This should produce one output directory per sample, including a
'contig.fa' file in an 'idba_out' subdirectory.	

### Compress Assemblies

Much of the data produced by the assembly process should be retainined but is
not likely to be required, so can be safely compressed. The script
`compress_assemblies.ax3` will created a bzipped tarball for each assembly
containing the files not likely to be required in future, while leaving the
contigs.fa and scaffolds.fa uncompressed. This should be run as:

`qsub compress_assemblies.ax3`

The resulting contents of the `${SCRATCH}/idba` directory should then be copied
back to codon.

### Contig Renaming

The contigs will have non-unique names between assemblies, so prior to doing
anything which may confuse the origins of these, they should be renamed to
include the sample identifier. This can be simply achieved using the
`tag_contigs` script on codon, which can be run directly rather than needing to
be queued. This requires a single parameter, 'in_dir', which is the path to the
directory of assemblies. Running this script will create a second copy of the
contigs named 'sample_id.contig.fa', with each contig id also prefixed with the
sample id.

`bin/tag_contigs --in_dir /data/florinash/assembly/IDBA_UD/`

## Second pass assembly

1) Extract unassembled reads
2) Reassemble - runtime = 36 hours, vmax=196Gb (real figures...)
3) COntig renaming
