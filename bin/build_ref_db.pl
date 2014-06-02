#!/usr/bin/perl

######################################################################
#
# $HeadURL$
# $Author$
# $Revision$
# $Date$
#
######################################################################

use warnings;
use strict;

use Getopt::Long;
use Net::FTP;
use File::Path qw(remove_tree);
use Archive::Tar;
use Bio::DB::Taxonomy;

{

    my $clean;
    GetOptions( 'clean' => \$clean );    #clean old taxdata and refresh...

    my $tax_host = 'ftp.ncbi.nlm.nih.gov';
    my $host     = 'ftp.ensemblgenomes.org';
    my $ref_data = '/data/florinash/reference_data/';
    my $tax_data = $ref_data . '/taxonomy';

    my $ftp = Net::FTP->new( $host, Passive => 1, Debug => 0 ) or die "Cannot connect to $host:$!";
    $ftp->login( "anonymous", 'bsshelp@imperial.ac.uk' ) or die "Cannot login: " . $ftp->message;
    $ftp->binary();
    $ftp->cwd("/pub/bacteria/current/fasta") or die "Cannot chdir to /pub/bacteria: " . $ftp->message;
    my $dir = $ftp->pwd();

    my $release = $1 if ( $dir =~ /release-([0-9]+)/ );
    die "Could not parse release version from $dir" unless ($release);

    my $tax_ftp = Net::FTP->new( $tax_host, Passive => 1, Debug => 0 ) or die "Cannot connect to $host:$!";
    $tax_ftp->login( "anonymous", 'bsshelp@imperial.ac.uk' ) or die "Cannot login: " . $ftp->message;
    $tax_ftp->binary();

    if ($clean) {
        remove_tree($tax_data) if ( -d $tax_data );
        mkdir "$tax_data" or die "Error creating $tax_data: $!";
        chdir "$tax_data" or die "Error chdiring to $tax_data: $!";
        $tax_ftp->get('/pub/taxonomy/taxdump.tar.gz');

        my $tar = Archive::Tar->new;
        $tar->read('taxdump.tar.gz');
        $tar->extract();
        unlink('taxdump.tar.gz') or die "Error removing taxdump.tar.gz: $!";
    }
    print "Creating taxonomy object...\n";
    my $taxonomy = Bio::DB::Taxonomy->new(
                                           -source    => 'flatfile',
                                           -directory => $tax_data,
                                           -nodesfile => "$tax_data/nodes.dmp",
                                           -namesfile => "$tax_data/names.dmp"
                                         );
    print "Done...\n";

    my $ensembl_db = $ref_data . "/ensembl_$release";
    #if ($clean) {
    #remove_tree($ensembl_db) if ( -d $ensembl_db );
    #mkdir $ensembl_db or die "Error creating $ensembl_db: $!";
    #chdir $ensembl_db or die "Error chdiring to $ensembl_db: $!";

    #my @collections = $ftp->ls();
    #foreach my $collection (@collections) {
    #    if ( $collection =~ /^bacteria_[0-9]+_collection/ ) {
    #        $ftp->cwd($collection) or die "Error chdiring to $collection: " . $ftp->message;
    #        my @genomes = $ftp->ls();
    #        foreach my $genome (@genomes) {
    #            $ftp->cwd("$genome/dna") or die "Error chdiring to $genome/dna:" . $ftp->message;
    #            my @files = $ftp->ls();
    #            foreach my $file (@files) {
    #                if ( $file =~ /\.dna\.genome\.fa\.gz/ ) {
    #                    print "Downloading $file...\n";
    #                    $ftp->get($file);
    #                }
    #            }
    #            $ftp->cwd("../../") or die "Error chdiring to ../..: " . $ftp->message;
    #        }
    #        $ftp->cwd("../") or die "Error chdiring to ..: " . $ftp->message;
    #    }
    #}
    #} else {
	chdir $ensembl_db or die "Error chdiring to $ensembl_db: $!";
    #}

    my $species = $ftp->get('/pub/bacteria/current/species_EnsemblBacteria.txt');
    open SPECIES, 'species_EnsemblBacteria.txt' or die "Error opening species_EnsemblBacteria.txt: $!";
    while ( my @fields = split( /\t/, <SPECIES> ) ) {
        next if ( $fields[0] =~ /^#/ );
        print "$fields[0]\t$fields[3]\n";
	my $taxon = $taxonomy->get_taxon(-taxonid=>$fields[3]);
	if (defined($taxon)) {
	    my $tree = Bio::Tree::Tree->new(-node => $taxon);
	    my @taxa = $tree->get_nodes;
	    my $species =  $tree->find_node(-rank => 'species');
	    print "Species = " . $species->scientific_name() . "\n";
	};
    }
    close SPECIES;


}
