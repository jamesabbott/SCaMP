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

use Net::FTP;

{
    my $host     = 'ftp.ensemblgenomes.org';
    my $ref_data = '/data/florinash/reference_data/db';

    my $ftp = Net::FTP->new( $host, Passive => 1, Debug => 0 ) or die "Cannot connect to $host:$!";
    $ftp->login( "anonymous", 'bsshelp@imperial.ac.uk' ) or die "Cannot login: " . $ftp->message;
    $ftp->cwd("/pub/bacteria/current/fasta") or die "Cannot chdir to /pub/bacteria: " . $ftp->message;

    my $dir = $ftp->pwd();

    my $release;
    if ( $dir =~ /release-([0-9]+)/ ) {
        $release = $1;
    }

    my $ensembl_db = $ref_data . "/ensembl_$release";
    if ( !-d $ensembl_db ) { mkdir $ensembl_db or die "Error creating $ensembl_db: $!" }
    chdir $ensembl_db or die "Error chdiring to $ensembl_db: $!";

    my @collections = $ftp->ls();
    foreach my $collection (@collections) {
        if ( $collection =~ /^bacteria_[0-9]+_collection/ ) {
            $ftp->cwd($collection) or die "Error chdiring to $collection: " . $ftp->message;
            my @genomes = $ftp->ls();
            foreach my $genome (@genomes) {
                $ftp->cwd("$genome/dna") or die "Error chdiring to $genome/dna:" . $ftp->message;
                my @files = $ftp->ls();
                foreach my $file (@files) {
                    if ( $file =~ /\.dna\.genome\.fa\.gz/ ) {
                        print "Downloading $file...\n";
                        $ftp->get($file);
                    }
                }
                $ftp->cwd("../../") or die "Error chdiring to ../..: " . $ftp->message;
            }
            $ftp->cwd("../") or die "Error chdiring to ..: " . $ftp->message;
        }
    }

}
