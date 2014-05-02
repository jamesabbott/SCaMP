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

use File::Find::Rule;
use Bio::SeqIO;

{
    my @files = File::Find::Rule->file()->name("*R1_val_1.fq.gz")->in("/data/florinash/trimmed");
    @files = map { $_->[0] }
      sort { $a->[1] cmp $b->[1] }
      map { [ $_, /trimmed\/([A-Z0-9_\-]+)/ ] } @files;

    foreach my $file (@files) {
	my $sample = $1 if($file=~/trimmed\/([A-Z0-9_\-]+)/);   
	my $count = `zcat $file \| echo \$((\`wc -l\`/4 ))`;
	print "$sample\t$count";
    }

}
