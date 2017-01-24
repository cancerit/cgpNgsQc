#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
#
# Author: Keiran Raine <cgpit@sanger.ac.uk>
#
# This file is part of cgpNgsQc.
#
# cgpNgsQc is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’."
########## LICENCE ##########


use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use Getopt::Long;
use File::Spec;
use File::Path qw(make_path);
use Pod::Usage qw(pod2usage);

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Sanger::CGP::NgsQc;
use Sanger::CGP::NgsQc::VerifyBamId;

my $options = &setup;

my $verify = Sanger::CGP::NgsQc::VerifyBamId->new($options->{'bam'});
$verify->set_ascat($options->{'ascat'}) if(exists $options->{'ascat'});
$verify->set_snps($options->{'snps'}) if(exists $options->{'snps'});
$verify->set_downsample($options->{'downsamp'});
$verify->filter_loci($options->{'out'});
$verify->run_verifyBam;

if(defined $options->{'json'}) {
  if($options->{'json'} eq '-') {
    print $verify->result_to_json,"\n";
  }
  else {
    open my $OUT, '>', $options->{'json'};
    print $OUT $verify->result_to_json,"\n" or die "Failed to write to $options->{json}: $!";
    close $OUT;
  }
}
exit 0;

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'b|bam=s' => \$opts{'bam'},
              'o|outdir=s' => \$opts{'out'},
              'a|ascat=s' => \$opts{'ascat'},
              's|snps=s' => \$opts{'snps'},
              'd|downsamp=i' => \$opts{'downsamp'},
              'j|json=s' => \$opts{'json'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if(defined $opts{'v'}) {
    print 'VERSION: '.Sanger::CGP::NgsQc->VERSION,"\n";
    exit 0;
  }

  # check the required params
  pod2usage(-message => qq{\nERROR: 'outdir' must be defined.\n}, -verbose => 1,  -output => \*STDERR) unless(defined $opts{'out'});
  pod2usage(-message => qq{\nERROR: 'bam' must be defined.\n}, -verbose => 1,  -output => \*STDERR) unless(defined $opts{'bam'});
  pod2usage(-message => qq{\nERROR: 'bam' must exist ($opts{bam}).\n}, -verbose => 1,  -output => \*STDERR) unless(-e $opts{'bam'});

  delete $opts{'ascat'} unless(defined $opts{'ascat'});
  pod2usage(-message => qq{\nERROR: 'ascat' must exist when specified ($opts{ascat}).\n}, -verbose => 1,  -output => \*STDERR) if(exists $opts{'ascat'} && !-e $opts{'ascat'});

  if(defined $opts{'snps'}) {
    pod2usage(-message => qq{\nERROR: 'snps' must exist when specified ($opts{snps}).\n}, -verbose => 1,  -output => \*STDERR) if(exists $opts{'snps'} && !-e $opts{'snps'});
  }
  else {
    delete $opts{'snps'};
  }

  $opts{'downsamp'} = 1 unless(defined $opts{'downsamp'});

  ## setup out location
  make_path($opts{'out'}) or die $! unless(-e $opts{'out'});

  return \%opts;
}

__END__

=head1 NAME

verifyBamHomChk.pl - Runs verify BAM

=head1 SYNOPSIS

verifyBamHomChk.pl [options]

  Required parameters:

    -outdir   -o  Directory for output

    -bam      -b  BAM file to process

  Optional parameters:

    -snps     -s  VCF file of SNP locations as described here:
                      http://genome.sph.umich.edu/wiki/VerifyBamID
                        ([b]gzip compressed works too)
                      defaults to included GRCh37 SNP6 loci.

    -downsamp -d  Subsample the SNPs at rate of 1 in X [1]

    -ascat    -a  Exclude LOH regions based on ASCAT segments file []

    -json     -j  Output summary as JSON string '-' for STDOUT.

  Other:
    -help     -h  Brief help message.
    -man      -m  Full documentation.
    -version  -v  Shows version.
