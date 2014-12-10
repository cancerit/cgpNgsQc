#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd.
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

use Getopt::Long;
use File::Spec;
use File::Path qw(make_path);
use Pod::Usage qw(pod2usage);

use Sanger::CGP::NgsQc;
use Sanger::CGP::NgsQc::CompareGenotypes;

my $options = &setup;

my $compare = Sanger::CGP::NgsQc::CompareGenotypes->new();
$compare->add_tumour_bams(@{$options->{'t_bams'}});
$compare->add_normal_bam($options->{'n_bam'});
$compare->compare($options->{'out'});

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
              'tb|tumour_bams=s@' => \$opts{'t_bams'},
              'nb|normal_bam=s' => \$opts{'n_bam'},
              'o|outdir=s' => \$opts{'out'},
              'j|json' => \$opts{'json'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if(defined $opts{'v'}) {
    print 'VERSION: '.Sanger::CGP::NgsQc->VERSION,"\n";
    exit 0;
  }

  # check the required params
  pod2usage(-message => qq{\nERROR: 'outdir' must be defined.\n}, -verbose => 1,  -output => \*STDERR) unless(defined $opts{'out'});
  pod2usage(-message => qq{\nERROR: 'tumour_bams' must be defined.\n}, -verbose => 1,  -output => \*STDERR) unless(defined $opts{'t_bams'});
  pod2usage(-message => qq{\nERROR: 'normal_bam' must be defined.\n}, -verbose => 1,  -output => \*STDERR) unless(defined $opts{'n_bam'});
  pod2usage(-message => qq{\nERROR: 'normal_bam' must exist ($opts{n_bam}).\n}, -verbose => 1,  -output => \*STDERR) unless(-e $opts{'n_bam'});

  ## setup out location
  make_path($opts{'out'}) or die $! unless(-e $opts{'out'});

  return \%opts;
}

__END__

=head1 NAME

compareBamGenotypes.pl - Compare a set of BAM files from the same donor.

=head1 SYNOPSIS

compareBamGenotypes.pl [options]

  Required parameters:

    -outdir       -o    Directory for output.
    -tumour_bams  -tb   Tumour BAM file(s) to process.
    -normal_bam   -nb   Normal BAM file to process.

  Optional parameters:
    -json         -j  Output summary as JSON string '-' for STDOUT.

  Other:
    -help         -h  Brief help message.
    -man          -m  Full documentation.
    -version      -v  Shows version.
