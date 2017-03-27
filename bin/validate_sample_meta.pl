#! /usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
#
# Author: Yaobo Xu <cgpit@sanger.ac.uk>
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
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use Carp;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader;

my $opts = setup_options();
my $obj = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader->new($opts);

if($obj->{validate_status}) {
  warn "\n--------\nSamples validated successfully!\n--------\n";
} else {
  warn "\n--------\nValidation failed, please check output:$opts->{'o'} for issues.\n--------\n";
}

exit;

sub setup_options {
  my %opts = (
            'f' => '',
            );
  GetOptions(
          'h|help' => \$opts{'h'},
          'm|man' => \$opts{'m'},
          'i|in=s' => \$opts{'i'},
          'o|out=s' => \$opts{'o'},
          'f|format:s' => \$opts{'f'},
          'v|version' => \$opts{'v'},
          ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if(defined $opts{'v'}) {
    print 'VERSION: '.Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader->VERSION,"\n";
    exit 0;
  }

  pod2usage(-message => "\nERROR: Option 'i|in' required and the file must exist .", -verbose => 0) if(!defined $opts{'i'} || ! -e $opts{'i'});
  pod2usage(-message => "\nERROR: Option 'o|out' required.", -verbose => 0) unless(defined $opts{'o'});
  pod2usage(-message => "\nERROR: Option 'f|format' must be be one of the 3 formats [tsv|xls|xlsx].", -verbose => 0) if($opts{'f'} !~ m/^xls$|^xlsx$|^tsv$/i);
  pod2usage(-message => "\nERROR: Input does not have .xls extension.", -verbose => 0) if($opts{'f'} =~ m/^xls$/i && $opts{'i'} !~ m/\.xls$/i);
  pod2usage(-message => "\nERROR: Input does not have .xlsx extension.", -verbose => 0) if($opts{'f'} =~ m/^xlsx$/i && $opts{'i'} !~ m/\.xlsx$/i);
  pod2usage(-message => "\nERROR: Output does not have .xls extension.", -verbose => 0) if($opts{'f'} =~ m/^xls$/i && $opts{'o'} !~ m/\.xls$/i);
  pod2usage(-message => "\nERROR: Output does not have .xlsx extension.", -verbose => 0) if($opts{'f'} =~ m/^xlsx$/i && $opts{'o'} !~ m/\.xlsx$/i);
  pod2usage(-message => "\nERROR: Input file is the same as output file.", -verbose => 0) if($opts{'i'} eq $opts{'o'});
  return \%opts;
}

__END__

=head1 NAME

validate_smaple_meta.pl - Validate sample meta data and coresponding bam files, upon successful validation, UUIDs will be assigned and md5sum of bam files will be produced.

=head1 SYNOPSIS

validate_smaple_meta.pl [options]

  Required parameters:
    -in          -i    Input file. First column must be DONOR_ID.
    -out         -o    Output file (It'll be in the format as specified with '-format'. Extension name will be added when missing).
    -format      -f    Input file format (xlsx, xls or tsv).

  Other:
    -help        -h    Brief help message.
    -man         -m    Full documentation.
    -version     -v    Shows version.

=head1 OPTIONS

=over 8

=item B<-input>

Input file name. The input can be in tsv, xlsx or xls format. First line will be taken as header (column names) if it starts with '#' or 'Donor_ID', otherwise it'll be considered as no header.

=item B<-out>

Output file name. It'll be in the format as specified with '-format'. Appropriate extension name will be added if it's missing. On a successful validation, if input format is 'xls' or 'xlsx', an tsv file will be generated as well.

=item B<-format>

Input file format [tsv|xlsx|xls].

=item B<-help>

Prints the help for this script

=item B<-man>

Prints the man page for this script

=back

=head1 DESCRIPTION

B<validate_smaple_meta.pl>

=cut
