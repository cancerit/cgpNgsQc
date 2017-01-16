#! /usr/bin/perl

# perl-5.16.3 -I /software/CGP/canpipe/live/lib/perl5 validate_smaple_meta.pl
use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use Carp;

use Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader;

my $opts = setup_options();
my $obj = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader->new($opts);

if($obj->validate_status) {
  warn "\n--------\nSamples validated successfully! Please supply the '.tsv' output.\n--------\n"; # hide the tsv/xls/xlsx here
} else {
  warn "\n--------\nValidation failed, please check output:$opts->{'o'} for issues.\n--------\n";
}

exit;

sub setup_options {
  my %opts = ('a' => 0,
            'f' => '',
            );
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'i|in=s' => \$opts{'i'},
          'o|out=s' => \$opts{'o'},
					'f|format:s' => \$opts{'f'},
          'a|check-all' => \$opts{'a'},
					) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  pod2usage(-message => "\nERROR: Option 'i|in' required and the file must exist .", -verbose => 0) if(!defined $opts{'i'} || ! -e $opts{'i'});
  pod2usage(-message => "\nERROR: Option 'o|out' required.", -verbose => 0) unless(defined $opts{'o'});
  pod2usage(-message => "\nERROR: Option 'f|format' must be be one of the 3 formats tsv|xls|xlsx, default to tsv.", -verbose => 0) if($opts->{'f'} !~ m/^xls$|^xlsx$|^tsv$/i);
  pod2usage(-message => "\nERROR: Input does not have .xls extension.", -verbose => 0) if($opts->{'f'} =~ m/^xls$/i && $opts->{'i'} !~ m/\.xls$/i);
  pod2usage(-message => "\nERROR: Input does not have .xlsx extension.", -verbose => 0) if($opts->{'f'} =~ m/^xlsx$/i && $opts->{'i'} !~ m/\.xlsx$/i);
  pod2usage(-message => "\nERROR: Output does not have .xls extension.", -verbose => 0) if($opts->{'f'} =~ m/^xls$/i && $opts->{'o'} !~ m/\.xls$/i);
  pod2usage(-message => "\nERROR: Output does not have .xlsx extension.", -verbose => 0) if($opts->{'f'} =~ m/^xlsx$/i && $opts->{'o'} !~ m/\.xlsx$/i);
  pod2usage(-message => "\nERROR: Input file is the same as output file.", -verbose => 0) if($opts->{'o'} eq $opts->{'o'});
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
    -format      -f    Input file format (xlsx, xls or tsv. Optional, default to tsv).
    -check-all   -a    Optional, check all reads in a bam instead of the first million reads. Use when some of RG IDs were not found in the first million reads.

  Other:
    -help        -h   Brief help message.
    -man         -m   Full documentation.

=head1 OPTIONS

=over 8

=item B<-input>

Input file name. The input can be in tsv, xlsx or xls format. First line will be taken as header (column names) if it starts with '#' or 'Donor_ID', otherwise it'll be considered as no header.

=item B<-out>

Output file name. It'll be in the format as specified with '-format'. Appropriate extension name will be added if it's missing. On a successful validation, if input format is 'xls' or 'xlsx', an tsv file will be generated as well.

=item B<-format>

Input file format [tsv|xlsx|xls]. Default to 'tsv'.

=item B<-check-all>

Check all reads in a bam to validate read groups. Optional, use when reads of a read group are NOT in the first million reads.

=item B<-help>

Prints the help for this script

=item B<-man>

Prints the man page for this script

=back

=head1 DESCRIPTION

B<validate_smaple_meta.pl>

=cut
