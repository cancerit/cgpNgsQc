#! /usr/bin/perl

# perl-5.16.3 -I /software/CGP/canpipe/live/lib/perl5 validate_smaple_meta.pl
use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use Carp;

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::WriteExcel; # to write xls output
use Excel::Writer::XLSX; # to write xlsx output

use Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader;

my $opts = setup_options();
my $input = $opts->{'i'};
my $output = $opts->{'o'};
my $check_all = 'n';
if (defined $opts->{'a'}) {
  $check_all = 'y';
}
#print "check_all setting: $check_all.\n";
my $format;

if (defined $opts->{'f'} && $opts->{'f'} =~ m/^xls$|^xlsx$|^tsv$/i) {
  $format = $opts->{'f'};
} elsif (! defined $opts->{'f'}) {
  $format = 'tsv';
} else {
  die("--------\nError: Does not recognise the specified file format, accept values are: 'tsv', 'xls' and 'xlsx'.\n--------\n");
}

my @in;
### read in input as array of lines (concatenated columns with tabs)
if ($format eq 'xlsx'|| $format eq 'xls') {
  my $book = ReadData ($input);
  my @rows = Spreadsheet::Read::rows($book->[1]);
  foreach my $row (@rows) {
    foreach my $cell (@$row) {
      $cell = '' if (! defined $cell);
      # print $cell."\n";
    }
    push @in, (join "\t", @$row);
  }
} else { # assume it's a tsv
  open (my $FH, '<', $input) or croak("--------\nError trying to open input file $input: $!\n--------\n");
  chomp(@in = <$FH>);
  close $FH;
}

my ($valid_results, $valid_status) = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples(\@in, $check_all);
# print "below is the out:\n";
# foreach my $line (@$valid_results) {print "$line\n";}

### output into the same format again on failed validation, an extra tsv if (input is excel) on a successful validation!
if ($format eq 'xlsx') {
  # rename output, if it does not ends with .xlsx
  if ($output !~ /\.xlsx$/) {
    $output .= '.xlsx';
  }
  my $workbook = Excel::Writer::XLSX->new($output);
  my $worksheet = $workbook->add_worksheet();
  my $row = 0;
  foreach my $line (@$valid_results) {
    my @array = split("\t", $line);
    $worksheet->write_row($row, 0, \@array);
    $row++;
  }
  if ($valid_status eq "green") {
    $output .= ".tsv";
    warn "\n--------\nSamples validated successfully! Please supply the tsv output file.\ntsv output: $output\n--------\n";
    open OUT, ">$output" or croak("--------\nCan not open the output: $output to write. $!\n--------\n");
    foreach my $line (@$valid_results) {print OUT "$line\n";}
    close OUT;
  }
} elsif ($format eq 'xls') {
  # rename output, if it does not ends with .xls
  if ($output !~ /\.xls$/) {
    $output .= '.xls';
  }
  my $workbook = Spreadsheet::WriteExcel->new($output);
  my $worksheet = $workbook->add_worksheet();
  my $row = 0;
  foreach my $line (@$valid_results) {
    my @array = split("\t", $line);
    $worksheet->write_row($row, 0, \@array);
    $row++;
  }
  if ($valid_status eq "green") {
    $output .= ".tsv";
    warn "\n--------\nSamples validated successfully! Please supply the tsv output file.\ntsv output: $output\n--------\n";
    open OUT, ">$output" or croak("--------\nCan not open the output: $output to write. $!\n--------\n");
    foreach my $line (@$valid_results) {print OUT "$line\n";}
    close OUT;
  }
} else { # output to tsv format
  if ($valid_status eq "green") {
    warn "\n--------\nSamples validated successfully! Please supply the tsv output file.\ntsv output: $output\n--------\n";
  }
  open OUT, ">$output" or croak("--------\nCan not open the output: $output to write. $!\n--------\n");
  foreach my $line (@$valid_results) {print OUT "$line\n";}
  close OUT;
}

exit;

sub setup_options {
  my %opts;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'i|in=s' => \$opts{'i'},
          'o|out=s' => \$opts{'o'},
					'f|format=s' => \$opts{'f'},
          'a|check-all' => \$opts{'a'},
					) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  pod2usage(-message => "Option 'i|in' required and the file must exist .", -verbose => 0) if(!defined $opts{'i'} || ! -e $opts{'i'});
  pod2usage(-message => "Option 'o|out' required.", -verbose => 0) unless(defined $opts{'o'});
  return \%opts;
}

__END__

=head1 NAME

validate_smaple_meta.pl - Validate sample meta data and coresponding bam files, upon successful validation, UUIDs will be assigned and md5sum of bam files will be produced.

=head1 SYNOPSIS

validate_smaple_meta.pl [options]

  Required parameters:
    -in          -i    Input file.
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
