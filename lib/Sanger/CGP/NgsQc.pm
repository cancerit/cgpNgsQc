package Sanger::CGP::NgsQc;

use strict;
use base 'Exporter';
use Bio::DB::Sam;

our $VERSION = '0.0.1';
our @EXPORT = qw($VERSION);

sub bam_sample_name {
  my $bam = shift;
  my @lines = split /\n/, Bio::DB::Sam->new(-bam => $bam)->header->text;
  my $sample;
  for(@lines) {
    if($_ =~ m/^\@RG.*\tSM:([^\t]+)/) {
      $sample = $1;
      last;
    }
  }
  die "Failed to determine sample from BAM header\n" unless(defined $sample);
  return $sample;
}

1;
