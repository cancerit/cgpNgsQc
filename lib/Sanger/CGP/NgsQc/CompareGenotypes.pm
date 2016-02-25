package Sanger::CGP::NgsQc::CompareGenotypes;

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


use Sanger::CGP::NgsQc;

use strict;
use warnings FATAL=>'all';
use Const::Fast qw(const);
use FindBin qw($Bin);
use File::ShareDir qw(dist_dir);
use File::Which qw(which);
use Capture::Tiny qw(capture);
use autodie qw(:all);
use JSON;
use Carp qw(croak);

const my $ALLELE_COUNT => q{%s --min-base-qual 20 --min-map-qual 35 --output-file /dev/stdout --loci-file %s -b %s};
const my $MIN_ALLELE_PROP => 0.21;
const my $MIN_READS_ALLELE_COUNT => 5;
const my $MIN_READS_XY => 3;

sub new {
  my $class = shift;
  my $self = {
              };
  bless $self, $class;
  return $self;
}

sub add_tumour_bams {
  my ($self, @bams) = @_;
  $self->{'tumour_bams'} = \@bams;
  1;
}

sub add_normal_bam {
  my ($self, $bam) = @_;
  $self->{'normal_bam'} = $bam;
  1;
}

sub result_to_json {
  my ($self, $top_key) = @_;
  my $final;
  if(defined $top_key) {
    $final = {$top_key => $self->{'summary'}};
  }
  else {
    $final = $self->{'summary'};
  }
  return encode_json($final);
}

sub compare {
  my ($self, $outdir) = @_;
  my $sm_normal = $self->genotype($outdir, $self->{'normal_bam'});
  $self->gender($outdir, $self->{'normal_bam'});
  my @sm_tumours;
  for my $tum_bam(@{$self->{'tumour_bams'}}) {
    my $sm_tumour = $self->genotype($outdir, $tum_bam);
    $self->gender($outdir, $tum_bam);
    $self->compare_samples($outdir, $sm_normal, $sm_tumour);
  }
  1;
}

sub compare_samples {
  my ($self, $outdir, $normal, $tumour) = @_;
  my @normal_calls_geno = @{$self->{'genotypes'}->{$normal}};
  my @tumour_calls_geno = @{$self->{'genotypes'}->{$tumour}};
  my $loci_count = scalar @normal_calls_geno;
  my ($informative, $matched) = (0,0);
  for(my $i=0; $i<$loci_count; $i++) {
    next if($normal_calls_geno[$i] eq q{.} || $tumour_calls_geno[$i] eq q{.});
    $informative++;
    $matched++ if($normal_calls_geno[$i] eq $tumour_calls_geno[$i]);
  }

  my @normal_calls_gender = @{$self->{'genders'}->{$normal}};
  my @tumour_calls_gender = @{$self->{'genders'}->{$tumour}};
  my $gender_loci_count = scalar @normal_calls_gender;
  my ($xy, $gender_match) = (0,0);
  for(my $i=0; $i<$gender_loci_count; $i++) {
    $gender_match++ if($normal_calls_gender[$i] eq $tumour_calls_gender[$i]);
    $xy++ if($normal_calls_gender[$i] eq 'XY');
  }
  my $gender = 'XX';
  $gender = 'XY' if($xy > 0); # any male specific loci present in normal results in XY call
  my %result_struct = ('sample' => $tumour,
                       'genotype' => {'frac_informative_genotype' =>  $informative / $loci_count,
                                      'frac_matched_genotype' => $matched / $informative},
                       'gender' => {'frac_match_gender' =>  $gender_match / $gender_loci_count,
                                    'gender' => $gender});
  push @{$self->{'summary'}{'tumours'}}, \%result_struct;

  # general
  $self->{'summary'}{'compared_against'} = $normal;
  $self->{'summary'}{'total_loci_genotype'} = $loci_count;
  $self->{'summary'}{'total_loci_gender'} = $gender_loci_count;

  my $result = "$outdir/${tumour}_vs_$normal.genotype.txt";
  open my $RES,'>',$result;
  print $RES "#GENERAL LOCI\n" or die "Failed to write line to $result: $!";
  for my $key(sort keys %{$result_struct{'genotype'}}) {
    print $RES "$key\t$result_struct{genotype}{$key}\n" or die "Failed to write line to $result: $!";
  }
  print $RES "genotype\t$tumour\t".join(q{,},@tumour_calls_geno)."\n" or die "Failed to write line to $result: $!";
  print $RES "genotype\t$normal\t".join(q{,},@normal_calls_geno)."\n" or die "Failed to write line to $result: $!";
  print $RES "#GENDER LOCI\n" or die "Failed to write line to $result: $!";
  for my $key(sort keys %{$result_struct{'gender'}}) {
    print $RES "$key\t$result_struct{gender}{$key}\n" or die "Failed to write line to $result: $!";
  }
  print $RES "gender\t$tumour\t".join(q{,},@tumour_calls_gender)."\n" or die "Failed to write line to $result: $!";
  print $RES "gender\t$normal\t".join(q{,},@normal_calls_gender)."\n" or die "Failed to write line to $result: $!";
  close $RES;
  1;
}

sub genotype {
  my ($self, $outdir, $bam) = @_;
  my $sample = Sanger::CGP::NgsQc::bam_sample_name($bam);
  my $command = sprintf $ALLELE_COUNT, alleleCounter(), default_genotype_loci(), $bam;
  warn "Running: $command\n";
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "An error occurred while executing:\n\t$command\nERROR: $stderr\n" if($exit);

  my @genotype;
  my @results = split /\n/, $stdout;

  my $augmented = "$outdir/$sample.full_genotype.tsv";
  open my $AUG, '>', $augmented;
  for my $line(@results) {
    if($line =~ m/^#/) {
      print $AUG "$line\tGENOTYPE\n" or die "Failed to write line to $augmented: $!";
      next;
    }
    my ($chr, $pos, $a, $c, $g, $t, $good) = split /\t/, $line;
    my $geno = _calculate_genotype_from_allele_count($a, $c, $g, $t, $good);
    print $AUG "$line\t$geno\n" or die "Failed to write line to $augmented: $!";
    push @genotype, $geno;
  }
  close $AUG;
  $self->{'genotypes'}->{$sample} = \@genotype;
  return $sample;
}

sub gender {
  my ($self, $outdir, $bam) = @_;
  my $sample = Sanger::CGP::NgsQc::bam_sample_name($bam);
  my $command = sprintf $ALLELE_COUNT, alleleCounter(), default_gender_loci(), $bam;
  warn "Running: $command\n";
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "An error occurred while executing:\n\t$command\nERROR: $stderr\n" if($exit);

  my @gender_chrs;
  my @results = split /\n/, $stdout;

  my $augmented = "$outdir/$sample.full_gender.tsv";
  open my $AUG, '>', $augmented;
  for my $line(@results) {
    if($line =~ m/^#/) {
      print $AUG "$line\tGENDER\n" or die "Failed to write line to $augmented: $!";
      next;
    }
    my $good = (split /\t/, $line)[-1];
    my $gender = _calculate_gender_from_allele_count($good);
    print $AUG "$line\t$gender\n" or die "Failed to write line to $augmented: $!";
    push @gender_chrs, $gender;
  }
  close $AUG;
  $self->{'genders'}->{$sample} = \@gender_chrs;
  return $sample;
}

sub _calculate_genotype_from_allele_count{
  my ($a_a,$a_c,$a_g,$a_t,$good) = @_;
  my $geno;
  return q{.} if($good < $MIN_READS_ALLELE_COUNT);

  my @counts;
  push @counts, ['A', $a_a] if($a_a/$good >= $MIN_ALLELE_PROP);
  push @counts, ['C', $a_c] if($a_c/$good >= $MIN_ALLELE_PROP);
  push @counts, ['G', $a_g] if($a_g/$good >= $MIN_ALLELE_PROP);
  push @counts, ['T', $a_t] if($a_t/$good >= $MIN_ALLELE_PROP);

  my $entries = scalar @counts;
  if($entries == 0) {
    $geno = q{.};
  }
  elsif($entries == 1) {
    $geno = $counts[0][0].$counts[0][0];
  }
  else {
    @counts = sort {$b->[1]<=>$a->[1]}  @counts; # reverse sorts by the counts
    $geno = join(q{}, sort {$a cmp $b} $counts[0][0], $counts[1][0]); # then sort the alleles into the string
  }
  die "Error calculating genotype from allele counts $a_a,$a_c,$a_g,$a_t,$good.\n" if((length $geno)>2 || (length $geno) == 0);
  return $geno;
}

sub _calculate_gender_from_allele_count{
  my $good = shift;
  my $gender = 'XX';
  $gender = 'XY' if($good > $MIN_READS_XY);
  return $gender;
}

sub alleleCounter {
  my $alleleCounter = "$Bin/alleleCounter";
  $alleleCounter = which('alleleCounter') unless(-e $alleleCounter);
  die "Unable to find alleleCounter executable\n" unless(-e $alleleCounter);
  return $alleleCounter;
}

sub default_genotype_loci {
  my $data_path .= share_dir().'/genotype/general.tsv';
  return $data_path;
}

sub default_gender_loci {
  my $data_path .= share_dir().'/genotype/gender.tsv';
  return $data_path;
}

sub share_dir {
  my $data_path = shift;
  unless(defined $data_path && -e $data_path) {
    $data_path = "$Bin/../share";
    $data_path = dist_dir('Sanger-CGP-NgsQc') unless(-e "$data_path/genotype");
  }
  return $data_path;
}

1;

