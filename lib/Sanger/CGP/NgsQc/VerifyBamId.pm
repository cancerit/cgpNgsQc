package Sanger::CGP::NgsQc::VerifyBamId;

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
use IO::Uncompress::Gunzip qw($GunzipError);

const my $MIN_MAP_Q => 10;

const my $VERIFY => qq{%s --noPhoneHome --precise --maxDepth 200 --minMapQ $MIN_MAP_Q --minQ 13 --maxQ 40 --grid 0.05 --ignoreOverlapPair --self --vcf %s --bam %s --out %s};

## flags set to discard:
# 4 - read unmapped (0x4)
# 256 - not primary alignment (0x100)
# 512 - read fails platform/vendor quality checks (0x200)
# 1024 - read is PCR or optical duplicate (0x400)
# 2048 - supplementary alignment (0x800)
const my $CRAMTOBAM => qq{%s view -F 3844 -L %s -q $MIN_MAP_Q -b %s | tee %s | samtools index - %s.bai};

sub new {
  my ($class, $bam) = @_;
  my $self = {'total_snps' => 0,
              'excluded_snps' => 0,
              'retained_snps' => 0,
              'downsamp' => 1,
              'workspace' => undef,
              };
  bless $self, $class;
  $self->set_bam($bam) if(defined $bam);
  return $self;
}

sub set_bam {
  my ($self, $bam) = @_;
  $self->_init($bam);
  1;
}

sub set_workspace {
  my ($self, $workspace) = @_;
  $self->{'workspace'} = $workspace;
  1;
}

sub set_downsample {
  my ($self, $downsamp) = @_;
  $self->{'downsamp'} = $downsamp;
  1;
}

sub set_ascat {
  my ($self, $ascat) = @_;
  $self->{'ascat'} = $ascat;
  1;
}

sub set_snps {
  my ($self, $snp_vcf) = @_;
  $self->{'snp_vcf'} = $snp_vcf;
  1;
}

sub get_snps {
  my $self = shift;
  unless(exists $self->{'snp_vcf'}) {
    $self->{'snp_vcf'} = &default_snp_vcf;
  }
  $self->{'snp_vcf'};
}

sub _init {
  my ($self, $bam) = @_;
  $self->{'sample'} = Sanger::CGP::NgsQc::bam_sample_name($bam);
  $self->{'bam'} = $bam;
  $self->{'total_snps'} = 0;
  $self->{'excluded_snps'} = 0;
  $self->{'retained_snps'} = 0;
  delete $self->{'filtered_snps'};
  1;
}

sub result_to_json {
  my ($self, $top_key) = @_;
  my $final;
  if(defined $top_key) {
    $final = {$top_key => $self->{'result'}};
  }
  else {
    $final = $self->{'result'};
  }
  return encode_json($final);
}

sub run_verifyBam {
  my $self = shift;
  my $snps = $self->{'filtered_snps'};
  my $bam = $self->{'bam'};
  my ($out_stub) = $snps =~ m/(.*)_snps[.]vcf$/;

  $bam = $self->subsamp_cram_to_bam() if($bam =~ m/\.cram$/);

  my $verifyBamID = "$Bin/verifyBamId";
  $verifyBamID = which('verifyBamId') unless(-e $verifyBamID);
  die "Unable to find verifyBamID executable\n" unless(-e $verifyBamID);
  my $command = sprintf $VERIFY, $verifyBamID, $snps, $bam, $out_stub;
  warn "Running: $command\n";
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "An error occurred while executing:\n\t$command\nERROR: $stderr\n" if($exit);

  if($self->{'bam'} =~ m/\.cram$/) {
    unlink $bam;
    unlink $bam.'.bai';
  }

  my $result = "$out_stub.selfSM";
  open my $RES, '<', $result;
  my $head = <$RES>;
  my ($this_sample, $snp_count, $reads, $avg_depth, $contam) = (split /\t/, <$RES>)[0,3,4,5,6];
  close $RES;

  my %rg_data;

  $result = "$out_stub.selfRG";
  open $RES, '<', $result;
  while(my $line = <$RES>) {
    next if($line =~ m/^#/);
    chomp $line;
    my ($rg, $rg_reads, $rg_avg_depth, $rg_contam) = (split /\t/, $line)[1,4,5,6];
    $rg_data{$rg} = { 'reads_used' => $rg_reads,
                      'avg_depth' => $rg_avg_depth,
                      'contamination' => $rg_contam,};
  }
  close $RES;

  my %contam_data = ($this_sample => {'snps_used' => $snp_count,
                                      'reads_used' => $reads,
                                      'avg_depth' => $avg_depth,
                                      'contamination' => $contam,
                                      'by_readgroup' => \%rg_data,
                                      }
                                    );

  $self->{'result'} = \%contam_data;
}

sub subsamp_cram_to_bam {
  my $self = shift;
  my $sample = $self->{'sample'};
  my $outdir = $self->{'workspace'};
  my $bed_locs = "$outdir/${sample}_snps.bed";
  my $sub_bam = "$outdir/${sample}_snps.bam";
  $self->vcf_to_bed($bed_locs);

  my $samtools = which('samtools');
  my $command = sprintf $CRAMTOBAM, $samtools, $bed_locs, $self->{'bam'}, $sub_bam, $sub_bam;
  warn "Running: $command\n";
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "An error occurred while executing:\n\t$command\nERROR: $stderr\n" if($exit);
  unlink $bed_locs;
  return $sub_bam;
}

sub vcf_to_bed {
  my ($self, $tmp_bed) = @_;

  open my $OUT, '>', $tmp_bed;
  open my $SNPS, '<', $self->{'filtered_snps'};
  while(<$SNPS>) {
    next if $_ =~ m/^#/;
    chomp $_;
    my ($chr, $pos) = (split /\t/, $_)[0,1];
    printf $OUT "%s\t%d\t%d\n", $chr, $pos-1, $pos;
  }
  close $SNPS;
  close $OUT;
  return 1;
}

sub get_loh_regions {
  my $self = shift;
  my %loh_by_chr;
  if(exists $self->{'ascat'}) {
    open my $SEG, '<', $self->{'ascat'};
    while(my $line = <$SEG>) {
      chomp $line;
      my (undef, $chr, $from, $to, $wt_total, $wt_minor, $mt_total, $mt_minor) = split q{,}, $line;
      next unless(defined $mt_minor);
      push @{$loh_by_chr{$chr}}, [$from, $to] if($mt_minor < 1); ## this value may not be correct
    }
    close $SEG;
  }
  return \%loh_by_chr;
}

sub filter_loci {
  my $self  = shift;
  my $loh_regions = $self->get_loh_regions;
  my $loci = $self->get_snps;
  warn "Original Loci: $loci\n";
  my $sample = $self->{'sample'};
  my $one_in_x = $self->{'downsamp'};
  my $outdir = $self->{'workspace'};
  make_path($outdir) unless(-d $outdir);
  my $filtered_vcf = "$outdir/${sample}_snps.vcf";

  open my $VCF, '>', $filtered_vcf;
  my $z = new IO::Uncompress::Gunzip $loci, {MultiStream => 1};
  my $header = <$z>;
  print $VCF $header or die $!;
  my $last_chr = q{.};
  my @ranges;
  LINE: while(my $line = <$z>) {
    $self->{'total_snps'}++;
    next if($. % $one_in_x != 0);
    next if($line =~ m/^#/); # just incase there are multiple comment lines
    my ($chr, $pos) = $line =~ m/^([^\t]+)\t([[:digit:]]+)/;
    if($chr ne $last_chr) {
      if(exists $loh_regions->{$chr}) {
        @ranges = @{$loh_regions->{$chr}};
      }
      else {
        @ranges = ();
      }
      $last_chr = $chr;
    }
    for my $range(@ranges) {
      if($pos >= $range->[0] && $pos <= $range->[1]) {
        $self->{'excluded_snps'}++;
        next LINE;
      }
    }
    $self->{'retained_snps'}++;
    print $VCF $line or die $!;
  }
  close $z;
  close $VCF;
  $self->{'filtered_snps'} = $filtered_vcf;
  1;
}

sub default_snp_vcf {
  my $data_path = shift; # only for use by test cases
  # default to location of running code, if not present
  # likely that module has been installed
  # so try installed area
  unless(defined $data_path && -e $data_path) {
    $data_path = "$Bin/../share";
    $data_path = dist_dir('Sanger-CGP-NgsQc') unless(-e "$data_path/SNP6");
  }
  $data_path .= '/SNP6/GRCh37.vcf.gz';
  return $data_path;
}

1;
