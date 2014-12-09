package Sanger::CGP::NgsQc::VerifyBamId;

use Sanger::CGP::NgsQc;

use strict;
use warnings FATAL=>'all';
use Const::Fast qw(const);
use FindBin qw($Bin);
use File::ShareDir qw(module_dir);
use File::Which qw(which);
use Capture::Tiny qw(capture);
use autodie qw(:all);
use JSON;
use IO::Uncompress::Gunzip qw($GunzipError);

const my $VERIFY => q{%s --precise --maxDepth 200 --minMapQ 10 --minQ 13 --maxQ 40 --grid 0.05 --ignoreOverlapPair --self --vcf %s --bam %s --out %s};

sub new {
  my ($class, $bam) = @_;
  my $self = {'total_snps' => 0,
              'excluded_snps' => 0,
              'retained_snps' => 0,
              'downsamp' => 1,
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

sub result_to_fh {
  my ($self, $fh) = @_;

  ...;

  1;
}

sub run_verifyBam {
  my $self = shift;
  my $snps = $self->{'filtered_snps'};
  my $bam = $self->{'bam'};
  my ($out_stub) = $snps =~ m/(.*)_snps[.]vcf$/;

  my $verifyBamID = "$Bin/verifyBamID";
  $verifyBamID = which('verifyBamID') unless(-e $verifyBamID);
  die "Unable to find verifyBamID executable\n" unless(-e $verifyBamID);
  my $command = sprintf $VERIFY, $verifyBamID, $snps, $bam, $out_stub;
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "An error occurred while executing:\n\t$command\nERROR: $stderr\n" if($exit);

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
  my ($self, $outdir) = @_;
  my $loh_regions = $self->get_loh_regions;
  my $loci = $self->get_snps;
  my $sample = $self->{'sample'};
  my $one_in_x = $self->{'downsamp'};
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
    $data_path = module_dir('Sanger::CGP::NgsQc::VerifyBamId') unless(-e "$data_path/SNP6");
  }
  $data_path .= '/SNP6/GRCh37.vcf.gz';
  return $data_path;
}

1;
