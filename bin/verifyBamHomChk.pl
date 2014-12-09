#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use FindBin qw($Bin);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use File::Which qw(which);
use File::Path qw(make_path);
use Bio::DB::Sam;
use Capture::Tiny qw(capture);

const my $VERIFY => q{%s --precise --maxDepth 200 --minMapQ 10 --minQ 13 --maxQ 40 --grid 0.05 --ignoreOverlapPair --self --vcf %s --bam %s --out %s};

my $options = &setup;

my $ascat_loh = get_loh_regions($options->{'loh_bed'});
my $sample = sample_name($options->{'bam'});
my $filtered_vcf = filter_loci($options->{'snps'}, $ascat_loh, $options->{'out'}, $sample, $options->{'downsamp'});
my ($reported_sample, $contamination) = run_verifyBam($filtered_vcf, $bam_in);
print join("\t", $reported_sample, $contamination),"\n";

sub run_verifyBam {
  my ($snps, $bam) = @_;
  my $verifyBamID = "$Bin/verifyBamID";
  $verifyBamID = which('verifyBamID') unless(-e $verifyBamID);
  die "Unable to find verifyBamID executable\n" unless(-e $verifyBamID);
  my ($out_stub) = $snps =~ m/(.*)_snps[.]vcf$/;
  my $command = sprintf $VERIFY, $verifyBamID, $snps, $bam, $out_stub;
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "An error occurred while executing:\n\t$command\nERROR: $stderr\n" if($exit);
  my $result = "$out_stub.selfSM";
  open my $RES, '<', $result;
  my $head = <$RES>;
  my ($this_sample, $contam) = (split /\t/, <$RES>)[0,6];
  close $RES;
  return ($this_sample, $contam);
}

sub sample_name {
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

sub get_loh_regions {
  my ($ascat) = @_;
  my %loh_by_chr;
  if(defined $ascat) {
    open my $SEG, '<', $ascat;
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
  my ($loci, $loh_regions, $outdir, $sample, $one_in_x) = @_;
  make_path($outdir) unless(-d $outdir);
  my $filtered_vcf = "$outdir/${sample}_snps.vcf";

  open my $VCF, '>', $filtered_vcf;
  my $z = new IO::Uncompress::Gunzip $loci;
  my $header = <$z>;
  print $VCF $header or die $!;
  my $last_chr = q{.};
  my @ranges;
  my ($excluded_snps, $total_snps, $included_snps) = (0,0,0);
  LINE: while(my $line = <$z>) {
    $total_snps++;
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
        $excluded_snps++;
        next LINE;
      }
    }
    $included_snps++;
    print $VCF $line or die $!;
  }
  close $z;
  close $VCF;
  warn "Total SNPs: $total_snps\n";
  warn "Excluded as HOM: $excluded_snps\n";
  warn "Retained SNPs: $included_snps\n";
  return $filtered_vcf;
}

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'b|bam=s' => \$opts{'bam'},
              'o|out_prefix=s' => \$opts{'out'},
              'l|loh_bed=s' => \$opts{'loh'},
              's|snps=s' => \$opts{'snps'},
              'd|downsamp=i' => \$opts{'downsamp'},
  ) or pod2usage(2);
}

__END__

=head1 NAME

verifyBamHomChk.pl - Runs verify BAM

=head1 SYNOPSIS

verifyBamHomChk.pl [options]

  Required parameters:

    -out_prefix  -o   Stub for output filenames

    -bam         -b   BAM file to process

  Optional parameters:

    -snps        -s   VCF file of SNP locations as described here:
                        http://genome.sph.umich.edu/wiki/VerifyBamID
                          ([b]gzip compressed works too)
                        defaults to included GRCh37 SNP6 loci.

    -downsamp    -d   Subsample the SNPs at rate of 1 in X [1]

    -loh_bed     -l   Exclude SNPs from these ranges []

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
