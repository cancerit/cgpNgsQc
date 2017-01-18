# /software/perl-5.16.3/bin/prove -v -I /software/CGP/canpipe/live/lib/perl5 -I lib t/QCvalidateSampleMeta.t
use strict;
use Test::More;
use Test::Fatal;
use Test::Exception;
use Const::Fast qw(const);
use File::Spec;
use FindBin qw($Bin);
use File::Temp;
use Data::Dumper;

use Bio::DB::HTS;

const my $MODULE => 'Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader';

my $test_data = "$Bin/TestData/sampleValidation";
my $test_xlsx_file = File::Spec->catfile($test_data,'test-meta.xlsx');
my $test_sample_1_file = File::Spec->catfile($test_data,'test_sample.1.bam');

subtest 'Initialisation checks' => sub {
  local $SIG{__WARN__} = sub {};
  my $tmp = File::Temp->new(TEMPLATE => '/tmp/sampleValidation_outputXXXXXX', SUFFIX => '.txt');
  use_ok($MODULE);
  new_ok($MODULE => [{ i => $test_xlsx_file, o => $tmp, f => 'xlsx', a => 0}]);
};

subtest 'bam_has_reads_more_than_threshold' => sub {
  local $SIG{__WARN__} = sub {};
  my (@got, @expect);
  @expect = (1, 100);
  @got = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::bam_has_reads_more_than_threshold($test_sample_1_file, 100);
  is_deeply(\@got, \@expect, 'if more than 100');
  @expect = (0, 222);
  @got = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::bam_has_reads_more_than_threshold($test_sample_1_file, 300);
  is_deeply(\@got, \@expect, 'if more than 300');
};

subtest 'validate_bam'  => sub {
  # local $SIG{__WARN__} = sub {};
  my $got;
  $got = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_bam('test', 'random', 'file should not exist please', 0, 'GRCh37d5', 100);
  is($got,'failed:invalid is_normal value; bam is not found', 'bad is_normal and if input exists');
  $got = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_bam('test', 'y', $test_sample_1_file, 0, 'GRCh37d5', 100);
  is($got,'failed:not all RG found, use --check-all', 'missing RG if requires check all');
  $got = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_bam('test', 'y', $test_sample_1_file, 1, 'GRCh37d5', 100);
  is($got,'pp-remapped', 'use check all to find all RG');
};

subtest 'validate_samples'  => sub {
  local $SIG{__WARN__} = sub {};
  my (@got, $test_input, @expect, @temp);

  # input with header
  $test_input = $test_xlsx_file;
  @expect = (
    "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\ttest_col_2ID2\trelative_file_path\tbam_state\tDonor_UUID\tTissue_UUID\tSample_UUID\tbam_md5sum",
    "fake donor 1\tfake tissue 1\tn\ttest\tt\ttest_sample.1.bam\tpp-remapped\t.*\td795d6494fea8527f2bd50d9accf4214",
    "fake donor 1\tfake tissue 1\ty\ttest 2\t\ttest_sample.2.bam\traw\t.*\t468436c355465c8440e5c7d26e169a85"
  );
  @temp = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 1, count_base_number => 100, genome_build => 'GRCh37d5'});
  @got = @{$temp[0]};
  is($got[0], $expect[0], 'input with header, good input, header');
  like($got[1], qr/$expect[1]$/, 'input with header, good input, output line');
  like($got[2], qr/$expect[2]$/, 'input with header, good input, output line, raw bam');
  is($temp[1], 1, 'input with header, good input, traffic_light');

  # not to check all
  @temp = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'});
  @expect = (
    "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\ttest_col_2ID2\trelative_file_path\tbam_state",
    "fake donor 1\tfake tissue 1\tn\ttest\tt\ttest_sample.1.bam\tfailed:not all RG found, use --check-all",
    "fake donor 1\tfake tissue 1\ty\ttest 2\t\ttest_sample.2.bam\traw"
  );
  is_deeply($temp[0], \@expect, 'input with header, good input, not to check all, output line');
  is($temp[1], 0, 'input with header, good input, not to check all, traffic_light');

  # input has no header
  $test_input = File::Spec->catfile($test_data,'test-meta_no_header.xlsx');
  @expect = (
    "Donor_ID\tTissue_ID\tis_normal (Yes/No, Y/N)\tSample_ID\trelative_file_path\tbam_state\tDonor_UUID\tTissue_UUID\tSample_UUID\tbam_md5sum",
    "fake donor 1\tfake tissue 1\tn\ttest\ttest_sample.1.bam\tpp-remapped\t.*\td795d6494fea8527f2bd50d9accf4214"
  );
  @temp = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 1, count_base_number => 100, genome_build => 'GRCh37d5'});
  @got = @{$temp[0]};
  is($got[0], $expect[0], 'input has no header, good input, header');
  like($got[1], qr/$expect[1]$/, 'input has no header, good input, output line');
  is($temp[1], 1, 'input has no header, good input, traffic_light');

  $test_input = File::Spec->catfile($test_data,'test-meta_invalid_isNormalColumn.xlsx');
  @expect = (
    "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\ttest_col_2ID2\trelative_file_path\tbam_state",
    "fake donor 1\tfake tissue 1\tZ\ttest\tt\ttest_sample.1.bam\tfailed:invalid is_normal value"
  );
  @temp = Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 1, count_base_number => 100, genome_build => 'GRCh37d5'});
  @got = @{$temp[0]};
  is($got[0], $expect[0], 'input with header, good input, bad is_normal, header');
  like($got[1], qr/$expect[1]$/, 'input with header, good input, bad is_normal, output line');
  is($temp[1], 0, 'input with header, good input, bad is_normal, output line, traffic_light');


  # test dies
  $test_input = File::Spec->catfile($test_data,'does no exist.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/Error: can not find input/, 'invalid input file';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'invalid input file, expecting to die';

  $test_input = File::Spec->catfile($test_data,'test-meta_no_sampleIDcolumn.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/Error: Sample_ID column is missing!/, 'no sample_id header';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'no sample_id header, expecting to die';

  $test_input = File::Spec->catfile($test_data,'test-meta_dup_donorIDcolumn.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/Error: Duplicated header Donor_ID found!/, 'dup donor_id header';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'dup donor_id header, expecting to die';

  $test_input = File::Spec->catfile($test_data,'test-meta_no_header_1columnShort.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/Error: Require 5 columns, 4 found!/, 'no header, a column short';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'no header, a column short, expecting to die';

  $test_input = File::Spec->catfile($test_data,'test-meta_dup_sampleID.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/Error: Duplicated Sample_ID found for /, 'dup sample id';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'dup sample id, expecting to die';

  $test_input = File::Spec->catfile($test_data,'test-meta_dup_bamPath.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/Error: Duplicated relative_file_path found: /, 'dup bam file';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'dup bam file, expecting to die';

  $test_input = File::Spec->catfile($test_data,'test-meta_invalid_bam.xlsx');
  throws_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } qr/$test_xlsx_file/, 'not a real bam';
  dies_ok { Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader::validate_samples({i => $test_input, f => 'xlsx', mod => 0, a => 0, count_base_number => 100, genome_build => 'GRCh37d5'}) } 'not a real bam, expecting to die';
};

done_testing();

__DATA__
