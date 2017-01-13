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

my $test_data = "$Bin/../share/sampleValidation";
my $test_xlsx_file = join('/',$test_data,'test-meta.xlsx');
my $test_sample_7_file = join('/',$test_data,'test_sample_7.bam');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  # local $SIG{__WARN__} = sub {};
  new_ok($MODULE => [{ i => $test_xlsx_file, o => 'test_output', f => 'xlsx', a => 0}]);
};

# subtest 'bam_has_reads_more_than_threshold' => sub {
#   local $SIG{__WARN__} = sub {};
#   my $got;
#   $got = $MODULE::bam_has_reads_more_than_threshold($test_sample_7_file, 1_000_000);
#   is($got,'yes','if more than 1,000,000');
#   $got = $MODULE::bam_has_reads_more_than_threshold($test_sample_7_file, 2_000_000);
#   is($got,1006738,'if more than 2,000,000');
# };
#
# subtest 'validate_bam'  => sub {
#   local $SIG{__WARN__} = sub {};
#   my $got;
#   $got = $MODULE::validate_bam('test', 'random', 'file should not exist please', 0);
#   is($got,'failed:invalid is_normal value; bam is not found;', 'bad is_normal and if input exists');
#   $got = $MODULE::validate_bam('test', 'y', $test_sample_7_file, 0);
#   is($got,'failed:not all RG found, use --check-all;', 'missing RG if requires check all');
#   $got = $MODULE::validate_bam('test', 'y', $test_sample_7_file, 1);
#   is($got,'pp-remapped', 'use check all to find all RG');
# };
#
# subtest 'validate_samples'  => sub {
#   local $SIG{__WARN__} = sub {};
#   my (@got, @test_input, @expect, @temp);
#
#   # input with header
#   @test_input = (
#     "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\trelative_file_path",
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file"
#   );
#   @expect = (
#     "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\trelative_file_path\tbam_state\tDonor_UUID\tTissue_UUID\tSample_UUID\tbam_md5sum",
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file\tpp-remapped\t.*\t5612e7ea0e88252900c377406af165d8"
#   );
#   @temp = $MODULE::validate_samples({i => \@test_input, mod => 0, a => 1});
#   @got = @{$temp[0]};
#   is($got[0], $expect[0], 'input with header, good input, header');
#   like($got[1], qr/$expect[1]$/, 'input with header, good input, output line');
#   is($temp[1], 1, 'input with header, good input, traffic_light');
#
#   # not to check all
#   @temp = $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0});
#   @expect = (
#     "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\trelative_file_path\tbam_state",
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file\tfailed:not all RG found, use --check-all;"
#   );
#   is_deeply($temp[0], \@expect, 'input with header, good input, not to check all, output line');
#   is($temp[1], 0, 'input with header, good input, not to check all, traffic_light');
#
#   # input has no header
#   @test_input = (
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file"
#   );
#   @expect = (
#     "Donor_ID\tTissue_ID\tis_normal (Yes/No, Y/N)\tSample_ID\trelative_file_path\tbam_state\tDonor_UUID\tTissue_UUID\tSample_UUID\tbam_md5sum",
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file\tpp-remapped\t.*\t5612e7ea0e88252900c377406af165d8"
#   );
#   @temp = $MODULE::validate_samples({i => \@test_input, mod => 0, a => 1});
#   @got = @{$temp[0]};
#   is($got[0], $expect[0], 'input has no header, good input, header');
#   like($got[1], qr/$expect[1]$/, 'input has no header, good input, output line');
#   is($temp[1], 1, 'input has no header, good input, traffic_light');
#
#   # test dies
#   @test_input = (
#     "#Tissue_ID\tis_normal\tSample_ID\trelative_file_path",
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file"
#   );
#   throws_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } qr/Error: Donor_ID column is missing!/, 'no donor_id header';
#   dies_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } 'no donor_id header, expecting to die';
#   @test_input = (
#     "#Donor_ID\tDonor_ID\tis_normal\tSample_ID\trelative_file_path",
#     "fake donor 1\tfake tissue\tY\ttest\t$test_sample_7_file"
#   );
#   throws_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } qr/Error: Duplicated header Donor_ID found!/, 'dup donor_id header';
#   dies_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } 'dup donor_id header, expecting to die';
#   @test_input = (
#     "fake donor 1\tfake tissue\ttest\t$test_sample_7_file"
#   );
#   throws_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } qr/Error: Require 5 columns, 4 found!/, 'no header, a column short';
#   dies_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } 'no header, a column short, expecting to die';
#   @test_input = (
#   "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\trelative_file_path",
#   "fake donor 1\tfake tissue\tY\ttest\tbam file 1",
#   "fake donor 1\tfake tissue\tY\ttest2\tbam file 2",
#   "fake donor 1\tfake tissue1\tN\ttest\t$test_sample_7_file"
#   );
#   throws_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } qr/Error: Duplicated Sample_ID found for /, 'dup sample id';
#   dies_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } 'dup sample id, expecting to die';
#   @test_input = (
#   "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\trelative_file_path",
#   "fake donor 1\tfake tissue\tY\ttest\tbam file 1",
#   "fake donor 1\tfake tissue\tY\ttest2\tbam file 2",
#   "fake donor 1\tfake tissue1\tN\ttest3\tbam file 1"
#   );
#   throws_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } qr/Error: Duplicated relative_file_path found: /, 'dup bam file';
#   dies_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } 'dup bam file, expecting to die';
#   @test_input = (
#   "#Donor_ID\tTissue_ID\tis_normal\tSample_ID\trelative_file_path",
#   "fake donor 1\tfake tissue\tY\ttest\t$test_xlsx_file"
#   );
#   throws_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } qr/$test_xlsx_file/, 'not a real bam';
#   dies_ok { $MODULE::validate_samples({i => \@test_input, mod => 0, a => 0}) } 'not a real bam, expecting to die';
# };

done_testing();

__DATA__
