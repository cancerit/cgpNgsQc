use strict;
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use File::Spec;
use FindBin qw($Bin);
use File::Temp qw/ :seekable /;
use Data::Dumper;

use Bio::DB::HTS;

const my $MODULE => 'Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader';

my $test_data = "$Bin/../share/sampleValidation";
my $test_xlsx_file = join('/',$test_data,'test-meta.xlsx');
my $test_sample_7_file = join('/',$test_data,'test_sample_7.bam');
my $test_sample_8_file = join('/',$test_data,'test_sample_8.bam');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

subtest 'Per RG checks' => sub {
  my $tmp;
  my $test_obj = _create_test_object();

  eval{
    $tmp = File::Temp->new(TEMPLATE => '/tmp/Bam_StatsXXXXXX', SUFFIX => '.bas.tmp');
    my $fname = $tmp->filename;
    $tmp->unlink_on_destroy( 1 );
    $test_obj->bas($tmp);
    close $tmp;

    my $got = tsv_to_data($fname);
    my $expected = tsv_to_data($test_bas_file);

    is_deeply($got,$expected, "Per RG checks compared to: $test_bas_file") or BAIL_OUT "Merged results will be invalid if per-RG data is";
    1;
  };if($@){
    BAIL_OUT $@;
  }
};

subtest 'Non-object funcions' => sub {

  subtest '_add_to_qplot' => sub {
    my $target = []; ##
    my $target2 = [];
    my $qual_ref1 = [1,2,3,4,5]; ## an array of qualities
    my $qual_ref2 = [1,3,4,5,4]; ## an array of qualities
    my $reverse_y = 1;
    my $reverse_n = 0;

    my $target_after_n = [
      [undef,1],[undef,undef,1],[undef,undef,undef,1],[undef,undef,undef,undef,1],[undef,undef,undef,undef,undef,1]
    ];

    my $target_after_n2 = [
      [undef,2],[undef,undef,1,1],[undef,undef,undef,1,1],[undef,undef,undef,undef,1,1],[undef,undef,undef,undef,1,1]
    ];

    my $target_after_y = [
      [undef,undef,undef,undef,undef,1],[undef,undef,undef,undef,1],[undef,undef,undef,1],[undef,undef,1],[undef,1]
    ];
    my $target_after_y2 = [
      [undef,undef,undef,undef,1,1],[undef,undef,undef,undef,1,1],[undef,undef,undef,1,1],[undef,undef,1,1],[undef,2]
    ];

    PCAP::Bam::Stats::_add_to_qplot($target, $qual_ref1, $reverse_n);

    is_deeply($target, $target_after_n, '_add_to_qplot');
    PCAP::Bam::Stats::_add_to_qplot($target, $qual_ref2, $reverse_n);
    is_deeply($target, $target_after_n2, '_add_to_qplot');

    PCAP::Bam::Stats::_add_to_qplot($target2, $qual_ref1, $reverse_y);
    is_deeply($target2, $target_after_y, '_add_to_qplot');
    PCAP::Bam::Stats::_add_to_qplot($target2, $qual_ref2, $reverse_y);
    is_deeply($target2, $target_after_y2, '_add_to_qplot');
  };

  subtest 'down_pop_quals' => sub {
    my $plot_vals = [
      [1,2,3,4],
      [5,undef,6,7,8]
    ]; ## some vals to process

    my $plot_vals_after = [
      [10,9,7,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
      [26,21,21,15,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ]; ## what the array_ref will look like after processing

    PCAP::Bam::Stats::down_pop_quals($plot_vals);

    is_deeply($plot_vals, $plot_vals_after, 'down_pop_quals');
  };
};

subtest 'Object funcions' => sub {

  subtest 'calc_frac_' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->calc_frac_properly_paired()),'0.6429','calc_frac_properly_paired');
    is(sprintf('%.4f',$test_obj->calc_frac_unmapped()),'0.2500','calc_frac_unmapped');
    is(sprintf('%.4f',$test_obj->calc_frac_duplicate_reads()),'0.1667','calc_frac_duplicate_reads');
  };

  subtest 'med_insert_size' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->med_insert_size(),100,'med_insert_size');
  };

  subtest 'mean_insert_size' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->mean_insert_size()),'262.8889','mean_insert_size');
  };

  subtest 'properly_mapped_ratio' => sub {
    my $test_obj = _create_test_object();
    is(sprintf('%.4f',$test_obj->properly_mapped_ratio()),'0.6429','properly_mapped_ratio');
  };

  subtest '_count_' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->count_total_reads(1),14,'count_total_reads read1');
    is($test_obj->count_total_reads(2),10,'count_total_reads read2');

    is($test_obj->count_properly_paired(),9,'count_properly_paired');

    is($test_obj->count_unmapped(1),3,'count_unmapped read1');
    is($test_obj->count_unmapped(2),3,'count_unmapped read2');

    is($test_obj->count_duplicate_reads(1),4,'count_duplicate_reads read1');
    is($test_obj->count_duplicate_reads(2),0,'count_duplicate_reads read2');

    is($test_obj->count_total_mapped_bases(1),270,'count_total_mapped_bases read1');
    is($test_obj->count_total_mapped_bases(2),192,'count_total_mapped_bases read2');

    is($test_obj->count_total_divergent_bases(1),42,'count_total_divergent_bases read1');
    is($test_obj->count_total_divergent_bases(2),30,'count_total_divergent_bases read2');
  };

  subtest 'cal_insert_size_sd' => sub {
     my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->cal_insert_size_sd()),'278.0310','cal_insert_size_sd');
  };

};

sub tsv_to_data {
  my $file = shift;
  my @data;
  open (my $fh, '<', $file) or fail "Could not open $file for reading: $!";
  while(my $line = <$fh>) {
    chomp $line;
    push @data, [split /\t/, $line];
  }
  close($fh);
  return \@data;
}

done_testing();

__DATA__
