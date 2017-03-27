package Sanger::CGP::NgsQc::ValidateSampleMetaAndBamHeader;

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
#
# Author: Yaobo Xu <cgpit@sanger.ac.uk>
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

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use Carp;
use Cwd;

use Sanger::CGP::NgsQc;

use File::ShareDir qw(dist_dir);
use List::Util qw(first);
use FindBin qw($Bin);
use Spreadsheet::Read;
use Spreadsheet::WriteExcel; # to write xls output
use Excel::Writer::XLSX; # to write xlsx output
use Bio::DB::HTS;
use File::Spec;
use Data::UUID;
use Digest::file qw(digest_file_hex);


const my $HEADER_DONOR_ID => 'Donor_ID';
const my $HEADER_TISSUE_ID => 'Tissue_ID';
const my $HEADER_IS_NORMAL => 'is_normal';
const my $HEADER_IS_NORMAL_FOR_DONOR => 'is_normal_for_donor';
const my $HEADER_SAMPLE_ID => 'Sample_ID';
const my $HEADER_BAM => 'relative_file_path';
const my $HEADER_DONOR_UUID => 'Donor_UUID';
const my $HEADER_TISSUE_UUID => 'Tissue_UUID';
const my $HEADER_SAMPLE_UUID => 'Sample_UUID';
const my $HEADER_BAM_MD5 => 'bam_md5sum';
const my $HEADER_BAM_STATE => 'bam_state';

const my %FAIL_CODES => (
  '010' => 'invalid is_normal value',
  '011' => 'invalid is_normal_for_donor value',
  '012' => 'no nominated normal_for_donor',
  '013' => 'nominated normal_for_donor is not normal',
  '014' => 'multiple nominated normal_for_donor',
  '015' => 'no normal sample for donor',
  '016' => 'no tumour sample for donor',
  '020' => 'bam is not found',                           # can not find the bam file;
  '021' => 'bam has too few reads',                      # bam has read fewer than a threshold;
  '030' => 'no RG line',                                 # bam header has no RG line;
  '031' => 'duplicated RG ID in header',                 # bam header has duplicatd RG IDs;
  '040' => 'no ID',                                      # Each RG line has to have an unique ID (of other RG lines);
  '060' => 'no LB',                                      # Each RG line has to have an LB tag;
  '070' => 'no PL',                                      # Each RG line has to have an PL tag;
  '080' => 'invalid PL',                                 # Valid PL values are: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and PACBIO. Error will be give for other values.
  '090' => 'no SM',                                      # Each RG line has to have an SM tag.
  '100' => 'SM not matching',                            # In one of the RG lines, SM value does not match to the Sample_ID in the excel file;
  '110' => 'read has no RG ID',                          # One or more reads have no RG ID.
  '120' => 'read has invalid RG ID',                     # One or more reads have RG ID not defined in the header.
  '140' => 'checked all, not all RG found',              # Checked all reads in the bam not all RG IDs were found.
);

#const my @VALID_PLATFORM => qw(CAPILLARY ILLUMINA LS454 SOLID HELICOS IONTORRENT ONT PACBIO);
const my @VALID_PLATFORM => qw(ILLUMINA); # can only process illinina data
const my $VALID_IS_NORMALS => qr/^yes$|^y$|^no$|^n$/i;
const my $VALID_IS_NORMALS_IS => qr/^yes$|^y$/i;
const my $VALID_IS_NORMAL_FOR_DONORS => qr/^yes$|^y$|^$/i;
const my $VALID_IS_NORMAL_FOR_DONORS_IS => qr/^yes$|^y$/i;
const my $FMT_ERR => qq{--------\nError: %s\n--------\n};

1;

sub new {
  my ($class, $options) = @_;
  my $self = { };
  bless $self, $class;
  $self->_init($options);
  return $self;
}

sub _init {
  my ($self, $options) = @_;
  croak sprintf $FMT_ERR, "can not find input file $options->{'i'}. $!" unless (-e $options->{'i'});
  open my $OUT, ">", $options->{'o'} or croak sprintf $FMT_ERR, "can not open the output: $options->{'o'} to write. $!";
  close $OUT;
  $options->{'mod'} = 1; # not in testing mode
  $options->{'count_base_number'} = 1_000_000;
  $options->{'genome_build'} = 'GRCh37d5' unless (exists $options->{'genome_build'}); #TODO may need a dict for this option
  my $validate_status = validate_samples($options);
  $self->{'validate_status'} = $validate_status;
  return 1;
}

sub input_to_array {
  my ($input, $format) =@_;
  my @in;
  ### read in input as array of lines (concatenated columns with tabs)
  if ($format eq 'xlsx'|| $format eq 'xls') {
    my $book = ReadData($input);
    my @rows = Spreadsheet::Read::rows($book->[1]);
    foreach my $row (@rows) {
      foreach my $cell (@$row) {
        $cell = '' if (! defined $cell);
        # print $cell."\n";
      }
      push @in, (join "\t", @$row);
    }
  } else { # assume it's a tsv
    open my $FH, '<', $input or confess sprintf $FMT_ERR, "trying to open input file $input: $!";
    chomp(@in = <$FH>);
    close $FH;
  }
  return \@in;
}

sub validate_samples {
  my ($opts) = @_;
  my $in_array;
  if (! -e $opts->{'i'}) {
    confess sprintf $FMT_ERR, "can not find input $opts->{'i'}";
  }
  $in_array = input_to_array($opts->{'i'}, $opts->{'f'});

  my @lines = @$in_array;

  my ($donor_index, $tissue_index, $is_normal_index, $is_normal_for_donor_index, $sample_index, $bam_index);
  my @out;

  ### check input table headers
  my $header_line = '';

  if ($lines[0] =~ m/^\#|^DONOR_ID/i) {
    warn "--------\nTake the first line as header line.\n";
    $header_line = $lines[0]; shift @lines;
    # assign column index
    my @colnames = split "\t", $header_line;
    my %header_checks;
    for my $i (0..$#colnames) {
      if ($colnames[$i] =~ m/$HEADER_DONOR_ID/i) {
        $donor_index = $i;
        if (! exists $header_checks{$HEADER_DONOR_ID}) {
          warn "Column $colnames[$i]:$donor_index -> $HEADER_DONOR_ID\n";
          $header_checks{$HEADER_DONOR_ID} = 1;
        } else {
          croak sprintf $FMT_ERR, "Duplicated header $HEADER_DONOR_ID found!";
        }
      } elsif ($colnames[$i] =~ m/$HEADER_TISSUE_ID/i) {
        $tissue_index = $i;
        warn "Column $colnames[$i]:$tissue_index -> $HEADER_TISSUE_ID\n";
        if (! exists $header_checks{$HEADER_TISSUE_ID}) {
          $header_checks{$HEADER_TISSUE_ID} = 1;
        } else {
          croak sprintf $FMT_ERR, "Duplicated header $HEADER_TISSUE_ID found!";
        }
      } elsif ($colnames[$i] =~ m/^$HEADER_IS_NORMAL$/i) {
        $is_normal_index = $i;
        warn "Column $colnames[$i]:$is_normal_index -> $HEADER_IS_NORMAL\n";
        if (! exists $header_checks{$HEADER_IS_NORMAL}) {
          $header_checks{$HEADER_IS_NORMAL} = 1;
        } else {
          croak sprintf $FMT_ERR, "Duplicated header $HEADER_IS_NORMAL found!";
        }
      } elsif ($colnames[$i] =~ m/^$HEADER_IS_NORMAL_FOR_DONOR$/i) {
        $is_normal_for_donor_index = $i;
        warn "Column $colnames[$i]:$is_normal_for_donor_index -> $HEADER_IS_NORMAL_FOR_DONOR\n";
        if (! exists $header_checks{$HEADER_IS_NORMAL_FOR_DONOR}) {
          $header_checks{$HEADER_IS_NORMAL_FOR_DONOR} = 1;
        } else {
          croak sprintf $FMT_ERR, "Duplicated header $HEADER_IS_NORMAL_FOR_DONOR found!";
        }
      } elsif ($colnames[$i] =~ m/$HEADER_SAMPLE_ID/i) {
        $sample_index = $i;
        warn "Column $colnames[$i]:$sample_index -> $HEADER_SAMPLE_ID\n";
        if (! exists $header_checks{$HEADER_SAMPLE_ID}) {
          $header_checks{$HEADER_SAMPLE_ID} = 1;
        } else {
          croak sprintf $FMT_ERR, "Duplicated header $HEADER_SAMPLE_ID found!";
        }
      } elsif ($colnames[$i] =~ m/$HEADER_BAM/i) {
        $bam_index = $i;
        warn "Column $colnames[$i]:$bam_index -> $HEADER_BAM\n";
        if (! exists $header_checks{$HEADER_BAM}) {
          $header_checks{$HEADER_BAM} = 1;
        } else {
          croak sprintf $FMT_ERR, "Duplicated header $HEADER_BAM found!";
        }
      }
    }

    if (!exists $header_checks{$HEADER_DONOR_ID}) {croak sprintf $FMT_ERR, "$HEADER_DONOR_ID column is missing!";};
    if (!exists $header_checks{$HEADER_TISSUE_ID}) {croak sprintf $FMT_ERR, "$HEADER_TISSUE_ID column is missing!";};
    if (!exists $header_checks{$HEADER_IS_NORMAL}) {croak sprintf $FMT_ERR, "$HEADER_IS_NORMAL column is missing!";};
    if (!exists $header_checks{$HEADER_IS_NORMAL_FOR_DONOR}) {croak sprintf $FMT_ERR, "$HEADER_IS_NORMAL_FOR_DONOR column is missing!";};
    if (!exists $header_checks{$HEADER_SAMPLE_ID}) {croak sprintf $FMT_ERR, "$HEADER_SAMPLE_ID column is missing!";};
    if (!exists $header_checks{$HEADER_BAM}) {croak sprintf $FMT_ERR, "$HEADER_BAM column is missing!";};

  } else {
    warn "Found no header, assume following:\n  column 1: $HEADER_DONOR_ID\n  column 2: $HEADER_TISSUE_ID\n  column 3: $HEADER_IS_NORMAL\n  column 4: $HEADER_IS_NORMAL_FOR_DONOR\n  column 5: $HEADER_SAMPLE_ID\n  column 6: $HEADER_BAM\n";
    my @temp = split "\t", $lines[0];
    if (scalar @temp >= 6) {
      ($donor_index, $tissue_index, $is_normal_index, $is_normal_for_donor_index, $sample_index, $bam_index) = (0, 1, 2, 3, 4, 5);
      $header_line = join("\t", ($HEADER_DONOR_ID, $HEADER_TISSUE_ID, $HEADER_IS_NORMAL, $HEADER_IS_NORMAL_FOR_DONOR, $HEADER_SAMPLE_ID, $HEADER_BAM))."\t" x (scalar @temp - 6);
    } else {
      croak sprintf $FMT_ERR, "Require 6 columns, ".scalar @temp." found!";
    }
  }
  warn "--------\n\n";

  my $ncol = scalar split("\t", $header_line);
  my ($in_volume,$in_directories,$in_file_name) = File::Spec->splitpath($opts->{'i'});
  ### store IDs and check if there's duplicated sample_id or relative_file_path
  my (%sample_id_checks, %bam_path_checks, @donor_ids, @tissue_ids, @is_normals, @is_normal_for_donors, @sample_ids, @bam_files);
  for my $i(0..$#lines) {
    my @temp = split "\t", $lines[$i];
    if (scalar @temp == $ncol) {
      my ($sample_id, $bam) = ($temp[$sample_index], $temp[$bam_index]);
      $bam =~ s/^\s+|\s+$//g; #trimming off leading and tailing spaces in file path
      $bam = File::Spec->catpath($in_volume, $in_directories, $bam);
      if (! exists $sample_id_checks{$sample_id}) {
        $sample_id_checks{$sample_id} = 1; push @sample_ids, $sample_id;
      } else {
        croak sprintf $FMT_ERR, "Duplicated $HEADER_SAMPLE_ID found for $sample_id!";
      }
      if (! exists $bam_path_checks{$bam}) {
        $bam_path_checks{$bam} = 1; push @bam_files, $bam;
      } else {
        croak sprintf $FMT_ERR, "Duplicated $HEADER_BAM found: $bam!";
      }
      push @is_normals, $temp[$is_normal_index];
      push @is_normal_for_donors, $temp[$is_normal_for_donor_index];
      push @donor_ids, $temp[$donor_index];
      push @tissue_ids, $temp[$tissue_index];
    } else {
      croak sprintf $FMT_ERR, "Expected $ncol columns, found ".scalar @temp." on line:$lines[$i]";
    }
  }

  ### validate bam headers
  my $traffic_light = 1; # traffic light for UUID generation and md5sum check
  my @bam_states;
  if (scalar @sample_ids == scalar @lines && scalar @sample_ids == scalar @bam_files && scalar @sample_ids == scalar @donor_ids
    && scalar @sample_ids == scalar @tissue_ids && scalar @sample_ids == scalar @is_normals && scalar @sample_ids == scalar @is_normal_for_donors) {
    for my $i(0..$#sample_ids) {
      warn "#--- processing sample: $sample_ids[$i]\n";
      my $bam_state = validate_bam($sample_ids[$i], $is_normals[$i], $is_normal_for_donors[$i], $bam_files[$i], $opts->{'genome_build'}, $opts->{'count_base_number'});
      push @bam_states, $bam_state;
      warn "done! ---#\n";
    }
    my $new_bam_states_ref = check_donor_states(\@donor_ids, \@is_normals, \@is_normal_for_donors, \@bam_states, );
    @bam_states = @$new_bam_states_ref;
    push @out, "$header_line\t$HEADER_BAM_STATE";
    for my $i(0..$#sample_ids) {
      if ($bam_states[$i] eq 'pp-remapped') { # if read group checks are ok
        warn "ovarall validate states for $sample_ids[$i]: pp_remapped_bam!\n";
      } elsif ($bam_states[$i] eq 'raw') {
        warn "ovarall validate states for $sample_ids[$i]: raw_bam!\n";
      } else {
        warn "ovarall validate states for $sample_ids[$i]: invalid bam! $bam_states[$i]\n";
      }
      push @out, "$lines[$i]\t$bam_states[$i]";
    }
  } else {
    confess sprintf $FMT_ERR, "script error 001: Different length of sample_ids and bams!"; # should never happen
  }

  for my $bam_state (@bam_states) {
    if ($bam_state !~ m/^raw$|^pp\-remapped$/) {
      $traffic_light = 0;
      last;
    }
  }

  ### generate UUIDs and md5sums
  if ($traffic_light == 1) {
    #start to produce UUIDs and md5sums
    @out = ();
    push @out, join("\t", ($header_line, $HEADER_BAM_STATE, $HEADER_DONOR_UUID, $HEADER_TISSUE_UUID, $HEADER_SAMPLE_UUID, $HEADER_BAM_MD5));
    my (%exist_donor_uuid, %exist_tissue_uuid);
    for (my $i = 0;$i < scalar @bam_files; $i++) {
      my ($donor_uuid, $tissue_uuid, $sample_uuid, $md5sum);
      if (! exists $exist_donor_uuid{$donor_ids[$i]}) {
        warn "--------\ngenerating donor UUID for: $donor_ids[$i]\n";
        $donor_uuid = get_uuid();
        $exist_donor_uuid{$donor_ids[$i]} = $donor_uuid;
      } else {
        warn "--------\nuse existing donor UUID for: $donor_ids[$i]\n";
        $donor_uuid = $exist_donor_uuid{$donor_ids[$i]};
      }

      if (! exists $exist_tissue_uuid{$tissue_ids[$i]}) {
        warn "generating tissue UUID for: $tissue_ids[$i]\n";
        $tissue_uuid = get_uuid();
        $exist_tissue_uuid{$tissue_ids[$i]} = $tissue_uuid;
      } else {
        warn "use existing tissue UUID for: $tissue_ids[$i]\n";
        $tissue_uuid = $exist_tissue_uuid{$tissue_ids[$i]};
      }
      warn "generating sample UUID for: $sample_ids[$i]\n";
      $sample_uuid = get_uuid();

      my $left = scalar @bam_files - $i - 1;
      warn "generating md5sum for: $bam_files[$i]. [$left bam(s) to go]\n";
      $md5sum = get_md5sum($bam_files[$i]);
      push @out, join("\t", ($lines[$i], $bam_states[$i], $donor_uuid, $tissue_uuid, $sample_uuid, $md5sum));
      warn "done!\n--------\n";
    }
  }

  if ($opts->{'mod'} == 1) {
    write_results(\@out, $opts->{'o'}, $opts->{'f'}, $traffic_light);
    return $traffic_light;
  } else { # return array ref when in test mode, instead of writing to file
    return \@out, $traffic_light;
  }
}

sub get_md5sum {
  my $file = shift;
  my $cksum = digest_file_hex($file, "MD5");
  if ($cksum) {
    return $cksum;
  } else {
    return 0;
  }
}

sub get_uuid {
  my $ug = Data::UUID->new;
  my $uuid = $ug->to_string($ug->create);
  return lc $uuid;
}

sub validate_bam {
  my ($sample_ID, $is_normal, $is_normal_for_donor, $align_file, $genome_build, $count_base_number) = @_; # $count_base_number is for testing the script
  my @bam_state;

  if ($is_normal !~ $VALID_IS_NORMALS) {
    push @bam_state, $FAIL_CODES{'010'};
  }

  if ($is_normal_for_donor !~ $VALID_IS_NORMAL_FOR_DONORS) {
    push @bam_state, $FAIL_CODES{'011'};
  }

  if (! -e $align_file) {
    warn "can not find the file: $align_file\n";
    return 'failed:'.join '; ', @bam_state, $FAIL_CODES{'020'};
  }

  warn "--------\nstart to validate alignment file: $align_file\n--------\n";
  my $header_sets = bam_header_to_hash($align_file);
  my %headers; # to handle duplicates!

  if (exists $header_sets->{RG}) {
    my %reported_error;
    for my $rg_line (@{$header_sets->{RG}}) {

      if (! exists $rg_line->{ID} && ! exists $reported_error{$FAIL_CODES{'040'}}) {
        warn "RG line has no ID.\n";
        $reported_error{$FAIL_CODES{'040'}} = 1;
        push @bam_state, $FAIL_CODES{'040'};
      } elsif (exists $rg_line->{ID}) {
        if (exists $headers{'RG'}{$rg_line->{ID}} && ! exists $reported_error{$FAIL_CODES{'031'}}) {
          warn "duplicated RG ID:$rg_line->{ID} in header.\n";
          $reported_error{$FAIL_CODES{'031'}} = 1;
          push @bam_state, $FAIL_CODES{'031'};
        } elsif (! exists $headers{'RG'}{$rg_line->{ID}}) {
          $headers{'RG'}{$rg_line->{ID}} = 0; # handle duplicates also = 1 when rg found in a read.
        }
      }

      if (! exists $rg_line->{PL} && ! exists $reported_error{$FAIL_CODES{'070'}}) {
        warn "RG line has no PL.\n";
        $reported_error{$FAIL_CODES{'070'}} = 1;
        push @bam_state, $FAIL_CODES{'070'};
      } elsif (! defined(first { $_ eq $rg_line->{PL} } @VALID_PLATFORM) && ! exists $reported_error{$FAIL_CODES{'080'}}) {
        warn "RG line has invalid PL.\n";
        $reported_error{$FAIL_CODES{'080'}} = 1;
        push @bam_state, $FAIL_CODES{'080'};
      }

      if (! exists $rg_line->{LB} && ! exists $reported_error{$FAIL_CODES{'060'}}) {
        warn "RG line has no LB\n";
        $reported_error{$FAIL_CODES{'060'}} = 1;
        push @bam_state, $FAIL_CODES{'060'};
      }

      if (! exists $rg_line->{SM} && ! exists $reported_error{$FAIL_CODES{'090'}}) {
        warn "RG line has no SM.\n";
        $reported_error{$FAIL_CODES{'090'}} = 1;
        push @bam_state, $FAIL_CODES{'090'};
      } elsif ($rg_line->{SM} ne $sample_ID && ! exists $reported_error{$FAIL_CODES{'100'}}) {
        warn "RG line has invalid SM.\n";
        $reported_error{$FAIL_CODES{'100'}} = 1;
        push @bam_state, $FAIL_CODES{'100'};
      }
    }
  } else {
    warn "no RG lines!\n";
    push @bam_state, $FAIL_CODES{'030'};
  }

  my ($enough_reads, $number_of_counted) = bam_has_reads_more_than_threshold($align_file, $count_base_number);
  if ($enough_reads == 0) {
    warn "Alignment file: $align_file has too few reads ($number_of_counted reads found).\n";
    push @bam_state, $FAIL_CODES{'021'};
  }

  my $match_pp_remapped_sq_pg = 0;

  if (scalar @bam_state == 0) { # check if it's a pp_remapped_bam
    my $valid_sq_hash_ref = sq_hash_from_bam_header($genome_build);
    my %valid_sq_hash = %$valid_sq_hash_ref;
    if(exists $header_sets->{SQ}) {
      for my $sq_line (@{$header_sets->{SQ}}) {
        if (exists $sq_line->{SN} && exists $sq_line->{LN} && exists $sq_line->{AS} && exists $sq_line->{M5} && exists $sq_line->{SP}) {
          if (! exists $headers{'SQ'}{$sq_line->{SN}} ) {
            $headers{'SQ'}{$sq_line->{SN}} = 1;
            if (exists $valid_sq_hash{$sq_line->{SN}} && $valid_sq_hash{$sq_line->{SN}}{LN} eq $sq_line->{LN} && $valid_sq_hash{$sq_line->{SN}}{AS} eq $sq_line->{AS} && $valid_sq_hash{$sq_line->{SN}}{M5} eq $sq_line->{M5} && $valid_sq_hash{$sq_line->{SN}}{SP} eq $sq_line->{SP}) {
              $match_pp_remapped_sq_pg = 1;
            } else {
              warn "sq line does not match to a pp_remapped_bam sq header line!\n";
              $match_pp_remapped_sq_pg = 0; last;
            }
          } else {
            warn "duplicated sq name found: $sq_line->{SN}!\n";
            $match_pp_remapped_sq_pg = 0; last;
          }
        } else {
          warn "SQ line has missing field!\n";
          $match_pp_remapped_sq_pg = 0; last;
        }
      }
      for (keys %valid_sq_hash) {
        if ( ! exists $headers{'SQ'}{$_} ) {
          warn "reference seq $_ is not found in the bam!\n";
          $match_pp_remapped_sq_pg = 0; last;
        }
      }
      for (keys %{$headers{'SQ'}}) {
        if ( ! exists $valid_sq_hash{$_} ) {
          warn "$_ is not a valid seq name in pp_remapped_bam requirements!\n";
          $match_pp_remapped_sq_pg = 0; last;
        }
      }
    }

    if(exists $header_sets->{PG} && $match_pp_remapped_sq_pg == 1 ) {
      my $valid_pg_hash_ref = valid_pp_remapped_pg_hash();
      my %valid_pg_hash = %$valid_pg_hash_ref;
      my @PG_bwa_lines;
      for my $pg_line (@{$header_sets->{PG}}) {
        if(exists $pg_line->{PN} && $pg_line->{PN} eq $valid_pg_hash{PN}) {
          push @PG_bwa_lines, $pg_line;
        }
      }

      if (scalar @PG_bwa_lines != 0) {
        for my $pg_line (@PG_bwa_lines) {
          if (exists $pg_line->{CL}) {
            my @words = split /\s+/, $pg_line->{CL};
            my %valid_words;
            my $count = 0;
            for my $i (0..$#words) {
              if ($words[$i] eq '-K' && ! exists $valid_words{'-K'}) {
                $valid_words{$words[$i]} = 1;
                if ($words[$i + 1] == $valid_pg_hash{CL}{'-K'}) {
                  $count += 1;
                }
              } else {
                if (! exists $valid_words{$words[$i]} && exists $valid_pg_hash{CL}{$words[$i]}) {
                  $valid_words{$words[$i]} = 1;
                  $count += 1;
                }
              }
            }
            if ($count == 6) {
              $match_pp_remapped_sq_pg = 1;
            } else {
              warn "bwa command line in head does not match pp_remapped_bam pg line!\n";
              $match_pp_remapped_sq_pg = 0; last;
            }
          } else {
            warn "PG CL is not found!\n";
            $match_pp_remapped_sq_pg = 0; last;
          }
        }
      } else {
        warn "No valid PG BWA header lines!\n";
        $match_pp_remapped_sq_pg = 0;
      }
    } else {
      warn "No valid PG header lines / invalid SQ lines\n";
      $match_pp_remapped_sq_pg = 0;
    }

    if ($match_pp_remapped_sq_pg == 1) {
      warn "this is a pp_remapped_bam!\n";
    } else {
      warn "this is a raw_bam!\n";
    }

    my $sam = Bio::DB::HTS->new(-bam => $align_file);
    my $bam = $sam->hts_file;
    my $header = $bam->header_read;
    my $processed_x = 0;
    my $start = time;
    my $processed_x_mill = 0;

    my %reported_error;

    while (my $a = $bam->read1($header)) {
      $processed_x += 1;

      my ($rg) = ($a->get_tag_values('RG'));
      if (defined $rg) {
        if (exists $headers{'RG'}{$rg}) {
          if ($headers{'RG'}{$rg} == 0) {
            $headers{'RG'}{$rg} = 1;
          }
        } else {
          if (! exists $reported_error{$FAIL_CODES{'120'}}) {
            warn "found a read with a RG tag not in header.\n";
            $reported_error{$FAIL_CODES{'120'}} = 1;
            push @bam_state, $FAIL_CODES{'120'};
          }
        }
      } else {
        if (! exists $reported_error{$FAIL_CODES{'110'}}) {
          warn "found a read with no RG tag.\n";
          $reported_error{$FAIL_CODES{'110'}} = 1;
          push @bam_state, $FAIL_CODES{'110'};
        }
      }

      if($processed_x % $count_base_number == 0) {
        $processed_x_mill++;
        my $end = time;
        my $elapsed = $end - $start;
        $start = $end;
        warn "$processed_x reads processed [time elapsed: ${elapsed}s].\n";
        last;
      }
    }


    my %not_found;
    for my $rg_id (keys %{$headers{'RG'}}) {
      if ($headers{'RG'}{$rg_id} == 0) {
        $not_found{$rg_id} = 0;
      }
    }

    if (scalar (keys %not_found) > 0) {
      warn "checked $count_base_number reads, reads with RG:".join(', ', keys %not_found)." not found, checking read groups of all reads, will take upto 2 hours, depending on the vloume of reads!\n";
      my $cmd = sprintf q@samtools view -F 80 '%s' | cut -f 12-@, $align_file;
      my %all_RG_IDs_from_reads;
      my ($pid, $process);
      $pid = open $process, q{-|}, $cmd;
      my $count = 0;
      while (my $tmp = <$process>) {
        chomp $tmp;
        $count += 1;
        my ($rgid) = $tmp =~ m/^.*RG\:Z\:([^\t]+)/;
        if ( ! exists $all_RG_IDs_from_reads{$rgid} ) {
          $all_RG_IDs_from_reads{$rgid} = 1;
          if (exists $not_found{$rgid}) {
            delete $not_found{$rgid};
          }
          if (scalar (keys %not_found) == 0) {
            warn "parsed $count * 2 reads, found reads of RG ID: $rgid, found all RG IDs.\n";
            last;
          } else {
            warn "parsed $count * 2 reads, found reads of RG ID: $rgid, still searching for ".join(', ', keys %not_found).".\n";
          }
        }
      }
      # close $process; # this close always return error with autodie, removing 'last' can get rid of the error as well
      {
        no autodie;
        unless (close($process)) {
          croak sprintf $FMT_ERR, "samtools - cut pipe close error: $!" if $!;
        }
      }

      for my $found_id (keys %all_RG_IDs_from_reads) {
        if (! exists $headers{'RG'}{$found_id}) {
          if (! exists $reported_error{$FAIL_CODES{'120'}}) {
            warn "found a read with a RG tag not in header.\n";
            $reported_error{$FAIL_CODES{'120'}} = 1;
            push @bam_state, $FAIL_CODES{'120'};
          }
        }
      }
      for my $header_id (keys %{$headers{'RG'}}) {
        if (! exists $all_RG_IDs_from_reads{$header_id}) {
          warn "$header_id is not in the reads.\n";
          warn "checked all reads, reads with RG:$header_id not found!\n";
          if (! exists $reported_error{$FAIL_CODES{'140'}}) {
            $reported_error{$FAIL_CODES{'140'}} = 1;
            push @bam_state, $FAIL_CODES{'140'};
          }
        }
      }
    } else {
      warn "processed $count_base_number reads in total, found all RG IDs in these reads.\n";
    }

  } else {
    warn "this is an invalid bam! Errors are: ".join('; ', @bam_state).".\n";
  }

  if (scalar @bam_state == 0 && $match_pp_remapped_sq_pg == 1) { # if read group checks are ok
    warn "bam_state conclusion: this is a pp_remapped_bam!\n";
    return "pp-remapped";
  } elsif (scalar @bam_state == 0 && $match_pp_remapped_sq_pg == 0) {
    warn "bam_state conclusion: this is a raw_bam!\n";
    return "raw";
  } else {
    warn "bam_state conclusion: this is an invalid bam! Errors are: ".join('; ', @bam_state).".\n";
    return 'failed:'.join '; ', @bam_state;
  }
}

sub valid_pp_remapped_pg_hash {
  my %valid_pp_remapped_pg_hash;
  $valid_pp_remapped_pg_hash{PN} = 'bwa';
  $valid_pp_remapped_pg_hash{CL} = {'/opt/wtsi-cgp/bin/bwa' => 1,  'mem' => 1, '-Y' => 1, '-K' => 100000000, '-p' => 1, '-R' => 1};
  return \%valid_pp_remapped_pg_hash;
}

sub sq_hash_from_bam_header {
  my $genome_build = shift;
  my %valid_sq_hash;
  my $data_path = "$Bin/../share";
  $data_path = dist_dir('Sanger-CGP-NgsQc') unless(-e "$data_path/bam_headers");
  $data_path .= "/bam_headers/$genome_build.header";
  my $headers = bam_header_to_hash($data_path);
  my @sq_lines;
  if (exists $headers->{SQ}) {
    @sq_lines = @{$headers->{SQ}};
  } else {
    confess sprintf $FMT_ERR, "no SQ line found in $data_path: $_. $!";
  }
  for my $sq (@sq_lines) {
    if (exists $sq->{SN} && exists $sq->{LN} && exists $sq->{AS} && exists $sq->{M5} && exists $sq->{SP}) {
      $valid_sq_hash{$sq->{SN}} = { LN => $sq->{LN},	AS => $sq->{AS},	M5 => $sq->{M5},	SP => $sq->{SP} };
    } else {
      confess sprintf $FMT_ERR, "invalid SQ line found in $data_path: $sq->{SN} : $_. $!";
    }
  }
  return \%valid_sq_hash;
}

sub check_donor_states {
  warn "\n#--- check donors\' normal/tumour/nominated_normal.. \n";
  my ($donor_ids_ref, $is_normals_ref, $is_normal_for_donors_ref, $bam_states_ref) = @_;
  my %donors;
  for my $i (0..scalar @$is_normals_ref - 1) {
    if (! exists $donors{${$donor_ids_ref}[$i]}) {
      $donors{${$donor_ids_ref}[$i]}{'nomal'} = 0;
      $donors{${$donor_ids_ref}[$i]}{'nominated'} = 0;
      $donors{${$donor_ids_ref}[$i]}{'tumour'} = 0;
    }
    if(${$is_normals_ref}[$i] =~ $VALID_IS_NORMALS) {
      if (${$is_normals_ref}[$i] =~ $VALID_IS_NORMALS_IS) {
        $donors{${$donor_ids_ref}[$i]}{'nomal'} += 1;
      } else {
        $donors{${$donor_ids_ref}[$i]}{'tumour'} += 1;
      }
    }
    if(${$is_normal_for_donors_ref}[$i] =~ $VALID_IS_NORMAL_FOR_DONORS && ${$is_normal_for_donors_ref}[$i] =~ $VALID_IS_NORMAL_FOR_DONORS_IS) {
      $donors{${$donor_ids_ref}[$i]}{'nominated'} += 1;
    }
  }
  for my $i (0..scalar @$is_normals_ref - 1) {
    my @bam_state;
    if (
    $donors{${$donor_ids_ref}[$i]}{'nomal'} == 0) {
      warn "found no normal for donor: ${$donor_ids_ref}[$i]\n";
      push @bam_state, $FAIL_CODES{'015'};
    }
    if (
    $donors{${$donor_ids_ref}[$i]}{'tumour'} == 0) {
      warn "found no tumour for donor: ${$donor_ids_ref}[$i]\n";
      push @bam_state, $FAIL_CODES{'016'};
    }
    if (
    $donors{${$donor_ids_ref}[$i]}{'nominated'} == 0) {
      warn "found no nominated_normal for donor: ${$donor_ids_ref}[$i]\n";
      push @bam_state, $FAIL_CODES{'012'};
    }
    if (
    $donors{${$donor_ids_ref}[$i]}{'nominated'} > 1) {
      warn "found more than 1 nominated_normal for donor: ${$donor_ids_ref}[$i]\n";
      push @bam_state, $FAIL_CODES{'014'};
    }
    if(${$is_normals_ref}[$i] !~ $VALID_IS_NORMALS_IS && ${$is_normal_for_donors_ref}[$i] =~ $VALID_IS_NORMAL_FOR_DONORS_IS) {
      warn "nominated_normal sample is not normal for donor: ${$donor_ids_ref}[$i]\n";
      push @bam_state, $FAIL_CODES{'013'};
    }
    if (scalar @bam_state > 0) {
      if (${$bam_states_ref}[$i] eq 'raw' || ${$bam_states_ref}[$i] eq 'pp-remapped') {
        ${$bam_states_ref}[$i] = 'failed:'.join '; ', @bam_state;
      } else {
        ${$bam_states_ref}[$i] .= ';'.join '; ', @bam_state;
      }
    }
  }
  return $bam_states_ref;
}

sub bam_header_to_hash { # note @CO line is not necessarily tab delimited.
  my $align_file = shift;
  my %header_sets;

  my $sam = Bio::DB::HTS->new(-bam => $align_file);
  my @lines = split "\n", $sam->header->text;

  for my $line (@lines) {
    chomp $line;
    my ($type, @bits) = split /\t/, $line;
    $type =~ s/^\@//;
    my %fields;
    for my $bit (@bits) {
      my ($ft, $val) = $bit =~ m/([^:]+):(.+)/;
      $fields{$ft}=$val;
      # warn "adding $type: $ft - $val\n";
    }
    push @{$header_sets{$type}}, \%fields;
  }
  return \%header_sets;
}

sub bam_has_reads_more_than_threshold {
  warn "Counting reads...\n";
  my ($align_file, $threshold) = @_;
  my $sam = Bio::DB::HTS->new(-bam => $align_file);
  my $bam = $sam->hts_file;
  my $header = $bam->header_read;
  my $processed_x = 0;
  my $start = time;
  my $processed_x_mill = 0;
  my $enough = 0;

  while (my $a = $bam->read1($header)) {
    $processed_x += 1;
    if ($processed_x >= $threshold) {
      $enough = 1;
      last;
    }
  }

  my $end = time;
  my $elapsed = $end - $start;
  if ($enough) {
    warn "bam has more than $threshold reads [time elapsed: ${elapsed}s].\n";
  } else {
    warn "bam has only $processed_x reads, which is fewer than required ($threshold). [time elapsed: ${elapsed}s].\n";
  }
  return $enough, $processed_x;
}

sub write_results {
  my ($out_array, $output, $format, $validate_status) = @_;
  if ($format eq 'xlsx') {
    # rename output, if it does not ends with .xlsx
    my $workbook = Excel::Writer::XLSX->new($output);
    my $worksheet = $workbook->add_worksheet();
    my $row = 0;
    foreach my $line (@$out_array) {
      my @array = split "\t", $line;
      $worksheet->write_row($row, 0, \@array);
      $row++;
    }
  } elsif ($format eq 'xls') {
    # rename output, if it does not ends with .xls
    my $workbook = Spreadsheet::WriteExcel->new($output);
    my $worksheet = $workbook->add_worksheet();
    my $row = 0;
    foreach my $line (@$out_array) {
      my @array = split "\t", $line;
      $worksheet->write_row($row, 0, \@array);
      $row++;
    }
  } else {
    write_to_file($out_array, $output);
  }

  # write an extra output file on sucessful validation if input format is xlsx or xls
  if ($validate_status && $format =~ m/^xls$|^xlsx$/i) {
    $output =~ s/\.xls$|\.xlsx$/\.tsv/i;
    write_to_file($out_array, $output);
  }
}

sub write_to_file {
  my ($array_ref, $output) = @_;
  open my $OUT, ">", $output or confess sprintf $FMT_ERR, "can not open the output: $output to write. $!";
  foreach my $line (@$array_ref) {
    chomp $line; print $OUT "$line\n";
  }
  close $OUT;
}


__DATA__

--- pp_remapped_bam_header ---
@SQ	SN:1	LN:249250621	AS:GRCh37d5	M5:1b22b98cdeb4a9304cb5d48026a85128	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:2	LN:243199373	AS:GRCh37d5	M5:a0d9851da00400dec1098a9255ac712e	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:3	LN:198022430	AS:GRCh37d5	M5:fdfd811849cc2fadebc929bb925902e5	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:4	LN:191154276	AS:GRCh37d5	M5:23dccd106897542ad87d2765d28a19a1	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:5	LN:180915260	AS:GRCh37d5	M5:0740173db9ffd264d728f32784845cd7	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:6	LN:171115067	AS:GRCh37d5	M5:1d3a93a248d92a729ee764823acbbc6b	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:7	LN:159138663	AS:GRCh37d5	M5:618366e953d6aaad97dbe4777c29375e	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:8	LN:146364022	AS:GRCh37d5	M5:96f514a9929e410c6651697bded59aec	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:9	LN:141213431	AS:GRCh37d5	M5:3e273117f15e0a400f01055d9f393768	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:10	LN:135534747	AS:GRCh37d5	M5:988c28e000e84c26d552359af1ea2e1d	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:11	LN:135006516	AS:GRCh37d5	M5:98c59049a2df285c76ffb1c6db8f8b96	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:12	LN:133851895	AS:GRCh37d5	M5:51851ac0e1a115847ad36449b0015864	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:13	LN:115169878	AS:GRCh37d5	M5:283f8d7892baa81b510a015719ca7b0b	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:14	LN:107349540	AS:GRCh37d5	M5:98f3cae32b2a2e9524bc19813927542e	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:15	LN:102531392	AS:GRCh37d5	M5:e5645a794a8238215b2cd77acb95a078	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:16	LN:90354753	AS:GRCh37d5	M5:fc9b1a7b42b97a864f56b348b06095e6	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:17	LN:81195210	AS:GRCh37d5	M5:351f64d4f4f9ddd45b35336ad97aa6de	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:18	LN:78077248	AS:GRCh37d5	M5:b15d4b2d29dde9d3e4f93d1d0f2cbc9c	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:19	LN:59128983	AS:GRCh37d5	M5:1aacd71f30db8e561810913e0b72636d	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:20	LN:63025520	AS:GRCh37d5	M5:0dec9660ec1efaaf33281c0d5ea2560f	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:21	LN:48129895	AS:GRCh37d5	M5:2979a6085bfe28e3ad6f552f361ed74d	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:22	LN:51304566	AS:GRCh37d5	M5:a718acaa6135fdca8357d5bfe94211dd	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:X	LN:155270560	AS:GRCh37d5	M5:7e0e2e580297b7764e31dbc80c2540dd	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:Y	LN:59373566	AS:GRCh37d5	M5:1fa3474750af0948bdf97d5a0ee52e51	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:MT	LN:16569	AS:GRCh37d5	M5:c68f52674c9fb33aef52dcf399755519	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000207.1	LN:4262	AS:GRCh37d5	M5:f3814841f1939d3ca19072d9e89f3fd7	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000226.1	LN:15008	AS:GRCh37d5	M5:1c1b2cd1fccbc0a99b6a447fa24d1504	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000229.1	LN:19913	AS:GRCh37d5	M5:d0f40ec87de311d8e715b52e4c7062e1	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000231.1	LN:27386	AS:GRCh37d5	M5:ba8882ce3a1efa2080e5d29b956568a4	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000210.1	LN:27682	AS:GRCh37d5	M5:851106a74238044126131ce2a8e5847c	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000239.1	LN:33824	AS:GRCh37d5	M5:99795f15702caec4fa1c4e15f8a29c07	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000235.1	LN:34474	AS:GRCh37d5	M5:118a25ca210cfbcdfb6c2ebb249f9680	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000201.1	LN:36148	AS:GRCh37d5	M5:dfb7e7ec60ffdcb85cb359ea28454ee9	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000247.1	LN:36422	AS:GRCh37d5	M5:7de00226bb7df1c57276ca6baabafd15	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000245.1	LN:36651	AS:GRCh37d5	M5:89bc61960f37d94abf0df2d481ada0ec	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000197.1	LN:37175	AS:GRCh37d5	M5:6f5efdd36643a9b8c8ccad6f2f1edc7b	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000203.1	LN:37498	AS:GRCh37d5	M5:96358c325fe0e70bee73436e8bb14dbd	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000246.1	LN:38154	AS:GRCh37d5	M5:e4afcd31912af9d9c2546acf1cb23af2	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000249.1	LN:38502	AS:GRCh37d5	M5:1d78abec37c15fe29a275eb08d5af236	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000196.1	LN:38914	AS:GRCh37d5	M5:d92206d1bb4c3b4019c43c0875c06dc0	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000248.1	LN:39786	AS:GRCh37d5	M5:5a8e43bec9be36c7b49c84d585107776	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000244.1	LN:39929	AS:GRCh37d5	M5:0996b4475f353ca98bacb756ac479140	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000238.1	LN:39939	AS:GRCh37d5	M5:131b1efc3270cc838686b54e7c34b17b	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000202.1	LN:40103	AS:GRCh37d5	M5:06cbf126247d89664a4faebad130fe9c	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000234.1	LN:40531	AS:GRCh37d5	M5:93f998536b61a56fd0ff47322a911d4b	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000232.1	LN:40652	AS:GRCh37d5	M5:3e06b6741061ad93a8587531307057d8	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000206.1	LN:41001	AS:GRCh37d5	M5:43f69e423533e948bfae5ce1d45bd3f1	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000240.1	LN:41933	AS:GRCh37d5	M5:445a86173da9f237d7bcf41c6cb8cc62	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000236.1	LN:41934	AS:GRCh37d5	M5:fdcd739913efa1fdc64b6c0cd7016779	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000241.1	LN:42152	AS:GRCh37d5	M5:ef4258cdc5a45c206cea8fc3e1d858cf	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000243.1	LN:43341	AS:GRCh37d5	M5:cc34279a7e353136741c9fce79bc4396	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000242.1	LN:43523	AS:GRCh37d5	M5:2f8694fc47576bc81b5fe9e7de0ba49e	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000230.1	LN:43691	AS:GRCh37d5	M5:b4eb71ee878d3706246b7c1dbef69299	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000237.1	LN:45867	AS:GRCh37d5	M5:e0c82e7751df73f4f6d0ed30cdc853c0	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000233.1	LN:45941	AS:GRCh37d5	M5:7fed60298a8d62ff808b74b6ce820001	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000204.1	LN:81310	AS:GRCh37d5	M5:efc49c871536fa8d79cb0a06fa739722	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000198.1	LN:90085	AS:GRCh37d5	M5:868e7784040da90d900d2d1b667a1383	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000208.1	LN:92689	AS:GRCh37d5	M5:aa81be49bf3fe63a79bdc6a6f279abf6	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000191.1	LN:106433	AS:GRCh37d5	M5:d75b436f50a8214ee9c2a51d30b2c2cc	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000227.1	LN:128374	AS:GRCh37d5	M5:a4aead23f8053f2655e468bcc6ecdceb	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000228.1	LN:129120	AS:GRCh37d5	M5:c5a17c97e2c1a0b6a9cc5a6b064b714f	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000214.1	LN:137718	AS:GRCh37d5	M5:46c2032c37f2ed899eb41c0473319a69	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000221.1	LN:155397	AS:GRCh37d5	M5:3238fb74ea87ae857f9c7508d315babb	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000209.1	LN:159169	AS:GRCh37d5	M5:f40598e2a5a6b26e84a3775e0d1e2c81	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000218.1	LN:161147	AS:GRCh37d5	M5:1d708b54644c26c7e01c2dad5426d38c	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000220.1	LN:161802	AS:GRCh37d5	M5:fc35de963c57bf7648429e6454f1c9db	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000213.1	LN:164239	AS:GRCh37d5	M5:9d424fdcc98866650b58f004080a992a	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000211.1	LN:166566	AS:GRCh37d5	M5:7daaa45c66b288847b9b32b964e623d3	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000199.1	LN:169874	AS:GRCh37d5	M5:569af3b73522fab4b40995ae4944e78e	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000217.1	LN:172149	AS:GRCh37d5	M5:6d243e18dea1945fb7f2517615b8f52e	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000216.1	LN:172294	AS:GRCh37d5	M5:642a232d91c486ac339263820aef7fe0	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000215.1	LN:172545	AS:GRCh37d5	M5:5eb3b418480ae67a997957c909375a73	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000205.1	LN:174588	AS:GRCh37d5	M5:d22441398d99caf673e9afb9a1908ec5	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000219.1	LN:179198	AS:GRCh37d5	M5:f977edd13bac459cb2ed4a5457dba1b3	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000224.1	LN:179693	AS:GRCh37d5	M5:d5b2fc04f6b41b212a4198a07f450e20	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000223.1	LN:180455	AS:GRCh37d5	M5:399dfa03bf32022ab52a846f7ca35b30	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000195.1	LN:182896	AS:GRCh37d5	M5:5d9ec007868d517e73543b005ba48535	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000212.1	LN:186858	AS:GRCh37d5	M5:563531689f3dbd691331fd6c5730a88b	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000222.1	LN:186861	AS:GRCh37d5	M5:6fe9abac455169f50470f5a6b01d0f59	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000200.1	LN:187035	AS:GRCh37d5	M5:75e4c8d17cd4addf3917d1703cacaf25	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000193.1	LN:189789	AS:GRCh37d5	M5:dbb6e8ece0b5de29da56601613007c2a	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000194.1	LN:191469	AS:GRCh37d5	M5:6ac8f815bf8e845bb3031b73f812c012	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000225.1	LN:211173	AS:GRCh37d5	M5:63945c3e6962f28ffd469719a747e73c	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:GL000192.1	LN:547496	AS:GRCh37d5	M5:325ba9e808f669dfeee210fdd7b470ac	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:NC_007605	LN:171823	AS:GRCh37d5	M5:6743bd63b3ff2b5b8985d8933c53290a	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
@SQ	SN:hs37d5	LN:35477943	AS:GRCh37d5	M5:5b6a4b3a81a2d3c134b7d14bf6ad39f1	SP:human	UR:file:///lustre/scratch112/sanger/kr2/PanCancerFinal/map_ref/genome.fa
