#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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


use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share dist => 'share';

WriteMakefile(
  NAME          => 'Sanger-CGP-NgsQc',
  LICENSE       => 'agpl_3', # http://search.cpan.org/~dagolden/CPAN-Meta-2.142690/lib/CPAN/Meta/Spec.pm#license
  VERSION_FROM  => 'lib/Sanger/CGP/NgsQc.pm',
  EXE_FILES     => [qw( bin/verifyBamHomChk.pl
                        bin/compareBamGenotypes.pl
                        bin/validate_sample_meta.pl
                         )],
  PREREQ_PM     => {
                     'Const::Fast' => 0.014,
                     'Bio::DB::HTS' => 1.11,
                     'Capture::Tiny' => 0.27,
                     'File::ShareDir' => 1.102,
                     'File::ShareDir::Install' => 0.09,
                     'JSON' => 0.00,
                     'JSON::XS' => 0.00,
                     'Data::UUID' => 0.00,
                     'Digest::MD5' => 0.00,
                     'Spreadsheet::Read' => 0.00,
                     'Spreadsheet::WriteExcel' => 0.00,
                     'Excel::Writer::XLSX' => 0.00,
                   }
);

package MY;
use File::ShareDir::Install qw( postamble );
