#!/bin/bash

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

SOURCE_VERIFYBAM="https://github.com/statgen/verifyBamID/releases/download/v1.1.2/verifyBamID.1.1.2"
SOURCE_SAMTOOLS="https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2"

get_file () {
  if hash curl 2>/dev/null; then
    curl -sSL -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

get_distro () {
  EXT=""
  DECOMP=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
    DECOMP="-j"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
    DECOMP="-z"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi

  get_file $1.$EXT $2
  mkdir -p $1
  tar --strip-components 1 -C $1 $DECOMP -xf $1.$EXT
}

if [[ ($# -ne 1 && $# -ne 2) ]] ; then
  echo "Please provide an installation path and optionally perl lib paths to allow, e.g."
  echo "  ./setup.sh /opt/myBundle"
  echo "OR all elements versioned:"
  echo "  ./setup.sh /opt/cgpVcf-X.X.X /opt/PCAP-X.X.X/lib/perl"
  exit 0
fi

INST_PATH=$1

if [[ $# -eq 2 ]] ; then
  CGP_PERLLIBS=$2
fi


CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5

# allows user to knowingly specify other PERL5LIB areas.
if [ -z ${CGP_PERLLIBS+x} ]; then
  export PERL5LIB="$PERLROOT"
else
  export PERL5LIB="$PERLROOT:$CGP_PERLLIBS"
fi

export PATH=$INST_PATH/bin:$PATH

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

## grab cpanm and stick in workspace, then do a self upgrade into bin:
get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install Bio::DB::HTS before proceeding"
  exit 1;
fi

CHK=`which alleleCounter`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install alleleCount before proceeding: https://github.com/cancerit/alleleCount/releases"
  exit 1;
fi

CHK=`which samtools`
if [[ "x$CHK" == "x" ]] ; then
  COMPILE="samtools"
fi

perlmods=( "File::ShareDir::Install" )

set -e
for i in "${perlmods[@]}" ; do
  echo "Installing build prerequisite $i..."
  $CPANM --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
done

cd $INIT_DIR

echo "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
  echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
$CPANM --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH/ --installdeps . < /dev/null

cd $SETUP_DIR
get_file "verifyBamId" $SOURCE_VERIFYBAM
cp $SETUP_DIR/verifyBamId $INST_PATH/bin/.
chmod +x $INST_PATH/bin/verifyBamId
cd $INIT_DIR

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building samtools ..."
  cd $SETUP_DIR
  rm -rf samtools
  get_distro "samtools" $SOURCE_SAMTOOLS
  cd samtools
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU all all-htslib
  make install all all-htslib
  cd $SETUP_DIR
  rm -f samtools.tar.bz2
else
  echo "samtools exists - will not install"
fi

echo "Installing cgpNgsQc ..."
cd $INIT_DIR
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
