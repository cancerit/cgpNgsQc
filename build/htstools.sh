#!/bin/bash

VER_HTSLIB="1.7"
VER_SAMTOOLS="1.7"
VER_BIODBHTS="2.9"

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

set -e

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH
cd $INST_PATH
INST_PATH=`pwd`
mkdir -p $INST_PATH/bin

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"
export PATH="$INST_PATH/bin:$PATH"
export LD_LIBRARY_PATH="$INST_PATH/lib:$PATH"

cd $INIT_DIR

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR
cd $SETUP_DIR
## HTSLIB (tar.bz2)
if [ ! -e $SETUP_DIR/htslib.success ]; then
  rm -rf htslib
  mkdir -p htslib
  curl -sSL --retry 10 https://github.com/samtools/htslib/releases/download/${VER_HTSLIB}/htslib-${VER_HTSLIB}.tar.bz2 > distro.tar.bz2
  tar --strip-components 1 -C htslib -jxf distro.tar.bz2
  cd htslib
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make clean
  make -j2
  make install
  cd $SETUP_DIR
  rm -rf distro.*
  touch $SETUP_DIR/htslib.success
fi

mkdir -p distro # reused
## SAMTOOLS (tar.bz2)
if [ ! -e $SETUP_DIR/samtools.success ]; then
  curl -sSL --retry 10 https://github.com/samtools/samtools/releases/download/${VER_SAMTOOLS}/samtools-${VER_SAMTOOLS}.tar.bz2 > distro.tar.bz2
  tar --strip-components 1 -C distro -xjf distro.tar.bz2
  cd distro
  ./configure --enable-plugins --enable-libcurl --with-htslib=$INST_PATH --prefix=$INST_PATH
  make clean
  make -j2 all
  make install
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/samtools.success
fi

## Bio::DB::HTS (tar.gz)
if [ ! -e $SETUP_DIR/Bio-DB-HTS.success ]; then
  ## add perl deps
  cpanm --no-wget --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH Module::Build
  cpanm --no-wget --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH Bio::Root::Version

  curl -sSL --retry 10 https://github.com/Ensembl/Bio-DB-HTS/archive/${VER_BIODBHTS}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -zxf distro.tar.gz
  cd distro
  perl Build.PL --install_base=$INST_PATH --htslib=$INST_PATH
  ./Build
  ./Build test
  ./Build install
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/Bio-DB-HTS.success
fi

# cleanup all junk
cd $INIT_DIR
rm -rf $SETUP_DIR

exit 0
