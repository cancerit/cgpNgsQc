#!/bin/bash

SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2"
SOURCE_SAMTOOLS="https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2"
SOURCE_BIOBDHTS="https://github.com/Ensembl/Bio-HTS/archive/2.9.tar.gz"

get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    echo "ERROR: zip archives are not supported by default, if pulling from github replace .zip with .tar.gz"
    exit 1
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  rm -f $1.$EXT
  if hash curl 2>/dev/null; then
    curl --retry 10 -sS -o $1.$EXT -L $2
  else
    echo "ERROR: curl not found"
    exit 1
  fi
}

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

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$CHK" != "x" ]] ; then
  echo "Bio::DB::HTS already installed: $CHK"
  exit 0
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

cd $INIT_DIR

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

echo -n "Get htslib ..."
if [ -e $SETUP_DIR/htslibGet.success ]; then
  echo " already staged ...";
else
  echo
  cd $SETUP_DIR
  get_distro "htslib" $SOURCE_HTSLIB
  touch $SETUP_DIR/htslibGet.success
fi

cd $INIT_DIR
echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo " previously installed ...";
else
  echo
  cd $SETUP_DIR
  mkdir -p htslib
  tar --strip-components 1 -C htslib -jxf htslib.tar.bz2
  cd htslib
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU
  make install
  cd $SETUP_DIR
  touch $SETUP_DIR/htslib.success
fi

export HTSLIB=$INST_PATH

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$CHK" == "x" ]] ; then
  echo -n "Building Bio::DB::HTS ..."
  if [ -e $SETUP_DIR/biohts.success ]; then
    echo " previously installed ...";
  else
    echo
    cd $SETUP_DIR
    rm -rf bioDbHts
    get_distro "bioDbHts" $SOURCE_BIOBDHTS
    mkdir -p bioDbHts
    tar --strip-components 1 -C bioDbHts -zxf bioDbHts.tar.gz
    cd bioDbHts
    perlmods=( "ExtUtils::CBuilder" "Module::Build~0.42" "Bio::Root::Version~1.006924")
    for i in "${perlmods[@]}" ; do
      cpanm --no-wget --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
    done
    perl Build.PL --htslib=$HTSLIB --install_base=$INST_PATH
    ./Build
    ./Build test
    ./Build install
    cd $SETUP_DIR
    rm -f bioDbHts.tar.gz
    touch $SETUP_DIR/biohts.success
  fi
else
  echo "Bio::DB::HTS already installed ..."
fi

SAM=`which samtools`
if [[ "x$SAM" != "x" ]] ; then
  echo "Building samtools ..."
  cd $SETUP_DIR
  rm -rf samtools
  get_distro "samtools" $SOURCE_SAMTOOLS
  mkdir -p samtools
  tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
  cd samtools
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU all all-htslib
  make install all all-htslib
  cd $SETUP_DIR
  rm -f samtools.tar.bz2
fi

# cleanup all junk
cd $INIT_DIR
rm -rf $SETUP_DIR

exit 0
