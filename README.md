# cgpNgsQc

Collection of code for checking NSG sequencing results.

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |

Contents:

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Usage](#usage)
- [Dependancies](#dependancies)
- [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
- [LICENCE](#licence)

<!-- /TOC -->

## Usage

Please see the [wiki](https://github.com/cancerit/cgpNgsQc/wiki) for details of usage.

## Dependancies

* [Bio::DB::HTS](http://search.cpan.org/~rishidev/Bio-DB-HTS)
* [Samtools](https://github.com/samtools/samtools/releases)
* [VerifyBamId](http://genome.sph.umich.edu/wiki/VerifyBamID)
* [alleleCount](https://github.com/cancerit/alleleCount)

`build/htstools.sh` will build `Bio::DB::HTS` and `samtools` but it is primarily included for
the travis build so only use this if you don't already have these installed.

## Docker, Singularity and Dockstore

There is a pre-built image containing this codebase on quay.io.

* [dockstore-cgpwgs][ds-cgpwgs-git]
  * Contains additional tools for WGS analysis.

This was primarily designed for use with dockstore.org but can be used as a normal container.

The docker image is know to work correctly after import into a singularity image.

## LICENCE

```
Copyright (c) 2014-2018 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of cgpNgsQc.

cgpNgsQc is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
```

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/cgpNgsQc
[travis-master]: https://travis-ci.org/cancerit/cgpNgsQc.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/cgpNgsQc.svg?branch=develop
