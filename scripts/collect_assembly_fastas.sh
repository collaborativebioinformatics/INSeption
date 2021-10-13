#!/bin/bash

sourcedir=$1
targetdir=$2
# Script to extract all .fasta clusters and collect them in a central place
mkdir -p $targetdir

find ./${sourcedir} -name "contigs.fasta*" \
| while read p
  do
  Basename=$(basename ${p%.fasta})
  Parentdir=$(basename $(dirname "$p"))
  cp $p ./${targetdir}/${Basename}_${Parentdir}.fasta
  COUNTER=$((COUNTER + 1))
  done

