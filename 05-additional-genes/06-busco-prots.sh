#!/usr/bin/env bash

mkdir -p out

busco -i ../data/Z_stapfii/Z_stapfii.pep.all.fa \
  -o ./out/ \
  -m proteins \
  -l embryophyta_odb10 \
  -c 8
