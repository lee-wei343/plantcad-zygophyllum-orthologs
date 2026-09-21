#!/usr/bin/env bash

mkdir -p ../out/busco

busco -i ../data/A_thaliana/Arabidopsis_thaliana.TAIR10.pep.all.fa \
  -o ../out/busco/busco_A_thaliana \
  -m proteins \
  -l embryophyta_odb10 \
  -c 8
