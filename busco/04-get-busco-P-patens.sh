#!/usr/bin/env bash

mkdir -p out/busco

busco -i ../data/P_patens/Physcomitrium_patens.Phypa_V3.dna.toplevel.fa \
  -o ../out/busco/busco_P_patens \
  -m genome \
  -l embryophyta_odb10 \
  -c 8
