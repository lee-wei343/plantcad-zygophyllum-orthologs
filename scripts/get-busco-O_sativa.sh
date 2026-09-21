#!/usr/bin/env bash

mkdir -p ../out/busco

busco -i ../data/O_sativa/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
  -o ../out/busco/busco_O_sativa \
  -m genome \
  -l embryophyta_odb10 \
  -c 8
