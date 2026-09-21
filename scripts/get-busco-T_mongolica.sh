#!/usr/bin/env bash

mkdir -p ../out/busco

busco -i ../data/T_mongolica/T_mongolica_genome.fa \
  -o ../out/busco/busco_T_mongolica \
  -m genome \
  -l embryophyta_odb10 \
  -c 8
