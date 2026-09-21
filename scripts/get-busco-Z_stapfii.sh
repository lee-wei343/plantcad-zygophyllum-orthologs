#!/usr/bin/env bash

mkdir -p ../out/busco

busco -i ../data/Z_stapfii/Z_stapfii_genome.fa \
  -o ../out/busco/busco_Z_stapfii \
  -m genome \
  -l embryophyta_odb10 \
  -c 8
