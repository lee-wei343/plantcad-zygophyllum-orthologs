#!/usr/bin/env bash

mkdir -p out

busco -i ../02-orthofinder/data/proteoms/primary_transcripts/Z_stapffii.fa \
  --out_path ./out \
  --out stapffii_proteom \
  --download_path ../01-busco/busco_downloads/ \
  -m proteins \
  -l embryophyta_odb10 \
  -c 8
