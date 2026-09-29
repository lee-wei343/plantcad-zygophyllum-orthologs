#!/usr/bin/env bash

# ../../data/T_mongolica/T_mongolica_genes.gff3.pep == ../data/proteoms/T_mongolica.fa

sed -i '/^>/! s/[.*]$//' ../data/proteoms/T_mongolica.fa
