#!/usr/bin/env bash

set -euo pipefail

grep -F -f ./tps_ids.txt ../02-orthofinder/orthofinder_longest_isoform/Results_Oct15/Orthologues/Orthologues_T_mongolica/T_mongolica__v__Z_stapffii.tsv |
cut -f1 |
sort -u > ogs.txt
