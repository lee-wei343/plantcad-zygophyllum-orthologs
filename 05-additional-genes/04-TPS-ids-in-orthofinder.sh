#!/usr/bin/env bash

set -euo pipefail

grep -f tps_ids.txt ../02-orthofinder/Results_Sep29/Orthologues/Orthologues_T_mongolica/T_mongolica__v__Z_stapffii.tsv
