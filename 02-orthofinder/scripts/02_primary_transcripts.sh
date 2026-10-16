#!/usr/bin/env bash

# Z_stapffii headers are BRAKER-style (>g21993.t1, >g21993.t2) without a gene= tag
# Because of this running the normal for non-stapffii
python3 ~/miniforge3/envs/orthofinder/bin/tools/primary_transcript.py ../data/proteoms/

# and with the argument last_dot for stapffii, taking the gene name
# without the .t1
python3 ~/miniforge3/envs/orthofinder/bin/tools/primary_transcript.py ../data/proteoms/Z_stapffii.fa last_dot
