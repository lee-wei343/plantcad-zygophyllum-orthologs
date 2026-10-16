#!/usr/bin/env bash

# Results go to ../orthofinder_longest_isoform/Results_<date>/; OrthoFinder refuses to start if -o already exists.
orthofinder -f ../data/proteoms/primary_transcripts/ -s ../data/rooted_tree.nwk -t 8 -o ../orthofinder_longest_isoform
