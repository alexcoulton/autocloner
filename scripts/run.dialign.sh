#!/bin/bash
JOBPATH=$1
export DIALIGN2_DIR='INSERT_DIALIGN2_PATH_HERE'

dialign2-2 -n -anc -fa -fn $JOBPATH/seq/extended/alignments/all.align.rev $JOBPATH/seq/extended/seqs/all.fa
