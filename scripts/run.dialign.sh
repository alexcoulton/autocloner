#!/bin/bash
JOBPATH=$1
export DIALIGN2_DIR='/home/ac14037/project.phd.main/autoclonernew/primer/pipeline/dialign2_dir'

dialign2-2 -n -anc -fa -fn $JOBPATH/seq/extended/alignments/all.align.rev $JOBPATH/seq/extended/seqs/all.fa
