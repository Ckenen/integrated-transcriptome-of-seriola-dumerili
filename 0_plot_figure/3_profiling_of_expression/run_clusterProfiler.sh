#!/bin/bash

set +u
source activate clusterProfiler

/data/chenzonggui/gaotishi/0_Plot_figure/3_profiling_of_expression/clusterProfiler.R $1 $2

conda deactivate