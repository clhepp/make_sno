#!/bin/bash

sindx=$SLURM_ARRAY_TASK_ID
#sindx=1

thisJob=$(sed -n ''${sindx}' p' jobDates.drv)

echo 'this job: ' ${thisJob}

/usr/cluster/matlab/2016a/bin/matlab -nosplash -nodesktop -nodisplay \
 -r "addpath('/home/chepplew/projects/sno/makeSNO'); make_IASI_lrCRIS_tsno_frm_rtp('${thisJob} '); exit;"

