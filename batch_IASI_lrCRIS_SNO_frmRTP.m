function batch_IASI_lrCRIS_SNO_frmRTP()
cd /home/chepplew/projects/sno/makeSNO/batchJobs

% Get the driver file of job dates - should be one calender month
 fh = fopen('jobDates.drv','r');
 fdates = textscan(fh,'%s');                 % cell array
 fclose(fh);
 njobs = numel(fdates{1});
 disp(['number of jobs: ' num2str(njobs)]);
 %%proc_airs_to_cris_mat(njob);

  %end

%%njobs = 2; clear sJobs;
%%for i=1:njobs sJobs{i} = fLst(i).name; end

fnbatch = './batch_IASI_lrCRIS_SNO_frmRTP.slurm';
FH = fopen(fnbatch,'w');
fprintf(FH,'#!/bin/bash\n\n');
fprintf(FH,'#SBATCH --job-name=ICsno\n');
  junk = '#SBATCH --output=ICsno-slurm-%N.%A.%a.out';
fprintf(FH,'%s\n',junk);
  junk = '#SBATCH --error=ICsno-slurm-%N.%A.%a.err';
fprintf(FH,'%s\n',junk);
fprintf(FH,'#SBATCH --partition=batch\n');
fprintf(FH,'#SBATCH --qos=short\n');
fprintf(FH,'#SBATCH --account=pi_strow\n');
%%fprintf(FH,'#SBATCH -N1\n');
fprintf(FH,'#SBATCH --mem-per-cpu=10000\n');
fprintf(FH,'#SBATCH --cpus-per-task 1\n');
fprintf(FH,'#SBATCH --array=1-%d\n\n',njobs);
fprintf(FH,'MATLAB=''/usr/cluster/matlab/2016a/bin/matlab''\n');
fprintf(FH,'MATOPTS='' -nodisplay -nojvm -nosplash''\n\n');
  junk = sprintf('srun ./batch_IASI_lrCRIS_SNO_frmRTP.sh');
fprintf(FH,'%s\n',junk);

fclose(FH);

[stat, resn] = system( ['chmod 755 ' fnbatch] );
disp([num2str(stat) '   ' resn]);

pause(3);

[stat, resn] = system( ['sbatch ' fnbatch] );
disp([num2str(stat) '   ' resn]);
