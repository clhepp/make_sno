function run_makeSNO_batch()
%
% This script is called by run_makeSNO_batch.slurm
% create simple list of dates and save as txt file, (ln 13)
% call the appropriate SNO make script for desired sensor pair (ln 17)
%

cd /home/chepplew/projects/sno/makeSNO/batchJobs
addpath /home/chepplew/projects/sno/makeSNO
addpath /home/chepplew/projects/sno/makeSNO/src

%slurmindex  = str2num(getenv('SLURM_ARRAY_TASK_ID'));
slurmindex = 1;

%[st, instr] = system(sprintf('sed -n "%dp" ./sno_date_list.txt | tr -d "\n"', slurmindex));
 fh = fopen('./jobDates.drv','r');
 fdates = textscan(fh,'%s');                 % cell array
 fclose(fh);
 instr = cell2mat(fdates{1}(slurmindex));
 disp(instr)
%%%%make_IASI_CRIS_sno(instr,'low')
make_AIRS_CRIS_sno(instr,'low');
%%%%%makeIasiAirsSno(1);

end
