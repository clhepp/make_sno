function [] = make_IASI_lrCRIS_SNO_frmRTP(req_date,model)
%
% function make_IASI_lrCRIS_SNO_frmRTP.m
%
% INPUTS:
%    re_date: date string 'YYYY/MM/DD'
%    model:   'ERA' or 'ECMWF' to select the model used for the sarta calcs in the RTP subsets.
%
% OUTPUTS: None
%
% RESULTS: A .mat data file of accummulated SNO data, consisting:
%        - structure: 'pre' with geolocation fields for each sensor:
%        (time, lat, lon, ifov, satzen, solzen, landfrac, time-delay, great-circle angle).
%        - arrays: CrIS Obs, CrIS Calcs, IASI Obs and IASI Calcs.
%        - and the RTP subset file names, their head structures, and separation criteria.
%
% Notes: After checking the existance of the daily clear subset RTP files, loads the
%  structures and selects centre-track FORs then compares each CrIS and IASI geolocation field 
%  to find those that pass the separation criteria. Since any one CrIS Obs can find multiple IASI
%  Obs (and vice-versa), only unique pairs are extracted.
%  The separation criteria are hard-wired in variables: maxDtim and maxDphi.
%
%  Can be run stand-alone or with batch control script: batch_IASI_lrCRIS_SNO_frmRTP.m and
%      shell script: batch_IASI_lrCRIS_SNO_frmRTP.sh  

cd /home/chepplew/projects/sno/makeSNO

addpath /asl/matlib/h4tools                % rtpread.m
addpath /asl/matlib/rtptools               % rtpread_12.m

prcnam = mfilename('fullpath');

% Hardwire separation criteria:
maxDtim = 1200.0;        % seconds
maxDphi = 0.18;          % deg. 0.07 = 7.8 km Earth radius: 6371 km. 1 deg = 111 km.

% Check date string
whos req_date; disp(req_date); fprintf('\n');
try 
   D = datenum(req_date,'YYYY/MM/DD');
catch
   error('Incorrect Date Format')
end

% Check model type
all_models = {'ERA','ECMWF'};
model      = upper(model);
if(~ismember(model, all_models))
  error('Invalid choice of model')
  return
end

% Get Date to process
% e.g. cyear = '2012';  cmon  = '05';  cday  = '01';
cyear = req_date(1:4); cmon = req_date(6:7); cday = req_date(9:10);

% Get loRes CCAST CrIS clear subset files
cpath  = '/asl/rtp/rtp_cris_ccast_lowres/clear/';
if(strcmp(model,'ERA'))
  cfnam  = [cpath cyear '/cris_lr_era_d' cyear cmon cday '_clear.rtp'];
elseif(strcmp(model,'ECMWF'))
  cfnam  = [cpath cyear '/cris_lr_ecmwf_d' cyear cmon cday '_isarta_clear.rtp'];  % uses iasi sarta
end

% Check file exists
[fid, errMsg] = fopen(cfnam);
if(isempty(errMsg))   %if(exist(cfnam))
  cfild = dir(cfnam);
else
  disp([cfnam ' ' errMsg]);
  error('Error, CrIS file does not exist');
  return
end
fclose(fid); clear fid errMsg;

% Get IASI clear subset files (NB .rtp_1 and .rtp_2 files)
dpath  = '/asl/rtp/rtp_iasi1/clear/';
if(strcmp(model,'ERA'))
  dfnam = [dpath cyear '/iasi1_era_d' cyear cmon cday '_clear.rtp'];
elseif(strcmp(model,'ECMWF'))
  dfnam  = [dpath cyear '/iasi1_ecmwf_d' cyear cmon cday '_clear.rtp'];
end

%Check file exists
[fid, errMsg] = fopen([dfnam '_1']);
if(isempty(errMsg))
  dfild = dir([dfnam '_*']);
else
  disp([dfnam ' ' errMsg]);
  error('Error, IASI file does not exist');
  return
end
fclose(fid); clear fid errMsg;

% Load the two files
fprintf(1,'loading  %s \n  and %s \n',cfnam,dfnam)
[chead chatt cprof cpatt] = rtpread(cfnam);
[dhead dhatt dprof dpatt] = rtpread_12(dfnam);

cnobs = size(cprof.robs1,2);
dnobs = size(dprof.robs1,2);
disp(['Total number obs in original subset files: ' num2str(cnobs) ' and ' num2str(dnobs)]);

%{
% review geographical distro.
  addpath /asl/matlib/aslutil
  addpath /asl/matlib/time
  figure(1);clf;simplemap(cprof.rlat,cprof.rlon,ones(1,cnobs));title('CrIS clear');
  figure(2);clf;simplemap(dprof.rlat,dprof.rlon,ones(1,dnobs));title('IASI clear');
% review temporal distro.
  figure(3);clf;plot([1:cnobs],utc2dnum(cprof.rtime),'.',[1:dnobs],utc2dnum(dprof.rtime),'.');  
%}
% -----------------------------------------------
% Options for subsetting - default is center FORs
% ------------------------------------------------
% Subset to centre track
ccntr   = find(cprof.xtrack == 14 | cprof.xtrack == 15 | cprof.xtrack == 16 | cprof.xtrack == 17);
dcntr   = find(dprof.xtrack == 14 | dprof.xtrack == 15 | dprof.xtrack == 16 | dprof.xtrack == 17);
cnumcn  = numel(ccntr);
dnumcn  = numel(dcntr);
disp(['There are ' num2str(cnumcn) ' CrIS center ' num2str(dnumcn) ' IASI center FORs'])

% Subset to tropical ocean
cindto = find(abs(cprof.rlat) <= 40 & cprof.landfrac == 0);
dindto = find(abs(dprof.rlat) <= 40 & dprof.landfrac == 0);
cnumto = numel(cindto);
dnumto = numel(dindto);
disp(['There are ' num2str(cnumto) ' CrIS T.O. and ' num2str(dnumto) ' IASI T.O.'])

% Hardwire which of these subsets to use (so that the sorting routine does not need editing):
icsub = ccntr;        % or cindto
idsub = dcntr;        % or dindto
ncsub = numel(icsub);
ndsub = numel(idsub);

% ------------------------------------------------------------------
% Find nearest obs & locations
% ------------------------------------------------------------------
% pre-compute position vectors - saves a heap of time
fprintf('computing position vectors\n');
P1 = [;;];                              % CRIS
P2 = [;;];                              % IASI
for ic = 1:ncsub                        % Number of centres.
  ix = icsub(ic);                       % ix = ic
  P1 = [P1; cos(cprof.rlat(ix)*pi/180.0) * cos(cprof.rlon(ix)*pi/180.0), ...
        cos(cprof.rlat(ix)*pi/180.0)*sin(cprof.rlon(ix)*pi/180.0), sin(cprof.rlat(ix)*pi/180.0)];
end
for ii = 1:ndsub                        % ii = no. centers. 
  iy = idsub(ii);                       % iy = ii
  P2 = [P2; cos(dprof.rlat(iy)*pi/180.0) * cos(dprof.rlon(iy)*pi/180.0), ...
        cos(dprof.rlat(iy)*pi/180.0)*sin(dprof.rlon(iy)*pi/180.0), sin(dprof.rlat(iy)*pi/180.0)];
end    

% -------------------------------------------------------
fprintf('Computing separations\n');
dist = 0.0; m = 0; k = 1;
pos  = [;]; tdiff = [];
tic
for i = 1:1:ncsub
  ix = icsub(i); 
  for j = 1:1:ndsub 
    iy = idsub(j); 
    Dtim = abs( dprof.rtime(iy) - cprof.rtime(ix) );
    if Dtim <= maxDtim       
      m = m+1;
      Phi = real(acos( sum(P1(i,:).*P2(j,:)) )*180/pi);
      if Phi <= maxDphi                 % phi = 0.18 deg (0.0031 rad) => 20 km.
	dist(k) = Phi;                  % dist;
	pos     = [pos;[i,j]];          % pos(i: CrIS index. j: IASI index)
	tdiff   = [tdiff, Dtim];
	k = k+1;
      end
    end
  end
  if(mod(i,1000) == 0 ) fprintf('.'); end
end
toc
dist = real(dist);                         % can get complex phi.
fprintf('\n');
fprintf(1,'Number of close pairs: %4i\n', size(pos,1));
%{
 % checks
 figure(3);clf;plot(dist,'.');title('IASI CrIS clr TSNO test dist (deg)');  
 figure(3);clf;plot(tdiff,'.');title('IASI CrIS TSNO test yield TDiff (s)');
 figure(4);clf;plot(dprof.rtime(idsub(pos(:,2))) - cprof.rtime(icsub(pos(:,1))),'o' );
%}
% ---------------------------------------------------------------------------
% This gives us duplicate hits - so need to select only unique pairs.
%  this is done by sequentially testing for uniquness on each sensor.
upos = [];
if(size(pos,1) > 4 )
  un1 = [;]; un2 = [;]; upos = [;];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];             % unique CrIS
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
  fprintf('Number unique SNOs %d %d\n',size(upos));
end
if(size(upos,1) > 2)
% Record time and space separations (pre reloading of files)
  pre = struct;
  pre.iSnoTim = dprof.rtime(idsub(upos(:,2)));    pre.cSnoTim = cprof.rtime(icsub(upos(:,1)));
  pre.iSnoLat = dprof.rlat(idsub(upos(:,2)));     pre.cSnoLat = cprof.rlat(icsub(upos(:,1)));
  pre.iSnoLon = dprof.rlon(idsub(upos(:,2)));     pre.cSnoLon = cprof.rlon(icsub(upos(:,1)));
  pre.utdiff  = pre.iSnoTim - pre.cSnoTim;    
  pre.ilnfr   = dprof.landfrac(idsub(upos(:,2))); pre.clnfr   = cprof.landfrac(icsub(upos(:,1)));
  pre.isolz   = dprof.solzen(idsub(upos(:,2)));   pre.csolz   = cprof.solzen(icsub(upos(:,1)));
  pre.isatz   = dprof.satzen(idsub(upos(:,2)));   pre.csatz   = cprof.satzen(icsub(upos(:,1)));
  pre.ifov    = dprof.ifov(idsub(upos(:,2)));     pre.cfov    = cprof.ifov(icsub(upos(:,1)));
  pre.upos    = upos;
  uP1    = P1(upos(:,1),:);      uP2 = P2(upos(:,2),:);    % P1: CRIS. P2: IASI
  for i = 1:size(upos,1)
     uPhi(i)   = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi );
  end
  pre.uPhi = uPhi;

  crobs  = cprof.robs1(:,icsub(upos(:,1)));
  iaobs  = dprof.robs1(:,idsub(upos(:,2)));
  crcalc = cprof.rcalc(:,icsub(upos(:,1)));
  iacalc = dprof.rcalc(:,idsub(upos(:,2)));
%{
figure(3);clf;simplemap(pre.iSnoLat,pre.iSnoLon,ones(1,26));hold on;simplemap(pre.cSnoLat,pre.cSnoLon,2*ones(1,26));
 title('2012.04.05 TSNO IASI CrIS 3.3hrs 100km 1342 SNOs')
 %saveas(gcf,'./figs/20120405_iasi_cris_lr_clear_tsno.png','png');
%}

  % Save the data
  savpath = ['/asl/s1/chepplew/data/sno/iasi_cris/LR/' cyear];
  if ~exist(savpath) disp(['mkdir ' savpath]); mkdir(savpath); end;
  if(strcmp(model,'ECMWF'))
    savfn   = [cyear cmon cday '_iasi_lrCris_clear_SNO_frm_isarta_ecmwf_RTP.mat'];
  elseif(strcmp(model,'ERA'))
    savfn   = [cyear cmon cday '_iasi_lrCris_clear_SNO_frm_era_RTP.mat'];
  end
  disp(['Saving data to: ' savfn]);
  savvars = {'pre','crobs','crcalc','iaobs','iacalc','chead','dhead',...
  'prcnam','cfnam','dfnam','maxDtim','maxDphi'};
  save([savpath '/' savfn],savvars{:});
end
