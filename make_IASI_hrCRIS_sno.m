function [] = make_IASI_hrCRIS_sno(req_date);
%
% function [] = make_IASI_hrCRIS_sno(par1)
%
% produce files of SNO using CRIS hi-res from CCAST processing at UMBC.ASL
%
%1. Sources & Inputs:
%   /asl/data/cris/ccast/sdr60_hr/YYYY/JJJ/SDR_d20130828_t0006589.mat
%   /asl/data/IASI/L1C/YYYY/MM/DD/IASI_xxx_1C_M01_20131020005656Z_20131020005959Z
%   
%2. Destination & Outputs:
%   /asl/s1/chepplew/projects/sno/iasi_cris/HR/
%
%Notes:  HM: L1a_err, in the SDR files is a 30 x nscan array, one value for each FOR:
%        1 = bad, 0 = OK.   
%

cd /home/chepplew/projects/sno/makeSNO/batchJobs;

addpath /asl/matlib/time                             % tai2dnum
addpath /home/chepplew/projects/sno/makeSNO
addpath /home/chepplew/projects/sno/makeSNO/src
%addpath /asl/rtp_prod/iasi/readers                   % readl1c_epsflip_all
%addpath /home/chepplew/myLib/matlib/aslutil          % utc2tai2000 used by readl1c_epsflip_all
%addpath /home/chepplew/myLib/matlib                  % tai2utc1958.m
%
savDir = '/asl/s1/chepplew/data/sno/iasi_cris/HR/';

% Criteria for separation and delay
maxDtim   = 0.014;                  % day. 600/86400 secs =  Mattime is decimal day (20 mins=0.0139)
maxDphi   = 0.18;                   % deg. 0.07 = 7.8 km Earth radius: 6371 km. 1 deg = 111 km.

% Check date string
whos req_date; disp(req_date); fprintf('\n');
try 
   D = datenum(req_date,'YYYY/MM/DD');
catch
   error('Incorrect Date Format')
end

% Get day to process & convert to required formats
% e.g. req_date = '2014/12/05';
strYr   = req_date(1:4);   strMn = req_date(6:7);    strDy = req_date(9:10);
numYr   = str2num(strYr);  numMn = str2num(strMn);   numDy = str2num(strDy);
junk    = sprintf('%4d/%02d/%02d',numYr-1,12,31);
jday    = datenum(req_date)-datenum(junk);  clear junk;           % needed for CRIS directory

% ************    Get CRIS_hr granule files  *****************
crDir = '/asl/data/cris/ccast/sdr60_hr/';
crDir = sprintf('%s%s/%02d/',crDir, strYr,jday);
crLst = dir(strcat(crDir, 'SDR_d', strYr, strMn, strDy,'_t*.mat'));
disp(['Found ' num2str(numel(crLst)) ' CrIS CCAST L1c Files']);

clat = []; clon = []; ctim = [];  cfov = [];  cxtrak = []; catrak = [];
for fn = 1:numel(crLst);
  load(strcat(crDir,crLst(fn).name));
  
  errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
  junk   = reshape(geo.FORTime(~L1a_err),[],1);       % <- (30 x 60) per granule. nu-secs since 1958
    %[yr mn dy hr] = tai2utc1958(junk*1.e-6);
    %mins = fix((hr - fix(hr))*60.0);  
    %secs = ( ((hr - fix(hr))*60.0 - mins)*60.0);
    %hr   = fix(hr);
    %tmpt = datenum(yr, mn, dy, hr, mins, secs);    clear junk;
  tmpt = iet2dnum(junk);    clear junk;
  clat = [clat; reshape(geo.Latitude(:,~L1a_err),[],1)]; 
  clon = [clon; reshape(geo.Longitude(:,~L1a_err),[],1)];         % <- (16200 x 1) per granule

  cradLW = reshape(rLW(:,:,~L1a_err),717,[],1);             % <- (717 x 16200) per gran
  cradMW = reshape(rMW(:,:,~L1a_err),869,[],1);             % <- (869 x 16200) per gran
  cradSW = reshape(rSW(:,:,~L1a_err),637,[],1);             % <- (637 x 16200) per gran

% generate FOV, along-track, cross-track indexes - only works with fixed array size
  if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
  [sz1 sz2 sz3] = size(geo.Latitude);
  if(sz1 ~= 9 | mod(sz2,30) ~= 0) fprintf('ERROR: granule size is wrong\n'); end
% granule size is okay:- proceed;
  count = 0;
  for m=1:sz3
    for i=1:sz2 
      for j=1:sz1
        tmpv(j, i, m) = j;          % <- nominal (9 x 30 x 60)
        tmpx(j, i, m) = i;
	tmpa(j, i, m) = m;
      end
    end
  end
  cfov   = [cfov; reshape(tmpv(:,~L1a_err),[],1)];     clear tmpv;
  cxtrak = [cxtrak; reshape(tmpx(:,~L1a_err),[],1)];   clear tmpx;
  catrak = [catrak; reshape(tmpa(:,~L1a_err),[],1)];   clear tmpa;
  
% Need to expand ctim to match every FOV so that sub-setting works
  junk   = repmat(tmpt,9,1);           clear tmpt;
  ctim   = [ctim; reshape(junk,[],1)]; clear junk;           % <- (16200 x 1) times for every FOV

  fprintf(1,'.');
end
fprintf(1,'\n');

% Subset to center tracks 15,16 in preparation for SNO processing.
  cntr = find(cxtrak == 15 | cxtrak == 16);
  disp(['Found ' num2str(numel(cntr)) ' CrIS center track FOVs'])
  cCnLat = clat(cntr);                         % <- [N x 1] arrays
  cCnLon = clon(cntr);
  cCnTim = ctim(cntr);
  cCnFov = cfov(cntr);
%{
 figure(1);clf(1);simplemap(cCnLat,cCnLon,cCnTim-cCnTim(1));
%}
% *******************   Get IASI granule files  *******************
fprintf('loading IASI data\n');
iadir = strcat('/asl/data/IASI/L1C/',req_date,'/');
iaLst = dir(strcat(iadir,'IASI_xxx_1C_M02*'));
disp(['Found ' num2str(numel(iaLst)) ' IASI L1C files']);

iafov = [;]; ialat = [;]; ialon = [;]; iatim = [;]; iasolzen = [;]; iascnlin = [;];
iaAtrk = []; iaXtrk = [];
for fn = 1:numel(iaLst);                        % fails at fn=423 IASI_xxx_1C_M02_20130827213858Z_20130827214154Z.gz corrupt file
 s = readl1c_epsflip_all(strcat(iadir,iaLst(fn).name)); 
 iafov    = [iafov; reshape(s.IASI_FOV',[],1)];            % [N x 4] -> [2640 x 1]
 ialat    = [ialat; reshape(s.Latitude',[],1)];
 ialon    = [ialon; reshape(s.Longitude',[],1)];
 iasolzen = [iasolzen; reshape(s.Solar_Zenith',[],1)];
 iascnlin = [iascnlin; reshape(s.Scan_Line',[],1)];        % [690 x 4] values 1:23 (4 times)
    junk  = reshape(s.Time2000',[],1);                     % 2000 epoch
 iatim    = [iatim; iasi2mattime(junk)];         % convert to matlab time
 fprintf(1,'.');
end
fprintf(1,'\n');
 ixtrk = [;]; iatrk = [;];
 for j = 1:fix(size(iascnlin,1)/120)              % retain row x col structure
  is = (j-1)*30 +1;  ie = j*30;
  iatrk(is:ie,[1 2 3 4]) = j;
  ixtrk(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';    
 end
%
ixtrk  = reshape(ixtrk',[],1);
iatrk  = reshape(iatrk',[],1);

% Subset center track FOVs
 icntr = find(ixtrk == 15 | ixtrk == 16);
 disp(['Found ' num2str(numel(icntr)) ' IASI center track FOVs'])
 iCnTim = iatim(icntr);                           % <- [N x 1] arrays
 iCnLat = ialat(icntr);
 iCnLon = ialon(icntr);
 iCnFov = iafov(icntr);
%{
 addpath /asl/matlib/aslutil
 figure(1);clf(1);simplemap(iCnLat,iCnLon,iCnTim-iCnTim(1));
%}

% --------------------------------------------------------------
%        pre-compute position vectors - saves a heap of time
% --------------------------------------------------------------
fprintf('computing position vectors\n');
P1 = [;;];                              % CRIS
P2 = [;;];                              % IASI
for ic = 1:length(cCnLat)               % Number of centres.
    P1 = [P1; cos(cCnLat(ic)*pi/180.0) * cos(cCnLon(ic)*pi/180.0), ...
         cos(cCnLat(ic)*pi/180.0)*sin(cCnLon(ic)*pi/180.0), sin(cCnLat(ic)*pi/180.0) ];
end
for ii = 1:length(iCnLat)
    P2 = [P2; cos(iCnLat(ii)*pi/180.0) * cos(iCnLon(ii)*pi/180.0), ...
         cos(iCnLat(ii)*pi/180.0)*sin(iCnLon(ii)*pi/180.0), sin(iCnLat(ii)*pi/180.0) ];
end    
%
% save(strcat(savDir,'sno_Iasi_hrCris_ln145_',strYr, strMn, strDy,'.mat'));

% Faster Algorithm 

dti1  = 2.5347e-06;        % <- (day) IASI time between adjacent FORs within Centre group
dti2  = 9.2593e-05;        % <- IASI time between each centre FOV group
dtc1  = 0.200/86400;       % <- CRIS time between adjacent FORs within center group
dtc2  = 7.800/86400;       % <- CRIS time between each center FOR group
wndta = fix(8*maxDtim/dti2)+1;    % <- no. IASI samples in time window set by maxDtim.

irsa  = max(1,find(iCnTim > cCnTim(1),1) - wndta);    % sample to start IASI
iren  = irsa + 2*wndta;

%{
fprintf('Finding nearest approach\n');
smPos = [;];  smTd = []; smPhi = 75.0; 
for jj = 5:18:length(cCnLat)
  for ii = irsa:8:length(iCnLat)
    cnPhi = real(acos( sum(P1(jj,:).*P2(ii,:)) )*180/pi);
    if(cnPhi < smPhi) 
      smPhi   = [smPhi,cnPhi]; 
      smPos   = [smPos;[jj,ii]];                    % pos(ii: AIRS index. jj: CRIS)
      smTd    = [smTd,(iCnTim(ii) - cCnTim(jj))];
    end
  end
  if(~mod(jj,365)) fprintf('.'); end
end
%}
%{ 
jsz = size(smPhi,2);
figure(1);clf(1);plot((1:jsz),smPhi,'k.',(1:jsz-1),smTd*24*60,'g.');
  legend('smPhi','smTd mins');
figure(1);scatter(smTd*24*60, smPhi(1:end-1));xlabel('delay mins');ylabel('distance deg');
  grid on;
% save('/asl/s1/chepplew/projects/sno/airs_cris/HR/20130827_snapshot.mat');
%}
%        Compute separations and save indexes when criteria are met
fprintf(1,'Computing separations\n');
dist = 0.0; m = 0; k = 1;
pos  = [;]; tdiff = [];
tic
for i = 1:1:length(cCnLat)
  %irst = max(1,find(iCnTim > cCnTim(i),1) - wndta);
  %irsp = min(irst + 2*wndta, length(iCnLat));
  for j = 1:1:length(iCnTim)  % irst:1:irsp
    Dtim = abs(iCnTim(j) - cCnTim(i));
    if Dtim <= maxDtim      
      m = m+1;
      Phi = real(acos( sum(P1(i,:).*P2(j,:)) )*180/pi);    % P1: CrIS, P2: IASI.
      if Phi <= maxDphi                  % phi = 0.18 deg (0.0031 rad) => 20 km.
	dist(k) = Phi;                   % dist;
	pos     = [pos;[i,j]];           % pos(i: CrIS index. j: IASI)
	tdiff   = [tdiff,(iCnTim(j) - cCnTim(i))];
	k = k+1;
      end
    end
  end
  if(~mod(i,1000)) fprintf(1,'.'); end
end
toc
dist = real(dist);                         % can get complex phi.
fprintf(1,'Found %i matchups\n',k);

% save('/asl/s1/chepplew/projects/sno/iasi_cris/HR/20130827_t007d15km.mat');

% This gives us duplicate hits - so need to select only unique pairs.
%  this is done by sequentially testing for uniquness on each sensor.
if(size(pos,1) > 4 )
  un1 = [;]; un2 = [;]; upos = [;];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
  fprintf('Number unique SNOs %d %d\n',size(upos));

% Record time and space separations (pre reloading of files)
  pre.iSnoTim   = iCnTim(upos(:,2));    pre.cSnoTim = cCnTim(upos(:,1));
  pre.iSnoLat   = iCnLat(upos(:,2));    pre.cSnoLat = cCnLat(upos(:,1));
  pre.iSnoLon   = iCnLon(upos(:,2));    pre.cSnoLon = cCnLon(upos(:,1));
  pre.utdiff    = pre.iSnoTim - pre.cSnoTim;    
  pre.cSnoFov   = cCnFov(upos(:,2));
%
  uP1    = P1(upos(:,1),:);      uP2 = P2(upos(:,2),:);    % P1: CRIS. P2: IASI
  for i = 1:size(upos,1)
     uPhi(i)   = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi );
  end
  pre.uPhi = uPhi;

% *********** SECTION 3 - Get Sensor Obs for the SNOs ***********
% ************  Load up CRIS SNO Obs data ************

  fprintf(1,'Loading CIRS SNO Obs\n');
  ps = 1; pe = 0; cCount = 0;  tcrLW = [;]; tcrMW = [;]; tcrSW = [;];
  tcatrk = [];    tcxtrk = []; tcfov = [];
  tctim  = [];     tclat = []; tclon = [];  tcszn = [];  tcsazn = []; tclnfr = [];
  for fn = 1:numel(crLst)
    load(strcat(crDir,crLst(fn).name));
    errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
    junk   = reshape(geo.FORTime(~L1a_err),[],1);       % <- (30 x 60) per granule. nu-secs since 1958
      %[yr mn dy hr] = tai2utc1958(junk*1.e-6);
      %mins = fix((hr - fix(hr))*60.0);  
      %secs = ( ((hr - fix(hr))*60.0 - mins)*60.0);
      %hr   = fix(hr);
      %tmpt = datenum(yr, mn, dy, hr, mins, secs);        clear junk;
  tmpt = iet2dnum(junk);    clear junk;
% Need to expand ctim to match every FOV so that sub-setting works
      junk = repmat(tmpt,9,1);           clear tmpt;
    ctim = reshape(junk,[],1);           clear junk;           % <- (16200 x 1) times for every FOV
    clat = reshape(geo.Latitude(:,~L1a_err),[],1); 
    clon = reshape(geo.Longitude(:,~L1a_err),[],1); 
    cszn = reshape(geo.SolarZenithAngle(:,~L1a_err),[],1);
    cradLW = reshape(rLW(:,:,~L1a_err),717,[],1);             % <- (717 x 16200) per gran
    cradMW = reshape(rMW(:,:,~L1a_err),869,[],1);             % <- (869 x 16200) per gran
    cradSW = reshape(rSW(:,:,~L1a_err),637,[],1);             % <- (637 x 16200) per gran

% generate FOV, along-track, cross-track indexes - only works with fixed array size
    if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
    [sz1 sz2 sz3] = size(geo.Latitude);
    if(sz1 ~= 9 | mod(sz2,30) ~= 0) fprintf('ERROR: granule size is wrong\n'); end
% granule size is okay:- proceed;
    for m=1:sz3
      for i=1:sz2 
        for j=1:sz1
          tmpv(j, i, m) = j;          % <- nominal (9 x 30 x 60)
          tmpx(j, i, m) = i;
  	  tmpa(j, i, m) = m;
        end
      end
    end
    cfov   = reshape(tmpv(:,~L1a_err),[],1);     clear tmpv;
    cxtrak = reshape(tmpx(:,~L1a_err),[],1);     clear tmpx;
    catrak = reshape(tmpa(:,~L1a_err),[],1);     clear tmpa;

% subset on center track (FORs 15 and 16)
    tCnInd   = find(cxtrak == 15 | cxtrak == 16);
    tCnLat   = clat(tCnInd);
    pe       = ps + numel(tCnLat) - 1;
    fprintf('%d %d ',ps,pe);
    junk   = ismember(upos(:,1),[ps:pe]);       cSams = upos(junk,1);  clear junk;
    if ( length(cSams >= 1) )
      fprintf('%d %d ',fn,length(cSams));
  % subset onto center track b4 using upos index.
      indx = find(cxtrak == 15 | cxtrak == 16);
      junk = cradLW(:,indx);      tcrLW  = [tcrLW, junk(:,cSams-ps+1)];  clear junk;
      junk = cradMW(:,indx);      tcrMW  = [tcrMW, junk(:,cSams-ps+1)];  clear junk;
      junk = cradSW(:,indx);      tcrSW  = [tcrSW, junk(:,cSams-ps+1)];  clear junk;
      junk = cfov(indx);          tcfov  = [tcfov; junk(cSams-ps+1)];    clear junk;
      junk = cxtrak(indx);        tcxtrk = [tcxtrk; junk(cSams-ps+1)];   clear junk;
      junk = catrak(indx);        tcatrk = [tcatrk; junk(cSams-ps+1)];   clear junk;
      junk = clat(indx);          tclat  = [tclat; junk(cSams-ps+1)];    clear junk;
      junk = clon(indx);          tclon  = [tclon; junk(cSams-ps+1)];    clear junk;
      junk = ctim(indx);          tctim  = [tctim; junk(cSams-ps+1)];    clear junk;
      junk = cszn(indx);          tcszn  = [tcszn; junk(cSams-ps+1)];    clear junk;
    end
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);  
  end      % <- end: for fn=1:numel(CrLst) 
%{
 figure(1);clf(1);simplemap(tclat, tclon, tctim)
%}
% ****************   Load up IASI SNO Obs data ************
  fprintf('Loading IASI SNO Obs\n');
  ps = 1; pe = 0; iCount = 0;  tirad = [;];
  tiatrk = []; tixtrk = []; tifov = [];  iarad = [;];  tirad = [;];
  titim  = [];  tilat = []; tilon = [];  tiszn = [];  tisazn = []; tilnfr = [];
  iafov  = [];  ialat = []; ialon =[];   iatim = [];   iaszn = []; iascln = []; 
  for fn = 1:numel(iaLst);                        % fails at fn=423 IASI_xxx_1C_M02_20130827213858Z_20130827214154Z.gz corrupt file
   s = readl1c_epsflip_all(strcat(iadir,iaLst(fn).name)); 
     junk  = permute(s.IASI_Radiances,[2 1 3]);            % [660 x 8461]
     junk2 = reshape(junk,[],8461,1);                      % [4 x 660 x 8461] per granule
   iarad    = junk2;  clear junk junk2;                    % [2640 x 8461] per gran
   iafov    = reshape(s.IASI_FOV',[],1);                   % [N x 4] -> [2640 x 1] arrays
   ialat    = reshape(s.Latitude',[],1);
   ialon    = reshape(s.Longitude',[],1);
     junk = reshape(s.Time2000',[],1);                     % 2000 epoch
   iatim    = iasi2mattime(junk);  clear junk;             % convert to matlab time
   iaszn    = reshape(s.Solar_Zenith',[],1);
   iascln   = reshape(s.Scan_Line',[],1);                  % [690 x 4] values 1:23 (4 times)
%   ixtrk = [;]; iatrk = [;]; junk = [;];
%   for j = 1:fix(size(iascln,1)/120)                        % retain row x col structure
%    is = (j-1)*30*4 +1;  ie = j*30*4;
%    iatrk(is:ie)  = [ones(1,30)*j; ones(1,30)*j; ones(1,30)*j; ones(1,30)*j];
%    ixtrk(is:ie)  = [1:30; 1:30; 1:30; 1:30];    
%   end
%   iatrk = iatrk';  ixtrk = ixtrk';
 ixtrk = [;]; iatrk = [;];
 for j = 1:fix(size(iascln,1)/120)              % retain row x col structure
  is = (j-1)*30 +1;  ie = j*30;
  iatrk(is:ie,[1 2 3 4]) = j;
  ixtrk(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';    
 end
 ixtrk  = reshape(ixtrk',[],1);  iatrk  = reshape(iatrk',[],1);
% subset on center track (FORs 15 and 16) to update counters
    tCnInd   = find(ixtrk == 15 | ixtrk == 16);
    tCnLat   = ialat(tCnInd);
    pe       = ps + numel(tCnLat) - 1;
    fprintf('%d %d ',ps,pe);
    junk   = ismember(upos(:,2),[ps:pe]);       iSams = upos(junk,2);  clear junk;
    if ( length(iSams) >= 1 )
      fprintf('%d %d ',fn,length(iSams));
  % subset onto center track b4 using upos index.
      indx = find(ixtrk == 15 | ixtrk == 16);
      junk = iarad(indx,:);        tirad  = [tirad; junk(iSams-ps+1,:)];  clear junk;
      junk = iafov(indx);          tifov  = [tifov; junk(iSams-ps+1)];    clear junk;
      junk = ixtrk(indx);          tixtrk = [tixtrk; junk(iSams-ps+1)];   clear junk;
      junk = iatrk(indx);          tiatrk = [tiatrk; junk(iSams-ps+1)];   clear junk;
      junk = ialat(indx);          tilat  = [tilat; junk(iSams-ps+1)];    clear junk;
      junk = ialon(indx);          tilon  = [tilon; junk(iSams-ps+1)];    clear junk;
      junk = iatim(indx);          titim  = [titim; junk(iSams-ps+1)];    clear junk;
      junk = iaszn(indx);          tiszn  = [tiszn; junk(iSams-ps+1)];    clear junk;
    end
    ps = pe + 1;
    iCount = iCount + length(iSams); fprintf('%d \n',iCount);  
  end      % <- end: for fn=1:numel(iaLst) 
%{
 figure(1);clf(1);simplemap(tilat, tilon, titim)
%}

% Saving data ****************
  %fprintf('Saving data\n');
  %save(strcat(savDir,'sno_iasi_hrcris_',strYr,strMn,strDy,'.mat'));    % 17Jul2014
  clear itime ilat ilon isolzen ri;
  clear ctime clat clon csolzen csatzen rc tdiff dist;
  itime = titim;  ilat = tilat; ilon = tilon;  isolzen = tiszn; ri = tirad; 
  ctime = tctim;  clat = tclat; clon = tclon;  csolzen = tcszn; cfov = tcfov;
    rc  = [tcrLW; tcrMW; tcrSW];  fc = [vLW; vMW; vSW];
    csatzen = tcsazn; dist = pre.uPhi'; tdiff = pre.utdiff;
  savVars = {'strYr','strMn','strDy', 'itime','ilat','ilon','isolzen','ri',...
            'ctime','clat','clon','csolzen','csatzen','rc','fc','tdiff',...
	    'dist','upos','maxDtim','maxDphi'};
  savFN = strcat(savDir,'sno_iasi_crisHr_asl_',strYr,strMn,strDy,'.mat')
  disp(['Saving file: ' savFN]);
  save(savFN,savVars{:});


end   % ********** END if(size(pos,1) > 1)

%{
load('~strow/Matlab/Iasi/iasi_f.mat');     % fiasi [8461 x 1]
np = size(rc,2);
bc = real(rad2bt(fc,rc));
bi = real(rad2bt(fiasi,ri'));
 ichn = find(fiasi > 900,1);     % 1022
 cchn = find(fc > 900,1);        % 404
 figure(3);clf(3);plot( [1:np],bc(cchn,:),'b.',[1:np],bi(ichn,:),'g.');
 figure(3);clf;scatter(dist,tdiff*24*60,'.');
 nanmean(bc(cchn,:),2) - nanmean(bi(ichn,:),2)  % 236.9370 - 236.9432
%}
