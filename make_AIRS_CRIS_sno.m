function [] = make_AIRS_CRIS_sno(req_date, cris_res)
%{
makeAIRS_hrCRIS_sno

to produce files of SNO using CRIS hi-res from CCAST processing at UMBC.ASL

1. Sources & Inputs:
   /asl/data/cris/ccast/sdr60_hr/YYYY/JJJ/SDR_d20130828_t0006589.mat
   /asl/data/airs/AIRIBRAD/YYYY/JJJ/AIRS.2013.08.28.239.L1B.AIRS_Rad.v5.0.22.0.G13343142628.hdf
   
2. Destination & Outputs:
   /asl/s1/chepplew/projects/sno/airs_cris/HR/   

Notes:  HM: L1a_err, in the SDR files is a 30 x nscan array, one value for each FOR:
        1 = bad, 0 = OK.   
%}

cd /home/chepplew/projects/sno/makeSNO;
addpath /asl/matlib/time                             % tai2dnum
addpath /asl/rtp_prod/airs/readers                   % readl1b_all
addpath /asl/rtp_prod/airs/utils                     % f_default_1lb.mat
warning('off');

% Check date string
whos req_date; disp(req_date); fprintf('\n');
try 
   D = datenum(req_date,'YYYY/MM/DD');
catch
   error('Incorrect Date Format')
   return
end

% Check CrIS Resolution
cris_res = upper(cris_res);
all_res  = {'HIGH','LOW'};
if(~ismember(cris_res,all_res))
  error('Invalid CrIS resolution. Choices are: high or low');
  return
end
  
if(strcmp(cris_res,'HIGH')) 
  savDir = '/asl/s1/chepplew/data/sno/airs_cris/HR/';
  npLW = 717;  npMW = 869;   npSW = 637;
end
if(strcmp(cris_res,'LOW'))
  savDir = '/asl/s1/chepplew/data/sno/airs_cris/LR/';
  npLW = 717;  npMW = 437;   npSW = 163;
end


% Criteria for separation and delay
maxDtim   = 0.007;                % day. 600/86400 secs =  Mattime is decimal day (20 mins=0.0139)
maxDphi   = 0.13;                   % deg. 0.07 = 7.8 km Earth radius: 6371 km. 1 deg = 111 km.

% Get day to process & convert to required formats
strYr   = req_date(1:4);   strMn = req_date(6:7);    strDy = req_date(9:10);
numYr   = str2num(strYr);  numMn = str2num(strMn);   numDy = str2num(strDy);
junk    = sprintf('%4d/%02d/%02d',numYr-1,12,31);
jday    = datenum(req_date)-datenum(junk);  clear junk;           % needed for CRIS directory

% ************    Get CRIS_hr daily files  *****************
if(strcmp(cris_res,'HIGH')) crDir = '/asl/data/cris/ccast/sdr60_hr/';  end
if(strcmp(cris_res,'LOW'))  crDir = '/asl/data/cris/ccast/sdr60/';  end
crDir = sprintf('%s%s/%03d/',crDir, strYr,jday);
crLst = dir(strcat(crDir, 'SDR_d', strYr, strMn, strDy,'_t*.mat'));
disp(['Found ' num2str(numel(crLst)) ' CrIS CCAST L1c Files']);

clat = []; clon = []; ctim = [];  cfov = [];  cxtrak = []; catrak = [];
for fn = 1:numel(crLst);
  load(strcat(crDir,crLst(fn).name));
  
  errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
   junk  = reshape(geo.FORTime(~L1a_err),[],1);       % <- (30 x 60) per granule. nu-secs since 1958
  tmpt   = iet2dnum(junk);    clear junk;
  clat   = [clat; reshape(geo.Latitude(:,~L1a_err),[],1)]; 
  clon   = [clon; reshape(geo.Longitude(:,~L1a_err),[],1)];         % <- (16200 x 1) per granule

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

  fprintf('.');
end
fprintf('\n');

% Subset to center tracks 15,16 in preparation for SNO processing.
  cntr = find(cxtrak == 15 | cxtrak == 16);
  disp(['Found ' num2str(numel(cntr)) ' CrIS center track FOVs'])
  cCnLat = clat(cntr);
  cCnLon = clon(cntr);
  cCnTim = ctim(cntr);
  cCnFov = cfov(cntr);

% *******************   Get AIRS daily files  *******************
fprintf('Loading AIRS geo\n'); 
arDir = strcat('/asl/data/airs/AIRIBRAD/',strYr,'/',sprintf('%03d',jday),'/');
arLst = dir(sprintf('%sAIRS.%4d.%02d.%02d*.hdf',arDir,numYr,numMn,numDy));
disp(['Found ' num2str(numel(arLst)) ' AIRIBRAD Files']);

alat = []; alon = []; atim = []; aAtrk = []; aXtrk = []; asolzn = [];
for fn = 1:numel(arLst);
 [eqXtai,f,gdata] = readl1b_all(strcat(arDir,arLst(fn).name));
 alat    = [alat; gdata.rlat'];
 alon    = [alon; gdata.rlon'];
 atim    = [atim; airs2dnum(gdata.rtime)'];               % convert to matlab time.
 aAtrk   = [aAtrk; gdata.atrack'];
 aXtrk   = [aXtrk; gdata.xtrack'];
 asolzn  = [asolzn; gdata.solzen'];
 fprintf('.');
end
fprintf('\n');

cntr = find(aXtrk == 43 | aXtrk == 44 | aXtrk == 45 | aXtrk == 46 | ...
            aXtrk == 47 | aXtrk == 48);
disp(['Found ' num2str(numel(cntr)) ' AIRS center track FOVs'])

aCnTim = atim(cntr);
aCnLat = alat(cntr);
aCnLon = alon(cntr);
aCnSzn = asolzn(cntr);

%  ****************** SECTION 3: Get indexes of SNOs  ********************
%            pre-compute position vectors - saves a heap of time *********

fprintf('computing position vectors\n');
P1 = [;];                % AIRS position
P2 = [;];                % CRIS position
for ii = 1:length(aCnLat)                                   % Number of Swath.
    P1 = [P1;[ cos(aCnLat(ii)*pi/180.0) * cos(aCnLon(ii)*pi/180.0), ...
         cos(aCnLat(ii)*pi/180.0)*sin(aCnLon(ii)*pi/180.0), sin(aCnLat(ii)*pi/180.0) ]];
end
for ii = 1:length(cCnLat)
    P2 = [P2;[ cos(cCnLat(ii)*pi/180.0) * cos(cCnLon(ii)*pi/180.0), ...
         cos(cCnLat(ii)*pi/180.0)*sin(cCnLon(ii)*pi/180.0), sin(cCnLat(ii)*pi/180.0) ]];
end    

% Faster Algorithm 

dta1  = 0.0225/86400;      % <- AIRS time between adjacent samples within Centre group
dta2  = 2.555/86400;       % <- AIRS time between each centre FOV group
dtc1  = 0.200/86400;       % <- CRIS time between adjacent FORs within center group
dtc2  = 7.800/86400;       % <- CRIS time between each center FOR group
wndta = fix(6*maxDtim/dta2)+1;    % <- no. AIRS samples in time window set by maxDtim.

arsa  = find(aCnTim > cCnTim(1),1) + wndta;    % sample to start AIRS
aren  = arsa + 2*wndta;


smPos = [;];  smTd = []; smPhi = 75.0; 
for jj = 5:18:length(cCnLat)
  for ii = arsa:6:length(aCnLat)
    cnPhi = real(acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi);
    if(cnPhi < smPhi) 
      smPhi   = [smPhi,cnPhi]; 
      smPos   = [smPos;[ii,jj]];                    % pos(ii: AIRS index. jj: CRIS)
      smTd    = [smTd,(aCnTim(ii) - cCnTim(jj))];
    end
  end
  if(~mod(jj,365)) fprintf('.'); end
end
 
% save('/asl/s1/chepplew/projects/sno/airs_cris/HR/20130827_snapshot.mat');

%        Compute separations and save indexes when criteria are met
fprintf('Computing separations\n');
dist = 0.0; m = 0; k = 1;
pos = [;]; tdiff = [];
tic
for jj = 1:1:length(cCnLat)
  arst = max(1,find(aCnTim > cCnTim(jj),1) - wndta);
  arsp = min(arst + 2*wndta, length(aCnLat));
  for ii = arst:1:arsp
    Dtim = abs(aCnTim(ii) - cCnTim(jj));
    if Dtim <= maxDtim   % 0.007   % maxDtim       
      m = m+1;
      % dist = distance(arCnLat(ja),arCnLon(ja),iaCnLat(ji),iaCnLon(ji)); % way too slow
      Phi = real(acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi);
      if Phi <= maxDphi                % phi = 0.18 deg (0.0031 rad) => 20 km.
	dist(k) = Phi;                % dist;
	pos     = [pos;[ii,jj]];      % pos(ii: AIRS index. jj: CRIS)
	tdiff   = [tdiff,(aCnTim(ii) - cCnTim(jj))];
	k = k+1;
      end
    end
  end
  if(~mod(jj,1000)) fprintf('.'); end
end
toc
dist = real(dist);                         % can get complex phi.
fprintf('\n');

% save('/asl/s1/chepplew/projects/sno/airs_cris/HR/20130827_t008d2.mat');

fprintf('Number matches found %d\n',size(pos,1));

% This gives us duplicate hits - so need to select only unique pairs.
%  this is done by sequentially testing for uniqueness on each sensor.
if(size(pos,1) > 10 )
  un1 = [;]; un2 = [;]; upos = [;];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
  fprintf('Number unique SNOs %d %d\n',size(upos));

% Record time and space separations (pre reloading of files)
  pre.aSnoTim   = aCnTim(upos(:,1));    pre.cSnoTim = cCnTim(upos(:,2));
  pre.aSnoLat   = aCnLat(upos(:,1));    pre.cSnoLat = cCnLat(upos(:,2));
  pre.aSnoLon   = aCnLon(upos(:,1));    pre.cSnoLon = cCnLon(upos(:,2));
  pre.utdiff    = pre.aSnoTim - pre.cSnoTim;    
  pre.cSnoFov   = cCnFov(upos(:,2));
%
  uP1    = P1(upos(:,1),:);      uP2 = P2(upos(:,2),:);    % P1: AIRS, P2: CRIS.
  for i = 1:size(upos,1)
     uPhi(i)   = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi );
  end
  pre.uPhi = uPhi;


% **** Section 4 Reload AIRS and CRIS files and save only required Obs.
%      ----------------------------------------------------------------
  fprintf('loading AIRS obs\n');
  ps = 1; pe = 0; aCount = 0;  tarad = [;]; taAtrk = []; taXtrk = []; tasozn = [];
  tasazn  = [];   talnfr = []; tatim = [];  talat  = []; talon  = [];
  for fn = 1:numel(arLst); % fn = 1;
    [eqXtai,f,gdata] = readl1b_all(strcat(arDir,arLst(fn).name));
    atim   = airs2dnum(gdata.rtime);               % convert to matlab time.
    tXtrk  = gdata.xtrack;                             % select centre track FOVs
    tCnInd = find(tXtrk == 43 | tXtrk == 44 | tXtrk == 45 | tXtrk == 46 | ...
                  tXtrk == 47 | tXtrk == 48);
    pe = ps + numel(gdata.rlat(tCnInd)) - 1;

  %  tests which data points we want to keep match to appropriate FOR
    fprintf('%d %d ',ps,pe);
    junk   = ismember(upos(:,1),[ps:pe]);       aSams = upos(junk,1);    clear junk;
    if ( length(aSams >= 1)  )
      fprintf('%d %d ',fn,length(aSams));
      % subset onto center track b4 using pos index.
      indx = find(gdata.xtrack == 43 | gdata.xtrack == 44 | gdata.xtrack == 45 | ...
                  gdata.xtrack == 46 | gdata.xtrack == 47 | gdata.xtrack == 48);
      junk = gdata.robs1(:,indx);  tarad  = [tarad, junk(:,aSams-ps+1) ]; clear junk;
      junk = gdata.xtrack(indx);   taXtrk = [taXtrk, junk(aSams-ps+1)];   clear junk;
      junk = gdata.atrack(indx);   taAtrk = [taAtrk, junk(aSams-ps+1)];   clear junk;
      junk = gdata.rlat(indx);     talat  = [talat, junk(aSams-ps+1)];    clear junk;
      junk = gdata.rlon(indx);     talon  = [talon, junk(aSams-ps+1)];    clear junk;
      junk = gdata.solzen(indx);   tasozn = [tasozn, junk(aSams-ps+1)];   clear junk;
      junk = atim(indx);           tatim  = [tatim, junk(aSams-ps+1)];    clear junk;
      junk = gdata.landfrac(indx); talnfr = [talnfr,junk(aSams-ps+1)];    clear junk;
      junk = gdata.satzen(indx);   tasazn = [tasazn,junk(aSams-ps+1)];    clear junk;
    end
    ps = pe+1;
    aCount = aCount + length(aSams);  fprintf('%d \n',aCount);
%   fprintf('.');
  end
  fprintf('\n');

% ************  Load up CRIS SNO Obs data ************

  fprintf('Loading CIRS SNO Obs\n');
  ps = 1; pe = 0; cCount = 0;  tcrLW = [;]; tcrMW = [;]; tcrSW = [;];
  tcatrk = [];    tcxtrk = []; tcfov = [];
  tctim  = [];     tclat = []; tclon = [];  tcsozn = []; tcsazn = []; tclnfr = [];
  snVara = [;];    cnObs = [;];
  for fn = 1:numel(crLst)
    load(strcat(crDir,crLst(fn).name));
    errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
    junk   = reshape(geo.FORTime(~L1a_err),[],1);       % <- (30 x 60) per granule. nu-secs since 1958
    tmpt   = iet2dnum(junk);    clear junk;
% Need to expand ctim to match every FOV so that sub-setting works
      junk = repmat(tmpt,9,1);           clear tmpt;
    ctim   = reshape(junk,[],1);         clear junk;           % <- (16200 x 1) times for every FOV
    clat   = reshape(geo.Latitude(:,~L1a_err),[],1); 
    clon   = reshape(geo.Longitude(:,~L1a_err),[],1);         % <- (16200 x 1) per granule
    cszn   = reshape(geo.SolarZenithAngle(:,~L1a_err),[],1);
    csaz   = reshape(geo.SatelliteZenithAngle(:,~L1a_err),[],1);

    cradLW = reshape(rLW(:,:,~L1a_err),npLW,[],1);             % <- (717 x 16200) per gran
    cradMW = reshape(rMW(:,:,~L1a_err),npMW,[],1);             % <- (869 x 16200) per gran
    cradSW = reshape(rSW(:,:,~L1a_err),npSW,[],1);             % <- (637 x 16200) per gran

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
    cfov     = reshape(tmpv(:,~L1a_err),[],1);     clear tmpv;
    cxtrak   = reshape(tmpx(:,~L1a_err),[],1);     clear tmpx;
    catrak   = reshape(tmpa(:,~L1a_err),[],1);     clear tmpa;

    tCnInd   = find(cxtrak == 15 | cxtrak == 16);
    tCnLat   = clat(tCnInd);
    pe       = ps + numel(tCnLat) - 1;
    fprintf('%d %d ',ps,pe);
    junk   = ismember(upos(:,2),[ps:pe]);       cSams = upos(junk,2);  clear junk;
    if ( length(cSams >= 1) )
      fprintf('%d %d ',fn,length(cSams));
  % subset onto center track b4 using pos index.
      indx = find(cxtrak == 15 | cxtrak == 16);
      junk = cradLW(:,indx);      tcrLW   = [tcrLW,  junk(:,cSams-ps+1)];  clear junk;
      junk = cradMW(:,indx);      tcrMW   = [tcrMW,  junk(:,cSams-ps+1)];  clear junk;
      junk = cradSW(:,indx);      tcrSW   = [tcrSW,  junk(:,cSams-ps+1)];  clear junk;
      junk = cfov(indx);          tcfov   = [tcfov;  junk(cSams-ps+1)];    clear junk;
      junk = cxtrak(indx);        tcxtrk  = [tcxtrk; junk(cSams-ps+1)];    clear junk;
      junk = catrak(indx);        tcatrk  = [tcatrk; junk(cSams-ps+1)];    clear junk;
      junk = clat(indx);          tclat   = [tclat;  junk(cSams-ps+1)];    clear junk;
      junk = clon(indx);          tclon   = [tclon;  junk(cSams-ps+1)];    clear junk;
      junk = ctim(indx);          tctim   = [tctim;  junk(cSams-ps+1)];    clear junk;
      junk = cszn(indx);          tcsozn  = [tcsozn; junk(cSams-ps+1)];    clear junk;
      junk = csaz(indx);          tcsazn  = [tcsazn; junk(cSams-ps+1)];    clear junk;
    end
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);  
  end      % <- end: for fn=1:numel(CrLst) 

% Saving data ****************
  clear atime alat alon asolzen asatzen ra;
  clear ctime clat clon cfov csolzen csatzen rc tdiff dist;
  atime = tatim'; alat = talat; alon = talon; asolzen = tasozn; asatzen = tasazn;   ra = tarad; 
  alandfrac = talnfr; 
  ctime = tctim;  clat = tclat; clon = tclon; csolzen = tcsozn; csatzen = tcsazn; cfov = tcfov;
  rc    = [tcrLW; tcrMW; tcrSW]; fc = [vLW; vMW; vSW];
  csatzen = tcsazn; dist = pre.uPhi'; tdiff = pre.utdiff;
  savFN   = strcat(savDir, strYr, '/sno_airs_cris_asl_',strYr,strMn,strDy,'.mat')
  savVars = {'atime','alat','alon','asolzen','asatzen','ra',...
             'ctime','clat','clon','csolzen','csatzen','cfov','rc',...
	     'fc','tdiff','dist','upos','alandfrac','maxDtim','maxDphi'};
  fprintf('Saving data to file: %s\n',savFN);
  save(savFN, savVars{:});


end

