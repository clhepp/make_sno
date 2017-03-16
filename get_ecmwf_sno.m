function [prof, head, pattr] = get_ecmwf_sno(prof, head, pattr);

% create SNO adjunt RTP file with ECMWF model fields
% Get ctim, clat, clon, from SNO set (ana_jpl_iasi_cris_sno.m)
% Note: the JPL SNO data has IASI time different to CrIS time

% N.B. prof.vectors must be [1 x nobs] (they come rotated form the SNO data).

cd /home/chepplew/projects/sno/makeSNO

addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /asl/rtp_prod2/grib
addpath /asl/matlib/time

% --------------------
% 1. Get the SNO data
% --------------------
srcpth = '/asl/s1/chepplew/data/sno/iasi_cris/JPL/';
srclst = dir([srcpth 'sno_iasi_cris20*.mat']);
disp(['Found ' num2str(numel(srclst)) ' SNO files']);
ip = 1;
   X=load([srcpth srclst(ip).name]);


% -------------------------
% set up the RTP structures
% -------------------------
hattr = {};
hattr{1} = {'header',    'instid',    'CrIS '};
hattr{2} = {'header'    'reader'    'ccast2rtp '};
hattr{3} = {'header'    'topo'    'usgs_deg10_dem '};

head = struct;
Y=load('/home/chepplew/projects/cris/cris_freq_2grd.mat'); fcris = Y.vchan;
Y=load('/home/chepplew/projects/iasi/f_iasi.mat');
head.ichan = Y.ichan_iasi;
head.vchan = Y.f_iasi;
clear Y;

pattr = {};
pattr{1} = {'profiles','model','NA '};

prof=struct;
%prof.rtime = X.ctim;
prof.rtime = dnum2tai(X.itime);
prof.rlat  = X.clat;
prof.rlon  = X.clon;

% ------------------------
% load in the ECMWF fields
% ------------------------
 [prof, head, pattr] = fill_ecmwf(prof, head, pattr);

% -------------------------
% Prep for running the RTA
% ------------------------
% build config struct
addpath /asl/rtp_prod2/emis
addpath /asl/rtp_prod2/util

cfg = struct;
cfg.model = 'ecmwf';

klayers_exec = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
if isfield(cfg, 'klayers_exec')
    klayers_exec = cfg.klayers_exec;
end

sarta_exec  = '/asl/packages/sartaV108/BinV201/sarta_iasi_may09_wcon_nte';
%sarta_exec  = '/asl/packages/sartaV108/BinV201/sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4';
if isfield(cfg, 'sarta_exec')
    sarta_exec = cfg.sarta_exec;
end

head2.pfields = 5;

% Add landfrac, etc.
[head, hattr, prof, pattr] = rtpadd_usgs_10dem(head,hattr,prof,pattr);

% Add Dan Zhou's emissivity and Masuda emis over ocean
% Dan Zhou's one-year climatology for land surface emissivity and
% standard routine for sea surface emissivity
% needs satzen
nobs = size(prof.rlat,2);
prof.satzen = zeros(1,nobs);
fprintf(1, '>>> Running rtp_ad_emis...');
[prof,pattr] = rtp_add_emis_single(prof,pattr);

% ----------------------------------------
%  Run KLAYERS and SARTA
% ----------------------------------------
 tmp = mktemp();
  outfiles = rtpwrite_12(tmp,head,hattr,prof,pattr);
  s1Path = '/tmp/';
  %disp(['tmp = ', tmp]);

  ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
  ofn_1 = [tmp '.kla_1'];  ofn_2 = [tmp '.kla_2'];
  ofn_3 = [tmp '.sar_1'];  ofn_4 = [tmp '.sar_2'];

  % run klayers on first half
  %unix([klayers_exec ' fin=' ifn_1 ' fout=' ofn_1 ' > ' s1Path '/klayers_stdout']);
  unix([klayers_exec ' fin=' ifn_1 ' fout=' ofn_1 ' > /dev/null']);

  % run sarta on first half
  %eval(['! ' sarta_exec ' fin=' ofn_1 ' fout=' ofn_3 ' > sartastdout1.txt']);
  eval(['! ' sarta_exec ' fin=' ofn_1 ' fout=' ofn_3 ' > /dev/null']);

  % run klayers on second half
  %unix([klayers_exec ' fin=' ifn_2 ' fout=' ofn_2 ' > ' s1Path '/klayers_stdout']);
  unix([klayers_exec ' fin=' ifn_2 ' fout=' ofn_2 ' > /dev/null']);

  % run sarta on second half
  %eval(['! ' sarta_exec ' fin=' ofn_2 ' fout=' ofn_4 ' > sartastdout1.txt']);
  eval(['! ' sarta_exec ' fin=' ofn_2 ' fout=' ofn_4 ' > /dev/null']);

  % read the results files back in
  cfin = [tmp '.sar'];

% $$$   [hd ha pd pa] = rtpread_12(cfin);
  [~,~,ptemp,~] = rtpread_12(cfin);
  prof.rcalc = ptemp.rcalc;
  clear ptemp;

  % silently delete temporary files:
  unlink(ifn_1); unlink(ifn_2); unlink(ofn_1); unlink(ofn_2); unlink(ofn_3);
  unlink(ofn_4); unlink(tmp);


%{
% Check for outliers in radiance units using CrIS Obs - calc
addpath /home/chepplew/myLib/matlib/math          % remove_6sigma
disp(['Removing outliers']);
crbias = X.rc - X.i2rc;
clear g;
for i=1:1317
   n  = remove_6sigma(crbias(i,:));
   nn = remove_6sigma(crbias(i,n));
   g(i).n = n(nn);
end

% Now find unique set of bad SNO samples
x = [];
[~, psz] = size(crbias);
for i=1:1317
   x = [x setdiff(1:psz,g(i).n)];
end
x  = unique(x);
ig = setdiff(1:psz,x);
% ------------------------------------------------

addpath /asl/packages/airs_decon/source         % hamm_app
junk     = single(hamm_app(double(prof.rcalc(:,ig))));
btic_h   = real(rad2bt(head.vchan,junk));
junk     = single(hamm_app(double(X.rc(:,ig))));
btco_h   = real(rad2bt(fcris,junk));
junk     = single(hamm_app(double(X.ri(:,ig))));
btio_h   = real(rad2bt(head.vchan,junk));
junk     = single(hamm_app(double(X.i2rc(:,ig))));
bti2co_h = real(rad2bt(fcris,junk));

bticm    = nanmean(btic_h,2);
btcom    = nanmean(btco_h,2);
btiom    = nanmean(btio_h,2);
bti2com  = nanmean(bti2co_h,2);
btics    = nanstd(btic_h,0,2);
btcos    = nanstd(btco_h,0,2);
btios    = nanstd(btio_h,0,2);
whos bticm btcom btiom btcos btios btics
figure(1);clf;plot(fcris,btcom,'-',fiasi,btiom,'-',fiasi,bticm,'-',fcris,bti2com,'-');
  xlim([640 1100]);
  xlabel('wavenumber cm^{-1}');ylabel('BT K');title('2012.04 JPL IC SNO 3840 pairs');
 legend('CrIS Obs mean','IASI Obs mean','Sarta.clear.mean','I2C mean','Location','South');
%}

