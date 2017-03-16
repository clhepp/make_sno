function [ds] = load_sno(type,cris_res)

% function [ds] = load_sno()
%
% INPUTS:
%     types are: airs and cris; iasi and cris; airs and iasi.
%     CrIS is either high-resolution or low-resolution. Options 'high', 'low'
% OUTPUTS
%     structure ds, with two components: 1. for sensor 1 (iasi or airs)
%                   2: for sensor 2 (cris or iasi)
%          and various fields:
%          lat
%          lon
%          time
%          fov
%          solzen
%          satzen
%          rad
% the following are only in ds(1)
%          tdiff
%          dist     great circle angle (deg)
%          

% Initialize
r_earth = 6371;        % radius of Earth (km)

% Determine type

all_types  = {'airs_cris','iasi_cris','airs_iasi'};
all_res    = {'HIGH','LOW'};
type       = lower(type);
cris_res   = upper(cris_res);

if(~ismember(type,all_types))
  error('Not a valid SNO type pair, options are: airs_cris, iasi_cris, airs_iasi')
  return
end
if(~ismember(cris_res, all_res))
  error('Not valid option for CrIS spectral resolution. Options are: HR or LR')
  return
end
if(strcmp(cris_res,'HIGH')) cres='HR'; end
if(strcmp(cris_res,'LOW'))  cres='LR'; end

ds(1) = struct('lat',[],'lon',[],'time',[],'solzen',[],'satzen',[],'fov',[],...
               'tdiff',[],'dist',[],'rad',single([]) );
ds(2) = struct('lat',[],'lon',[],'time',[],'solzen',[],'satzen',[],'fov',[],...
               'tdiff',[],'dist',[],'rad',single([]) );

% --------------- test and load up IASI CRIS SNOs -----------------------------
if(strcmp(type,'iasi_cris') & strcmp(cris_res,'HIGH'))
  srcdir = ['/asl/s1/chepplew/data/sno/' type '/' cris_res '/'];
  fnlst = dir(strcat(srcdir,'sno_iasi_crisHr_asl_2016*.mat'));
elseif(strcmp(type,'iasi_cris') & strcmp(cris_res,'LOW'))
  srcdir = ['/asl/s1/chepplew/data/sno/' type '/' cris_res '/'];
  fnlst = dir(strcat(srcdir,'sno_iasi_crisHr_asl_2016*.mat'));
end
if(strcmp(type,'iasi_cris'))
  for fn = 1:numel(fnlst)
    load(strcat(srcdir, fnlst(fn).name));
    ds(1).lat    = [ds(1).lat, ilat'];
    ds(1).lon    = [ds(1).lon, ilon'];
    ds(1).time   = [ds(1).time, itime'];
    ds(1).fov    = [ds(1).fov, ifov'];
    ds(1).solzen = [ds(1).solzen, isolzen'];
    ds(1).satzen = [ds(1).satzen, isatzen'];
    ds(1).rad    = [ds(1).rad; ri];
    ds(2).lat    = [ds(2).lat, clat'];
    ds(2).lon    = [ds(2).lon, clon'];
    ds(2).time   = [ds(2).time, ctime'];
    ds(2).fov    = [ds(2).fov, cfov'];
    ds(2).solzen = [ds(2).solzen, csolzen'];
    ds(2).satzen = [ds(2).satzen, csatzen'];
    ds(2).rad    = [ds(2).rad; rc];             %<- rc' to be fixed in SNO maker
   
    ds(1).tdiff  = [ds(1).tdiff, tdiff'];
    %dist         = r_earth*tan(ds(1).dist *3.14159/180;     convert frm angle to geometric dist (km)
    ds(1).dist   = [ds(1).dist, dist'];
    fprintf(1,'.');
  end
end

% --------------- test and load up AIRS CRIS SNOs -----------------------------
if(strcmp(type,'airs_cris') & strcmp(cris_res,'HIGH'))
  srcdir = ['/asl/s1/chepplew/data/sno/' type '/' cres '/2016/'];
  fnlst = dir(strcat(srcdir,'sno_airs_cris_asl_2016*.mat'));
elseif(strcmp(type,'airs_cris') & strcmp(cris_res,'LOW'))
  srcdir = ['/asl/s1/chepplew/data/sno/' type '/' cres '/2016/'];
  fnlst = dir(strcat(srcdir,'sno_airs_cris_asl_2016*.mat'));
end
if(strcmp(type,'airs_cris'))
  for fn = 1:numel(fnlst)
    load(strcat(srcdir, fnlst(fn).name));
    ds(1).lat    = [ds(1).lat, alat];
    ds(1).lon    = [ds(1).lon, alon];
    ds(1).time   = [ds(1).time, atime];
    ds(1).fov    = [];
    ds(1).solzen = [ds(1).solzen, asolzen];
    ds(1).satzen = [ds(1).satzen, asatzen];
    ds(1).rad    = [ds(1).rad; ra];
    ds(2).lat    = [ds(2).lat, clat'];
    ds(2).lon    = [ds(2).lon, clon'];
    ds(2).time   = [ds(2).time, ctime'];
    ds(2).fov    = [ds(2).fov, cfov'];
    ds(2).solzen = [ds(2).solzen, csolzen'];
    ds(2).satzen = [ds(2).satzen, csatzen'];
    ds(2).rad    = [ds(2).rad; rc];             %<- rc' to be fixed in SNO maker
   
    ds(1).tdiff  = [ds(1).tdiff, tdiff'];
    %dist         = r_earth*tan(ds(1).dist *3.14159/180;     convert frm angle to geometric dist (km)
    ds(1).dist   = [ds(1).dist, dist'];
    fprintf(1,'.');
  end
end
