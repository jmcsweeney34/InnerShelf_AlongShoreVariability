close all 
clear all 

buoydir='/Volumes/InnerShelf1/NOAA_Buoy_Data/';

%% Download text data 

fname=[buoydir 'SantaMaria_46011/46011_2017_hourlyMETdata.rtf'];
% formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';

[M.yr, M.mo, M.dy, M.hr, M.mn, M.wdir, M.wspd, M.gst, M.wvht, M.dpd, ...
    M.apd, M.mwd, M.pres, M.atmp, M.wtmp, M.dewp, M.vis, M.tide]=...
    textread(fname,formatSpec,'commentstyle','shell');

press=str2double(M.pres(11:end));
yr=str2double(M.yr(11:end-1)); 
mo=str2double(M.mo(11:end));
dy=str2double(M.dy(11:end));
hr=str2double(M.hr(11:end));
mn=str2double(M.mn(11:end));
wvht=str2double(M.wvht(11:end));
wspd=str2double(M.wspd(11:end));
wdir=str2double(M.wdir(11:end));

dn_atm=datenum(yr,mo,dy,hr,mn,0);

wind_e=cos(deg2rad(wdir)).*(wspd);
wind_n=sin(deg2rad(wdir)).*(wspd);

clear press 
%% 

dn=datenum(2017,9,1):4/24:datenum(2017,11,5);
wnorth=interp1(dn_atm,pl33tn(wind_n,1,33),dn);
weast=interp1(dn_atm,pl33tn(wind_e,1,33),dn);

close all
figure('position',[         159         377        1045         394]);
subplot(3,1,1)
quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,24),pl33tn(wind_n,1,24),'autoscale','off')
subplot(3,1,2)
quiver(dn,dn.*0,weast,wnorth);
subplot(3,1,3)
quiver(dn,dn.*0,weast,wnorth,'autoscale','off');





