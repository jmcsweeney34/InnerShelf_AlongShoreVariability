% Written by Jack McSweeney 
% June 6, 2019 

% Trying to show definitively that there are two tides - work for
% manuscript describing the internal tides 


close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

%mooring data 
moordir='/Volumes/InnerShelf1/JackProcessing/Level2/';
veldir='/Volumes/InnerShelf1/JackProcessing/ADCPS/';
cleandatadir='/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/';

% coastline
load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalcoast.mat');

%bathy
% bathy=load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalBathy.mat');
bathy=load('/Volumes/InnerShelf1/Moorings/gathered_grids.mat');

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/JPOmanuscript/CrossSections/';

% Mooring locations 
load('/Volumes/InnerShelf1/Moorings/MoorLocs_AlongshoreVar.mat');

%% Lat Long to XY Coordinates  
[coast.x, coast.y]=lltoxy_PtSal(coast.lat,coast.lon);
coast.patchx=[10*1000;  coast.x;  10*1000]; 
coast.patchy=[-31259.1251413617; coast.y; 38034.4551997054]; 
[bathy.x, bathy.y]=lltoxy_PtSal(bathy.G.g200m.lat,bathy.G.g200m.lon);
[x, y]=lltoxy_PtSal(lat,lon);


%% Load Data 
% These are all extrapolated, cleaned up products with 1 min resolution 

cleandatadir='/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/cleanup/';

MS100=load([cleandatadir 'MS100.mat']); MS100.dn=MS100.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/ms100_alpha_beta_Cw.mat','alpha','Cw');
MS100.alpha=tempo.alpha; MS100.Cw=tempo.Cw;  clear tempo;

OC50=load([cleandatadir 'OC50.mat']); OC50.dn=OC50.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc50_alpha_beta_Cw.mat','alpha','Cw');
OC50.alpha=tempo.alpha; OC50.Cw=tempo.Cw; clear tempo;

OC40N=load([cleandatadir 'OC40N.mat']); OC40N.dn=OC40N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc40n_alpha_beta_Cw.mat','alpha','Cw');
OC40N.alpha=tempo.alpha; OC40N.Cw=tempo.Cw; clear tempo;

OC40S=load([cleandatadir 'OC40S.mat']); OC40S.dn=OC40S.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc40s_alpha_beta_Cw.mat','alpha','Cw');
OC40S.alpha=tempo.alpha; OC40S.Cw=tempo.Cw; clear tempo;

OC32N=load([cleandatadir 'OC32N.mat']); OC32N.dn=OC32N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc32n_alpha_beta_Cw.mat','alpha','Cw');
OC32N.alpha=tempo.alpha; OC32N.Cw=tempo.Cw; clear tempo;

OC32S=load([cleandatadir 'OC32S.mat']); OC32S.dn=OC32S.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc32s_alpha_beta_Cw.mat','alpha','Cw');
OC32S.alpha=tempo.alpha; OC32S.Cw=tempo.Cw;clear tempo;

OC25NA=load([cleandatadir 'OC25NA.mat']); OC25NA.dn=OC25NA.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25na_alpha_beta_Cw.mat','alpha','Cw');
OC25NA.alpha=tempo.alpha; OC25NA.Cw=tempo.Cw; clear tempo;

OC25NB=load([cleandatadir 'OC25NB.mat']); OC25NB.dn=OC25NB.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25nb_alpha_beta_Cw.mat','alpha','Cw');
OC25NB.alpha=tempo.alpha; OC25NB.Cw=tempo.Cw; clear tempo;

OC25M=load([cleandatadir 'OC25M.mat']); OC25M.dn=OC25M.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25m_alpha_beta_Cw.mat','alpha','Cw');
OC25M.alpha=tempo.alpha; OC25M.Cw=tempo.Cw; clear tempo;

OC25SB=load([cleandatadir 'OC25SB.mat']); OC25SB.dn=OC25SB.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25sb_alpha_beta_Cw.mat','alpha','Cw');
OC25SB.alpha=tempo.alpha; OC25SB.Cw=tempo.Cw; clear tempo;

OC25SA=load([cleandatadir 'OC25SA.mat']); OC25SA.dn=OC25SA.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25sa_alpha_beta_Cw.mat','alpha','Cw');
OC25SA.alpha=tempo.alpha; OC25SA.Cw=tempo.Cw; clear tempo;

OC17N=load([cleandatadir 'OC17N.mat']); OC17N.dn=OC17N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc17n_alpha_beta_Cw.mat','alpha','Cw');
OC17N.alpha=tempo.alpha; OC17N.Cw=tempo.Cw;  clear tempo;

OC17S=load([cleandatadir 'OC17S.mat']); OC17S.dn=OC17S.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc17s_alpha_beta_Cw.mat','alpha','Cw');
OC17S.alpha=tempo.alpha; OC17S.Cw=tempo.Cw; clear tempo;

OC10N=load([cleandatadir 'OC10N.mat']); OC10N.dn=OC10N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc10n_alpha_beta_Cw.mat','alpha','Cw');
OC10N.alpha=tempo.alpha; OC10N.Cw=tempo.Cw; clear tempo;

STR3B=load([cleandatadir 'STR3B.mat']); STR3B.dn=STR3B.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/str3b_alpha_beta_Cw.mat','alpha','Cw');
STR3B.alpha=tempo.alpha; STR3B.Cw=tempo.Cw; clear tempo;

NRL50N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL50N.mat');
PS50=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS50.mat');
VB50N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB50N.mat');
VB50S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB50S.mat');
NRL50S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL50S.mat');
PS40N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS40N.mat');
PS40M=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS40M.mat');
PS40S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS40S.mat');
NRL35N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL35N.mat');
PS35M=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS35M.mat');
NRL35S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL35S.mat');
PS30N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS30N.mat');
PS30M=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS30M.mat');
% NEED TO CLEAN THIS FILE
PS30S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS30S.mat');

VB30N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB30N.mat');
VB30S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB30S.mat');
VB25N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB25N.mat');
VB25S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB25S.mat');
NRL20N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL20N.mat');
NRL20S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL20S.mat');
NW2=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NW2.mat');
RW3=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/RW3.mat');
BW1=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/BW1.mat');
PW5=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PW5.mat');
STR4F=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/STR4F.mat');
STR5B=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/STR5B.mat');
STR6B=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/STR6B.mat');

%% Let's low pass filter the velocities 

filt=16; %33; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 

MS100.east_res=MS100.Usig-repmat(nanmean(MS100.Usig),length(MS100.sig),1);
MS100.east_sd16=MS100.east_res-pl66tn(MS100.east_res,dt,filt)';

OC50.east_res=OC50.Usig-repmat(nanmean(OC50.Usig),length(OC50.sig),1);
OC50.east_sd16=OC50.east_res-pl66tn(OC50.east_res,dt,filt)';

NRL50N.east_res=NRL50N.Usig-repmat(nanmean(NRL50N.Usig),size(NRL50N.Usig,1),1);
NRL50N.east_sd16=NRL50N.east_res-pl66tn(NRL50N.east_res,dt,filt)';

PS40M.east_res=PS40M.Usig-repmat(nanmean(PS40M.Usig),size(PS40M.Usig,1),1);
PS40M.east_sd16=PS40M.east_res-pl66tn(PS40M.east_res,dt,filt)';

PS50.east_res=PS50.Usig-repmat(nanmean(PS50.Usig),size(PS50.Usig,1),1);
PS50.east_sd16=PS50.east_res-pl66tn(PS50.east_res,dt,filt)';

VB50N.east_res=VB50N.Usig-repmat(nanmean(VB50N.Usig),size(VB50N.Usig,1),1);
VB50N.east_sd16=VB50N.east_res-pl66tn(VB50N.east_res,dt,filt)';

VB50S.east_res=VB50S.Usig-repmat(nanmean(VB50S.Usig),size(VB50S.Usig,1),1);
VB50S.east_sd16=VB50S.east_res-pl66tn(VB50S.east_res,dt,filt)';

NRL50S.east_res=NRL50S.Usig-repmat(nanmean(NRL50S.Usig),size(NRL50S.Usig,1),1);
NRL50S.east_sd16=NRL50S.east_res-pl66tn(NRL50S.east_res,dt,filt)'; 

OC25NA.east_res=OC25NA.Usig-repmat(nanmean(OC25NA.Usig),length(OC25NA.sig),1);
OC25NA.east_sd16=OC25NA.east_res-pl66tn(OC25NA.east_res,dt,filt)';

OC25NB.east_res=OC25NB.Usig-repmat(nanmean(OC25NB.Usig),length(OC25NB.sig),1);
OC25NB.east_sd16=OC25NB.east_res-pl66tn(OC25NB.east_res,dt,filt)';

OC25M.east_res=OC25M.Usig-repmat(nanmean(OC25M.Usig),length(OC25M.sig),1);
OC25M.east_sd16=OC25M.east_res-pl66tn(OC25M.east_res,dt,filt)';

OC25SB.east_res=OC25SB.Usig-repmat(nanmean(OC25SB.Usig),length(OC25SB.sig),1);
OC25SB.east_sd16=OC25SB.east_res-pl66tn(OC25SB.east_res,dt,filt)';

OC25SA.east_res=OC25SA.Usig-repmat(nanmean(OC25SA.Usig),length(OC25SA.sig),1);
OC25SA.east_sd16=OC25SA.east_res-pl66tn(OC25SA.east_res,dt,filt)';

OC10N.east_res=OC10N.Usig-repmat(nanmean(OC10N.Usig),length(OC10N.sig),1);
OC10N.east_sd16=OC10N.east_res-pl66tn(OC10N.east_res,dt,filt)';

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL20N.sig = (0:dsig:1)' ; % sigma 
NRL20N.east_res=NRL20N.Usig-repmat(nanmean(NRL20N.Usig),length(NRL20N.sig),1);
NRL20N.east_sd16=NRL20N.east_res-pl66tn(NRL20N.east_res,dt,filt)';

PS30M.sig = (0:dsig:1)' ; % sigma 
PS30M.east_res=PS30M.Usig-repmat(nanmean(PS30M.Usig),length(PS30M.sig),1);
PS30M.east_sd16=PS30M.east_res-pl66tn(PS30M.east_res,dt,filt)';

PS30S.sig = (0:dsig:1)' ; % sigma 
PS30S.east_res=PS30S.Usig-repmat(nanmean(PS30S.Usig),length(PS30S.sig),1);
PS30S.east_sd16=PS30S.east_res-pl66tn(PS30S.east_res,dt,filt)';

VB25N.sig = (0:dsig:1)' ; % sigma 
VB25N.east_res=VB25N.Usig-repmat(nanmean(VB25N.Usig),length(VB25N.sig),1);
VB25N.east_sd16=VB25N.east_res-pl66tn(VB25N.east_res,dt,filt)';

VB30N.sig = (0:dsig:1)' ; % sigma 
VB30N.east_res=VB30N.Usig-repmat(nanmean(VB30N.Usig),length(VB30N.sig),1);
VB30N.east_sd16=VB30N.east_res-pl66tn(VB30N.east_res,dt,filt)';

VB30S.sig = (0:dsig:1)' ; % sigma 
VB30S.east_res=VB30S.Usig-repmat(nanmean(VB30S.Usig),length(VB30S.sig),1);
VB30S.east_sd16=VB30S.east_res-pl66tn(VB30S.east_res,dt,filt)';

NRL20S.sig = (0:dsig:1)' ; % sigma 
NRL20S.east_res=NRL20S.Usig-repmat(nanmean(NRL20S.Usig),length(NRL20S.sig),1);
NRL20S.east_sd16=NRL20S.east_res-pl66tn(NRL20S.east_res,dt,filt)';

NW2.east_res=NW2.Usig-repmat(nanmean(NW2.Usig),length(NW2.sig),1);
NW2.east_sd16=NW2.east_res-pl66tn(NW2.east_res,dt,filt)';

BW1.east_res=BW1.Usig-repmat(nanmean(BW1.Usig),length(BW1.sig),1);
BW1.east_sd16=BW1.east_res-pl66tn(BW1.east_res,dt,filt)';

PW5.east_res=PW5.Usig-repmat(nanmean(PW5.Usig),length(PW5.sig),1);
PW5.east_sd16=PW5.east_res-pl66tn(PW5.east_res,dt,filt)';

STR5B.east_res=STR5B.Usig-repmat(nanmean(STR5B.Usig),length(STR5B.sig),1);
STR5B.east_sd16=STR5B.east_res-pl66tn(STR5B.east_res,dt,filt)';

STR6B.east_res=STR6B.Usig-repmat(nanmean(STR6B.Usig),length(STR6B.sig),1);
STR6B.east_sd16=STR6B.east_res-pl66tn(STR6B.east_res,dt,filt)';

%% Sept 11 - 50m 
close all

dn1=datenum(2017,9,11,9,0,0);
dn2=datenum(2017,9,12,12,0,0);
dt=2*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 

map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(6,1,[0.01 0.1],[0.08 0.05],[0.05 0.08]);


axes(hax(1))
i=findnearest_JACK(OC50.dn,dn1); j=findnearest_JACK(OC50.dn,dn2); 
pcolorjw(OC50.dn_sig(:,i:j),OC50.H(:,i:j),OC50.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.93   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(OC50.dn_sig(:,i:j),OC50.H(:,i:j),OC50.Tsig(:,i:j),[1:1:33],'k');
contour(OC50.dn_sig(:,i:j),OC50.H(:,i:j),OC50.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 52); xaxis(dn1,dn2);
title([ datestr(dn1,'mmm dd') ' - ' datestr(dn2,'dd')])

axes(hax(2))
i=findnearest_JACK(NRL50N.dn,dn1); j=findnearest_JACK(NRL50N.dn,dn2); 
pcolorjw(NRL50N.dn_sig(:,i:j),NRL50N.H(:,i:j),NRL50N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NRL50N.dn_sig(:,i:j),NRL50N.H(:,i:j),NRL50N.Tsig(:,i:j),[1:1:33],'k');
contour(NRL50N.dn_sig(:,i:j),NRL50N.H(:,i:j),NRL50N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 52); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(PS50.dn,dn1); j=findnearest_JACK(PS50.dn,dn2); 
pcolorjw(PS50.dn_sig(:,i:j),PS50.H(:,i:j),PS50.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(PS50.dn_sig(:,i:j),PS50.H(:,i:j),PS50.Tsig(:,i:j),[1:1:33],'k');
contour(PS50.dn_sig(:,i:j),PS50.H(:,i:j),PS50.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 52); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);

axes(hax(4))
i=findnearest_JACK(VB50N.dn,dn1); j=findnearest_JACK(VB50N.dn,dn2); 
pcolorjw(VB50N.dn_sig(:,i:j),VB50N.H(:,i:j),VB50N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB50N.dn_sig(:,i:j),VB50N.H(:,i:j),VB50N.Tsig(:,i:j),[1:1:33],'k');
contour(VB50N.dn_sig(:,i:j),VB50N.H(:,i:j),VB50N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 52); xaxis(dn1,dn2);

axes(hax(5))
i=findnearest_JACK(VB50S.dn,dn1); j=findnearest_JACK(VB50S.dn,dn2); 
pcolorjw(VB50S.dn_sig(:,i:j),VB50S.H(:,i:j),VB50S.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB50S.dn_sig(:,i:j),VB50S.H(:,i:j),VB50S.Tsig(:,i:j),[1:1:33],'k');
contour(VB50S.dn_sig(:,i:j),VB50S.H(:,i:j),VB50S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 52); xaxis(dn1,dn2);

axes(hax(6))
i=findnearest_JACK(NRL50S.dn,dn1); j=findnearest_JACK(NRL50S.dn,dn2); 
pcolorjw(NRL50S.dn_sig(:,i:j),NRL50S.H(:,i:j),NRL50S.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NRL50S.dn_sig(:,i:j),NRL50S.H(:,i:j),NRL50S.Tsig(:,i:j),[1:1:33],'k');
contour(NRL50S.dn_sig(:,i:j),NRL50S.H(:,i:j),NRL50S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on');
yaxis(1, 52); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir '50m' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%% Sept 11 - 25m 
close all

dn1=datenum(2017,9,11,18,0,0);
dn2=datenum(2017,9,12,17,0,0);
dt=2*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(5,1,[0.01 0.1],[0.08 0.05],[0.05 0.08]);


axes(hax(1))
i=findnearest_JACK(OC25NA.dn,dn1); j=findnearest_JACK(OC25NA.dn,dn2); 
pcolorjw(OC25NA.dn_sig(:,i:j),OC25NA.H(:,i:j),OC25NA.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(OC25NA.dn_sig(:,i:j),OC25NA.H(:,i:j),OC25NA.Tsig(:,i:j),[1:1:33],'k');
contour(OC25NA.dn_sig(:,i:j),OC25NA.H(:,i:j),OC25NA.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 26); xaxis(dn1,dn2);
% title([ datestr(dn1,'mmm dd') ' - ' datestr(dn2,'dd')])

axes(hax(2))
i=findnearest_JACK(OC25NB.dn,dn1); j=findnearest_JACK(OC25NB.dn,dn2); 
pcolorjw(OC25NB.dn_sig(:,i:j),OC25NB.H(:,i:j),OC25NB.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC25NB.dn_sig(:,i:j),OC25NB.H(:,i:j),OC25NB.Tsig(:,i:j),[1:1:33],'k');
contour(OC25NB.dn_sig(:,i:j),OC25NB.H(:,i:j),OC25NB.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 26); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(OC25M.dn,dn1); j=findnearest_JACK(OC25M.dn,dn2); 
pcolorjw(OC25M.dn_sig(:,i:j),OC25M.H(:,i:j),OC25M.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC25M.dn_sig(:,i:j),OC25M.H(:,i:j),OC25M.Tsig(:,i:j),[1:1:33],'k');
contour(OC25M.dn_sig(:,i:j),OC25M.H(:,i:j),OC25M.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 26); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);

axes(hax(4))
i=findnearest_JACK(OC25SB.dn,dn1); j=findnearest_JACK(OC25SB.dn,dn2); 
pcolorjw(OC25SB.dn_sig(:,i:j),OC25SB.H(:,i:j),OC25SB.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC25SB.dn_sig(:,i:j),OC25SB.H(:,i:j),OC25SB.Tsig(:,i:j),[1:1:33],'k');
contour(OC25SB.dn_sig(:,i:j),OC25SB.H(:,i:j),OC25SB.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on');
yaxis(1, 26); xaxis(dn1,dn2);

axes(hax(5))
i=findnearest_JACK(OC25SA.dn,dn1); j=findnearest_JACK(OC25SA.dn,dn2); 
pcolorjw(OC25SA.dn_sig(:,i:j),OC25SA.H(:,i:j),OC25SA.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC25SA.dn_sig(:,i:j),OC25SA.H(:,i:j),OC25SA.Tsig(:,i:j),[1:1:33],'k');
contour(OC25SA.dn_sig(:,i:j),OC25SA.H(:,i:j),OC25SA.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on');
yaxis(1, 26); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'OC25' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%%  Cross Section at Northern Most 
close all

dn1=datenum(2017,9,11,6,0,0);
dn2=datenum(2017,9,12,19,0,0);
dt=4*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(7,1,[0.02 0.1],[0.07 0.05],[0.05 0.09]);


axes(hax(1))
i=findnearest_JACK(MS100.dn,dn1); j=findnearest_JACK(MS100.dn,dn2-14/24); 
pcolorjw(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[1:1:33],'k');
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',25:25:100);
yaxis(1, 102); xaxis(dn1,dn2);
title([ datestr(dn1,'dd') ' - ' datestr(dn2,'dd') ' September 2017'])

axes(hax(2))
i=findnearest_JACK(OC50.dn,dn1+6/24); j=findnearest_JACK(OC50.dn,dn2-6/24); 
pcolorjw(OC50.dn_sig(:,i:j),OC50.H(:,i:j),OC50.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC50.dn_sig(:,i:j),OC50.H(:,i:j),OC50.Tsig(:,i:j),[1:1:33],'k');
contour(OC50.dn_sig(:,i:j),OC50.H(:,i:j),OC50.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',15:15:50);
yaxis(1, 51); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(OC40N.dn,dn1+8/24); j=findnearest_JACK(OC40N.dn,dn2-4.5/24); 
hold on 
contour(OC40N.dn_sig(:,i:j),OC40N.H(:,i:j),OC40N.Tsig(:,i:j),[1:1:33],'k');
contour(OC40N.dn_sig(:,i:j),OC40N.H(:,i:j),OC40N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:40,'yticklabel',10:10:40);
yaxis(1, 41); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);

axes(hax(4))
i=findnearest_JACK(OC32N.dn,dn1+9.5/25); j=findnearest_JACK(OC32N.dn,dn2-3/24); 
hold on 
contour(OC32N.dn_sig(:,i:j),OC32N.H(:,i:j),OC32N.Tsig(:,i:j),[1:1:33],'k');
contour(OC32N.dn_sig(:,i:j),OC32N.H(:,i:j),OC32N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:30,'yticklabel',10:10:30);
yaxis(1, 34); xaxis(dn1,dn2);

axes(hax(5))
i=findnearest_JACK(OC25NB.dn,dn1+12/24); j=findnearest_JACK(OC25NB.dn,dn2-1.5/24); 
pcolorjw(OC25NB.dn_sig(:,i:j),OC25NB.H(:,i:j),OC25NB.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC25NB.dn_sig(:,i:j),OC25NB.H(:,i:j),OC25NB.Tsig(:,i:j),[1:1:33],'k');
contour(OC25NB.dn_sig(:,i:j),OC25NB.H(:,i:j),OC25NB.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:10:25,'yticklabel',5:10:25);
yaxis(1, 26); xaxis(dn1,dn2);

axes(hax(6))
i=findnearest_JACK(OC17N.dn,dn1+13/24); j=findnearest_JACK(OC17N.dn,dn2); 
hold on 
contour(OC17N.dn_sig(:,i:j),OC17N.H(:,i:j),OC17N.Tsig(:,i:j),[1:1:33],'k');
contour(OC17N.dn_sig(:,i:j),OC17N.H(:,i:j),OC17N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:25,'yticklabel',5:5:25);
yaxis(1, 18); xaxis(dn1,dn2);

axes(hax(7))
i=findnearest_JACK(OC10N.dn,dn1+14/24); j=findnearest_JACK(OC10N.dn,dn2); 
pcolorjw(OC10N.dn_sig(:,i:j),OC10N.H(:,i:j),OC10N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(OC10N.dn_sig(:,i:j),OC10N.H(:,i:j),OC10N.Tsig(:,i:j),[1:1:33],'k');
contour(OC10N.dn_sig(:,i:j),OC10N.H(:,i:j),OC10N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on','ytick',2:2:10,'yticklabel',2:2:10);
yaxis(1, 10); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CrossSection_North' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%%  Cross Section2, North of Pt Sal  
close all

dn1=datenum(2017,9,11,6,0,0);
dn2=datenum(2017,9,12,19,0,0);
dt=4*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(5,1,[0.02 0.1],[0.07 0.05],[0.05 0.09]);


axes(hax(1))
i=findnearest_JACK(MS100.dn,dn1); j=findnearest_JACK(MS100.dn,dn2-14/24); 
pcolorjw(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[1:1:33],'k');
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',25:25:100);
yaxis(1, 102); xaxis(dn1,dn2);
title([ datestr(dn1,'dd') ' - ' datestr(dn2,'dd') ' September 2017'])

axes(hax(2))
i=findnearest_JACK(NRL50N.dn,dn1+6/24); j=findnearest_JACK(NRL50N.dn,dn2-6/24); 
pcolorjw(NRL50N.dn_sig(:,i:j),NRL50N.H(:,i:j),NRL50N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NRL50N.dn_sig(:,i:j),NRL50N.H(:,i:j),NRL50N.Tsig(:,i:j),[1:1:33],'k');
contour(NRL50N.dn_sig(:,i:j),NRL50N.H(:,i:j),NRL50N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',15:15:50);
yaxis(1, 51); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(NRL35N.dn,dn1+8/24); j=findnearest_JACK(NRL35N.dn,dn2-4.5/24); 
hold on 
contour(NRL35N.dn_sig(:,i:j),NRL35N.H(:,i:j),NRL35N.Tsig(:,i:j),[1:1:33],'k');
contour(NRL35N.dn_sig(:,i:j),NRL35N.H(:,i:j),NRL35N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:40,'yticklabel',5:5:40);
yaxis(1, 36); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);

axes(hax(4))
i=findnearest_JACK(NRL20N.dn,dn1+9.5/25); j=findnearest_JACK(NRL20N.dn,dn2-3/24); 
pcolorjw(NRL20N.dn_sig(:,i:j),NRL20N.H(:,i:j),NRL20N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NRL20N.dn_sig(:,i:j),NRL20N.H(:,i:j),NRL20N.Tsig(:,i:j),[1:1:33],'k');
contour(NRL20N.dn_sig(:,i:j),NRL20N.H(:,i:j),NRL20N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:30,'yticklabel',5:5:30);
yaxis(1, 21); xaxis(dn1,dn2);

axes(hax(5))
i=findnearest_JACK(NW2.dn,dn1+14/24); j=findnearest_JACK(NW2.dn,dn2); 
pcolorjw(NW2.dn_sig(:,i:j),NW2.H(:,i:j),NW2.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NW2.dn_sig(:,i:j),NW2.H(:,i:j),NW2.Tsig(:,i:j),[1:1:33],'k');
contour(NW2.dn_sig(:,i:j),NW2.H(:,i:j),NW2.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on','ytick',2:2:12,'yticklabel',2:2:12);
yaxis(1, 13); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CrossSection2_NPtSal' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%%  Cross Section 3,  Pt Sal  
close all

dn1=datenum(2017,9,11,6,0,0);
dn2=datenum(2017,9,12,19,0,0);
dt=4*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(6,1,[0.02 0.1],[0.07 0.05],[0.05 0.09]);


axes(hax(1))
i=findnearest_JACK(MS100.dn,dn1); j=findnearest_JACK(MS100.dn,dn2-14/24); 
pcolorjw(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[1:1:33],'k');
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',25:25:100);
yaxis(1, 102); xaxis(dn1,dn2);
title([ datestr(dn1,'dd') ' - ' datestr(dn2,'dd') ' September 2017'])

axes(hax(2))
i=findnearest_JACK(PS50.dn,dn1+6/24); j=findnearest_JACK(PS50.dn,dn2-9/24); 
pcolorjw(PS50.dn_sig(:,i:j),PS50.H(:,i:j),PS50.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(PS50.dn_sig(:,i:j),PS50.H(:,i:j),PS50.Tsig(:,i:j),[1:1:33],'k');
contour(PS50.dn_sig(:,i:j),PS50.H(:,i:j),PS50.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',15:15:50);
yaxis(1, 51); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(PS40M.dn,dn1+6/24); j=findnearest_JACK(PS40M.dn,dn2-6/24); 
pcolorjw(PS40M.dn_sig(:,i:j),PS40M.H(:,i:j),PS40M.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(PS40M.dn_sig(:,i:j),PS40M.H(:,i:j),PS40M.Tsig(:,i:j),[1:1:33],'k');
contour(PS40M.dn_sig(:,i:j),PS40M.H(:,i:j),PS40M.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:40);
yaxis(1, 41); xaxis(dn1,dn2);

axes(hax(4))
i=findnearest_JACK(PS35M.dn,dn1+8/24); j=findnearest_JACK(PS35M.dn,dn2-5/24); 
hold on 
contour(PS35M.dn_sig(:,i:j),PS35M.H(:,i:j),PS35M.Tsig(:,i:j),[1:1:33],'k');
contour(PS35M.dn_sig(:,i:j),PS35M.H(:,i:j),PS35M.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:40,'yticklabel',5:5:40);
yaxis(1, 36); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);

axes(hax(5))
i=findnearest_JACK(PS30M.dn,dn1+9.5/25); j=findnearest_JACK(PS30M.dn,dn2-3/24); 
pcolorjw(PS30M.dn_sig(:,i:j),PS30M.H(:,i:j),PS30M.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(PS30M.dn_sig(:,i:j),PS30M.H(:,i:j),PS30M.Tsig(:,i:j),[1:1:33],'k');
contour(PS30M.dn_sig(:,i:j),PS30M.H(:,i:j),PS30M.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:30,'yticklabel',10:10:30);
yaxis(1, 31); xaxis(dn1,dn2);

axes(hax(6))
i=findnearest_JACK(BW1.dn,dn1+10/24); j=findnearest_JACK(BW1.dn,dn2); 
pcolorjw(BW1.dn_sig(:,i:j),BW1.H(:,i:j),BW1.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(BW1.dn_sig(:,i:j),BW1.H(:,i:j),BW1.Tsig(:,i:j),[1:1:33],'k');
contour(BW1.dn_sig(:,i:j),BW1.H(:,i:j),BW1.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on','ytick',2:2:12,'yticklabel',2:2:12);
yaxis(1, 13); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CrossSection3_PtSal' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%%  Cross Section 4,  South Pt Sal  
close all

dn1=datenum(2017,9,11,6,0,0);
dn2=datenum(2017,9,12,19,0,0);
dt=4*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(5,1,[0.02 0.1],[0.07 0.05],[0.05 0.09]);


axes(hax(1))
i=findnearest_JACK(MS100.dn,dn1); j=findnearest_JACK(MS100.dn,dn2-14/24); 
pcolorjw(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[1:1:33],'k');
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',25:25:100);
yaxis(1, 102); xaxis(dn1,dn2);
title([ datestr(dn1,'dd') ' - ' datestr(dn2,'dd') ' September 2017'])

axes(hax(2))
i=findnearest_JACK(VB50N.dn,dn1+3/24); j=findnearest_JACK(VB50N.dn,dn2-9/24); 
pcolorjw(VB50N.dn_sig(:,i:j),VB50N.H(:,i:j),VB50N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB50N.dn_sig(:,i:j),VB50N.H(:,i:j),VB50N.Tsig(:,i:j),[1:1:33],'k');
contour(VB50N.dn_sig(:,i:j),VB50N.H(:,i:j),VB50N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',15:15:50);
yaxis(1, 52); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(VB30N.dn,dn1+6/24); j=findnearest_JACK(VB30N.dn,dn2-6/24); 
pcolorjw(VB30N.dn_sig(:,i:j),VB30N.H(:,i:j),VB30N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB30N.dn_sig(:,i:j),VB30N.H(:,i:j),VB30N.Tsig(:,i:j),[1:1:33],'k');
contour(VB30N.dn_sig(:,i:j),VB30N.H(:,i:j),VB30N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:40);
yaxis(1, 31); xaxis(dn1,dn2);

axes(hax(4))
i=findnearest_JACK(VB25N.dn,dn1+8/24); j=findnearest_JACK(VB25N.dn,dn2-4.5/24); 
pcolorjw(VB25N.dn_sig(:,i:j),VB25N.H(:,i:j),VB25N.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB25N.dn_sig(:,i:j),VB25N.H(:,i:j),VB25N.Tsig(:,i:j),[1:1:33],'k');
contour(VB25N.dn_sig(:,i:j),VB25N.H(:,i:j),VB25N.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:40,'yticklabel',5:5:40);
yaxis(1, 26); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);


axes(hax(5))
i=findnearest_JACK(STR4F.dn,dn1+10/24); j=findnearest_JACK(STR4F.dn,dn2); 
hold on 
contour(STR4F.dn_sig(:,i:j),STR4F.H(:,i:j),STR4F.Tsig(:,i:j),[1:0.5:33],'r');
contour(STR4F.dn_sig(:,i:j),STR4F.H(:,i:j),STR4F.Tsig(:,i:j),[1:1:33],'k');
contour(STR4F.dn_sig(:,i:j),STR4F.H(:,i:j),STR4F.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on','ytick',2:2:12,'yticklabel',2:2:12);
yaxis(1, 10.75); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CrossSection4_SPtSal' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%%  Cross Section 5,  VB  
close all

dn1=datenum(2017,9,11,6,0,0);
dn2=datenum(2017,9,12,19,0,0);
dt=4*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(5,1,[0.02 0.1],[0.07 0.05],[0.05 0.09]);


axes(hax(1))
i=findnearest_JACK(MS100.dn,dn1); j=findnearest_JACK(MS100.dn,dn2-14/24); 
pcolorjw(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[1:1:33],'k');
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',25:25:100);
yaxis(1, 102); xaxis(dn1,dn2);
title([ datestr(dn1,'dd') ' - ' datestr(dn2,'dd') ' September 2017'])

axes(hax(2))
i=findnearest_JACK(VB50S.dn,dn1+3/24); j=findnearest_JACK(VB50S.dn,dn2-9/24); 
pcolorjw(VB50S.dn_sig(:,i:j),VB50S.H(:,i:j),VB50S.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB50S.dn_sig(:,i:j),VB50S.H(:,i:j),VB50S.Tsig(:,i:j),[1:1:33],'k');
contour(VB50S.dn_sig(:,i:j),VB50S.H(:,i:j),VB50S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',15:15:50);
yaxis(1, 52); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(VB30S.dn,dn1+6/24); j=findnearest_JACK(VB30S.dn,dn2-6/24); 
pcolorjw(VB30S.dn_sig(:,i:j),VB30S.H(:,i:j),VB30S.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(VB30S.dn_sig(:,i:j),VB30S.H(:,i:j),VB30S.Tsig(:,i:j),[1:1:33],'k');
contour(VB30S.dn_sig(:,i:j),VB30S.H(:,i:j),VB30S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:40);
yaxis(1, 31); xaxis(dn1,dn2);

axes(hax(4))
i=findnearest_JACK(VB25S.dn,dn1+8/24); j=findnearest_JACK(VB25S.dn,dn2-3/24); 
hold on 
contour(VB25S.dn_sig(:,i:j),VB25S.H(:,i:j),VB25S.Tsig(:,i:j),[1:1:33],'k');
contour(VB25S.dn_sig(:,i:j),VB25S.H(:,i:j),VB25S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:40,'yticklabel',5:5:40);
yaxis(1, 26); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);


axes(hax(5))
i=findnearest_JACK(STR5B.dn,dn1+10/24); j=findnearest_JACK(STR5B.dn,dn2); 
pcolorjw(STR5B.dn_sig(:,i:j),STR5B.H(:,i:j),STR5B.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(STR5B.dn_sig(:,i:j),STR5B.H(:,i:j),STR5B.Tsig(:,i:j),[1:1:33],'k');
contour(STR5B.dn_sig(:,i:j),STR5B.H(:,i:j),STR5B.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on','ytick',2:2:12,'yticklabel',2:2:12);
yaxis(1, 10.75); xaxis(dn1,dn2);


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CrossSection5_VB' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end

%%  Cross Section 6,  South  
close all

dn1=datenum(2017,9,11,6,0,0);
dn2=datenum(2017,9,12,19,0,0);
dt=4*1800/86400;
labels = datestr(dn1:dt:dn2,'HH:MM'); % extract
labels(1:2:end,:) = nan; % remove every other one

vmin=-0.25;
vmax=0.25;
isobold=15; 


map=flipud(cbrewer('div','RdBu',40));
figure('position',[20    55   747   750]);
hax=tight_subplot(5,1,[0.02 0.1],[0.07 0.05],[0.05 0.09]);


axes(hax(1))
i=findnearest_JACK(MS100.dn,dn1); j=findnearest_JACK(MS100.dn,dn2-14/24); 
pcolorjw(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
cc=colorbar;
set(cc,'position',[0.92964467750426   0.08  0.0229521012122679 0.87],...
    'ticks',[vmin:0.05:vmax]);
hold on 
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[1:1:33],'k');
contour(MS100.dn_sig(:,i:j),MS100.H(:,i:j),MS100.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',25:25:100);
yaxis(1, 102); xaxis(dn1,dn2);
title([ datestr(dn1,'dd') ' - ' datestr(dn2,'dd') ' September 2017'])

axes(hax(2))
i=findnearest_JACK(NRL50S.dn,dn1+3/24); j=findnearest_JACK(NRL50S.dn,dn2-9/24); 
pcolorjw(NRL50S.dn_sig(:,i:j),NRL50S.H(:,i:j),NRL50S.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NRL50S.dn_sig(:,i:j),NRL50S.H(:,i:j),NRL50S.Tsig(:,i:j),[1:1:33],'k');
contour(NRL50S.dn_sig(:,i:j),NRL50S.H(:,i:j),NRL50S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',15:15:50);
yaxis(1, 52); xaxis(dn1,dn2);

axes(hax(3))
i=findnearest_JACK(NRL35S.dn,dn1+6/24); j=findnearest_JACK(NRL35S.dn,dn2-6/24); 
hold on 
contour(NRL35S.dn_sig(:,i:j),NRL35S.H(:,i:j),NRL35S.Tsig(:,i:j),[1:1:33],'k');
contour(NRL35S.dn_sig(:,i:j),NRL35S.H(:,i:j),NRL35S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',10:10:40,'yticklabels',10:10:40);
yaxis(1, 36); xaxis(dn1,dn2);

axes(hax(4))
i=findnearest_JACK(NRL20S.dn,dn1+8/24); j=findnearest_JACK(NRL20S.dn,dn2-3/24); 
pcolorjw(NRL20S.dn_sig(:,i:j),NRL20S.H(:,i:j),NRL20S.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(NRL20S.dn_sig(:,i:j),NRL20S.H(:,i:j),NRL20S.Tsig(:,i:j),[1:1:33],'k');
contour(NRL20S.dn_sig(:,i:j),NRL20S.H(:,i:j),NRL20S.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',[],'layer','top','box','on','ytick',5:5:40,'yticklabel',5:5:40);
yaxis(1, 21); xaxis(dn1,dn2);
% t1=ylabel('Meters above Bed (m)');
% set(t1,'position',[736949.349448122          5.10490236747556          1.00000000000001]);


axes(hax(5))
i=findnearest_JACK(STR6B.dn,dn1+10/24); j=findnearest_JACK(STR6B.dn,dn2); 
pcolorjw(STR6B.dn_sig(:,i:j),STR6B.H(:,i:j),STR6B.east_sd16(:,i:j))
colormap(map)
caxis([vmin vmax])
hold on 
contour(STR6B.dn_sig(:,i:j),STR6B.H(:,i:j),STR6B.Tsig(:,i:j),[1:1:33],'k');
contour(STR6B.dn_sig(:,i:j),STR6B.H(:,i:j),STR6B.Tsig(:,i:j),[isobold isobold],'k','linewidth',2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on','ytick',2:2:12,'yticklabel',2:2:12);
yaxis(1, 10.75); xaxis(dn1,dn2);

if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CrossSection6_South' datestr(dn1,'mmmdd') '-' datestr(dn2,'dd')]);
end


