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
anadir= '/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/JPOmanuscript/StratifcationMaps/';

% Mooring locations 
load('/Volumes/InnerShelf1/Moorings/MoorLocs_AlongshoreVar.mat');

%% Lat Long to XY Coordinates  
[coast.x, coast.y]=lltoxy_PtSal(coast.lat,coast.lon);
coast.patchx=[10*1000;  coast.x;  10*1000]; 
coast.patchy=[-31259.1251413617; coast.y; 38034.4551997054]; 
[bathy.x, bathy.y]=lltoxy_PtSal(bathy.G.g200m.lat,bathy.G.g200m.lon);
[x, y]=lltoxy_PtSal(lat,lon);


%% Download Temp and Velocity Data - 
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


%% Compare timeseries of alpha at 50 isobaths 
close all 

dn1=datenum(2017,9,8,0,0,0);
dn2=datenum(2017,11,3,0,0,0);

dt=2;
labels = datestr(dn1:dt:dn2,'mm/dd'); % extract
labels(2:2:end,:) = nan; % remove every other one

figure('position',[46         442        1323         334]);
hold on 
plot([dn1 dn2],[0 0],'--k','linewidth',0.5);
h1=plot(OC50.DNrho,OC50.alpha,'color','k','linewidth',2);
h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.alpha,'color','r','linewidth',2);
h3=plot(PS50.subtidal.DNrho,PS50.subtidal.alpha,'color','g','linewidth',2);
h4=plot(VB50N.subtidal.DNrho,VB50N.subtidal.alpha,'color','b','linewidth',2);
h5=plot(VB50S.subtidal.DNrho,VB50S.subtidal.alpha,'color','m','linewidth',2);
h6=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.alpha,'color','c','linewidth',2);

legend([h1 h2 h3 h4 h5 h6],names(2:7));
ylabel('\alpha');

xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on');


if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'alphaAt50m']);
end

%% Compare timeseries of alpha/co at 50 isobaths 
close all 

dn1=datenum(2017,9,8,0,0,0);
dn2=datenum(2017,11,3,0,0,0);

dt=2;
labels = datestr(dn1:dt:dn2,'mm/dd'); % extract
labels(2:2:end,:) = nan; % remove every other one

figure('position',[46         442        1323         334]);
hold on 
plot([dn1 dn2],[0 0],'--k','linewidth',0.5);
h1=plot(OC50.DNrho,OC50.alpha./OC50.Cw,'color','k','linewidth',2);
h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.alpha./NRL50N.subtidal.Cw,'color','r','linewidth',2);
h3=plot(PS50.subtidal.DNrho,PS50.subtidal.alpha./PS50.subtidal.Cw,'color','g','linewidth',2);
h4=plot(VB50N.subtidal.DNrho,VB50N.subtidal.alpha./VB50N.subtidal.Cw,'color','b','linewidth',2);
h5=plot(VB50S.subtidal.DNrho,VB50S.subtidal.alpha./VB50S.subtidal.Cw,'color','m','linewidth',2);
h6=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.alpha./NRL50S.subtidal.Cw,'color','c','linewidth',2);

legend([h1 h2 h3 h4 h5 h6],names(2:7));
ylabel('\alpha/c_o');

xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on');

if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'alpha_DivByCoAt50mB']);
end
    

%% Compare timeseries of alpha at 50 isobaths 
close all 

dn1=datenum(2017,9,8,0,0,0);
dn2=datenum(2017,11,3,0,0,0);

dt=2;
labels = datestr(dn1:dt:dn2,'mm/dd'); % extract
labels(2:2:end,:) = nan; % remove every other one

figure('position',[46         442        1323         334]);
hold on 
% plot([dn1 dn2],[0 0],'--k','linewidth',0.5);
h1=plot(OC50.DNrho,OC50.Cw,'color','k','linewidth',2);
h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.Cw,'color','r','linewidth',2);
h3=plot(PS50.subtidal.DNrho,PS50.subtidal.Cw,'color','g','linewidth',2);
h4=plot(VB50N.subtidal.DNrho,VB50N.subtidal.Cw,'color','b','linewidth',2);
h5=plot(VB50S.subtidal.DNrho,VB50S.subtidal.Cw,'color','m','linewidth',2);
h6=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.Cw,'color','c','linewidth',2);

legend([h1 h2 h3 h4 h5 h6],names(2:7));
ylabel('c_o');

xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',14,'xtick',dn1:dt:dn2,...
    'xticklabel',labels,'layer','top','box','on');

if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'CoAt50mB']);
end
    
%% Plot Alpha Map

cutoff=0.0416666667442769; % must be within 1 hr 

    dt=6;
    
    dn1=datenum(2017,9,11,12,0,0);
    alpha=lon.*nan;
    
    
    ii=findnearest_JACK(MS100.DNrho,dn1); jj=find(ismember(names,'MS100'));
    if abs(MS100.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=MS100.alpha(ii);
    end
    
    ii=findnearest_JACK(OC50.DNrho,dn1); jj=find(ismember(names,'OC50'));
    if abs(OC50.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=OC50.alpha(ii);
    end
    
    ii=findnearest_JACK(OC40N.DNrho,dn1); jj=find(ismember(names,'OC40N'));
    if abs(OC40N.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=OC40N.alpha(ii);
    end
    
    ii=findnearest_JACK(OC40S.DNrho,dn1); jj=find(ismember(names,'OC40S'));
    if abs(OC40S.DNrho(ii)-dn1)< cutoff
        alpha(jj)=OC40S.alpha(ii);
    end
    
    ii=findnearest_JACK(OC32N.DNrho,dn1); jj=find(ismember(names,'OC32N'));
    if abs(OC32N.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=OC32N.alpha(ii);
    end
    
    ii=findnearest_JACK(OC32S.DNrho,dn1); jj=find(ismember(names,'OC32S'));
    if abs(OC32S.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=OC32S.alpha(ii);
    end
    
    ii=findnearest_JACK(OC25NA.DNrho,dn1); jj=find(ismember(names,'OC25NA'));
    if abs(OC25NA.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=OC25NA.alpha(ii);
    end
    
    ii=findnearest_JACK(OC25NB.DNrho,dn1); jj=find(ismember(names,'OC25NB'));
    if abs(OC25NB.DNrho(ii)-dn1)< cutoff 
        alpha(jj)=OC25NB.alpha(ii);
    end
    
    ii=findnearest_JACK(OC25M.DNrho,dn1); jj=find(ismember(names,'OC25M'));
    if abs(OC25M.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=OC25M.alpha(ii); 
    end
    
    ii=findnearest_JACK(OC25SB.DNrho,dn1); jj=find(ismember(names,'OC25SB'));
    if abs(OC25SB.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=OC25SB.alpha(ii);
    end
    
    ii=findnearest_JACK(OC25SA.DNrho,dn1); jj=find(ismember(names,'OC25SA'));
    if abs(OC25SA.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=OC25SA.alpha(ii);
    end
    
    ii=findnearest_JACK(OC17N.DNrho,dn1); jj=find(ismember(names,'OC17N'));
    if abs(OC17N.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=OC17N.alpha(ii);
    end
    
    ii=findnearest_JACK(OC17S.DNrho,dn1); jj=find(ismember(names,'OC17S'));
    if abs(OC17S.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=OC17S.alpha(ii);
    end
    
    ii=findnearest_JACK(OC10N.DNrho,dn1); jj=find(ismember(names,'OC10N'));
    if abs(OC10N.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=OC10N.alpha(ii);
    end
    
    ii=findnearest_JACK(STR3B.DNrho,dn1); jj=find(ismember(names,'STR3B'));
    if abs(STR3B.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=STR3B.alpha(ii);
    end
    
    ii=findnearest_JACK(NRL50N.subtidal.DNrho,dn1); jj=find(ismember(names,'NRL50N'));
    if abs(NRL50N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=NRL50N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(NRL50S.subtidal.DNrho,dn1); jj=find(ismember(names,'NRL50S'));
    if abs(NRL50S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=NRL50S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS50.subtidal.DNrho,dn1); jj=find(ismember(names,'PS50'));
    if abs(PS50.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS50.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(VB50N.subtidal.DNrho,dn1); jj=find(ismember(names,'VB50N'));
    if abs(VB50N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=VB50N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(VB50S.subtidal.DNrho,dn1); jj=find(ismember(names,'VB50S'));
    if abs(VB50S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=VB50S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS40N.subtidal.DNrho,dn1); jj=find(ismember(names,'PS40N'));
    if abs(PS40N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS40N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS40M.subtidal.DNrho,dn1); jj=find(ismember(names,'PS40M'));
    if abs(PS40M.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS40M.subtidal.alpha(ii);
    end
    ii=findnearest_JACK(PS40S.subtidal.DNrho,dn1); jj=find(ismember(names,'PS40S'));
    if abs(PS40S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS40S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(NRL35N.subtidal.DNrho,dn1); jj=find(ismember(names,'NRL35N'));
    if abs(NRL35N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=NRL35N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS35M.subtidal.DNrho,dn1); jj=find(ismember(names,'PS35M'));
    if abs(PS35M.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS35M.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(NRL35S.subtidal.DNrho,dn1); jj=find(ismember(names,'NRL35S'));
    if abs(NRL35S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=NRL35S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS30N.subtidal.DNrho,dn1); jj=find(ismember(names,'PS30N'));
    if abs(PS30N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS30N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS30M.subtidal.DNrho,dn1); jj=find(ismember(names,'PS30M'));
    if abs(PS30M.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=PS30M.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(PS30S.subtidal.DNrho,dn1); jj=find(ismember(names,'PS30S'));
    if abs(PS30S.subtidal.DNrho(ii)-dn1)< cutoff
        alpha(jj)=PS30S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(VB30N.subtidal.DNrho,dn1); jj=find(ismember(names,'VB30N'));
    if abs(VB30N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=VB30N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(VB30S.subtidal.DNrho,dn1); jj=find(ismember(names,'VB30S'));
    if abs(VB30S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=VB30S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(VB25N.subtidal.DNrho,dn1); jj=find(ismember(names,'VB25N'));
    if abs(VB25N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=VB25N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(VB25S.subtidal.DNrho,dn1); jj=find(ismember(names,'VB25S'));
    if abs(VB25S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=VB25S.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(NRL20N.subtidal.DNrho,dn1); jj=find(ismember(names,'NRL20N'));
    if abs(NRL20N.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=NRL20N.subtidal.alpha(ii);
    end
    
    ii=findnearest_JACK(NRL20S.subtidal.DNrho,dn1); jj=find(ismember(names,'NRL20S'));
    if abs(NRL20S.subtidal.DNrho(ii)-dn1)< cutoff  
        alpha(jj)=NRL20S.subtidal.alpha(ii);
    end
    
    % Plot alpha at the various locations
    
    close all;
    load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/AlphaColormap.mat');
    map=flipud(map);
    
    figure('position',[440   100   479   698]);
    hold on
    patch(coast.patchx/1000,coast.patchy/1000,[1 1 1].*.8,'HandleVisibility','on');
    plot(coast.x/1000,coast.y/1000,'-','color',[1 1 1].*0.5,'linewidth',2);
    [c,h] = contour(bathy.x/1000,bathy.y/1000,-bathy.G.g200m.h,[-120 -100:10:10],...
        'color',[1 1 1].*0.8,'HandleVisibility','off','labelspacing',375);
    % plot(lon,lat,'.k','markersize',15);
    
    h1=scatter(x/1000,y/1000,1000,alpha,'.');
    cc=colorbar; colormap(map); caxis([-0.015, 0.015]);
    set(cc,'position',[0.709812108559493         0.535816618911175        0.0626304801670197         0.282234957020056]);
    cl=xlabel(cc,'\alpha');
    set(cl,'position',[-0.782666714986165      -0.00142130072346799                         0],...
        'fontweight','bold','fontsize',18);
    yaxis(-20,15);xaxis(-15,10);
    aspect_jack(gca,lat(1));
    title(datestr(dn1,'mmm dd HH:MM'));
    set(gca,'fontsize',14,'fontweight','bold');
    
    if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'alphaMap_' datestr(dn1,'mmmdd_HHMM')]);
    end
% end

