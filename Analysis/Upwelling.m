% Written by Jack McSweeney 
% July 23, 2019

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/JPOmanuscript/UpwellingAnalysis/';

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

%% Load tide data - port of San Luis 
tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');

%% Download wind data 
load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SantaMaria_46011/SantaMariaWind.mat');
dn_atm2=datenum(2017,9,1):2/24:datenum(2017,11,5);
wnorth=interp1(dn_atm,pl33tn(wind_n,1,36),dn_atm2);
weast=interp1(dn_atm,pl33tn(wind_e,1,36),dn_atm2);

% 36 hr filter because of milton 2009 
%% Plot Northward component of the winds

close all
figure('position',[440   458   700   340]); 
dn1=datenum(2017,9,6);
dn2=datenum(2017,11,2);
dt=7;

subplot(2,1,1);
quiver(dn_atm2,dn_atm2.*0,weast,wnorth,'k'); %,'autoscale','off');%,'maxheadsize',0)
% quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,33),pl33tn(wind_n,1,33),'k','autoscale','off','maxheadsize',0)
xaxis(dn1,dn2);
yaxis(-2,2);
hold on
plot(dn_atm2,wnorth/7.5,'r','linewidth',2)
ylabel(['Low-passed Wind speed' char(10) '(>36hrs, m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',-2:1:2,'yticklabels',[-2:1:2]*7.5);
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');

subplot(2,1,2);
plot([dn1 dn2],[0 0],'k','linewidth',1.5);
hold on 
plot(dn_atm,wind_n,'k','linewidth',2)
plot(dn_atm2,wnorth,'r','linewidth',2)
ylabel(['Wind speed' char(10) '(m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',-15:7.5:15);
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');
xaxis(dn1,dn2);
yaxis(-15,15);

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'windtimeseries']);
end 

%% lowpassed vs raw winds NE

close all; 
figure;
hold on; 
plot([-17.5 17.5],[0 0],'k','linewidth',1.5);
plot([0 0],[-17.5 17.5],'k','linewidth',1.5);
h1=plot(wind_e,wind_n,'.r','markersize',12);
h2=plot(weast,wnorth,'.k','markersize',12);
xaxis(-17.5,17.5); yaxis(-17.5,17.5);
grid on; hh=gca; hh.XMinorGrid='on'; 
xlabel('Eastward winds');
ylabel('Northward winds');
legend([h1 h2],'raw wind data','low-passed wind');
set(gca,'fontweight','bold','fontsize',14);

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'wind_ne']);
end 

%%  Principle component analysis of winds 
[ma,mi,deg]=pca(wnorth,weast);

RR=weast+sqrt(-1)*wnorth;
RR_rot=RR*exp(sqrt(-1)*-deg);
wmajor=-real(RR_rot);
wminor=-imag(RR_rot);

close all;
figure; 
hold on
plot([-17.5 17.5],[0 0],'k','linewidth',1.5);
plot([0 0],[-17.5 17.5],'k','linewidth',1.5);
h1=plot(weast,wnorth,'.k','markersize',12);
h2=plot(wminor,wmajor,'.','color',hex2rgb('#95E06C'),'markersize',12);
xaxis(-17.5,17.5); yaxis(-17.5,17.5);
grid on; hh=gca; hh.XMinorGrid='on'; 
xlabel('Eastward winds');
ylabel('Northward winds');
set(gca,'fontweight','bold','fontsize',14);
legend([h1 h2],'low-passed wind','rotated low-passed wind');

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'winds_rotated']);
end 

%% 
close all
figure('position',[173   612   897   181]); 
plot([dn1 dn2],[0 0],'k','linewidth',1.5);
hold on 
h1=plot(dn_atm,wind_n,'k','linewidth',2)
h2=plot(dn_atm2,wnorth,'r','linewidth',2)
h3=plot(dn_atm2,wmajor,'-','color',hex2rgb('#95E06C'),'linewidth',2);

ylabel(['Wind speed' char(10) '(m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',-15:7.5:15);
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');
xaxis(dn1,dn2);
yaxis(-15,15);
legend([h1 h2 h3],'raw wind','low-passed wind','rotated low-passed wind');

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'windtimeseries_rotated']);
end 

%% EOF analysis 
[F,A,D,~] = calc_eof(wmajor,'nan',0);

L=D*D';
percent_variance=100*diag(L*L')/trace(L*L');


% come back to this 9/26 -- need wind data from more than one place
%% Let's just pick from hand for now 

dwdt=interp1([dn_atm2(1:end-1)+dn_atm2(2:end)]/2,diff(wmajor),dn_atm2);

upwell=wmajor>-3;


dn1=datenum(2017,9,6);
dn2=datenum(2017,11,2);
dt=7; 

close all
figure('position',[173   612   897   181]); 
subplot(2,1,1)

plot([dn1 dn2],[0 0],'k','linewidth',1.5);
hold on 
h3=plot(dn_atm2,wmajor,'-','color',hex2rgb('#95E06C'),'linewidth',2);
h3=plot(dn_atm2(upwell),wmajor(upwell),'.r');
ylabel(['Wind speed' char(10) '(m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',-15:7.5:15);
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');
xaxis(dn1,dn2);
yaxis(-15,15);

subplot(2,1,2)
plot([dn1 dn2],[0 0],'k','linewidth',1.5);
hold on 
h3=plot(dn_atm2,dwdt,'-r','linewidth',2);

ylabel(['Wind speed' char(10) '(m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top');
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');
xaxis(dn1,dn2);

%% Plot Northward component of the winds

close all
figure('position',[440   609   889   189]); 
dn1=datenum(2017,9,6);
dn2=datenum(2017,11,2);
dt=7;

quiver(dn_atm2,dn_atm2.*0,weast,wnorth,'k'); %,'autoscale','off');%,'maxheadsize',0)
% quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,33),pl33tn(wind_n,1,33),'k','autoscale','off','maxheadsize',0)
xaxis(dn1,dn2);
yaxis(-2,2);
hold on
plot(dn_atm2,wmajor/7.5,'r','linewidth',2)
plot(dn_atm2(upwell),wmajor(upwell)/7.5,'.g','markersize',13)
ylabel(['Low-passed Wind speed' char(10) '(>36hrs, m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',-2:1:2,'yticklabels',[-2:1:2]*7.5);
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'upwellingtimes']);
end 

%%
dn=datenum(2017,9,6,0,0,0):(1/60)/24:datenum(2017,11,3,0,0,0);
% % datestr(dn_atm2(upwell))
uw1=dn>datenum(2017,9,1,10,0,0) & dn<datenum(2017,9,7,4,0,0);
uw2=dn>datenum(2017,9,11,6,0,0) & dn<datenum(2017,9,14,20,0,0);
uw3=dn>datenum(2017,9,26,16,0,0) & dn<datenum(2017,9,28,8,0,0);
uw4=dn>datenum(2017,10,4,0,0,0) & dn<datenum(2017,10,4,6,0,0);
uw5=dn>datenum(2017,10,9,6,0,0) & dn<datenum(2017,10,11,6,0,0);
uw6=dn>datenum(2017,10,15,6,0,0) & dn<datenum(2017,10,17,10,0,0);
uw7=dn>datenum(2017,10,23,20,0,0) & dn<datenum(2017,10,27,10,0,0);
uw8=dn>datenum(2017,10,30,4,0,0) & dn<datenum(2017,11,5,0,0,0);
uw=logical(uw1+uw2+uw3+uw4+uw5+uw6+uw7+uw8);

if 1
close all; 
figure('position',[440   651   914   147]); 
plot(dn_atm2(upwell),wmajor(upwell),'.k')
hold on
plot(dn(uw),dn(uw)*0,'.r')

set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top');
grid on; hh=gca; hh.XMinorGrid='on';
ax1=gca;
set(ax1,'color','none');
xaxis(dn1,dn2);
    if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'checkngUWindices']);
    end
end
%% Let's interpolate data to dn to make things easier 

east.oc50=interp2(OC50.dn,OC50.sig,OC50.Usig,dn,OC50.sig);
north.oc50=interp2(OC50.dn,OC50.sig,OC50.Vsig,dn,OC50.sig);
temp.oc50=interp2(OC50.dn,OC50.sig,OC50.Tsig,dn,OC50.sig);
H.oc50=interp2(OC50.dn,OC50.sig,OC50.H,dn,OC50.sig);

Nz = 50 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL50N.sig = (0:dsig:1)' ; % sigma  c
clear Nz sig 
east.nrl50n=interp2(NRL50N.dn,NRL50N.sig,NRL50N.Usig,dn,NRL50N.sig);
north.nrl50n=interp2(NRL50N.dn,NRL50N.sig,NRL50N.Vsig,dn,NRL50N.sig);
temp.nrl50n=interp2(NRL50N.dn,NRL50N.sig,NRL50N.Tsig,dn,NRL50N.sig);
H.nrl50n=interp2(NRL50N.dn,NRL50N.sig,NRL50N.H,dn,NRL50N.sig);

Nz = 25 ;  dsig = 1/Nz ; PS50.sig = (0:dsig:1)' ; clear Nz sig 
east.ps50=interp2(PS50.dn,PS50.sig,PS50.Usig,dn,PS50.sig);
north.ps50=interp2(PS50.dn,PS50.sig,PS50.Vsig,dn,PS50.sig);
temp.ps50=interp2(PS50.dn,PS50.sig,PS50.Tsig,dn,PS50.sig);
H.ps50=interp2(PS50.dn,PS50.sig,PS50.H,dn,PS50.sig);

east.vb50n=interp2(VB50N.dn,VB50N.sig,VB50N.Usig,dn,VB50N.sig);
north.vb50n=interp2(VB50N.dn,VB50N.sig,VB50N.Vsig,dn,VB50N.sig);
temp.vb50n=interp2(VB50N.dn,VB50N.sig,VB50N.Tsig,dn,VB50N.sig);
H.vb50n=interp2(VB50N.dn,VB50N.sig,VB50N.H,dn,VB50N.sig);

east.vb50s=interp2(VB50S.dn,VB50S.sig,VB50S.Usig,dn,VB50S.sig);
north.vb50s=interp2(VB50S.dn,VB50S.sig,VB50S.Vsig,dn,VB50S.sig);
temp.vb50s=interp2(VB50S.dn,VB50S.sig,VB50S.Tsig,dn,VB50S.sig);
H.vb50s=interp2(VB50S.dn,VB50S.sig,VB50S.H,dn,VB50S.sig);

Nz = 50;  dsig = 1/Nz ; NRL50S.sig = (0:dsig:1)' ; clear Nz sig 
east.nrl50s=interp2(NRL50S.dn,NRL50S.sig,NRL50S.Usig,dn,NRL50S.sig);
north.nrl50s=interp2(NRL50S.dn,NRL50S.sig,NRL50S.Vsig,dn,NRL50S.sig);
temp.nrl50s=interp2(NRL50S.dn,NRL50S.sig,NRL50S.Tsig,dn,NRL50S.sig);
H.nrl50s=interp2(NRL50S.dn,NRL50S.sig,NRL50S.H,dn,NRL50S.sig);

% check that the interpolations look right 
if 0
    moor='NRL50S'; %OC50, NRL50N, PS50, VB50N, VB50S, NRL50S   
    m=eval(moor);
    
    mm=eval(['H.' lower(moor)]);
    mmm=eval(['east.' lower(moor)]);
    mmmm=eval(['temp.' lower(moor)]);
    
    dn1=datenum(2017,9,11,0,0,0); dn2=datenum(2017,9,12,0,0,0);

    close all; figure('position',[25         389        1388         357]);
    subplot(2,1,1)
    i=findnearest_JACK(m.dn,dn1);
    ii=findnearest_JACK(m.dn,dn2);
    pcolorjw(m.dn_sig(:,i:ii),m.H(:,i:ii),m.Usig(:,i:ii));
    hold on 
    contour(m.dn_sig(:,i:ii),m.H(:,i:ii),m.Tsig(:,i:ii),[1:1:33],'k');
    caxis([-.4,.4]);
    xaxis(dn1,dn2);
    title([moor ' raw']);
    
    subplot(2,1,2)
    i=findnearest_JACK(dn,dn1);
    ii=findnearest_JACK(dn,dn2);    
    pcolorjw(repmat(dn(i:ii),size(mm,1),1),mm(:,i:ii),mmm(:,i:ii));
    hold on 
    contour(repmat(dn(i:ii),size(mm,1),1),mm(:,i:ii),mmmm(:,i:ii),[1:1:33],'k');
    caxis([-.4,.4]);
    xaxis(dn1,dn2);
    title([moor ' interp']);

end 

%% depth average currents during upwelling; 
% not super helpful 

close all; figure;
plot([-0.4 0.4],[0 0],'-k','linewidth',2); hold on;
plot([0 0],[-0.4 0.4],'-k','linewidth',2);
plot(nanmean(east.oc50(:,uw)),nanmean(north.oc50(:,uw)),'.k')
plot(nanmean(east.nrl50n(:,uw)),nanmean(north.nrl50n(:,uw)),'.r')
plot(nanmean(east.ps50(:,uw)),nanmean(north.ps50(:,uw)),'.g')
plot(nanmean(east.vb50n(:,uw)),nanmean(north.vb50n(:,uw)),'.m')
plot(nanmean(east.vb50s(:,uw)),nanmean(north.vb50s(:,uw)),'.c')
plot(nanmean(east.nrl50s(:,uw)),nanmean(north.nrl50s(:,uw)),'.y')
xaxis(-.4,0.4);yaxis(-0.4,0.4);
grid on;








