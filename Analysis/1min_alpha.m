% Written by Jack McSweeney 
% June 6, 2019 

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


%% Download Temp and Velocity Data - 
% These are all extrapolated, cleaned up products with 1 min resolution 

cleandatadir='/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/cleanup/';


OC50=load([cleandatadir 'OC50.mat']); OC50.dn=OC50.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc50_alpha_beta_Cw.mat','alpha','Cw');
OC50.alpha=tempo.alpha; OC50.Cw=tempo.Cw; clear tempo;

NRL50N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL50N.mat');
PS50=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS50.mat');
VB50N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB50N.mat');
VB50S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB50S.mat');
NRL50S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL50S.mat');

%% Calculate instantaneous alpha and co - OC50

OC50.alphai=OC50.dn.*nan; 
OC50.coi=OC50.dn.*nan; 

for ind=1:length(OC50.dn);
    om=2*pi/(12.4206*60*60);
    
    R = OC50.Rho_sig(:,ind);
    H0=OC50.H(:,ind);
    
    if ~isnan(nanmean(R))
        %  sort densities to obtain background profile
        RS = R;
        
        %  dither to get rid of common values
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        [rho,mm] = sort(RS,'descend') ;
        
        %  Calculate the wave speed
        dz = nanmean(diff(H0)) ;
        drhodz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(rho),H0)/dz;
        N2 = abs(-9.81*drhodz)/nanmean(rho) ;
        n2i=fillmissing(N2,'nearest');
        
        [c,ww] = nmodes(H0,n2i,om) ; % had to modify nmodes here
        [~,nn] = max(c) ;
        OC50.coi(ind) = c(nn) ;
        phi = ww(:,nn) ;
        phiz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(phi),H0)/dz;
        
        %  now calculate KdV parameters for the stratification
        OC50.alphai(ind) = -(3/2)*OC50.coi(ind)*nansum(phiz.^3)/nansum(phiz.^2) ;
    end
end

%% Reisolate the "prearrival" alpha as in paper 1 
load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/JPO_Manuscript/ReidentifyBores/GriddedBoresALL.mat','arrival');

for i =1:length(arrival.oc50)
    if arrival.oc50(i)>0
        ii=findnearest_JACK(OC50.dn,arrival.oc50(i)-0.5/24);
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)-1/60/24);
        OC50.pa_dn(i)=(OC50.dn(ii)+OC50.dn(jj))/2;
        mab=OC50.H(:,i);
        dz = nanmean(diff(mab)) ;
        r=nanmean(OC50.Rho_sig(:,ii:jj),2);
        drhodz=interp1([mab(1:end-1)+mab(2:end)]/2,diff(r),mab)/dz;
        N2 = abs(-9.8*drhodz)/nanmean(r) ;
        n2i=fillmissing(N2,'nearest');
        [c,w] = nmodes(mab,n2i,2*pi/(12.4206*60*60)) ; % had to modify nmodes here
        [vl,nn] = max(c) ;
        OC50.pa_Cw(i) = c(nn) ;
        phi = w(:,nn) ;
        phiz=interp1([mab(1:end-1)+mab(2:end)]/2,diff(phi),mab)/dz;
        OC50.pa_alpha(i) = -(3/2)*OC50.pa_Cw(i)*nansum(phiz.^3)/nansum(phiz.^2) ;
        
    else
        OC50.pa_Cw(i)=nan;
        OC50.pa_alpha(i)=nan; 
    end
end 

%% Check this alpha calc - OC50 
zoom='yes'; %yes or no 

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);
        dt=0.5;
end

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

close all; 
figure('position',[ 1   333        1309         459]); 
hax=tight_subplot(2,1,[0.07 0.1],[0.08 0.04],[0.05 0.02]);
axes(hax(1));
plot([dn1 dn2],[0 0],'--k'); hold on 
plot(OC50.dn,OC50.alphai,'-','color',0.7*[1 1 1],'linewidth',1);
plot(OC50.DNrho,OC50.alpha,'k','linewidth',3);
plot(OC50.pa_dn,OC50.pa_alpha,'r.','markersize',18);
% plot(arrival.oc50,OC50.pa_alpha,'g.','markersize',18);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',[],...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('\alpha (s^-^1)');

axes(hax(2));
plot(OC50.dn,OC50.coi,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(OC50.DNrho,OC50.Cw,'k','linewidth',3);
plot(OC50.pa_dn,OC50.pa_Cw,'r.','markersize',18);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',labels,...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('c_0 (ms^-^1)');

if 0
    set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'OC50alpha_co_comparison']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'OC50alpha_co_compZOOM']);
    end
end


%% Calculate instantaneous alpha and co - NRL50N

NRL50N.alphai=NRL50N.dn.*nan; 
NRL50N.coi=NRL50N.dn.*nan; 

for ind=1:length(NRL50N.dn);
    om=2*pi/(12.4206*60*60);
    
    R = NRL50N.Rho_sig(:,ind);
    H0=NRL50N.H(:,ind);
    
    if ~isnan(nanmean(R))
        %  sort densities to obtain background profile
        RS = R;
        
        %  dither to get rid of common values
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        [rho,mm] = sort(RS,'descend') ;
        
        %  Calculate the wave speed
        dz = nanmean(diff(H0)) ;
        drhodz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(rho),H0)/dz;
        N2 = abs(-9.81*drhodz)/nanmean(rho) ;
        n2i=fillmissing(N2,'nearest');
        
        [c,ww] = nmodes(H0,n2i,om) ; % had to modify nmodes here
        [~,nn] = max(c) ;
        NRL50N.coi(ind) = c(nn) ;
        phi = ww(:,nn) ;
        phiz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(phi),H0)/dz;
        
        %  now calculate KdV parameters for the stratification
        NRL50N.alphai(ind) = -(3/2)*NRL50N.coi(ind)*nansum(phiz.^3)/nansum(phiz.^2) ;
    end
end
    disp('calculations done!')

%% Check this alpha calc - NRL50N 

zoom='yes'; %yes or no 

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);
        dt=0.5;
end

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

close all; 
figure('position',[ 1   333        1309         459]); 
hax=tight_subplot(2,1,[0.07 0.1],[0.08 0.04],[0.05 0.02]);
axes(hax(1));
plot(NRL50N.dn,NRL50N.alphai,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.alpha,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',[],...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('\alpha (s^-^1)');
title('NRL50N');

axes(hax(2));
plot(NRL50N.dn,NRL50N.coi,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.Cw,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',labels,...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('c_0 (ms^-^1)');

if 0
    set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'NRL50Nalpha_co_comparison']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'NRL50Nalpha_co_compZOOM']);
    end
end

%% Calculate instantaneous alpha and co - PS50

PS50.alphai=PS50.dn.*nan; 
PS50.coi=PS50.dn.*nan; 

for ind=1:length(PS50.dn);
    om=2*pi/(12.4206*60*60);
    
    R = PS50.Rho_sig(:,ind);
    H0=PS50.H(:,ind);
    
    if ~isnan(nanmean(R))
        %  sort densities to obtain background profile
        RS = R;
        
        %  dither to get rid of common values
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        [rho,mm] = sort(RS,'descend') ;
        
        %  Calculate the wave speed
        dz = nanmean(diff(H0)) ;
        drhodz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(rho),H0)/dz;
        N2 = abs(-9.81*drhodz)/nanmean(rho) ;
        n2i=fillmissing(N2,'nearest');
        
        [c,ww] = nmodes(H0,n2i,om) ; % had to modify nmodes here
        [~,nn] = max(c) ;
        PS50.coi(ind) = c(nn) ;
        phi = ww(:,nn) ;
        phiz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(phi),H0)/dz;
        
        %  now calculate KdV parameters for the stratification
        PS50.alphai(ind) = -(3/2)*PS50.coi(ind)*nansum(phiz.^3)/nansum(phiz.^2) ;
    end
end
    disp('PS50 calculations done!')
    
%% Check this alpha calc - PS50 

zoom='no'; %yes or no 

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);
        dt=0.5;
end

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

close all; 
figure('position',[ 1   333        1309         459]); 
hax=tight_subplot(2,1,[0.07 0.1],[0.08 0.04],[0.05 0.02]);
axes(hax(1));
plot(PS50.dn,PS50.alphai,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(PS50.subtidal.DNrho,PS50.subtidal.alpha,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',[],...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('\alpha (s^-^1)');
title('PS50');

axes(hax(2));
plot(PS50.dn,PS50.coi,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(PS50.subtidal.DNrho,PS50.subtidal.Cw,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',labels,...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('c_0 (ms^-^1)');

if 0
    set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'PS50_alpha_co_comparison']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'PS50_alpha_co_compZOOM']);
    end
end

%% Calculate instantaneous alpha and co - VB50N

VB50N.alphai=VB50N.dn.*nan; 
VB50N.coi=VB50N.dn.*nan; 

for ind=1:length(VB50N.dn);
    om=2*pi/(12.4206*60*60);
    
    R = VB50N.Rho_sig(:,ind);
    H0=VB50N.H(:,ind);
    
    if ~isnan(nanmean(R))
        %  sort densities to obtain background profile
        RS = R;
        
        %  dither to get rid of common values
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        [rho,mm] = sort(RS,'descend') ;
        
        %  Calculate the wave speed
        dz = nanmean(diff(H0)) ;
        drhodz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(rho),H0)/dz;
        N2 = abs(-9.81*drhodz)/nanmean(rho) ;
        n2i=fillmissing(N2,'nearest');
        
        [c,ww] = nmodes(H0,n2i,om) ; % had to modify nmodes here
        [~,nn] = max(c) ;
        VB50N.coi(ind) = c(nn) ;
        phi = ww(:,nn) ;
        phiz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(phi),H0)/dz;
        
        %  now calculate KdV parameters for the stratification
        VB50N.alphai(ind) = -(3/2)*VB50N.coi(ind)*nansum(phiz.^3)/nansum(phiz.^2) ;
    end
end
    disp('VB50N calculations done!')
    
%% Check this alpha calc - VB50N 

zoom='no'; %yes or no 

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);
        dt=0.5;
end

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

close all; 
figure('position',[ 1   333        1309         459]); 
hax=tight_subplot(2,1,[0.07 0.1],[0.08 0.04],[0.05 0.02]);
axes(hax(1));
plot(VB50N.dn,VB50N.alphai,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(VB50N.subtidal.DNrho,VB50N.subtidal.alpha,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',[],...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('\alpha (s^-^1)');
title('VB50N');

axes(hax(2));
plot(VB50N.dn,VB50N.coi,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(VB50N.subtidal.DNrho,VB50N.subtidal.Cw,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',labels,...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('c_0 (ms^-^1)');

if 0
    set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'VB50N_alpha_co_comparison']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'VB50N_alpha_co_compZOOM']);
    end
end

%% Calculate instantaneous alpha and co - VB50S

VB50S.alphai=VB50S.dn.*nan; 
VB50S.coi=VB50S.dn.*nan; 

for ind=1:length(VB50S.dn);
    om=2*pi/(12.4206*60*60);
    
    R = VB50S.Rho_sig(:,ind);
    H0=VB50S.H(:,ind);
    
    if ~isnan(nanmean(R))
        %  sort densities to obtain background profile
        RS = R;
        
        %  dither to get rid of common values
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        [rho,mm] = sort(RS,'descend') ;
        
        %  Calculate the wave speed
        dz = nanmean(diff(H0)) ;
        drhodz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(rho),H0)/dz;
        N2 = abs(-9.81*drhodz)/nanmean(rho) ;
        n2i=fillmissing(N2,'nearest');
        
        [c,ww] = nmodes(H0,n2i,om) ; % had to modify nmodes here
        [~,nn] = max(c) ;
        VB50S.coi(ind) = c(nn) ;
        phi = ww(:,nn) ;
        phiz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(phi),H0)/dz;
        
        %  now calculate KdV parameters for the stratification
        VB50S.alphai(ind) = -(3/2)*VB50S.coi(ind)*nansum(phiz.^3)/nansum(phiz.^2) ;
    end
end
    disp('VB50S calculations done!')
    
%% Check this alpha calc - VB50S 

zoom='no'; %yes or no 

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);
        dt=0.5;
end

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

close all; 
figure('position',[ 1   333        1309         459]); 
hax=tight_subplot(2,1,[0.07 0.1],[0.08 0.04],[0.05 0.02]);
axes(hax(1));
plot(VB50S.dn,VB50S.alphai,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(VB50S.subtidal.DNrho,VB50S.subtidal.alpha,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',[],...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('\alpha (s^-^1)');
title('VB50S');

axes(hax(2));
plot(VB50S.dn,VB50S.coi,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(VB50S.subtidal.DNrho,VB50S.subtidal.Cw,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',labels,...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('c_0 (ms^-^1)');

if 0
    set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'VB50S_alpha_co_comparison']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'VB50S_alpha_co_compZOOM']);
    end
end

%% Calculate instantaneous alpha and co - NRL50S

NRL50S.alphai=NRL50S.dn.*nan; 
NRL50S.coi=NRL50S.dn.*nan; 

for ind=1:length(NRL50S.dn);
    om=2*pi/(12.4206*60*60);
    
    R = NRL50S.Rho_sig(:,ind);
    H0=NRL50S.H(:,ind);
    
    if ~isnan(nanmean(R))
        %  sort densities to obtain background profile
        RS = R;
        
        %  dither to get rid of common values
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        RS = RS + 1e-5*randn(size(RS)) ;
        [rho,mm] = sort(RS,'descend') ;
        
        %  Calculate the wave speed
        dz = nanmean(diff(H0)) ;
        drhodz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(rho),H0)/dz;
        N2 = abs(-9.81*drhodz)/nanmean(rho) ;
        n2i=fillmissing(N2,'nearest');
        
        [c,ww] = nmodes(H0,n2i,om) ; % had to modify nmodes here
        [~,nn] = max(c) ;
        NRL50S.coi(ind) = c(nn) ;
        phi = ww(:,nn) ;
        phiz=interp1([H0(1:end-1)+H0(2:end)]/2,diff(phi),H0)/dz;
        
        %  now calculate KdV parameters for the stratification
        NRL50S.alphai(ind) = -(3/2)*NRL50S.coi(ind)*nansum(phiz.^3)/nansum(phiz.^2) ;
    end
end
    disp('NRL50S calculations done!')
    
%% Check this alpha calc - NRL50S 

zoom='yes'; %yes or no 

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);
        dt=0.5;
end

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

close all; 
figure('position',[ 1   333        1309         459]); 
hax=tight_subplot(2,1,[0.07 0.1],[0.08 0.04],[0.05 0.02]);
axes(hax(1));
plot(NRL50S.dn,NRL50S.alphai,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.alpha,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',[],...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('\alpha (s^-^1)');
title('NRL50S');

axes(hax(2));
plot(NRL50S.dn,NRL50S.coi,'-','color',0.7*[1 1 1],'linewidth',1);
hold on
plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.Cw,'k','linewidth',3);
xaxis(dn1,dn2); 
set(gca,'fontweight','bold','xtick',dn1:dt:dn2,'xticklabels',labels,...
    'TickDir', 'both','TickLength',[0.005 0.005]);
ylabel('c_0 (ms^-^1)');

if 0
set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'NRL50S_alpha_co_comparison']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'NRL50S_alpha_co_compZOOM']);
    end
end

%% Filter the Velocities for the Kinetic Energy Calculations  

filt=33; %33; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 
filt1=16; % the filter cutoff; 16 hours here 
filt3= 3/60; % the filter cutoff; 3 min here

OC50.north_res=OC50.Vsig-repmat(nanmean(OC50.Vsig),length(OC50.sig),1);
OC50.east_res=OC50.Usig-repmat(nanmean(OC50.Usig),length(OC50.sig),1);
tempo=OC50.north_res-pl66tn(OC50.north_res,dt,filt1)';
OC50.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC50.east_res-pl66tn(OC50.east_res,dt,filt1)';
OC50.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC50.Wsig-pl66tn(OC50.Wsig,dt,filt1)';
OC50.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC50 velocity filtering');

Nz = 50 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL50N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
NRL50N.north_res=NRL50N.Vsig-repmat(nanmean(NRL50N.Vsig),length(NRL50N.sig),1);
NRL50N.east_res=NRL50N.Usig-repmat(nanmean(NRL50N.Usig),length(NRL50N.sig),1);
tempo=NRL50N.north_res-pl66tn(NRL50N.north_res,dt,filt1)';
NRL50N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50N.east_res-pl66tn(NRL50N.east_res,dt,filt1)';
NRL50N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50N.Wsig-pl66tn(NRL50N.Wsig,dt,filt1)';
NRL50N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished NRL50N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS50.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS50.north_res=PS50.Vsig-repmat(nanmean(PS50.Vsig),length(PS50.sig),1);
PS50.east_res=PS50.Usig-repmat(nanmean(PS50.Usig),length(PS50.sig),1);
tempo=PS50.north_res-pl66tn(PS50.north_res,dt,filt1)';
PS50.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS50.east_res-pl66tn(PS50.east_res,dt,filt1)';
PS50.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS50.Wsig-pl66tn(PS50.Wsig,dt,filt1)';
PS50.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS50 velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB50N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB50N.north_res=VB50N.Vsig-repmat(nanmean(VB50N.Vsig),length(VB50N.sig),1);
VB50N.east_res=VB50N.Usig-repmat(nanmean(VB50N.Usig),length(VB50N.sig),1);
tempo=VB50N.north_res-pl66tn(VB50N.north_res,dt,filt1)';
VB50N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50N.east_res-pl66tn(VB50N.east_res,dt,filt1)';
VB50N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50N.Wsig-pl66tn(VB50N.Wsig,dt,filt1)';
VB50N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB50N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB50S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB50S.north_res=VB50S.Vsig-repmat(nanmean(VB50S.Vsig),length(VB50S.sig),1);
VB50S.east_res=VB50S.Usig-repmat(nanmean(VB50S.Usig),length(VB50S.sig),1);
tempo=VB50S.north_res-pl66tn(VB50S.north_res,dt,filt1)';
VB50S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50S.east_res-pl66tn(VB50S.east_res,dt,filt1)';
VB50S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50S.Wsig-pl66tn(VB50S.Wsig,dt,filt1)';
VB50S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB50S velocity filtering');

Nz = 50 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL50S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
NRL50S.north_res=NRL50S.Vsig-repmat(nanmean(NRL50S.Vsig),length(NRL50S.sig),1);
NRL50S.east_res=NRL50S.Usig-repmat(nanmean(NRL50S.Usig),length(NRL50S.sig),1);
tempo=NRL50S.north_res-pl66tn(NRL50S.north_res,dt,filt1)';
NRL50S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50S.east_res-pl66tn(NRL50S.east_res,dt,filt1)';
NRL50S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50S.Wsig-pl66tn(NRL50S.Wsig,dt,filt1)';
NRL50S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished NRL50S velocity filtering');

%% Calculate Kinetic Energy Denisty

OC50.rho=nanmean(nanmean(OC50.Rho_sig));
OC50.ke_sd=nanmean(0.5*OC50.rho*(OC50.east_sd.^2+OC50.north_sd.^2+OC50.vert_sd.^2));

NRL50N.rho=nanmean(nanmean(NRL50N.Rho_sig));
NRL50N.ke_sd=nanmean(0.5*NRL50N.rho*(NRL50N.east_sd.^2+NRL50N.north_sd.^2+NRL50N.vert_sd.^2));

PS50.rho=nanmean(nanmean(PS50.Rho_sig));
PS50.ke_sd=nanmean(0.5*PS50.rho*(PS50.east_sd.^2+PS50.north_sd.^2+PS50.vert_sd.^2));

VB50N.rho=nanmean(nanmean(VB50N.Rho_sig));
VB50N.ke_sd=nanmean(0.5*VB50N.rho*(VB50N.east_sd.^2+VB50N.north_sd.^2+VB50N.vert_sd.^2));

VB50S.rho=nanmean(nanmean(VB50S.Rho_sig));
VB50S.ke_sd=nanmean(0.5*VB50S.rho*(VB50S.east_sd.^2+VB50S.north_sd.^2+VB50S.vert_sd.^2));

NRL50S.rho=nanmean(nanmean(NRL50S.Rho_sig));
NRL50S.ke_sd=nanmean(0.5*NRL50S.rho*(NRL50S.east_sd.^2+NRL50S.north_sd.^2+NRL50S.vert_sd.^2));
%% Save the 50m mooring data that we really want  

% make this a little lighter 
try 
NRL50N=rmfield(NRL50N,{'u','v','w','t','mab'});
PS50=rmfield(PS50,{'u','v','w','t','mab'});
VB50N=rmfield(VB50N,{'u','v','w','t','mab'});
VB50S=rmfield(VB50S,{'u','v','w','t','mab'});
NRL50S=rmfield(NRL50S,{'u','v','w','t','mab'});
catch 
end 

if 0
    save([anadir '50mMoorings.mat'],'OC50','NRL50N','PS50','VB50N','VB50S','NRL50S','-v7.3');
    disp(['Saved ' anadir '50mMoorings.mat']);
end 

%% Toss bad ke data 
OC50.ke_sd(OC50.ke_sd>45)=nan;
NRL50N.ke_sd(NRL50N.ke_sd>45)=nan;
PS50.ke_sd(PS50.ke_sd>45)=nan;
VB50N.ke_sd(VB50N.ke_sd>45)=nan;
VB50S.ke_sd(VB50S.ke_sd>45)=nan;
NRL50S.ke_sd(NRL50S.ke_sd>45)=nan;

%% ke averaged over 2 hours of each bore at OC50 
% (post bore =pb) 

for i =1:length(arrival.oc50)
    if arrival.oc50(i)>0
        ii=findnearest_JACK(OC50.dn,arrival.oc50(i));
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)+2/24);
        OC50.pb_dn(i)=(OC50.dn(ii)+OC50.dn(jj))/2;
        OC50.pb_ke(i)=nanmean(OC50.ke_sd(ii:jj));  
        
        ii=findnearest_JACK(OC50.dn,arrival.oc50(i));
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)+20/24/60);
        OC50.pb_dn20min(i)=(OC50.dn(ii)+OC50.dn(jj))/2;
        OC50.pb_ke20min(i)=nanmean(OC50.ke_sd(ii:jj));  
    else
        OC50.pb_ke(i)=nan;
    end
end 

%% OC KEsd - time series with the bore arrivals and the bore-averaged KE 
zoom='no';

switch zoom
    case 'no'
        dn1=datenum(2017,9,7);
        dn2=datenum(2017,11,2);
        dt=1;
        labels=datestr(dn1:dt:dn2,'mm/dd');
        labels(2:2:end,:) = nan; % remove every other one
    case 'yes'
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,15);
        dt=0.25;
        labels=datestr(dn1:dt:dn2,'mm/dd');
        labels(2:4:end,:) = nan; % remove every other one
        labels(3:4:end,:) = nan; % remove every other one
        labels(4:4:end,:) = nan; % remove every other one
end

close all; 
figure('position',[ 33         532        1350         260]); 
hold on;
plot(OC50.dn,OC50.ke_sd,'-','color',0.8*[1 1 1],'linewidth',1.25);
plot(OC50.dn,pl33tn(OC50.ke_sd,1/60,1),'-k','linewidth',3);
for i=1:length(arrival.oc50)
    plot([arrival.oc50(i) arrival.oc50(i)],[0 35],'-','color',hex2rgb('#46B5CE'),'linewidth',1.5);
end
plot(OC50.pb_dn ,OC50.pb_ke,'.r','markersize',20);

xaxis(dn1,dn2); yaxis(0,35);
ylabel('KE_s_d (J m^-^3)');
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels,'TickDir', 'both','TickLength',[0.005 0.005],...
    'fontweight','bold','fontsize',14,...
    'position',[0.0451851851851852         0.127307692307692         0.931111111111111         0.797692307692308]);
box on; grid on  

if 0
    set(gcf,'PaperPositionMode','auto');
    switch zoom
        case 'no'
            print(gcf,'-dpng','-r300',[anadir 'OC50_KEtimeseries']);
        case 'yes'
            print(gcf,'-dpng','-r300',[anadir 'OC50_KEtimeseriesZOOM']);
    end 
end


%% Is the KE of a bore correlated to pre-arrival alpha?? 
% doesn't look like it..  

[xfun, S]=polyfit(OC50.pa_alpha,OC50.pb_ke20min,1);
x0=-0.025:0.005:0.025;
[fit, delta]=polyval(xfun,x0,S);
[R ~]=corrcoef(OC50.pa_alpha,OC50.pb_ke20min);
R2= R(1,2).^2;

close all 
figure; hold on; 
plot([0 0],[0 15],'--k');
% plot(OC50.pa_alpha,OC50.pb_ke,'k.','markersize',15);
% ylabel(['Bore KE_s_d' char(10) '(avg 2 hours after bore passage)']);

plot(OC50.pa_alpha,OC50.pb_ke20min,'k.','markersize',15);
% plot(x0,fit,'-r')
ylabel(['Bore KE_s_d' char(10) '(avg 30 min after bore passage)']);
xaxis(-0.015,0.025);
xlabel('pre-arrival \alpha');
box on; grid on ;

if 0
    set(gcf,'PaperPositionMode','auto');    
%     print(gcf,'-dpng','-r300',[anadir 'OC50_alpha_ke_2hrs']);
    print(gcf,'-dpng','-r300',[anadir 'OC50_alpha_ke_30min']);
end


% corrcoef(OC50.pa_alpha,OC50.pb_ke)

%% look at change in alpha as bore comes - OC50 only  

for i =1:length(arrival.oc50)
    if arrival.oc50(i)>0
        ii=findnearest_JACK(OC50.dn,arrival.oc50(i)-0.5/24);
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)-1/60/24);
        prebore=nanmean(OC50.alphai(ii:jj));

        ii=findnearest_JACK(OC50.dn,arrival.oc50(i)+1/60/24);
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)+0.5/24);
        postbore=nanmean(OC50.alphai(ii:jj));
        OC50.dalpha(i)=prebore-postbore; 
        
    else
        OC50.dalpha(i)=nan;
    end 
        
end 
 
%% Is the change of alpha across the bore relate to the pre bore alpha ? 
% yes

[xfun, S]=polyfit(OC50.pa_alpha,OC50.dalpha,1);
x0=-0.025:0.005:0.025;
[fit, delta]=polyval(xfun,x0,S);
[R ~]=corrcoef(OC50.pa_alpha,OC50.dalpha);
R2= R(1,2).^2;

close all 
figure; hold on; 
plot([0 0],[-0.005 0.025],'--k');
plot(OC50.pa_alpha,OC50.dalpha,'k.','markersize',15);
plot(x0,fit,'-r')
ylabel('d\alpha/dt across a bore (s^-^1)');
xaxis(-0.015,0.025);
xlabel('pre-arrival \alpha');
box on; grid on ;
set(gca,'fontweight','bold','fontsize',14);
annotation(gcf,'textbox',[0.1625 0.80257 0.054464 0.060714],...
    'String',{['R^2=' num2str(round(R2,3))]},'LineStyle','none',...
    'fontweight','bold','fontsize',14);


if 0
    set(gcf,'PaperPositionMode','auto');    
    print(gcf,'-dpng','-r300',[anadir 'OC50_dalpha_paalpha']);
end


%% Is the KE of a bore correlated to the change of alpha across the bore? 
% maybe

[xfun, S]=polyfit(OC50.dalpha,OC50.pb_ke20min,1);
x0=-0.025:0.005:0.025;
[fit, delta]=polyval(xfun,x0,S);
[R ~]=corrcoef(OC50.dalpha,OC50.pb_ke20min);
R2= R(1,2).^2;

ind=OC50.pa_alpha>0.002;

close all 
figure; hold on; 
plot([0 0],[0 30],'--k');
plot(OC50.dalpha,OC50.pb_ke20min,'k.','markersize',15);
% plot(OC50.dalpha(ind),OC50.pb_ke20min(ind),'.','color',hex2rgb('#9EE273'),'markersize',15);
plot(x0,fit,'-r','linewidth',1.5)
xlabel('d\alpha/dt across a bore (s^-^1)');
xaxis(-0.005,0.022);
yaxis(0,25);
ylabel(['Bore KE_s_d' char(10) '(avg 2 hrs after bore passage)']);
box on; grid on ;
set(gca,'fontweight','bold','fontsize',14);
annotation(gcf,'textbox',[0.135 0.80257 0.054464 0.060714],...
    'String',{['R^2=' num2str(round(R2,3))]},'LineStyle','none',...
    'fontweight','bold','fontsize',14);


if 0
    set(gcf,'PaperPositionMode','auto');    
    print(gcf,'-dpng','-r300',[anadir 'OC50_dalpha_ke20mins']);
end

%% %% Is the KE of a bore correlated to the change of alpha across the bore? 
% colored by pre-arrival alpha 

[xfun, S]=polyfit(OC50.dalpha,OC50.pb_ke20min,1);
x0=-0.025:0.005:0.025;
[fit, delta]=polyval(xfun,x0,S);
[R ~]=corrcoef(OC50.dalpha,OC50.pb_ke20min);
R2= R(1,2).^2;

map=(cbrewer('div','Spectral',11));
close all; figure; hold on 
plot([0 0],[0 30],'--k');
scatter(OC50.dalpha,OC50.pb_ke20min,350,OC50.pa_alpha,'.');
cc=colorbar; colormap(map);
ylabel(cc,'\alpha (s^-^1)');
caxis([-0.01 0.01]);
xlabel('d\alpha/dt across a bore (s^-^1)');
xaxis(-0.005,0.022);
yaxis(0,25);
ylabel(['Bore KE_s_d' char(10) '(avg 2 hrs after bore passage)']);
box on; grid on ;
set(gca,'fontweight','bold','fontsize',14);

annotation(gcf,'textbox',[0.135 0.80257 0.054464 0.060714],...
    'String',{['R^2=' num2str(round(R2,3))]},'LineStyle','none',...
    'fontweight','bold','fontsize',14);

if 0
    set(gcf,'PaperPositionMode','auto');    
    print(gcf,'-dpng','-r300',[anadir 'OC50_dalpha_ke20mins_alpha']);
end

%% same figure - just plotted a different way 
% conclusion - still not there

map=(cbrewer('seq','PuBuGn',11));
close all 
scatter(OC50.pa_alpha,OC50.dalpha,350,OC50.pb_ke20min,'.');
colorbar; colormap(map)
caxis([0 14]);
 
%% look at change in alpha as bore comes - OC50 only  

for i =1:length(arrival.oc50)
    if arrival.oc50(i)>0
        ii=findnearest_JACK(OC50.dn,arrival.oc50(i)-0.5/24);
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)-1/60/24);
        OC50.stratpa(i)=nanmean(OC50.Tsig(end,ii:jj)-OC50.Tsig(1,ii:jj));
        
        ii=findnearest_JACK(OC50.dn,arrival.oc50(i)+1/60/24);
        jj=findnearest_JACK(OC50.dn,arrival.oc50(i)+0.5/24);
        OC50.stratpb(i)=nanmean(OC50.Tsig(end,ii:jj)-OC50.Tsig(1,ii:jj));

        OC50.dstrat(i)=OC50.stratpa(i)-OC50.stratpb(i);

    else
        OC50.stratpa(i)=nan;
        OC50.stratpb(i)=nan;
        OC50.dstrat(i)=nan;
    end 
        
end 

%% 
[R ~]=corrcoef(OC50.pa_Cw,OC50.pb_ke20min);
R2= R(1,2).^2;

close all; 
% plot(OC50.stratpa,OC50.pa_alpha,'.','markersize',20)
% plot(OC50.dstrat,OC50.dalpha,'.','markersize',20)
plot(OC50.pa_Cw,OC50.pb_ke20min,'.','markersize',20)

%% 
map=flipud(cbrewer('div','Spectral',11));

close all; 
scatter(OC50.dalpha,OC50.pb_ke20min,350,OC50.dstrat,'.');
colorbar; colormap(map);
caxis([-0.6 0.6]);