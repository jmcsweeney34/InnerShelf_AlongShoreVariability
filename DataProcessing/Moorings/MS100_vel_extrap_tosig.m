% Written by Jack McSweeney 
% September 21, 2018 

% In this script we extrapolate the velocities to the surface and bottom, 
% interpolate the velocites onto sigma coordinates, and then (both in sigma
% and z coordinates) remove the vertical mean and then bandpass filter to 
% subtidal (>33hrs), semidiurnal (10-16 hrs), and high frequency (<2hrs).    

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

%mooring data 
moordir='/Volumes/InnerShelf1/JackProcessing/Level2/';
veldir='/Volumes/InnerShelf1/JackProcessing/ADCPS/';

% coastline
load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalcoast.mat');

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/';

%% Download Data 

ms100=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/MS100-A/dpl1/mat/adcp_clean.mat');
ms100b=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/MS100-A/dpl2/mat/adcp_clean.mat');


%% Check that north is true north 
% % % according to the readme, its in true north

%% Nan out bad data 
ms100.A.east_vel(:,58471:end)=nan;
ms100.A.north_vel(:,58471:end)=nan;
ms100.A.vert_vel(:,58471:end)=nan;

ms100b.A.east_vel(:,58434:end)=nan;
ms100b.A.north_vel(:,58434:end)=nan;
ms100b.A.vert_vel(:,58434:end)=nan;

%% Transform to uniform grid in z 
dn=datenum(2017,9,8,0,0,0):1/(24*60):datenum(2017,11,3,0,0,0);
mab=(0:3:102)'; 

u=interp2([ms100.A.time ms100b.A.time],ms100.A.mab,[ms100.A.east_vel ms100b.A.east_vel],dn,mab);
v=interp2([ms100.A.time ms100b.A.time],ms100.A.mab,[ms100.A.north_vel ms100b.A.north_vel],dn,mab);
w=interp2([ms100.A.time ms100b.A.time],ms100.A.mab,[ms100.A.vert_vel ms100b.A.vert_vel],dn,mab);

east(1:size(u,1),1:size(u,2))=nan;
north(1:size(v,1),1:size(v,2))=nan;
vert(1:size(v,1),1:size(v,2))=nan;

%% Extrapolate currents to surface and bottom - EAST 

for i = 1:length(dn)
    if rem(i,1000) == 0
        disp([num2str(i) ' of ' num2str(length(dn))]);
    end
    %Extrapolate to surface
    utmp=u(:,i);
    try
        utmpi=utmp(~isnan(utmp));
        uscale=min(utmpi)-0.3:0.025:max(utmpi)+0.3;
        
        [k, ~]=findnearest_JACK(utmp,utmpi(end-2:end));
        [x, xi]=polyfit(utmp(k),mab(k),1);
        z=polyval(x,uscale);
        
        [kk, ~]=findnearest_JACK(utmp,utmpi(end-4:end));
        [x2, x2i]=polyfit(utmp(kk),mab(kk),3);
        z2=polyval(x2,uscale);
        
        utmpii=utmp;
        
        diff1=abs(utmpi(end)-interp1(z,uscale,mab(mab>max(mab(k)))));
        diff2= abs(utmpi(end)-interp1(z2,uscale,mab(mab>max(mab(k)))));
        
        if  nanmean(diff1) < nanmean(diff2)
            utmpii(mab>max(mab(k)))=interp1(z,uscale,mab(mab>max(mab(k))));
        elseif nanmean(diff1) > nanmean(diff2)
            utmpii(mab>max(mab(k)))=interp1(z2,uscale,mab(mab>max(mab(k))));
        end
        
        % clunky way to impose zero gradient at the top
        utmpii(end)=utmpii(end-1);
        
        %Extrapolate to bottom
        [kk, ~]=findnearest_JACK(utmp,utmpi(1:1+2));
        [xx, xxi]=polyfit(utmp(kk),mab(kk),1);
        zz=polyval(xx,uscale);
        utmpii(mab<min(mab(kk)))=interp1(zz,uscale,mab(mab<min(mab(kk))));
        % clunky way to impose zero gradient at the bottom
        utmpii(1)=utmpii(2);
        east(:,i)=utmpii;
        
%         clear umtp utmpi k x z utmpii kk xx zz
        
    catch
        east(:,i)=nan;
    end
    if i==length(dn)
        disp('Done with East!')
    end 
end

%% Plotting to check 
% % % % For problem solving the extrapolations  - need to uncomment clear
% % % % statement above
if 0
close all 
figure('position',[20    94   593   711]);
ii=max((findnearest_JACK(ms100.A.time,dn(i))));
subplot(1,2,1);
plot(ms100.A.east_vel(:,ii),ms100.A.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
hold on 
plot(u(:,i),mab,'.-k','markersize',40);
plot(utmp(k),mab(k),'.r','markersize',25);
plot(uscale,z,'.m','markersize',10);
plot(uscale,z2,'.b','markersize',10);
plot(utmpii,mab,'.g','markersize',14);
plot(utmp(kk),mab(kk),'.r','markersize',25);
plot(uscale,zz,'.m','markersize',10);
yaxis(mab(1),mab(end));


subplot(1,2,2);
plot(ms100.A.east_vel(:,ii),ms100.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
hold on 
plot(utmpii,mab,'.k','markersize',15);
yaxis(mab(1),mab(end));
leg=legend('original','extrapolated and interpolated','location','northoutside');
set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end 

% For an end check on extrapolations 
if 1
    i=9286;
    if dn(i)<datenum(2017,10,5,19,0,0)
        ii=max((findnearest_JACK(ms100.A.time,dn(i))));
    else
        ii=max((findnearest_JACK(ms100b.A.time,dn(i))));
    end
    
    close all
    figure('position',[440   116   325   682]);
    if dn(i)<datenum(2017,10,5,19,0,0)
        plot(ms100.A.east_vel(:,ii),ms100.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    else
        plot(ms100b.A.east_vel(:,ii),ms100b.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    end
    hold on
    plot(east(:,i),mab,'.k','markersize',15);
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'position',[0.277692317782113         0.941348973607038         0.492307692307692        0.0388563049853372]);
    yaxis(mab(1),mab(end));
    xaxis(nanmin(east(:,i)-0.05),nanmax(east(:,i)+0.05));
end


%% Extrapolate currents to surface and bottom - North 

for i = 1:length(dn)
    if rem(i,1000) == 0
        disp([num2str(i) ' of ' num2str(length(dn))]);
    end
    %Extrapolate to surface
    vtmp=v(:,i);
    try
        vtmpi=vtmp(~isnan(vtmp));
        vscale=min(vtmpi)-0.3:0.025:max(vtmpi)+0.3;
        
        [k, ~]=findnearest_JACK(vtmp,vtmpi(end-2:end));
        [x, xi]=polyfit(vtmp(k),mab(k),1);
        z=polyval(x,vscale);
        
        [kk, ~]=findnearest_JACK(vtmp,vtmpi(end-4:end));
        [x2, x2i]=polyfit(vtmp(kk),mab(kk),3);
        z2=polyval(x2,vscale);
        
        vtmpii=vtmp;
        
        diff1=abs(vtmpi(end)-interp1(z,vscale,mab(mab>max(mab(k)))));
        diff2= abs(vtmpi(end)-interp1(z2,vscale,mab(mab>max(mab(k)))));
        
        if  nanmean(diff1) < nanmean(diff2)
            vtmpii(mab>max(mab(k)))=interp1(z,vscale,mab(mab>max(mab(k))));
        elseif nanmean(diff1) > nanmean(diff2)
            vtmpii(mab>max(mab(k)))=interp1(z2,vscale,mab(mab>max(mab(k))));
        end
        
        % clunky way to impose zero gradient at the top
        vtmpii(end)=vtmpii(end-1);
        
        %Extrapolate to bottom
        [kk, ~]=findnearest_JACK(vtmp,vtmpi(1:1+2));
        [xx, xxi]=polyfit(vtmp(kk),mab(kk),1);
        zz=polyval(xx,vscale);
        vtmpii(mab<min(mab(kk)))=interp1(zz,vscale,mab(mab<min(mab(kk))));
        % clunky way to impose zero gradient at the bottom
        vtmpii(1)=vtmpii(2);
        north(:,i)=vtmpii;
        
%         clear vtmp vtmpi k x z vtmpii kk xx zz
        
    catch
        north(:,i)=nan;
    end
    if i==length(dn)
        disp('Done with North!')
    end 
end

%% Plotting to check North extrapolations
% % % % For problem solving the extrapolations  - need to uncomment clear
% % % % statement above
if 0
    close all
    figure('position',[20    94   593   711]);
    ii=max((findnearest_JACK(ms100.A.time,dn(i))));
    subplot(1,2,1);
    plot(ms100.A.north_vel(:,ii),ms100.A.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
    hold on
    plot(v(:,i),mab,'.-k','markersize',40);
    plot(vtmp(k),mab(k),'.r','markersize',25);
    plot(vscale,z,'.m','markersize',10);
    plot(vscale,z2,'.b','markersize',10);
    plot(vtmpii,mab,'.g','markersize',14);
    plot(vtmp(kk),mab(kk),'.r','markersize',25);
    plot(vscale,zz,'.m','markersize',10);
    yaxis(mab(1),mab(end));
    
    subplot(1,2,2);
    plot(ms100.A.north_vel(:,ii),ms100.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    hold on
    plot(vtmpii,mab,'.k','markersize',15);
    yaxis(mab(1),mab(end));
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end

% For an end check on extrapolations
if 1
    i=11111;
    if dn(i)<datenum(2017,10,5,19,0,0)
        ii=max((findnearest_JACK(ms100.A.time,dn(i))));
    else
        ii=max((findnearest_JACK(ms100b.A.time,dn(i))));
    end
    
    close all
    figure('position',[440   116   325   682]);
    if dn(i)<datenum(2017,10,5,19,0,0)
        plot(ms100.A.north_vel(:,ii),ms100.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    else
        plot(ms100b.A.north_vel(:,ii),ms100b.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    end
    hold on
    plot(north(:,i),mab,'.k','markersize',15);
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'position',[0.277692317782113         0.941348973607038         0.492307692307692        0.0388563049853372]);
    yaxis(mab(1),mab(end));
    xaxis(nanmin(north(:,i)-0.05),nanmax(north(:,i)+0.05));
end


%% Extrapolate currents to surface and bottom - Vertical 
for i = 1:length(dn)
    if rem(i,1000) == 0
        disp([num2str(i) ' of ' num2str(length(dn))]);
    end
    %Extrapolate to surface
    wtmp=w(:,i);
    try
        wtmpi=wtmp(~isnan(wtmp));
        wscale=min(wtmpi)-0.3:0.025:max(wtmpi)+0.3;
        
        [k, ~]=findnearest_JACK(wtmp,wtmpi(end-2:end));
        [x, xi]=polyfit(wtmp(k),mab(k),1);
        z=polyval(x,wscale);
        
        [kk, ~]=findnearest_JACK(wtmp,wtmpi(end-4:end));
        [x2, x2i]=polyfit(wtmp(kk),mab(kk),3);
        z2=polyval(x2,wscale);
        
        wtmpii=wtmp;
        
        diff1=abs(wtmpi(end)-interp1(z,wscale,mab(mab>max(mab(k)))));
        diff2= abs(wtmpi(end)-interp1(z2,wscale,mab(mab>max(mab(k)))));
        
        if  nanmean(diff1) < nanmean(diff2)
            wtmpii(mab>max(mab(k)))=interp1(z,wscale,mab(mab>max(mab(k))));
        elseif nanmean(diff1) > nanmean(diff2)
            wtmpii(mab>max(mab(k)))=interp1(z2,wscale,mab(mab>max(mab(k))));
        end
        
        % clunky way to impose zero gradient at the top
        wtmpii(end)=wtmpii(end-1);
        
        %Extrapolate to bottom
        [kk, ~]=findnearest_JACK(wtmp,wtmpi(1:1+2));
        [xx, xxi]=polyfit(wtmp(kk),mab(kk),1);
        zz=polyval(xx,wscale);
        wtmpii(mab<min(mab(kk)))=interp1(zz,wscale,mab(mab<min(mab(kk))));
        % clunky way to impose zero gradient at the bottom
        wtmpii(1)=wtmpii(2);
        vert(:,i)=wtmpii;
        
%         clear wtmp wtmpi k x z wtmpii kk xx zz
        
    catch
        vert(:,i)=nan;
    end
    if i==length(dn)
        disp('Done with Vertical!')
    end 
end


%% Plotting to check Vertical extrapolations
% % % % For problem solving the extrapolations  - need to uncomment clear
% % % % statement above
if 0
    close all
    figure('position',[20    94   593   711]);
    ii=max((findnearest_JACK(ms100.A.time,dn(i))));
    subplot(1,2,1);
    plot(ms100.A.vert_vel(:,ii),ms100.A.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
    hold on
    plot(w(:,i),mab,'.-k','markersize',40);
    plot(wtmp(k),mab(k),'.r','markersize',25);
    plot(wscale,z,'.m','markersize',10);
    plot(wscale,z2,'.b','markersize',10);
    plot(wtmpii,mab,'.g','markersize',14);
    plot(wtmp(kk),mab(kk),'.r','markersize',25);
    plot(wscale,zz,'.m','markersize',10);
    yaxis(mab(1),mab(end));
    
    subplot(1,2,2);
    plot(ms100.A.vert_vel(:,ii),ms100.A.mab,'.-','color',[0.8 0.8 0.8],'markersize',15);
    hold on
    plot(wtmpii,mab,'.-k','markersize',15);
    yaxis(mab(1),mab(end));
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end

% For an end check on extrapolations
if 1
    i=11111;
    if dn(i)<datenum(2017,10,5,19,0,0)
        ii=max((findnearest_JACK(ms100.A.time,dn(i))));
    else
        ii=max((findnearest_JACK(ms100b.A.time,dn(i))));
    end
    
    close all
    figure('position',[440   116   325   682]);
    if dn(i)<datenum(2017,10,5,19,0,0)
        plot(ms100.A.vert_vel(:,ii),ms100.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    else
        plot(ms100b.A.vert_vel(:,ii),ms100b.A.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    end
    hold on
    plot(vert(:,i),mab,'.k','markersize',15);
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'position',[0.277692317782113         0.941348973607038         0.492307692307692        0.0388563049853372]);
    yaxis(mab(1),mab(end));
    xaxis(nanmin(north(:,i)-0.05),nanmax(north(:,i)+0.05));
end

%% Get rid of bad interpolations by high pass filter - get rid of everything less than 5 mins
east_good=pl33tn(east,1/60,5/60)'; 
EAST=fillmissing(east_good,'linear');

north_good=pl33tn(north,1/60,5/60)'; 
NORTH=fillmissing(north_good,'linear');

vert_good=pl33tn(vert,1/60,5/60)'; 
VERT=fillmissing(vert_good,'linear');

%%
% Plot to check this 
if 1
    var='e'; % e, n or v 
    if strcmp(var,'e')
        map=flipud(cbrewer('div','RdBu',40));
    elseif  strcmp(var,'n')
        map=flipud(cbrewer('div','PuOr',40));
    elseif  strcmp(var,'v')
        map=flipud(cbrewer('div','BrBG',40));
    end
    
%         dn1=datenum(2017,9,14);
%     dn2=datenum(2017,9,15);
    dn1=736947.684997976;
    dn2=736947.788231128; 
%     dn1=736952.373;
%   dn2=736952.52;
    
    
%     close all
    figure('position',[ 57          79        1272         726]);
    subplot(5,1,1)
    if strcmp(var,'e')
        pcolorjw(dn,mab,u);
    elseif  strcmp(var,'n')
        pcolorjw(dn,mab,v);
    elseif  strcmp(var,'v')
        pcolorjw(dn,mab,w);    
    end
    colormap(map); colorbar;
    caxis([-0.5 0.5]);
    xaxis(dn1,dn2); yaxis(mab(1),mab(end));
    set(gca,'layer','top');
    
    subplot(5,1,2)
    if strcmp(var,'e')
        pcolorjw(dn,mab,east);
    elseif  strcmp(var,'n')
        pcolorjw(dn,mab,north);
    elseif  strcmp(var,'v')
        pcolorjw(dn,mab,vert);    
    end
    colormap(map); colorbar;
    caxis([-0.5 0.5]);
    xaxis(dn1,dn2);yaxis(mab(1),mab(end));
    set(gca,'layer','top');

    subplot(5,1,3)
    if strcmp(var,'e')
        pcolorjw(dn,mab,east_good);
    elseif  strcmp(var,'n')
        pcolorjw(dn,mab,north_good);
    elseif  strcmp(var,'v')
        pcolorjw(dn,mab,vert_good);    
    end
    colormap(map); colorbar;
    caxis([-0.5 0.5]);
    xaxis(dn1,dn2);yaxis(mab(1),mab(end));
    set(gca,'layer','top');
    
    subplot(5,1,4)
    if strcmp(var,'e')
        pcolorjw(dn,mab,east-east_good);
    elseif  strcmp(var,'n')
        pcolorjw(dn,mab,north-north_good);
    elseif  strcmp(var,'v')
        pcolorjw(dn,mab,vert-vert_good);    
    end
    colormap(map); colorbar;
    caxis([-0.5 0.5]);
    xaxis(dn1,dn2);yaxis(mab(1),mab(end));
    set(gca,'layer','top');
    
    subplot(5,1,5)
    if strcmp(var,'e')
        pcolorjw(dn,mab,EAST);
    elseif  strcmp(var,'n')
        pcolorjw(dn,mab,NORTH);
    elseif  strcmp(var,'v')
        pcolorjw(dn,mab,VERT);    
    end
    colormap(map); colorbar;
    caxis([-0.5 0.5]);
    xaxis(dn1,dn2);yaxis(mab(1),mab(end)); 
    set(gca,'layer','top');
end


%% Clear variables except those needed 

clearvars -except dn mab EAST NORTH VERT ms100 ms100b anadir


%% Load BT tide to help with sigma coordinates 
tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');

% Here's the concern - lacking good pressure data;  gonna have to be
% clever
if 0
    close all
    figure('position',[1         507        1440         292]);
    plot(tide.dn,tide.wl+nanmean(ms100.A.depth),'r','linewidth',2);
    hold on
    plot(ms100.A.time,pl33tn(ms100.A.depth,0.4/60,1),'k','linewidth',2);
    plot(ms100b.A.time,pl33tn(ms100b.A.depth,0.4/60,1),'k','linewidth',2);
    datetick('x','keeplimits','keepticks');
    legend('BT tide based on depth mean','adcp depth');
    if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir  'ms100DepthIssue']);
    end 
end

% decision is that D will be based on the tide data from San Luis 

%% Let's get the data into the same sigma coordinate system

Nz = 25 ;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig/2:dsig:1)' ; % sigma 


%  for the bottom depth, filter to keep variability with periods of 1 hour
%  and longer

% D=interp1(tide.dn,tide.wl+nanmean(ms100.A.depth),dn);
D=interp1(tide.dn,tide.wl+88,dn);
D=pl33tn(D,1/60,1)';
D=fillmissing(D,'linear');
H=sig*D;

Usig = NaN*ones(length(sig),length(dn)) ;
Vsig = NaN*ones(length(sig),length(dn)) ;
Wsig = NaN*ones(length(sig),length(dn)) ;

%  OK now let's interpolate onto the sigma grid.  
for i = 1:length(dn) 
    if mod(i,1000) == 0
        disp([num2str(i) '/' num2str(length(dn))]) ;
    end
    
    for ii=1:length(sig)
        if ii==1
            kk=find(mab<H(ii,i));
        else
            kk=find(mab>H(ii-1,i)& mab<H(ii,i));
        end
        Usig(ii,i)=nanmean(EAST(kk,i));
        Vsig(ii,i)=nanmean(NORTH(kk,i));
        Wsig(ii,i)=nanmean(VERT(kk,i));
    end 
end
    
%% Plot the sigma and z against each other to show its good 
if 1
    map=flipud(cbrewer('div','RdBu',40));
    close all
    dn1=datenum(2017,9,14);
    dn2=datenum(2017,9,16);
    %     dn1=736947.684997976;
    %     dn2=736947.788231128;
    %     dn1=736952.373;
    %   dn2=736952.52;
    figure;
    subplot(2,1,1);
    pcolorjw(dn,mab,EAST);
    colormap(map); colorbar; caxis([-0.5 0.5]);
    xaxis(dn1,dn2);
    subplot(2,1,2);
    pcolorjw(repmat(dn,length(sig),1),H,Usig);
    colormap(map); colorbar; caxis([-0.5 0.5]);
    xaxis(dn1,dn2);
end

%% 
close all; 
figure('position',[440   188   309   610]); 
i=60300;
plot([0 0],[0 max(H(:,i))+1],'--k','linewidth',1) ;
hold on 
try
    if dn(i)<datenum(2017,10,5,19,0,0)
        ii=max((findnearest_JACK(ms100.A.time,dn(i))));
        h1=plot(ms100.A.east_vel(:,ii),ms100.A.mab,'.-','color',[0.5 0.5 0.5],'linewidth',2,'markersize',40);
    else
        ii=max((findnearest_JACK(ms100b.A.time,dn(i))));
        h1=plot(ms100b.A.east_vel(:,ii),ms100b.A.mab,'.-','color',[0.5 0.5 0.5],'linewidth',2,'markersize',40);
    end
    h2=plot(EAST(:,i),mab,'.-k','linewidth',2,'markersize',20);
    h3=plot(Usig(:,i),H(:,i),'.-r','linewidth',2,'markersize',20);
    legend([h1 h2 h3],'raw','z','sig'); 
catch
end

yaxis(0, max(H(:,i))+1);
xaxis(min(Usig(:,i))-0.1,max(Usig(:,i))+0.1);

%% Clear variables except those needed 
dn_sig=repmat(dn,length(sig),1);

clearvars -except dn mab EAST NORTH VERT H sig Usig Vsig Wsig anadir dn_sig D

if 0
    save([anadir 'ms100vel_zandsig.mat']);
end 






