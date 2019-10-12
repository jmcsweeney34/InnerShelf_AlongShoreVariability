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
anadir= '/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/FixingVelocities/';

%% Download Data 

if 0
    save([anadir 'MS100ExtrapolatedVel.mat']);
end 

%% Remove the vertical mean from the velocity data 
ms100.dn=dn;
ms100.dn_sig=repmat(dn,length(sig),1);
ms100.mab=mab;
ms100.sig=sig;
ms100.H=H;

ms100.north_z.total=NORTH;
ms100.north_sig.total=Vsig;
ms100.east_z.total=EAST;
ms100.east_sig.total=Usig;
ms100.vert_z.total=VERT;
ms100.vert_sig.total=Wsig;

ms100.north_z.depthmean=nanmean(NORTH);
ms100.north_sig.depthmean=nanmean(Vsig);
ms100.east_z.depthmean=nanmean(EAST);
ms100.east_sig.depthmean=nanmean(Usig);

% residual is total minus depth mean 
ms100.north_z.res=NORTH-repmat(ms100.north_z.depthmean,length(ms100.mab),1);
ms100.north_sig.res=Vsig-repmat(ms100.north_sig.depthmean,length(ms100.sig),1);
ms100.east_z.res=EAST-repmat(ms100.east_z.depthmean,length(ms100.mab),1);
ms100.east_sig.res=Usig-repmat(ms100.east_sig.depthmean,length(ms100.sig),1);


%% Lowpass the velocity data - isolate subtidal  (>33hrs)
% will use the residual for horizontal currents and the total for vertical 
filt=33; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 40 seconds for MS100 

ms100.north_z.subtidal=pl66tn(ms100.north_z.res,dt,filt)';
ms100.north_sig.subtidal=pl66tn(ms100.north_sig.res,dt,filt)';

ms100.east_z.subtidal=pl66tn(ms100.east_z.res,dt,filt)';
ms100.east_sig.subtidal=pl66tn(ms100.east_sig.res,dt,filt)';

ms100.vert_z.subtidal=pl66tn(ms100.vert_z.total,dt,filt)';
ms100.vert_sig.subtidal=pl66tn(ms100.vert_sig.total,dt,filt)';

%% Lowpass the velocity data - isolate diurnal (between 16 and 33)
% will use the residual for horizontal currents and the total for vertical 
filt=16; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 40 seconds for MS100 

ms100.north_z.diurnal=pl66tn(ms100.north_z.res-ms100.north_z.subtidal,dt,filt)';
ms100.north_sig.diurnal=pl66tn(ms100.north_sig.res-ms100.north_sig.subtidal,dt,filt)';

ms100.east_z.diurnal=pl66tn(ms100.east_z.res-ms100.east_z.subtidal,dt,filt)';
ms100.east_sig.diurnal=pl66tn(ms100.east_sig.res-ms100.east_sig.subtidal,dt,filt)';

ms100.vert_z.diurnal=pl66tn(ms100.vert_z.total-ms100.vert_z.subtidal,dt,filt)';
ms100.vert_sig.diurnal=pl66tn(ms100.vert_sig.total-ms100.vert_sig.subtidal,dt,filt)';

%% Lowpass the velocity data - isolate semidiurnal(between 10 and 16 )
% will use the residual for horizontal currents and the total for vertical 
filt=10; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 40 seconds for MS100 

ms100.north_z.semidiurnal=pl66tn(ms100.north_z.res-ms100.north_z.diurnal,dt,filt)';
ms100.north_sig.semidiurnal=pl66tn(ms100.north_sig.res-ms100.north_sig.diurnal,dt,filt)';

ms100.east_z.semidiurnal=pl66tn(ms100.east_z.res-ms100.east_z.diurnal,dt,filt)';
ms100.east_sig.semidiurnal=pl66tn(ms100.east_sig.res-ms100.east_sig.diurnal,dt,filt)';

ms100.vert_z.semidiurnal=pl66tn(ms100.vert_z.total-ms100.vert_z.diurnal,dt,filt)';
ms100.vert_sig.semidiurnal=pl66tn(ms100.vert_sig.total-ms100.vert_sig.diurnal,dt,filt)';

%% Lowpass the velocity data - isolate highfreq (less than 10 hrs )
ms100.north_z.highfreq=ms100.north_z.res-ms100.north_z.subtidal-ms100.north_z.diurnal-ms100.north_z.semidiurnal;
ms100.north_sig.highfreq=ms100.north_sig.res-ms100.north_sig.subtidal-ms100.north_sig.diurnal-ms100.north_sig.semidiurnal;

ms100.east_z.highfreq=ms100.east_z.res-ms100.east_z.subtidal-ms100.east_z.diurnal-ms100.east_z.semidiurnal;
ms100.east_sig.highfreq=ms100.east_sig.res-ms100.east_sig.subtidal-ms100.east_sig.diurnal-ms100.east_sig.semidiurnal;

ms100.vert_z.highfreq=ms100.vert_z.total-ms100.vert_z.subtidal-ms100.vert_z.diurnal-ms100.vert_z.semidiurnal;
ms100.vert_sig.highfreq=ms100.vert_sig.total-ms100.vert_sig.subtidal-ms100.vert_sig.diurnal-ms100.vert_sig.semidiurnal;

%% Sanity check - make sure everything adds up 
if 0
    
    figure('position',[ 44         385        1028         408]) ;
    subplot(3,1,1);
    map=flipud(cbrewer('div','PuOr',40));
    pcolorjw(ms100.dn,ms100.mab,ms100.north_z.total-...
        ms100.north_z.depthmean-...
        ms100.north_z.subtidal-...
        ms100.north_z.diurnal-...
        ms100.north_z.semidiurnal-...
        ms100.north_z.highfreq);
    caxis([-0.005 0.005])
    colormap(map);
    colorbar;
    datetick('x','keeplimits','keepticks')
    
    subplot(3,1,2);
    map=flipud(cbrewer('div','RdBu',40));
    pcolorjw(ms100.dn,ms100.mab,ms100.east_z.total-...
        ms100.east_z.depthmean-...
        ms100.east_z.subtidal-...
        ms100.east_z.diurnal-...
        ms100.east_z.semidiurnal-...
        ms100.east_z.highfreq);
    caxis([-0.005 0.005])
    colormap(map);
    colorbar;
    datetick('x','keeplimits','keepticks')
    
    subplot(3,1,3);
    map=flipud(cbrewer('div','BrBG',40));
    pcolorjw(ms100.dn,ms100.mab,ms100.vert_z.total-...
        ms100.vert_z.subtidal-...
        ms100.vert_z.diurnal-...
        ms100.vert_z.semidiurnal-...
        ms100.vert_z.highfreq);
    caxis([-0.005 0.005])
    colormap(map);
    colorbar;
    datetick('x','keeplimits','keepticks')
end

%% Lowpass the velocity data - isolate iw (less than 2 hrs )
% will use the residual for horizontal currents and the total for vertical 
filt=2; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 40 seconds for MS100 

ms100.north_z.iw=ms100.north_z.highfreq-pl66tn(ms100.north_z.highfreq,dt,filt)';
ms100.north_sig.iw=ms100.north_sig.highfreq-pl66tn(ms100.north_sig.highfreq,dt,filt)';

ms100.east_z.iw=ms100.east_z.highfreq-pl66tn(ms100.east_z.highfreq,dt,filt)';
ms100.east_sig.iw=ms100.east_sig.highfreq-pl66tn(ms100.east_sig.highfreq,dt,filt)';

ms100.vert_z.iw=ms100.vert_z.highfreq-pl66tn(ms100.vert_z.highfreq,dt,filt)';
ms100.vert_sig.iw=ms100.vert_sig.highfreq-pl66tn(ms100.vert_sig.highfreq,dt,filt)';



%% Let's save the filtered currents in both z and sig 
if 0
    save([anadir 'ms100filteredvel_zandsig.mat'],'ms100');
end 

%% Plot the NORTH total, residual, subtidal> 33 , diurnal (16-33), semidurnal (10-16), and high freq (<10) 
if 1
    close all
    save='n' ; % y or n
    
    dn1=datenum(2017,9,8);
    dn2=datenum(2017,11,3);
    dt=4;
    
    mapb=flipud(cbrewer('div','PuOr',40));
    figure('position',[35 69  1309  736]);
    hax=tight_subplot(6,1,[0.08 0.1],[0.09 0.07],[0.04 0.07]);

    axes(hax(1))        
    pcolorjw(ms100.dn, ms100.mab,ms100.north_z.total); hold on 
    colorbar('peer',hax(1),'Position',[0.93812 0.086957 0.023682 0.84239])
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title(['MS100' char(10) 'Total North Velocity (m/s)']);
    xaxis(dn1,dn2);

    axes(hax(2))
    pcolorjw(ms100.dn, ms100.mab,ms100.north_z.res); hold on
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title(['North Residual Velocity (m/s)']);
    xaxis(dn1,dn2);
    
    axes(hax(3))
    pcolorjw(ms100.dn, ms100.mab,ms100.north_z.subtidal); hold on 
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    xaxis(dn1,dn2);
    title('North Subtidal Velocity (m/s, >33hrs)');
    
    axes(hax(4))
    pcolorjw(ms100.dn, ms100.mab,ms100.north_z.diurnal); hold on 
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    xaxis(dn1,dn2);
    title('North Diurnal Velocity (m/s, 16-33 hrs)');
    
    axes(hax(5))
    pcolorjw(ms100.dn, ms100.mab,ms100.north_z.semidiurnal); hold on 
    xaxis(dn1,dn2);
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title('North Semidiurnal Velocity (m/s, 10-16 hrs)');
    colormap(mapb)
    
    axes(hax(6))
    pcolorjw(ms100.dn, ms100.mab,ms100.north_z.highfreq); hold on 
    xaxis(dn1,dn2);
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title('North High Frequency Velocity (m/s, < 10 hrs)');
    colormap(mapb)
    
    if strcmp(save,'y')
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'MS100NorthVelocity']);
    end
    
end



%% Plot the EAST total, residual, subtidal> 33 , diurnal (16-33), semidurnal (10-16), and high freq (<10) 
if 1
    close all
    save='n' ; % y or n
    
    dn1=datenum(2017,9,8);
    dn2=datenum(2017,11,3);
    dt=4;
    
    mapb=flipud(cbrewer('div','RdBu',40));
    figure('position',[35 69  1309  736]);
    hax=tight_subplot(6,1,[0.08 0.1],[0.09 0.07],[0.04 0.07]);

    axes(hax(1))        
    pcolorjw(ms100.dn, ms100.mab,ms100.east_z.total); hold on 
    colorbar('peer',hax(1),'Position',[0.95187 0.086957 0.023682 0.84239])
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title(['MS100' char(10) 'Total East Velocity (m/s)']);
    xaxis(dn1,dn2);

    axes(hax(2))
    pcolorjw(ms100.dn, ms100.mab,ms100.east_z.res); hold on
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title(['East Residual Velocity (m/s)']);
    xaxis(dn1,dn2);
    
    axes(hax(3))
    pcolorjw(ms100.dn, ms100.mab,ms100.east_z.subtidal); hold on 
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    xaxis(dn1,dn2);
    title('East Subtidal Velocity (m/s, >33hrs)');
    
    axes(hax(4))
    pcolorjw(ms100.dn, ms100.mab,ms100.east_z.diurnal); hold on 
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    xaxis(dn1,dn2);
    title('East Diurnal Velocity (m/s, 16-33 hrs)');
    
    axes(hax(5))
    pcolorjw(ms100.dn, ms100.mab,ms100.east_z.semidiurnal); hold on 
    xaxis(dn1,dn2);
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title('East Semidiurnal Velocity (m/s, 10-16 hrs)');
    colormap(mapb)
    
    axes(hax(6))
    pcolorjw(ms100.dn, ms100.mab,ms100.east_z.highfreq); hold on 
    xaxis(dn1,dn2);
    caxis([-0.5, 0.5])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title('East High Frequency Velocity (m/s, < 10 hrs)');
    colormap(mapb)
    
    if strcmp(save,'y')
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'MS100EastVelocity']);
    end
    
end


%% Plot the Vertical total, subtidal> 33 , diurnal (16-33), semidurnal (10-16), and high freq (<10) 
if 1
    close all
    save='n' ; % y or n
    
    dn1=datenum(2017,9,8);
    dn2=datenum(2017,11,3);
    dt=4;
    
    mapb=flipud(cbrewer('div','BrBG',40));
    figure('position',[35 69  1309  736]);
    hax=tight_subplot(5,1,[0.08 0.1],[0.09 0.07],[0.04 0.07]);

    axes(hax(1))        
    pcolorjw(ms100.dn, ms100.mab,ms100.vert_z.total); hold on 
    colorbar('peer',hax(1),'Position',[0.95187 0.086957 0.023682 0.84239])
    caxis([-0.1, 0.1])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title(['MS100' char(10) 'Total Vertical Velocity (m/s)']);
    xaxis(dn1,dn2);

    
    axes(hax(2))
    pcolorjw(ms100.dn, ms100.mab,ms100.vert_z.subtidal); hold on 
    caxis([-0.1, 0.1])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    xaxis(dn1,dn2);
    title('Vertical Subtidal Velocity (m/s, >33hrs)');
    
    axes(hax(3))
    pcolorjw(ms100.dn, ms100.mab,ms100.vert_z.diurnal); hold on 
    caxis([-0.1, 0.1])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    xaxis(dn1,dn2);
    title('Vertical Diurnal Velocity (m/s, 16-33 hrs)');
    
    axes(hax(4))
    pcolorjw(ms100.dn, ms100.mab,ms100.vert_z.semidiurnal); hold on 
    xaxis(dn1,dn2);
    caxis([-0.1, 0.1])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title('Vertical Semidiurnal Velocity (m/s, 10-16 hrs)');
    colormap(mapb)
    
    axes(hax(5))
    pcolorjw(ms100.dn, ms100.mab,ms100.vert_z.highfreq); hold on 
    xaxis(dn1,dn2);
    caxis([-0.1, 0.1])
    set(gca,'fontweight','bold','fontsize',14,'layer','top','xtick',dn1:dt:dn2,'xticklabel',datestr(dn1:dt:dn2,'mm/dd'));
    title('Vertical High Frequency Velocity (m/s, < 10 hrs)');
    colormap(mapb)
    
    if strcmp(save,'y')
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'MS100VertVelocity']);
    end
    
end

%% 
close all 
figure;
plot([-0.5 0.5],[0 0],'--k','linewidth',1); hold on 
plot([0 0],[-0.5 0.5],'--k','linewidth',1);
plot(nanmean(ms100.east_z.subtidal),nanmean(ms100.north_z.subtidal),'.k');
plot(nanmean(ms100.east_z.diurnal),nanmean(ms100.north_z.diurnal),'.g');
plot(nanmean(ms100.east_z.semidiurnal),nanmean(ms100.north_z.semidiurnal),'.r');





%% 

close all 
figure;
% quiver(ms100.A.time,ms100.A.time.*0,nanmean(ms100.east_lp),nanmean(ms100.north_hp),'autoscale','off','maxheadsize',0)
quiver(dn,dn.*0,nanmean(ms100.east_z.subtidal),nanmean(ms100.north_z.subtidal),'autoscale','off','maxheadsize',0)


%% 





