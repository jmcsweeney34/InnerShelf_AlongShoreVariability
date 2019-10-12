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
oc40s=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/OC40S-A/Level 1/deploy1/UVWtot_1min.mat');
oc40sb=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/OC40S-A/Level 1/deploy2/OC40S-A.mat');


%% Check that north is true north 
% % % as of 12/8/18  I'm not sure yet -- NEED to come back to this 


%% Remove bad data 

perbad=15;
oc40s.A.z=meshgrid(oc40s.A.zab,oc40s.A.datenum)';
oc40s.A.maskjm=double(oc40s.A.z<(40-(40*perbad/100)));
oc40s.A.maskjm(oc40s.A.maskjm==0)=nan;
oc40s.A.maskjm(:,oc40s.A.datenum<datenum(2017,9,9,0,0,0))=nan;

oc40sb.A.z=meshgrid(oc40sb.A.zab,oc40sb.A.datenum)';
oc40sb.A.maskjm=double(oc40sb.A.z<(40-(40*perbad/100)));
oc40sb.A.maskjm(oc40sb.A.maskjm==0)=nan;
oc40sb.A.maskjm(:,oc40sb.A.datenum<datenum(2017,9,9,0,0,0))=nan;


oc40s.u=oc40s.A.u.*oc40s.A.maskjm; oc40sb.u=oc40sb.A.u.*oc40sb.A.maskjm;
oc40s.v=oc40s.A.v.*oc40s.A.maskjm; oc40sb.v=oc40sb.A.v.*oc40sb.A.maskjm;
oc40s.w5=oc40s.A.w5.*oc40s.A.maskjm; oc40sb.w5=oc40sb.A.w5.*oc40sb.A.maskjm;

oc40s.u(:,41133:end)=nan; 
oc40sb.u(:,1:1007)=nan; oc40sb.u(:,34191:end)=nan;

oc40s.v(:,41133:end)=nan; 
oc40sb.v(:,1:1007)=nan; oc40sb.v(:,34191:end)=nan;

oc40s.w5(:,41133:end)=nan; 
oc40sb.w5(:,1:1007)=nan; oc40sb.w5(:,34191:end)=nan;

%% Transform to uniform grid in z 
dn=datenum(2017,9,8,0,0,0):1/(24*60):datenum(2017,11,1,0,0,0);
mab=(0:0.5:43)'; 

u=interp2([oc40s.A.datenum oc40sb.A.datenum],oc40s.A.zab,[oc40s.u oc40sb.u],dn,mab);
v=interp2([oc40s.A.datenum oc40sb.A.datenum],oc40s.A.zab,[oc40s.v oc40sb.v],dn,mab);
w=interp2([oc40s.A.datenum oc40sb.A.datenum],oc40s.A.zab,[oc40s.w5 oc40sb.w5],dn,mab);

%% Let's get the data into a sigma coordinate system
east(1:size(u,1),1:size(u,2))=nan;
north(1:size(v,1),1:size(v,2))=nan;
vert(1:size(v,1),1:size(v,2))=nan;

%% Extrapolate currents to surface and bottom - EAST 

for i = 1:length(dn)
    if rem(i,1000) == 0
        disp([num2str(i) ' of ' num2str(length(dn))]);
    end
    %Extrapolate to surface
%     utmp=u(:,i);
    utmp=sgolayfilt(double(u(:,i)),1,5);
    try
        utmpi=utmp(~isnan(utmp));
        uscale=min(utmpi)-0.3:0.025:max(utmpi)+0.3;
        warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
        
        [k, ~]=findnearest_JACK(utmp,utmpi(end-3:end));
        [x, xi]=polyfit(utmp(k),mab(k),1);
        z=polyval(x,uscale);

        [kk, ~]=findnearest_JACK(utmp,utmpi(end-6:end));
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
        [kk, ~]=findnearest_JACK(utmp,utmpi(1:1+3));
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
ii=max((findnearest_JACK(oc40s.A.datenum,dn(i))));
subplot(1,2,1);
plot(oc40s.u(:,ii),oc40s.A.zab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
plot(oc40n.u(:,ii),oc40n.mab,'.-','color',[0.8 0.8 0.8],'markersize',40);
hold on 
plot(utmpii,mab,'.-k','markersize',15);
yaxis(mab(1),mab(end));
leg=legend('original','extrapolated and interpolated','location','northoutside');
set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end 

% For an end check on extrapolations 
if 1
    i=1713;
        ii=max((findnearest_JACK(oc40s.A.datenum,dn(i))));
        
    close all
    figure('position',[440   116   325   682]);
        plot(oc40s.u(:,ii),oc40s.A.zab,'.','color',[0.8 0.8 0.8],'markersize',40);
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
%     vtmp=v(:,i);
    vtmp=sgolayfilt(double(v(:,i)),1,5);

    try
        vtmpi=vtmp(~isnan(vtmp));
        vscale=min(vtmpi)-0.3:0.025:max(vtmpi)+0.3;

        [k, ~]=findnearest_JACK(vtmp,vtmpi(end-3:end));
        [x, xi]=polyfit(vtmp(k),mab(k),1);
        z=polyval(x,vscale);
        
        [kk, ~]=findnearest_JACK(vtmp,vtmpi(end-6:end));
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
%         
        % clunky way to impose zero gradient at the top
        vtmpii(end)=vtmpii(end-1);
        
        %Extrapolate to bottom
        [kk, ~]=findnearest_JACK(vtmp,vtmpi(1:1+3));
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
    ii=max((findnearest_JACK(oc40n.dn,dn(i))));
    subplot(1,2,1);
    plot(oc40n.v(:,ii),oc40n.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
    plot(oc40n.v(:,ii),oc40n.mab,'.','color',[0.8 0.8 0.8],'markersize',40);
    hold on
    plot(vtmpii,mab,'.k','markersize',15);
    yaxis(mab(1),mab(end));
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end

% For an end check on extrapolations
if 1
    i=1713;
    ii=max((findnearest_JACK(oc40s.A.datenum,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(oc40s.v(:,ii),oc40s.A.zab,'.','color',[0.8 0.8 0.8],'markersize',40);
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
%     wtmp=w(:,i);
        wtmp=sgolayfilt(double(w(:,i)),1,5);

    try
        wtmpi=wtmp(~isnan(wtmp));
        wscale=min(wtmpi)-0.3:0.025:max(wtmpi)+0.3;
        
        [k, ~]=findnearest_JACK(wtmp,wtmpi(end-3:end));
        [x, xi]=polyfit(wtmp(k),mab(k),1);
        z=polyval(x,wscale);
        
        [kk, ~]=findnearest_JACK(wtmp,wtmpi(end-6:end));
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
        [kk, ~]=findnearest_JACK(wtmp,wtmpi(1:1+3));
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
    ii=max((findnearest_JACK(oc40n.dn,dn(i))));
    subplot(1,2,1);
    plot(oc40n.w5(:,ii),oc40n.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
    plot(oc40n.w5(:,ii),oc40n.mab,'.-','color',[0.8 0.8 0.8],'markersize',15);
    hold on
    plot(wtmpii,mab,'.-k','markersize',15);
    yaxis(mab(1),mab(end));
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end

% For an end check on extrapolations
if 0
    i=11111;
    ii=max((findnearest_JACK(oc40s.A.datenum,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(oc40s.A.w5(:,ii),oc40s.A.zab,'.','color',[0.8 0.8 0.8],'markersize',40);
    hold on
    plot(vert(:,i),mab,'.k','markersize',15);
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'position',[0.277692317782113         0.941348973607038         0.492307692307692        0.0388563049853372]);
    yaxis(mab(1),mab(end));
    xaxis(nanmin(north(:,i)-0.05),nanmax(north(:,i)+0.05));
end

%% Showing why we need 5 min filter 
dn1=datenum(2017,9,10,0,0,00);
dn2=datenum(2017,9,10,12,0,00);

close all 
figure('position',[ 1         336        1415         468]); 
map=flipud(cbrewer('div','RdYlBu',40));
subplot(2,1,1);
pcolorjw(dn,mab,u); colormap(gca,map); colorbar; caxis([-0.5 0.5]);xaxis(dn1,dn2);
subplot(2,1,2);
pcolorjw(dn,mab,east); colormap(gca,map); colorbar; caxis([-0.5 0.5]);xaxis(dn1,dn2);


%% Get rid of bad interpolations by high pass filter - get rid of everything less than 5 mins
east_good=pl33tn(east,1/60,5/60)'; 
EAST=fillmissing(east_good,'nearest');

north_good=pl33tn(north,1/60,5/60)'; 
NORTH=fillmissing(north_good,'nearest');

vert_good=pl33tn(vert,1/60,5/60)'; 
VERT=fillmissing(vert_good,'nearest');

%% 
dn1=datenum(2017,9,10,0,0,00);
dn2=datenum(2017,9,10,12,0,00);

map=flipud(cbrewer('div','RdYlBu',40));
% close all 
figure; 
subplot(2,1,1); pcolorjw(dn,mab,u);caxis([-0.5 0.5]);xaxis(dn1,dn2);
subplot(2,1,2); pcolorjw(dn,mab,EAST);caxis([-0.5 0.5]);xaxis(dn1,dn2);
colormap(map)

%%  Nan out data above surface 

Nz = 25;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig:dsig:1)' ; % sigma 


%  for the bottom depth, filter to keep variability with periods of 1 hour
%  and longer
tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');

D=interp1(tide.dn,pl33tn(tide.wl+40,6/60,1),dn);
%  do a harmonic fit to estimate depth when we don't have data
[TS,XO]=t_tide(D,'interval',1/60,'output','none') ;
nn = find(~isnan(D)) ;
n1 = nn(1) ;
D(1:n1) = XO(1:n1) + D(n1)-XO(n1) ;
n1 = nn(end) ;
D(n1:end) = XO(n1:end) + D(n1)-XO(n1) ;
D=fillmissing(D,'linear');
H=sig*D;

%% NAN out the velocity data above surface 
 
z=meshgrid(mab,dn)';
surf=meshgrid(D,mab);

maskjm=double(z<(surf));
maskjm(maskjm==0)=nan;

EAST(isnan(maskjm))=nan;
NORTH(isnan(maskjm))=nan;
VERT(isnan(maskjm))=nan;


%% Let's get the data into the same sigma coordinate system

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
    dn2=datenum(2017,9,20);
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




%% Clear variables except those needed 
dn_sig=repmat(dn,length(sig),1);

clearvars -except dn mab EAST NORTH VERT H sig Usig Vsig Wsig anadir dn_sig D

if 0
    save([anadir 'oc40svel_zandsig.mat']);
end 






