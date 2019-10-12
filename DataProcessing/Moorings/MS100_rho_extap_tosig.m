% Written by Jack McSweeney 
% September 21, 2018 

% In this script we extrapolate the temp to the surface and bottom, 
% interpolate the onto sigma coordinates, and then sort ro get background 
% stratification  

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

%mooring data 
moordir='/Volumes/InnerShelf1/JackProcessing/Level2/';

% coastline
load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalcoast.mat');

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/';

%% Download Data 
MS100=load([moordir 'MS100-T_30s.mat']);


%% Transform to uniform grid in z 
dn=datenum(2017,9,8,0,0,0):1/(24*60):datenum(2017,11,3,0,0,0);
mab=(0:3:102)'; 

Temp=interp2(MS100.dn,MS100.mab,MS100.temp,dn,mab);
temp(1:length(mab),1:length(dn))=nan;

%% Extrapolate currents to surface and bottom - EAST 

for i = 1:length(dn)
    if rem(i,1000) == 0
        disp([num2str(i) ' of ' num2str(length(dn))]);
    end
    %Extrapolate to surface
    utmp=Temp(:,i);
    try
        utmpi=utmp(~isnan(utmp));
        uscale=min(utmpi)-2:0.25:max(utmpi)+2;
        
        [k, ~]=findnearest_JACK(utmp,utmpi(end-2:end));
        [x, xi]=polyfit(utmp(k),mab(k),1);
        z=polyval(x,uscale);
        
%         [kk, ~]=findnearest_JACK(utmp,utmpi(end-4:end));
%         [x2, x2i]=polyfit(utmp(kk),mab(kk),3);
%         z2=polyval(x2,uscale);
        
        utmpii=utmp;
        
%         diff1=abs(utmpi(end)-interp1(z,uscale,mab(mab>max(mab(k)))));
%         diff2= abs(utmpi(end)-interp1(z2,uscale,mab(mab>max(mab(k)))));
%         
%         if  nanmean(diff1) < nanmean(diff2)
            utmpii(mab>max(mab(k)))=interp1(z,uscale,mab(mab>max(mab(k))));
%         elseif nanmean(diff1) > nanmean(diff2)
%             utmpii(mab>max(mab(k)))=interp1(z2,uscale,mab(mab>max(mab(k))));
%         end
%         
        % clunky way to impose zero gradient at the top
        utmpii(end)=utmpii(end-1);
        
        %Extrapolate to bottom
        [kk, ~]=findnearest_JACK(utmp,utmpi(1:1+2));
        [xx, xxi]=polyfit(utmp(kk),mab(kk),1);
        zz=polyval(xx,uscale);
        utmpii(mab<min(mab(kk)))=interp1(zz,uscale,mab(mab<min(mab(kk))));
        % clunky way to impose zero gradient at the bottom
        utmpii(1)=utmpii(2);
        temp(:,i)=utmpii;
        
%         clear umtp utmpi k x z utmpii kk xx zz
        
    catch
        temp(:,i)=nan;
    end
    if i==length(dn)
        disp('Done with Temp!')
    end 
end

%% Plotting to check 
% % % % For problem solving the extrapolations  - need to uncomment clear
% % % % statement above
if 0
close all 
figure('position',[20    94   593   711]);
ii=max((findnearest_JACK(MS100.dn,dn(i))));
subplot(1,2,1);
plot(MS100.temp(:,ii),MS100.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
hold on 
plot(temp(:,i),mab,'.-k','markersize',40);
plot(utmp(k),mab(k),'.r','markersize',25);
plot(uscale,z,'.m','markersize',10);
plot(uscale,z2,'.b','markersize',10);
plot(utmpii,mab,'.g','markersize',14);
plot(utmp(kk),mab(kk),'.r','markersize',25);
plot(uscale,zz,'.m','markersize',10);
yaxis(mab(1),mab(end));


subplot(1,2,2);
plot(MS100.temp(:,ii),MS100.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
hold on 
plot(utmpii,mab,'.k','markersize',15);
yaxis(mab(1),mab(end));
leg=legend('original','extrapolated and interpolated','location','northoutside');
set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end 

% For an end check on extrapolations 
if 1
    i=79455;
        ii=max((findnearest_JACK(MS100.dn,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(MS100.temp(:,ii),MS100.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
    hold on
    plot(temp(:,i),mab,'.k','markersize',15);
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'position',[0.277692317782113         0.941348973607038         0.492307692307692        0.0388563049853372]);
    yaxis(mab(1),mab(end));
    xaxis(nanmin(temp(:,i)-1),nanmax(temp(:,i)+1));
end

%% 

close all 
figure; 
subplot(2,1,1);
pcolorjw(MS100.dn,MS100.mab,MS100.temp); colorbar 
xaxis(dn(1),dn(end));
yaxis(mab(1),mab(end));

subplot(2,1,2);
pcolorjw(dn,mab,temp); colorbar 
xaxis(dn(1),dn(end));
yaxis(mab(1),mab(end));



%% Get rid of bad data
 
ii=findnearest_JACK(dn,737000.91825156); % 01-Nov-2017 22:02:16

temp(:,ii:end)=nan;



%% Clear variables except those needed 

clearvars -except dn mab MS100 temp anadir


%% Load BT tide to help with sigma coordinates 
tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');
ms100=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/MS100-A/dpl1/mat/adcp_clean.mat');


%% Let's get the data into the same sigma coordinate system

Nz = 25 ;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig/2:dsig:1)' ; % sigma 


%  for the bottom depth, filter to keep variability with periods of 1 hour
%  and longer

D=interp1(tide.dn,tide.wl+nanmean(ms100.A.depth),dn);
D=fillmissing(D,'linear');
D=pl33tn(D,1/60,1)';
H=sig*D;

Tsig = NaN*ones(length(sig),length(dn)) ;


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
        Tsig(ii,i)=nanmean(temp(kk,i));
    end
end
%     Tsig=fillmissing(Tsig,'linear');
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
    pcolorjw(dn,mab,temp);
    colormap(map); colorbar; caxis([10 23]);
    xaxis(dn1,dn2);
    subplot(2,1,2);
    pcolorjw(repmat(dn,length(sig),1),H,Tsig);
    colormap(map); colorbar; caxis([10 23]);
    xaxis(dn1,dn2);
end

%% Calculate Salinity 

sal=nanmean(MS100.sal); % 33.4589434431925
Psig = ones(size(sal))*D - H ; % pressure 
P = ones(size(sal))*D - mab ; % pressure 

rho=sw_dens(0*temp+sal,temp,P);
Rho_sig=sw_dens(0*Tsig+sal,Tsig,Psig) ;

dn_sig=repmat(dn,length(sig),1);



%% Clear variables except those needed 

clearvars -except dn MS100 mab temp sal rho H sig Tsig Rho_sig dn_sig anadir D

if 0
    save([anadir 'MS100TempRho_zandsig.mat']);
end 



