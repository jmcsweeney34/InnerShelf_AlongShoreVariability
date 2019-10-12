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
veldir='/Volumes/InnerShelf1/JackProcessing/ADCPS/';

% coastline
load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalcoast.mat');

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/';

%% Download Data 
OC25NA=load([moordir 'OC25NA-T_30s.mat']);


%% Transform to uniform grid in z 
dn=datenum(2017,9,6,0,0,0):1/(24*60):datenum(2017,11,3,0,0,0);
mab=(0:0.5:26.5)'; 

Temp=interp2(OC25NA.dn,OC25NA.mab,OC25NA.temp,dn,mab);
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
        kk=unique(kk);
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

% a few missing bits near the surface.  Not sure why, but this is a quick
% fix 
temp=fillmissing(temp,'nearest');

% one bad nan for some reason;  fixing manually
for i=63371
    for ii=1:length(mab)
        tmp=fillmissing(temp(ii,i-1:i+1),'linear');
        temp(ii,i)=tmp(2);
    end
end

%% Plotting to check 
% % % % For problem solving the extrapolations  - need to uncomment clear
% % % % statement above
if 0
close all 
figure('position',[20    94   593   711]);
ii=max((findnearest_JACK(OC25NA.dn,dn(i))));
subplot(1,2,1);
plot(OC25NA.temp(:,ii),OC25NA.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
plot(OC25NA.temp(:,ii),OC25NA.mab,'.','color',[0.8 0.8 0.8],'markersize',15);
hold on 
plot(utmpii,mab,'.k','markersize',15);
yaxis(mab(1),mab(end));
leg=legend('original','extrapolated and interpolated','location','northoutside');
set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end 

% For an end check on extrapolations 
if 1
    i=60000;
    ii=max((findnearest_JACK(OC25NA.dn,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(OC25NA.temp(:,ii),OC25NA.mab,'.','color',[0.8 0.8 0.8],'markersize',40);
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
pcolorjw(OC25NA.dn,OC25NA.mab,OC25NA.temp); colorbar 
xaxis(dn(1),dn(end));
yaxis(mab(1),mab(end));

subplot(2,1,2);
pcolorjw(dn,mab,temp); colorbar 
xaxis(dn(1),dn(end));
yaxis(mab(1),mab(end));


%% Clear variables except those needed 

clearvars -except dn mab OC25NA temp anadir veldir


%% Load BT tide to help with sigma coordinates 
oc25na=load('/Volumes/InnerShelf1/JackProcessing/ADCPS/OC25NA_1min_deploy1.mat');


%% Let's get the data into the same sigma coordinate system

Nz = 25 ;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig/2:dsig:1)' ; % sigma 


%  for the bottom depth, filter to keep variability with periods of 1 hour
%  and longer
D=interp1(oc25na.dn,oc25na.H,dn);
%  do a harmonic fit to estimate depth when we don't have data
[TS,XO]=t_tide(D,'interval',1/60,'output','none') ;
nn = find(~isnan(D)) ;
n1 = nn(1) ;
D(1:n1) = XO(1:n1) + D(n1)-XO(n1) ;
n1 = nn(end) ;
D(n1:end) = XO(n1:end) + D(n1)-XO(n1) ;
D=pl33tn(D,1/60,1)';

% H=sig*D;
for i=1:length(D)
    H(:,i)=0:D(i)/24:D(i);   
end

Tsig = NaN*ones(length(sig),length(dn)) ;

%  OK now let's interpolate onto the sigma grid.  
for i = 1:length(dn)
    if mod(i,1000) == 0
        disp([num2str(i) '/' num2str(length(dn))]) ;
    end
    
    for ii=1:length(sig)
        if ii==1
            kk=find(mab<H(ii,i));
            if isempty(kk)
                kk=1;
            end 
        else
            kk=find(mab>H(ii-1,i)& mab<H(ii,i));
        end
        Tsig(ii,i)=nanmean(temp(kk,i));
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
    pcolorjw(dn,mab,temp);
    colormap(map); colorbar; caxis([10 23]);
    xaxis(dn1,dn2);
    subplot(2,1,2);
    pcolorjw(repmat(dn,length(sig),1),H,Tsig);
    colormap(map); colorbar; caxis([10 23]);
    xaxis(dn1,dn2);
end


%% Calculate Salinity 
sal= 33.4188991428255;
Psig = ones(size(sal))*D - H ; % pressure 
P = ones(size(sal))*D - mab ; % pressure 

rho=sw_dens(0*temp+sal,temp,P);
Rho_sig=sw_dens(0*Tsig+sal,Tsig,Psig) ;

dn_sig=repmat(dn,length(sig),1);

%% NAN out the velcity data above surface 
 
z=meshgrid(mab,dn)';
surf=meshgrid(D,mab);

maskjm=double(z<(surf));
maskjm(maskjm==0)=nan;

temp(isnan(maskjm))=nan;
rho(isnan(maskjm))=nan;

%% Clear variables except those needed 

clearvars -except dn OC25NA mab temp sal rho H sig Tsig Rho_sig dn_sig anadir D

if 0
    save([anadir 'OC25NA_TempRho_zandsig.mat']);
end 


