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
OC32N=load([moordir 'OC32N-T_30s.mat']);


%% Transform to uniform grid in z 
dn=datenum(2017,9,6,0,0,0):1/(24*60):datenum(2017,11,3,0,0,0);
mab=(0:0.5:34)'; 

Temp=interp2(OC32N.dn,OC32N.mab,OC32N.temp,dn,mab);
temp(1:length(mab),1:length(dn))=nan;

%% Extrapolate Temperature to surface and bottom 

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

% % one bad nan for some reason;  fixing manually
% for i=63371
%     for ii=1:length(mab)
%         tmp=fillmissing(temp(ii,i-1:i+1),'linear');
%         temp(ii,i)=tmp(2);
%     end
% end

%% 

close all 
figure; 
subplot(2,1,1);
pcolorjw(OC32N.dn,OC32N.mab,OC32N.temp); colorbar 
xaxis(dn(1),dn(end));
yaxis(mab(1),mab(end));

subplot(2,1,2);
pcolorjw(dn,mab,temp); colorbar 
xaxis(dn(1),dn(end));
yaxis(mab(1),mab(end));

%% Nan out bad data 

temp(:,1:1375)=nan;

%% Get rid of data above surfave 

tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');

Nz = 25 ;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig:dsig:1)' ; % sigma 

D=interp1(tide.dn,tide.wl+32,dn);
D=fillmissing(D,'linear');
D=pl33tn(D,1/60,1)';
% H=sig*D;

for i=1:length(D)
    H(:,i)=0:D(i)/24:D(i);   
end

%% NAN out the velcity data above surface 
 
z=meshgrid(mab,dn)';
surf=meshgrid(D,mab);

maskjm=double(z<(surf));
maskjm(maskjm==0)=nan;

temp(isnan(maskjm))=nan;

%% Let's get the data into the same sigma coordinate system


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

%% Calculate Salinity 
sal= 33.4188991428255;
Psig = ones(size(sal))*D - H ; % pressure 
P = ones(size(sal))*D - mab ; % pressure 

rho=sw_dens(0*temp+sal,temp,P);
Rho_sig=sw_dens(0*Tsig+sal,Tsig,Psig) ;

dn_sig=repmat(dn,length(sig),1);
mab_sig=mab*D;

%% Clear variables except those needed 

clearvars -except dn OC32N mab temp sal rho sig Tsig Rho_sig dn_sig anadir H D

if 0
    save([anadir 'OC32N_TempRho_zandsig.mat']);
end 


