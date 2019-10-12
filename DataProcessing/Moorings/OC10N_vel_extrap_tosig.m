% Written by Jack McSweeney 
% October 5, 2018 
% 
% In this script we extrapolate the velocities to the surface and bottom, 
% interpolate the velocites onto sigma coordinates


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
oc10=load('/Volumes/InnerShelf1/Moorings/LanderADCPs/OC10ADCP/OC10N_A_ADCP_proc.mat');


%% Check that north is true north 
% % Jim rotated the currents according to the tidal and subtidal currents


%% Remove bad data  
% % Bad data was nanned out in processing 

%% Filter data from 15s to 60s 
dnmin=oc10.dn(1):60/(60*60*24):oc10.dn(end);
u_tmp=pl66tn(oc10.u,15/3600,60/3600)'; %1 min filter, 15 second data 
u=interp2(meshgrid(oc10.dn,oc10.z),meshgrid(oc10.z,oc10.dn)',u_tmp,meshgrid(dnmin,oc10.z),meshgrid(oc10.z,dnmin)');

v_tmp=pl66tn(oc10.v,15/3600,60/3600)'; %1 min filter, 15 second data 
v=interp2(meshgrid(oc10.dn,oc10.z),meshgrid(oc10.z,oc10.dn)',v_tmp,meshgrid(dnmin,oc10.z),meshgrid(oc10.z,dnmin)');

w_tmp=pl66tn(oc10.w,15/3600,60/3600)'; %1 min filter, 15 second data 
w=interp2(meshgrid(oc10.dn,oc10.z),meshgrid(oc10.z,oc10.dn)',w_tmp,meshgrid(dnmin,oc10.z),meshgrid(oc10.z,dnmin)');


%% Transform to uniform grid in z 
dn=datenum(2017,9,9,0,0,0):1/(24*60):datenum(2017,10,11,0,0,0);
mab=(0:0.25:12.5)'; 

u=interp2(dnmin,oc10.z,u,dn,mab);
v=interp2(dnmin,oc10.z,v,dn,mab);
w=interp2(dnmin,oc10.z,w,dn,mab);

east(1:size(u,1),1:size(u,2))=nan;
north(1:size(v,1),1:size(v,2))=nan;
vert(1:size(v,1),1:size(v,2))=nan;
%% Prep for extrapolation and sigma coordinates 
Nz = 25 ;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig/2:dsig:1)' ; % sigma 


%  for the bottom depth, filter to keep variability with periods of 1 hour
%  and longer

D=interp1(oc10.dn,oc10.D,dn);


% D=interp1([oc25na.dn oc25na_b.dn],[oc25na.H oc25na_b.H-.2],dn);
[TS,XO]=t_tide(D,'interval',1/60,'output','none') ;
nn = find(~isnan(D)) ;
n1 = nn(1) ;
D(1:n1) = XO(1:n1) + D(n1)-XO(n1) ;
n1 = nn(end) ;
D(n1:end) = XO(n1:end) + D(n1)-XO(n1) ;

%fill the gap (figured out indices manually 
D=pl33tn(D,1/60,1)';

% H=sig*D;
for i=1:length(D)
    H(:,i)=0:D(i)/24:D(i);   
end
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

%% NAN out the velcity data above surface 
 
z=meshgrid(mab,dn)';
surf=meshgrid(D,mab);

maskjm=double(z<(surf));
maskjm(maskjm==0)=nan;

EAST(isnan(maskjm))=nan;
NORTH(isnan(maskjm))=nan;
VERT(isnan(maskjm))=nan;

%% Showingfinal product
dn1=datenum(2017,9,10,0,0,00);
dn2=datenum(2017,9,10,12,0,00);

close all 
figure('position',[ 1         336        1415         468]); 
map=flipud(cbrewer('div','RdYlBu',40));
subplot(2,1,1);
pcolorjw(dn,mab,u); colormap(gca,map); colorbar; caxis([-0.5 0.5]);xaxis(dn1,dn2);
subplot(2,1,2);
pcolorjw(dn,mab,EAST); colormap(gca,map); colorbar; caxis([-0.5 0.5]);xaxis(dn1,dn2);


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
            if isempty(kk)
                kk=1;
            end 
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
%     close all
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
clearvars -except dn mab EAST NORTH VERT H sig Usig Vsig Wsig anadir dn_sig readme D

if 0
    save([anadir 'oc10N_vel_zandsig.mat']);
end 


