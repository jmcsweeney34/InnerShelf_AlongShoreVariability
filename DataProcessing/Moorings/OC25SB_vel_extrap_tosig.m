% Written by Jack McSweeney 
% October 2, 2018 
% 
% In this script we extrapolate the velocities to the surface and bottom, 
% interpolate the velocites onto sigma coordinates

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

% coastline
load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalcoast.mat');

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/';

%% Download Data 
oc25sb=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/OC25SB-A/Level2/OC25SB-A.mat');

%% Check that north is true north 
% % % need to check 1/31/19

%% Transform to uniform grid in z 
dn=datenum(2017,9,7,0,0,0):1/(24*60):datenum(2017,11,3,0,0,0);
mab=(0:0.5:26.5)'; 

u=interp2(oc25sb.A.datenum,oc25sb.A.zab,oc25sb.A.u.*oc25sb.A.mask,dn,mab);
v=interp2(oc25sb.A.datenum,oc25sb.A.zab,oc25sb.A.v.*oc25sb.A.mask,dn,mab);
w=interp2(oc25sb.A.datenum,oc25sb.A.zab,oc25sb.A.w5.*oc25sb.A.mask,dn,mab);



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
ii=max((findnearest_JACK(oc25nb.datenum,dn(i))));
subplot(1,2,1);
plot(oc25nb.u(:,ii),oc25nb.zab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
plot(oc25nb.u(:,ii),oc25nb.zab,'.-','color',[0.8 0.8 0.8],'markersize',40);
hold on 
plot(utmpii,mab,'.-k','markersize',15);
yaxis(mab(1),mab(end));
leg=legend('original','extrapolated and interpolated','location','northoutside');
set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end 

% For an end check on extrapolations 
if 1
    i=70013;
    ii=max((findnearest_JACK(oc25sb.A.datenum,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(oc25sb.A.u(:,ii),oc25sb.A.zab,'.','color',[0.8 0.8 0.8],'markersize',40);
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
    ii=max((findnearest_JACK(oc25na.dn,dn(i))));
    subplot(1,2,1);
    plot(oc25na.v(:,ii),oc25na.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
    plot(oc25na.v(:,ii),oc25na.mab,'.','color',[0.8 0.8 0.8],'markersize',40);
    hold on
    plot(vtmpii,mab,'.k','markersize',15);
    yaxis(mab(1),mab(end));
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end

% For an end check on extrapolations
if 0
    i=5713;
    ii=max((findnearest_JACK(oc25nb.datenum,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(oc25nb.v(:,ii),oc25nb.zab,'.','color',[0.8 0.8 0.8],'markersize',40);
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
    ii=max((findnearest_JACK(oc25na.dn,dn(i))));
    subplot(1,2,1);
    plot(oc25na.w(:,ii),oc25na.mab,'.','color',[0.5 0.5 0.5],'markersize',15);
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
    plot(oc25na.w(:,ii),oc25na.mab,'.-','color',[0.8 0.8 0.8],'markersize',15);
    hold on
    plot(wtmpii,mab,'.-k','markersize',15);
    yaxis(mab(1),mab(end));
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'Position',[0.60076 0.94937 0.26981 0.037271]);
end

% For an end check on extrapolations
if 1
    i=11111;
    ii=max((findnearest_JACK(oc25nb.datenum,dn(i))));
    
    close all
    figure('position',[440   116   325   682]);
    plot(oc25nb.w(:,ii),oc25nb.zab,'.','color',[0.8 0.8 0.8],'markersize',40);
    hold on
    plot(vert(:,i),mab,'.k','markersize',15);
    leg=legend('original','extrapolated and interpolated','location','northoutside');
    set(leg,'position',[0.277692317782113         0.941348973607038         0.492307692307692        0.0388563049853372]);
    yaxis(mab(1),mab(end));
    xaxis(nanmin(north(:,i)-0.05),nanmax(north(:,i)+0.05));
end

%% Showing why we need 5 min filter 
dn1=datenum(2017,10,22,0,0,00);
dn2=datenum(2017,10,22,12,0,00);

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
dn1=datenum(2017,10,22,0,0,00);
dn2=datenum(2017,10,22,12,0,00);

map=flipud(cbrewer('div','RdYlBu',40));
% close all 
figure; 
subplot(2,1,1); pcolorjw(dn,mab,u);caxis([-0.5 0.5]);xaxis(dn1,dn2);
subplot(2,1,2); pcolorjw(dn,mab,EAST);caxis([-0.5 0.5]);xaxis(dn1,dn2);
colormap(map)

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

clearvars -except dn mab EAST NORTH VERT oc25sb anadir
%     save([anadir 'oc25nb_vel_z.mat']);

%% Let's get the data into the same sigma coordinate system

%%%% Also need to fix this - its not quite right 

Nz = 25 ;  % number of vertical levels
dsig = 1/Nz ; % delta sigma
sig = (dsig/2:dsig:1)' ; % sigma 


%  for the bottom depth, filter to keep variability with periods of 1 hour
%  and longer

tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');
D=interp1(tide.dn,pl33tn(tide.wl+25,6/60,1),dn);
D=fillmissing(D,'linear');

[TS,XO]=t_tide(D,'interval',1/60,'output','none') ;
nn = find(~isnan(D)) ;
n1 = nn(1) ;
D(1:n1) = XO(1:n1) + D(n1)-XO(n1) ;
n1 = nn(end) ;
D(n1:end) = XO(n1:end) + D(n1)-XO(n1) ;

for i=1:length(D)
    H(:,i)=0:D(i)/24:D(i);   
end

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
    dn1=datenum(2017,10,24);
    dn2=datenum(2017,10,26);
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
i=10300;
plot([0 0],[0 max(H(:,i))+1],'--k','linewidth',1) ;
hold on
try
    ii=max((findnearest_JACK(oc25nb.datenum,dn(i))));
    h1=plot(oc25nb.u(:,ii),oc25nb.zab,'.-','color',[0.5 0.5 0.5],'linewidth',2,'markersize',40);
    h2=plot(EAST(:,i),mab,'.-k','linewidth',2,'markersize',20);
    h3=plot(Usig(:,i),H(:,i),'.-r','linewidth',2,'markersize',20);
    legend([h1 h2 h3],'raw','z','sig');
catch
end

yaxis(0, max(H(:,i))+1);
xaxis(min(Usig(:,i))-0.1,max(Usig(:,i))+0.1);


%% NAN out the velcity data above surface 
 
z=meshgrid(mab,dn)';
surf=meshgrid(D,mab);

maskjm=double(z<(surf));
maskjm(maskjm==0)=nan;

EAST(isnan(maskjm))=nan;
NORTH(isnan(maskjm))=nan;
VERT(isnan(maskjm))=nan;

%% Clear variables except those needed 
dn_sig=repmat(dn,length(sig),1);
clearvars -except dn mab EAST NORTH VERT H sig Usig Vsig Wsig anadir dn_sig readme

if 0
    save([anadir 'oc25sb_vel_zandsig.mat']);
end 






