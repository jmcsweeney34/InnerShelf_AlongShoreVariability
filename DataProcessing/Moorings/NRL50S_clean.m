% Written by Jack McSweeney 
% June 5, 2019 

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/';

%% Download data 

NRL50S=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/NRL50S/NRL50S_T.mat');
nrl50s=load('/Volumes/InnerShelf1/Moorings/NPS_SIO_Moorings/NRL50S/NRL50S_1min.mat');

%% Fix NRL temp - 50S  

NRL50S.dn=datenum(2017,9,6):60/(24*60*60):datenum(2017,11,2);
NRL50S.mab=flipud(NRL50S.T.zab);

NRL50S.temp(1:length(NRL50S.mab),1:length(NRL50S.dn))=nan;
NRL50S.temp(1,:)=interp1(datenum(NRL50S.T.datenum.z2),NRL50S.T.temp.z2,NRL50S.dn);
NRL50S.temp(2,:)=interp1(datenum(NRL50S.T.datenum.z3),NRL50S.T.temp.z3,NRL50S.dn);
NRL50S.temp(3,:)=interp1(datenum(NRL50S.T.datenum.z6),NRL50S.T.temp.z6,NRL50S.dn);
NRL50S.temp(4,:)=interp1(datenum(NRL50S.T.datenum.z9),NRL50S.T.temp.z9,NRL50S.dn);
NRL50S.temp(5,:)=interp1(datenum(NRL50S.T.datenum.z12),NRL50S.T.temp.z12,NRL50S.dn);
NRL50S.temp(6,:)=interp1(datenum(NRL50S.T.datenum.z15),NRL50S.T.temp.z15,NRL50S.dn);
NRL50S.temp(7,:)=interp1(datenum(NRL50S.T.datenum.z18),NRL50S.T.temp.z18,NRL50S.dn);
NRL50S.temp(8,:)=interp1(datenum(NRL50S.T.datenum.z21),NRL50S.T.temp.z21,NRL50S.dn);
NRL50S.temp(9,:)=interp1(datenum(NRL50S.T.datenum.z24),NRL50S.T.temp.z24,NRL50S.dn);
NRL50S.temp(10,:)=interp1(datenum(NRL50S.T.datenum.z27),NRL50S.T.temp.z27,NRL50S.dn);
NRL50S.temp(11,:)=interp1(datenum(NRL50S.T.datenum.z30),NRL50S.T.temp.z30,NRL50S.dn);
NRL50S.temp(12,:)=interp1(datenum(NRL50S.T.datenum.z33),NRL50S.T.temp.z33,NRL50S.dn);
NRL50S.temp(13,:)=interp1(datenum(NRL50S.T.datenum.z36),NRL50S.T.temp.z36,NRL50S.dn);
NRL50S.temp(14,:)=interp1(datenum(NRL50S.T.datenum.z39),NRL50S.T.temp.z39,NRL50S.dn);
NRL50S.temp(15,:)=interp1(datenum(NRL50S.T.datenum.z42),NRL50S.T.temp.z42,NRL50S.dn);
NRL50S.temp(16,:)=interp1(datenum(NRL50S.T.datenum.z50),NRL50S.T.temp.z50,NRL50S.dn);

% % clever way to identify good data; but want bad data haha 
% i=(NRL50S.dn>datenum(2017,9,9,18,1,0)).*(NRL50S.dn<datenum(2017,10,31,3,0,0));
% i(i==0)=nan;

i=find(NRL50S.dn<datenum(2017,9,10,19,55,0) | NRL50S.dn>datenum(2017,10,29,19,46,0));

NRL50S.dn(i)=[];
NRL50S.temp(:,i)=[];


%% Get rid of bad ADCP data - NRL50S 
i=find(nrl50s.A.datenum<datenum(2017,9,10,20,5,0) | nrl50s.A.datenum>datenum(2017,10,29,19,40,0));

nrl50s.A.datenum(i)=[];
nrl50s.A.u(:,i)=[];
nrl50s.A.v(:,i)=[];
nrl50s.A.w5(:,i)=[];

perbad=10;
i=find(nrl50s.A.zab>(max(nrl50s.A.zab)-max(nrl50s.A.zab)*perbad/100));
nrl50s.A.u(i,:)=nan;
nrl50s.A.v(i,:)=nan;
nrl50s.A.w5(i,:)=nan;

%% %% Transform to uniform grid in z 

dn=datenum(2017,9,10,0,0,0):1/(24*60):datenum(2017,10,30,0,0,0);
mab=(0:0.5:52)'; 

u=interp2(nrl50s.A.datenum,nrl50s.A.zab,nrl50s.A.u,dn,mab);
v=interp2(nrl50s.A.datenum,nrl50s.A.zab,nrl50s.A.v,dn,mab);
w=interp2(nrl50s.A.datenum,nrl50s.A.zab,nrl50s.A.w5,dn,mab);
t=interp2(NRL50S.dn,NRL50S.mab,NRL50S.temp,dn,mab);

%%  Set up sigma coordinates 

Nz = 50 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
sig = (0:dsig:1)' ; % sigma 

tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');
D=interp1(tide.dn,tide.wl+50,dn);
D=pl33tn(D,1/60,2)';
D=fillmissing(D,'linear');
H=sig*D;

dn_sig=repmat(dn,length(sig),1);
Usig(1:length(sig),1:length(dn))=nan;
Vsig(1:length(sig),1:length(dn))=nan;
Wsig(1:length(sig),1:length(dn))=nan;
Tsig(1:length(sig),1:length(dn))=nan;

%% Put the data in sigma coordinates 

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
        Tsig(ii,i)=nanmean(t(kk,i));
        Usig(ii,i)=nanmean(u(kk,i));
        Vsig(ii,i)=nanmean(v(kk,i));
        Wsig(ii,i)=nanmean(w(kk,i));
     end
end

%% Extrapolate to surface and bed 


for i =1:length(dn)
    if mod(i,1000) == 0
        disp([num2str(i) '/' num2str(length(dn))]) ;
    end
    
    utmp=Usig(:,i); 
    vtmp=Vsig(:,i); 
    wtmp=Wsig(:,i);
    ttmp=Tsig(:,i);  
    ztmp=H(:,i);
    
    % surface
    ks=find(~isnan(utmp),1,'last');
    wtmp(end-1:end)=0;
    wtmp(ks-3:end)=fillmissing(wtmp(ks-3:end),'pchip');
    vtmp(ks+1:end)=vtmp(ks);
    utmp(ks+1:end)=utmp(ks);

    if isnan(ttmp(ks))
        if ~isnan(nanmean(ttmp))
        ttmpi=fillmissing(ttmp(ks-5:ks),'linear');
        [xx,~]=polyfit(ttmpi,ztmp(ks-5:ks),2);
        tscale=min(ttmp)-2:0.25:max(ttmp)+2;
        zz=polyval(xx,tscale);
        ttmp(ks+1:end-1)=interp1(zz,tscale,ztmp(ks+1:end-1));
        ttmp(end)=ttmp(end-1);
        end 
    end 

    % bed
    kb=find(~isnan(utmp),1);
    wtmp(1:2)=0;
    wtmp(1:kb+3)=fillmissing(wtmp(1:kb+3),'pchip');
    vtmp(1:kb-1)=vtmp(kb);
    utmp(1:kb-1)=utmp(kb);

    kb=find(~isnan(ttmp),1);
    [xx, ~]=polyfit(ttmp(kb:kb+3),ztmp(kb:kb+3),1);
    tscale=min(ttmp)-1:0.25:max(ttmp)+1;
    zz=polyval(xx,tscale);
    ttmp(2:kb-1)=interp1(zz,tscale,ztmp(2:kb-1));
    ttmp(1)=ttmp(2);
%     
    
    % plot to show the extrapolations 
    if 0
        close all;
        figure('position',[125   396   845   399]);
        
        subplot(141)
        plot(t(:,i),mab,'.k','markersize',20);
        hold on
        plot(Tsig(:,i),ztmp,'.-r','markersize',15,'linewidth',2);
        plot(ttmp,ztmp,'.g','markersize',10);
        yaxis(0,51); title('temp');
        
        subplot(142)
        plot(u(:,i),mab,'.k','markersize',20);
        hold on
        plot(Usig(:,i),ztmp,'.-r','markersize',15,'linewidth',2);
        plot(utmp,ztmp,'.g','markersize',10);
        yaxis(0,51); title('u');
        
        subplot(143)
        plot(v(:,i),mab,'.k','markersize',20);
        hold on
        plot(Vsig(:,i),ztmp,'.-r','markersize',15,'linewidth',2);
        plot(vtmp,ztmp,'.g','markersize',10);
        yaxis(0,51); title('v');
        
        subplot(144)
        plot(w(:,i),mab,'.k','markersize',20);
        hold on
        plot(Wsig(:,i),ztmp,'.-r','markersize',15,'linewidth',2);
        plot(wtmp,ztmp,'.g','markersize',10);
        yaxis(0,51); title('w');
    end
    
    Usig(:,i)=fillmissing(utmp,'linear');
    Vsig(:,i)=fillmissing(vtmp,'linear');
    Wsig(:,i)=fillmissing(wtmp,'linear');
    Tsig(:,i)=fillmissing(ttmp,'linear');
    
end


%% figures to check the extrapolation 

% Temp 
if 0
    close all
    figure('position',[ 66         380        1139         416]);
    subplot(211);
    pcolorjw(dn,mab,t);
    colorbar; caxis([7 18]);
    
    subplot(212);
    pcolorjw(dn_sig,H,Tsig);
    colorbar; caxis([7 18]);
end

% u 
if 1
    map=flipud(cbrewer('div','RdBu',40));
    close all
    figure('position',[ 66         380        1139         416]);
    subplot(211);
    pcolorjw(dn,mab,u);
    colorbar; caxis([-0.5 0.5]); colormap(map);
    
    subplot(212);
    pcolorjw(dn_sig,H,Usig);
    colorbar; caxis([-0.5 0.5]); colormap(map);
end

% v
if 0
    map=flipud(cbrewer('div','PuOr',40));
    close all
    figure('position',[ 66         380        1139         416]);
    subplot(211);
    pcolorjw(dn,mab,v);
    colorbar; caxis([-0.5 0.5]); colormap(map);
    
    subplot(212);
    pcolorjw(dn_sig,H,Vsig);
    colorbar; caxis([-0.5 0.5]); colormap(map);
end

% w
if 0
    map=flipud(cbrewer('div','BrBG',40));
    close all
    figure('position',[ 66         380        1139         416]);
    subplot(211);
    pcolorjw(dn,mab,w);
    colorbar; caxis([-0.1 0.1]); colormap(map);
    
    subplot(212);
    pcolorjw(dn_sig,H,Wsig);
    colorbar; caxis([-0.1 0.1]); colormap(map);
end

%%  Estimate the subtidal field 

H0 = nanmean(nanmax(H));
g = 9.81 ;
f = (2*2*pi/(24*60*60))*sin(34.9457*pi/180) ;

DDN = 2*12.4206/24 ;
DNskip = 1/4 ;
minlength = 1400 ;
zmin = 3 ;
zmax = 45 ;

sal=33.4188991428255;
Psig = ones(size(sal))*nanmax(H) - H ; % pressure 
P = ones(size(sal))*nanmax(H) - mab ; % pressure 

rho=sw_dens(0*t+sal,t,P);
Rho_sig=sw_dens(0*Tsig+sal,Tsig,Psig) ;

om=2*pi/(12.4206*60*60);
dw=om/1000;

DN1 = dn(1) ;
DN2 = DN1+DDN ;
idx = 0 ;

NZ = 100 ;
dZ = H0/NZ ;
Z = (dZ/2:dZ:H0)' ;

% figure(1); set(figure(1),'position',[9   318   687   487]);
% figure(2); set(figure(2),'position',[ 715   320   560   420]);
while DN2<=dn(end)
    idx = idx+1 ;
    nn = find((dn>=DN1).*(dn<=DN2)) ;
    dntmp = dn(nn) ;
    first(idx)=DN1;
    disp(datestr(DN1));
    
    DNrho(idx)=nan;
    Rho(1:length(Z),idx)=nan;
    Cw(idx)=nan;
    Cwf(idx)=nan;
    Cg(idx)=nan;
    alpha(idx)=nan;
    beta(idx)=nan;
    Zpyc(idx)=nan;
    Rhopyc(idx)=nan;
    
    if ~isempty(nn)
        R = Rho_sig(:,nn) ;
        if (sum(sum(isnan(R)))==0)&&(length(nn)>minlength)
            DNrho(idx) = nanmean(dntmp) ; 
            %  sort densities to obtain background profile
            RS = [] ;
            for iz = 1:size(R,1)
                RS = [RS R(iz,:)] ;
            end
            %  I need to dither to get rid of common values
            RS = RS + 1e-5*randn(size(RS)) ;
            RS = RS + 1e-5*randn(size(RS)) ;
            RS = RS + 1e-5*randn(size(RS)) ;
            RS = RS + 1e-5*randn(size(RS)) ;
            [RS,mm] = sort(RS,'descend') ;
            NRS = length(RS) ;
            dz = H0/NRS ;
            ZRS = (dz/2:dz:H0)' ;
            Rho(:,idx) = interp1(ZRS,RS,Z) ;   
            
            %  Calculate the wave speed
            r = Rho(:,idx) ;
            dz = nanmean(diff(Z)) ;
%             N2 = abs(-g*dfdz(r,dz))/nanmean(r) ;
            drhodz=interp1([Z(1:end-1)+Z(2:end)]/2,diff(r),Z)/dz;
            N2 = abs(-g*drhodz)/nanmean(r) ;
            z = (dz:dz:H0-dz)' ;
            n2 = interp1(Z,N2,z,'linear','extrap') ;
            nn0 = find((z>=zmin).*(z<=zmax)) ;
            ztmp = z(nn0) ;
            rtmp = r(nn0) ;
            ntmp = sqrt(n2(nn0)) ;
            [Npyc(idx),nn0] = max(ntmp) ;
            Zpyc(idx) = ztmp(nn0) ;
            Rhopyc(idx) = rtmp(nn0) ;
            
            n2i=fillmissing(n2,'nearest');
            [c,ww] = nmodes(z,n2i,om) ; % had to modify nmodes here
            [vl,nn] = max(c) ;
            Cw(idx) = c(nn) ;
            phi = ww(:,nn) ;
            phiz=interp1([z(1:end-1)+z(2:end)]/2,diff(phi),z)/dz;
%             phiz = dfdz(phi,dz) ;
            [c,ww] = nmodesf(z,n2i,om,f) ;
            [vl,nn] = max(c) ;
            Cwf(idx) = c(nn) ;
            k1(idx)=om/Cwf(idx);
            N2b(idx) = nanmean(n2) ;
            
            %Estimate group speed with rotation 
            [cp_1,~] = nmodesf(z,n2i,om+dw,f) ; % had to modify nmodes here
            [~,nnb] = max(cp_1);
            Cp_1(idx)=cp_1(nnb);
            kp(idx)=(om+dw)/Cp_1(idx);
            
            [cp_2,~] = nmodesf(z,n2i,om-dw,f) ; % had to modify nmodes here
            [~,nnb] = max(cp_2);
            Cp_2(idx)=cp_2(nnb);
            km(idx)=(om-dw)/Cp_2(idx);
            
            Cg_1(idx)=dw./(kp(idx)-k1(idx));
            Cg_2(idx)=dw./(k1(idx)-km(idx));
            Cg(idx)=0.5*(Cg_1(idx)+Cg_2(idx));
            
            %  now calculate KdV parameters for the stratification
            alpha(idx) = -(3/2)*Cw(idx)*nansum(phiz.^3)/nansum(phiz.^2) ;
            beta(idx) = (1/2)*Cw(idx)*nansum(phi.^2)/nansum(phiz.^2) ;
        end       
    end
    
    if 0
        try
            jet_JACK
            figure(1);
            clf
            pcolorjw(R);    hold on
            colorbar
            title([datestr(DN1) 'Rho Me'])
            caxis([min(min(R)) max(max(R))]); colormap(map)
            
      
        catch
        end
        
        try
            figure(2);
            clf
            subplot(1,3,1);
            plot(RS,ZRS,'k','linewidth',2)
            title('Sorted Rho profile')
            
            subplot(1,3,2);
            plot(phi,z,'k','linewidth',2)
            title([datestr(DN1) char(10) ' phi'])
            
            subplot(1,3,3);
            plot(phiz,z,'k','linewidth',2)
            title('phiz')
            annotation(figure(2),'textbox',[0.75993 0.84048 0.12936 0.064286],...
                'String',['Alpha' char(10) num2str(alpha(idx))],'color','r',...
                'LineStyle','none','FitBoxToText','off','fontweight','bold');
            
            
            annotation(figure(2),'textbox',[0.69743 0.17857 0.12936 0.064286],...
                'String',['Beta' char(10) num2str(beta(idx))],'color','r',...
                'LineStyle','none','FitBoxToText','off','fontweight','bold');
            
        catch
            
        end
    end
    
    DN1 = DN1+DNskip;
    DN2 = DN1+DDN ;
    
%     pause
end

ii=find(~isnan(DNrho));
subtidal.DNrho=DNrho(ii);
subtidal.Rho= Rho(:,ii); 
subtidal.Cw=Cw(ii); 
subtidal.Cwf=Cwf(ii);
subtidal.Cg =Cg(ii);
subtidal.Z=Z;
subtidal.alpha=alpha(ii); 
subtidal.beta=beta(ii);

% pcolorjw(subtidal.DNrho,subtidal.Z,subtidal.Rho)

%% 
if 0
    save([anadir 'NRL50S.mat'],'Usig','Vsig','Wsig','Tsig','Rho_sig','Psig',...
        'dn_sig','sig','H','u','v','w','t','dn','mab','subtidal','-v7.3');
    disp('NRL50S.mat saved!');
end


