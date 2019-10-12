% Written by Jack McSweeney 
% September 21, 2018 

% Trying to figure out why my estimate of OC50 is different than Jim's 

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/';


%% Download mine and jim's 

load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/MS100TempRho_zandsig.mat');

%% Prep 

H0 = nanmean(D) ;
g = 9.81 ;
f = (2*2*pi/(24*60*60))*sin(34.9457*pi/180) ;

DDN = 2*12.4206/24 ;
DNskip = 1/4 ;
minlength = 1400 ;
zmin = 3 ;
zmax = 95 ;



%%  Compute background (sorted) stratification 
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
            [c,w] = nmodes(z,n2i,2*pi/(12.4206*60*60)) ; % had to modify nmodes here
            [vl,nn] = max(c) ;
            Cw(idx) = c(nn) ;
            phi = w(:,nn) ;
            phiz=interp1([z(1:end-1)+z(2:end)]/2,diff(phi),z)/dz;
%             phiz = dfdz(phi,dz) ;
            [c,w] = nmodesf(z,n2i,2*pi/(12.4206*60*60),f) ;
            [vl,nn] = max(c) ;
            Cwf(idx) = c(nn) ;
            
            N2b(idx) = nanmean(n2) ;
            
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


%% 

ii=find(DNrho==0);
DNrho(ii)=nan;
alpha(ii)=nan;
beta(ii)=nan;
Cw(ii)=nan;
Cwf(ii)=nan;
N2b(ii)=nan;
Npyc(ii)=nan;
Zpyc(ii)=nan;
Rhopyc(ii)=nan;
Rho(:,ii)=nan;


%% Cw, Alpha, Beta, N2  - Compare to jims

X=load('/Volumes/InnerShelf1/JimProcessing/OC50/OC50_BackgroundStrat_01.mat');

k=find(diff(X.DNrho)>1);


close all;
figure('position',[185   103   709   702]);
subplot(4,1,1);
h1=plot(DNrho,Cw,'-k','linewidth',2); hold on 
h3=plot(X.DNrho(1:k),X.Cw(1:k),'-g','linewidth',2); 
plot(X.DNrho(k+1:end),X.Cw(k+1:end),'-g','linewidth',2);
datetick('x','keeplimits','keepticks');
title ('Cw');
legend([h1 h3],'ms100','0C50');

subplot(4,1,2);
plot([datenum(2017,9,5) datenum(2017,11,3)],[0 0],'--k','linewidth',1) ;hold on 
h1=plot(DNrho,alpha,'-k','linewidth',2); 
h3=plot(X.DNrho(1:k),X.alpha(1:k),'-g','linewidth',2); 
plot(X.DNrho(k+1:end),X.alpha(k+1:end),'-g','linewidth',2);
datetick('x','keeplimits','keepticks');
title ('Alpha');

subplot(4,1,3);
h1=plot(DNrho,beta,'-k','linewidth',2); hold on 
h3=plot(X.DNrho(1:k),X.beta(1:k),'-g','linewidth',2); 
plot(X.DNrho(k+1:end),X.beta(k+1:end),'-g','linewidth',2);
datetick('x','keeplimits','keepticks');
title ('Beta');

subplot(4,1,4)
h1=plot(DNrho,N2b,'-k','linewidth',2); hold on 
h3=plot(X.DNrho(1:k),X.N2b(1:k),'-g','linewidth',2); 
plot(X.DNrho(k+1:end),X.N2b(k+1:end),'-g','linewidth',2);
datetick('x','keeplimits','keepticks');
title ('N^2');

%% 
jet_JACK
close all 
figure;
pcolorjw(repmat(DNrho,length(Z),1),repmat(Z,1,length(DNrho)),Rho);
hold on 
contour(DNrho(~isnan(DNrho)),Z,Rho(:,~isnan(DNrho)),[1019:.1:1026],'k') ;
title('Jim');
datetick('x','keeplimits','keepticks')
colorbar
colormap(flipud(map))




%% 
clearvars -except DNrho Rho N2b Cw Cwf Z Zpyc Npyc alpha beta anadir Rhopyc

if 0
    save([anadir 'ms100_alpha_beta_Cw.mat']);
end 







