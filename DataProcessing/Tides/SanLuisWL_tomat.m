clear 

buoydir='/Volumes/InnerShelf1/NOAA_Buoy_Data/';

k=csvimport([buoydir 'SanLuis_9412110/CO-OPS__9412110__wl_msl_Sept17.csv']);

dn_tide=k(2:end,1);
wlev_tide=k(2:end,2);
T(1:length(dn_tide),1)=nan;
WL(1:length(dn_tide),1)=nan;
for i=1:length(dn_tide)
    T(i)=datenum(cell2mat(dn_tide(i)));
    if ~isempty((cell2mat(wlev_tide(i))))
        WL(i)=(cell2mat(wlev_tide(i)));
    end
end
tide1.dn=T;
tide1.wl=WL;

kk=find(tide1.wl==0);
tide1.wl(kk)=nan;

clear k fname kk wlev_tide T WL i dn_tide

k=csvimport([buoydir 'SanLuis_9412110/CO-OPS__9412110__wl_msl_Oct2017.csv']);

dn_tide=k(2:end,1);
wlev_tide=k(2:end,2);
T(1:length(dn_tide),1)=nan;
WL(1:length(dn_tide),1)=nan;
for i=1:length(dn_tide)
    T(i)=datenum(cell2mat(dn_tide(i)));
    if ~isempty((cell2mat(wlev_tide(i))))
        WL(i)=(cell2mat(wlev_tide(i)));
    end
end
tide2.dn=T;
tide2.wl=WL;

kk=find(tide2.wl==0);
tide2.wl(kk)=nan;

clear k fname kk wlev_tide T WL i dn_tide

k=csvimport([buoydir 'SanLuis_9412110/CO-OPS__9412110__wl_msl_Nov2017.csv']);

dn_tide=k(2:end,1);
wlev_tide=k(2:end,2);
T(1:length(dn_tide),1)=nan;
WL(1:length(dn_tide),1)=nan;
for i=1:length(dn_tide)
    T(i)=datenum(cell2mat(dn_tide(i)));
    if ~isempty((cell2mat(wlev_tide(i))))
        WL(i)=(cell2mat(wlev_tide(i)));
    end
end
tide3.dn=T;
tide3.wl=WL;

kk=find(tide3.wl==0);
tide3.wl(kk)=nan;

clear k fname kk wlev_tide T WL i dn_tide

%%
dn=[tide1.dn; tide2.dn; tide3.dn];
wl=[tide1.wl; tide2.wl; tide3.wl];


save([buoydir 'SanLuis_9412110/SanLuisWaterLevel_msl.mat'],'dn','wl');