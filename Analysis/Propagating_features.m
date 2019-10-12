% Written by Jack McSweeney 
% July 23, 2019

close all
clear 

addpath(genpath('/Volumes/InnerShelf1/MatlabCode/'));

% analysis directory
anadir= '/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/JPOmanuscript/PropagatingFeatures/';

%% Load Data 
% These are all extrapolated, cleaned up products with 1 min resolution 

cleandatadir='/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/cleanup/';

MS100=load([cleandatadir 'MS100.mat']); MS100.dn=MS100.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/ms100_alpha_beta_Cw.mat','alpha','Cw');
MS100.alpha=tempo.alpha; MS100.Cw=tempo.Cw;  clear tempo;

OC50=load([cleandatadir 'OC50.mat']); OC50.dn=OC50.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc50_alpha_beta_Cw.mat','alpha','Cw');
OC50.alpha=tempo.alpha; OC50.Cw=tempo.Cw; clear tempo;

OC40N=load([cleandatadir 'OC40N.mat']); OC40N.dn=OC40N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc40n_alpha_beta_Cw.mat','alpha','Cw');
OC40N.alpha=tempo.alpha; OC40N.Cw=tempo.Cw; clear tempo;

OC40S=load([cleandatadir 'OC40S.mat']); OC40S.dn=OC40S.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc40s_alpha_beta_Cw.mat','alpha','Cw');
OC40S.alpha=tempo.alpha; OC40S.Cw=tempo.Cw; clear tempo;

OC32N=load([cleandatadir 'OC32N.mat']); OC32N.dn=OC32N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc32n_alpha_beta_Cw.mat','alpha','Cw');
OC32N.alpha=tempo.alpha; OC32N.Cw=tempo.Cw; clear tempo;

OC32S=load([cleandatadir 'OC32S.mat']); OC32S.dn=OC32S.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc32s_alpha_beta_Cw.mat','alpha','Cw');
OC32S.alpha=tempo.alpha; OC32S.Cw=tempo.Cw;clear tempo;

OC25NA=load([cleandatadir 'OC25NA.mat']); OC25NA.dn=OC25NA.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25na_alpha_beta_Cw.mat','alpha','Cw');
OC25NA.alpha=tempo.alpha; OC25NA.Cw=tempo.Cw; clear tempo;

OC25NB=load([cleandatadir 'OC25NB.mat']); OC25NB.dn=OC25NB.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25nb_alpha_beta_Cw.mat','alpha','Cw');
OC25NB.alpha=tempo.alpha; OC25NB.Cw=tempo.Cw; clear tempo;

OC25M=load([cleandatadir 'OC25M.mat']); OC25M.dn=OC25M.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25m_alpha_beta_Cw.mat','alpha','Cw');
OC25M.alpha=tempo.alpha; OC25M.Cw=tempo.Cw; clear tempo;

OC25SB=load([cleandatadir 'OC25SB.mat']); OC25SB.dn=OC25SB.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25sb_alpha_beta_Cw.mat','alpha','Cw');
OC25SB.alpha=tempo.alpha; OC25SB.Cw=tempo.Cw; clear tempo;

OC25SA=load([cleandatadir 'OC25SA.mat']); OC25SA.dn=OC25SA.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc25sa_alpha_beta_Cw.mat','alpha','Cw');
OC25SA.alpha=tempo.alpha; OC25SA.Cw=tempo.Cw; clear tempo;

OC17N=load([cleandatadir 'OC17N.mat']); OC17N.dn=OC17N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc17n_alpha_beta_Cw.mat','alpha','Cw');
OC17N.alpha=tempo.alpha; OC17N.Cw=tempo.Cw;  clear tempo;

OC17S=load([cleandatadir 'OC17S.mat']); OC17S.dn=OC17S.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc17s_alpha_beta_Cw.mat','alpha','Cw');
OC17S.alpha=tempo.alpha; OC17S.Cw=tempo.Cw; clear tempo;

OC10N=load([cleandatadir 'OC10N.mat']); OC10N.dn=OC10N.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/oc10n_alpha_beta_Cw.mat','alpha','Cw');
OC10N.alpha=tempo.alpha; OC10N.Cw=tempo.Cw; clear tempo;

STR3B=load([cleandatadir 'STR3B.mat']); STR3B.dn=STR3B.dn_sig(1,:);
tempo=load('/Volumes/InnerShelf1/JackAnalysis/Moorings/ShowTwoTides/AlphaBeta/str3b_alpha_beta_Cw.mat','alpha','Cw');
STR3B.alpha=tempo.alpha; STR3B.Cw=tempo.Cw; clear tempo;

NRL50N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL50N.mat');
PS50=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS50.mat');
VB50N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB50N.mat');
VB50S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB50S.mat');
NRL50S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL50S.mat');
PS40N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS40N.mat');
PS40M=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS40M.mat');
PS40S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS40S.mat');
NRL35N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL35N.mat');
PS35M=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS35M.mat');
NRL35S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL35S.mat');
PS30N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS30N.mat');
PS30M=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS30M.mat');
% NEED TO CLEAN THIS FILE
PS30S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PS30S.mat');

VB30N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB30N.mat');
VB30S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB30S.mat');
VB25N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB25N.mat');
VB25S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/VB25S.mat');
NRL20N=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL20N.mat');
NRL20S=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NRL20S.mat');

NW2=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/NW2.mat');
RW3=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/RW3.mat');
BW1=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/BW1.mat');
PW5=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/PW5.mat');
STR4F=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/STR4F.mat');
STR5B=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/STR5B.mat');
STR6B=load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/CleanUpMoorings/STR6B.mat');

%% Load tide data - port of San Luis 
tide=load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SanLuis_9412110/SanLuisWaterLevel_msl.mat');

%% Download wind data 
load('/Volumes/InnerShelf1/NOAA_Buoy_Data/SantaMaria_46011/SantaMariaWind.mat');
dn_atm2=datenum(2017,9,1):2/24:datenum(2017,11,5);
wnorth=interp1(dn_atm,pl33tn(wind_n,1,33),dn_atm2);
weast=interp1(dn_atm,pl33tn(wind_e,1,33),dn_atm2);
%% Filter the Velocities for the Kinetic Energy Calculations  

filt=33; %33; % the filter cutoff; 33 hours here 
dt= 1/60; % the sample period; 
filt1=16; % the filter cutoff; 16 hours here 
filt2=1; % the filter cutoff; 1 hours here 
filt3= 3/60; % the filter cutoff; 1 hours here

MS100.north_res=MS100.Vsig-repmat(nanmean(MS100.Vsig),length(MS100.sig),1);
MS100.east_res=MS100.Usig-repmat(nanmean(MS100.Usig),length(MS100.sig),1);
% MS100.north_subtidal=pl66tn(MS100.north_res,dt,filt)';
% MS100.east_subtidal=pl66tn(MS100.east_res,dt,filt)';
% MS100.vert_subtidal=pl66tn(MS100.Wsig,dt,filt)';
tempo=MS100.north_res-pl66tn(MS100.north_res,dt,filt1)';
MS100.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=MS100.east_res-pl66tn(MS100.east_res,dt,filt1)';
MS100.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=MS100.Wsig-pl66tn(MS100.Wsig,dt,filt1)';
MS100.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished MS100 velocity filtering');

OC50.north_res=OC50.Vsig-repmat(nanmean(OC50.Vsig),length(OC50.sig),1);
OC50.east_res=OC50.Usig-repmat(nanmean(OC50.Usig),length(OC50.sig),1);
% OC50.north_subtidal=pl66tn(OC50.north_res,dt,filt)';
% OC50.east_subtidal=pl66tn(OC50.east_res,dt,filt)';
% OC50.vert_subtidal=pl66tn(OC50.Wsig,dt,filt)';
tempo=OC50.north_res-pl66tn(OC50.north_res,dt,filt1)';
OC50.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC50.east_res-pl66tn(OC50.east_res,dt,filt1)';
OC50.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC50.Wsig-pl66tn(OC50.Wsig,dt,filt1)';
OC50.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC50 velocity filtering');

% % If I want to separate out high frequency and semidiurnal 
% tempo=OC50.north_res-pl66tn(OC50.north_res,dt,filt1)';
% OC50.north_semidiurnal=pl66tn(tempo,dt,filt2)'; clear tempo
% tempo=OC50.north_res-pl66tn(OC50.north_res,dt,filt2)';
% OC50.north_hf=pl66tn(tempo,dt,filt3)'; clear tempo
% tempo=OC50.east_res-pl66tn(OC50.east_res,dt,filt1)';
% OC50.east_semidiurnal=pl66tn(tempo,dt,filt2)'; clear tempo
% tempo=OC50.east_res-pl66tn(OC50.east_res,dt,filt2)';
% OC50.east_hf=pl66tn(tempo,dt,filt3)'; clear tempo
% tempo=OC50.Wsig-pl66tn(OC50.Wsig,dt,filt1)';
% OC50.vert_semidiurnal=pl66tn(tempo,dt,filt2)'; clear tempo
% tempo=OC50.Wsig-pl66tn(OC50.Wsig,dt,filt2)';
% OC50.vert_hf=pl66tn(tempo,dt,filt3)'; clear tempo;

Nz = 50 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL50N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
NRL50N.north_res=NRL50N.Vsig-repmat(nanmean(NRL50N.Vsig),length(NRL50N.sig),1);
NRL50N.east_res=NRL50N.Usig-repmat(nanmean(NRL50N.Usig),length(NRL50N.sig),1);
% NRL50N.north_subtidal=pl66tn(NRL50N.north_res,dt,filt)';
% NRL50N.east_subtidal=pl66tn(NRL50N.east_res,dt,filt)';
% NRL50N.vert_subtidal=pl66tn(NRL50N.Wsig,dt,filt)';
tempo=NRL50N.north_res-pl66tn(NRL50N.north_res,dt,filt1)';
NRL50N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50N.east_res-pl66tn(NRL50N.east_res,dt,filt1)';
NRL50N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50N.Wsig-pl66tn(NRL50N.Wsig,dt,filt1)';
NRL50N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished NRL50N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS50.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS50.north_res=PS50.Vsig-repmat(nanmean(PS50.Vsig),length(PS50.sig),1);
PS50.east_res=PS50.Usig-repmat(nanmean(PS50.Usig),length(PS50.sig),1);
% PS50.north_subtidal=pl66tn(PS50.north_res,dt,filt)';
% PS50.east_subtidal=pl66tn(PS50.east_res,dt,filt)';
% PS50.vert_subtidal=pl66tn(PS50.Wsig,dt,filt)';
tempo=PS50.north_res-pl66tn(PS50.north_res,dt,filt1)';
PS50.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS50.east_res-pl66tn(PS50.east_res,dt,filt1)';
PS50.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS50.Wsig-pl66tn(PS50.Wsig,dt,filt1)';
PS50.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS50 velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB50N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB50N.north_res=VB50N.Vsig-repmat(nanmean(VB50N.Vsig),length(VB50N.sig),1);
VB50N.east_res=VB50N.Usig-repmat(nanmean(VB50N.Usig),length(VB50N.sig),1);
% VB50N.north_subtidal=pl66tn(VB50N.north_res,dt,filt)';
% VB50N.east_subtidal=pl66tn(VB50N.east_res,dt,filt)';
% VB50N.vert_subtidal=pl66tn(VB50N.Wsig,dt,filt)';
tempo=VB50N.north_res-pl66tn(VB50N.north_res,dt,filt1)';
VB50N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50N.east_res-pl66tn(VB50N.east_res,dt,filt1)';
VB50N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50N.Wsig-pl66tn(VB50N.Wsig,dt,filt1)';
VB50N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB50N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB50S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB50S.north_res=VB50S.Vsig-repmat(nanmean(VB50S.Vsig),length(VB50S.sig),1);
VB50S.east_res=VB50S.Usig-repmat(nanmean(VB50S.Usig),length(VB50S.sig),1);
% VB50S.north_subtidal=pl66tn(VB50S.north_res,dt,filt)';
% VB50S.east_subtidal=pl66tn(VB50S.east_res,dt,filt)';
% VB50S.vert_subtidal=pl66tn(VB50S.Wsig,dt,filt)';
tempo=VB50S.north_res-pl66tn(VB50S.north_res,dt,filt1)';
VB50S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50S.east_res-pl66tn(VB50S.east_res,dt,filt1)';
VB50S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB50S.Wsig-pl66tn(VB50S.Wsig,dt,filt1)';
VB50S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB50S velocity filtering');

Nz = 50 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL50S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
NRL50S.north_res=NRL50S.Vsig-repmat(nanmean(NRL50S.Vsig),length(NRL50S.sig),1);
NRL50S.east_res=NRL50S.Usig-repmat(nanmean(NRL50S.Usig),length(NRL50S.sig),1);
% NRL50S.north_subtidal=pl66tn(NRL50S.north_res,dt,filt)';
% NRL50S.east_subtidal=pl66tn(NRL50S.east_res,dt,filt)';
% NRL50S.vert_subtidal=pl66tn(NRL50S.Wsig,dt,filt)';
tempo=NRL50S.north_res-pl66tn(NRL50S.north_res,dt,filt1)';
NRL50S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50S.east_res-pl66tn(NRL50S.east_res,dt,filt1)';
NRL50S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL50S.Wsig-pl66tn(NRL50S.Wsig,dt,filt1)';
NRL50S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished NRL50S velocity filtering');

OC40N.north_res=OC40N.Vsig-repmat(nanmean(OC40N.Vsig),length(OC40N.sig),1);
OC40N.east_res=OC40N.Usig-repmat(nanmean(OC40N.Usig),length(OC40N.sig),1);
% OC40N.north_subtidal=pl66tn(OC40N.north_res,dt,filt)';
% OC40N.east_subtidal=pl66tn(OC40N.east_res,dt,filt)';
% OC40N.vert_subtidal=pl66tn(OC40N.Wsig,dt,filt)';
tempo=OC40N.north_res-pl66tn(OC40N.north_res,dt,filt1)';
OC40N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC40N.east_res-pl66tn(OC40N.east_res,dt,filt1)';
OC40N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC40N.Wsig-pl66tn(OC40N.Wsig,dt,filt1)';
OC40N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC40N velocity filtering');

OC40S.north_res=OC40S.Vsig-repmat(nanmean(OC40S.Vsig),length(OC40S.sig),1);
OC40S.east_res=OC40S.Usig-repmat(nanmean(OC40S.Usig),length(OC40S.sig),1);
% OC40S.north_subtidal=pl66tn(OC40S.north_res,dt,filt)';
% OC40S.east_subtidal=pl66tn(OC40S.east_res,dt,filt)';
% OC40S.vert_subtidal=pl66tn(OC40S.Wsig,dt,filt)';
tempo=OC40S.north_res-pl66tn(OC40S.north_res,dt,filt1)';
OC40S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC40S.east_res-pl66tn(OC40S.east_res,dt,filt1)';
OC40S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC40S.Wsig-pl66tn(OC40S.Wsig,dt,filt1)';
OC40S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC40S velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS40N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS40N.north_res=PS40N.Vsig-repmat(nanmean(PS40N.Vsig),length(PS40N.sig),1);
PS40N.east_res=PS40N.Usig-repmat(nanmean(PS40N.Usig),length(PS40N.sig),1);
% PS40N.north_subtidal=pl66tn(PS40N.north_res,dt,filt)';
% PS40N.east_subtidal=pl66tn(PS40N.east_res,dt,filt)';
% PS40N.vert_subtidal=pl66tn(PS40N.Wsig,dt,filt)';
tempo=PS40N.north_res-pl66tn(PS40N.north_res,dt,filt1)';
PS40N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS40N.east_res-pl66tn(PS40N.east_res,dt,filt1)';
PS40N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS40N.Wsig-pl66tn(PS40N.Wsig,dt,filt1)';
PS40N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS40N velocity filtering');


Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS40M.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS40M.north_res=PS40M.Vsig-repmat(nanmean(PS40M.Vsig),length(PS40M.sig),1);
PS40M.east_res=PS40M.Usig-repmat(nanmean(PS40M.Usig),length(PS40M.sig),1);
% PS40M.north_subtidal=pl66tn(PS40M.north_res,dt,filt)';
% PS40M.east_subtidal=pl66tn(PS40M.east_res,dt,filt)';
% PS40M.vert_subtidal=pl66tn(PS40M.Wsig,dt,filt)';
tempo=PS40M.north_res-pl66tn(PS40M.north_res,dt,filt1)';
PS40M.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS40M.east_res-pl66tn(PS40M.east_res,dt,filt1)';
PS40M.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS40M.Wsig-pl66tn(PS40M.Wsig,dt,filt1)';
PS40M.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS40M velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS40S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS40S.north_res=PS40S.Vsig-repmat(nanmean(PS40S.Vsig),length(PS40S.sig),1);
PS40S.east_res=PS40S.Usig-repmat(nanmean(PS40S.Usig),length(PS40S.sig),1);
% PS40S.north_subtidal=pl66tn(PS40S.north_res,dt,filt)';
% PS40S.east_subtidal=pl66tn(PS40S.east_res,dt,filt)';
% PS40S.vert_subtidal=pl66tn(PS40S.Wsig,dt,filt)';
tempo=PS40S.north_res-pl66tn(PS40S.north_res,dt,filt1)';
PS40S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS40S.east_res-pl66tn(PS40S.east_res,dt,filt1)';
PS40S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS40S.Wsig-pl66tn(PS40S.Wsig,dt,filt1)';
PS40S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS40S velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS30M.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS30M.north_res=PS30M.Vsig-repmat(nanmean(PS30M.Vsig),length(PS30M.sig),1);
PS30M.east_res=PS30M.Usig-repmat(nanmean(PS30M.Usig),length(PS30M.sig),1);
% PS30M.north_subtidal=pl66tn(PS30M.north_res,dt,filt)';
% PS30M.east_subtidal=pl66tn(PS30M.east_res,dt,filt)';
% PS30M.vert_subtidal=pl66tn(PS30M.Wsig,dt,filt)';
tempo=PS30M.north_res-pl66tn(PS30M.north_res,dt,filt1)';
PS30M.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS30M.east_res-pl66tn(PS30M.east_res,dt,filt1)';
PS30M.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS30M.Wsig-pl66tn(PS30M.Wsig,dt,filt1)';
PS30M.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS30M velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
PS30S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
PS30S.north_res=PS30S.Vsig-repmat(nanmean(PS30S.Vsig),length(PS30S.sig),1);
PS30S.east_res=PS30S.Usig-repmat(nanmean(PS30S.Usig),length(PS30S.sig),1);
% PS30S.north_subtidal=pl66tn(PS30S.north_res,dt,filt)';
% PS30S.east_subtidal=pl66tn(PS30S.east_res,dt,filt)';
% PS30S.vert_subtidal=pl66tn(PS30S.Wsig,dt,filt)';
tempo=PS30S.north_res-pl66tn(PS30S.north_res,dt,filt1)';
PS30S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS30S.east_res-pl66tn(PS30S.east_res,dt,filt1)';
PS30S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PS30S.Wsig-pl66tn(PS30S.Wsig,dt,filt1)';
PS30S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PS30S velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB30N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB30N.north_res=VB30N.Vsig-repmat(nanmean(VB30N.Vsig),length(VB30N.sig),1);
VB30N.east_res=VB30N.Usig-repmat(nanmean(VB30N.Usig),length(VB30N.sig),1);
% VB30N.north_subtidal=pl66tn(VB30N.north_res,dt,filt)';
% VB30N.east_subtidal=pl66tn(VB30N.east_res,dt,filt)';
% VB30N.vert_subtidal=pl66tn(VB30N.Wsig,dt,filt)';
tempo=VB30N.north_res-pl66tn(VB30N.north_res,dt,filt1)';
VB30N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB30N.east_res-pl66tn(VB30N.east_res,dt,filt1)';
VB30N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB30N.Wsig-pl66tn(VB30N.Wsig,dt,filt1)';
VB30N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB30N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB30S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB30S.north_res=VB30S.Vsig-repmat(nanmean(VB30S.Vsig),length(VB30S.sig),1);
VB30S.east_res=VB30S.Usig-repmat(nanmean(VB30S.Usig),length(VB30S.sig),1);
% VB30S.north_subtidal=pl66tn(VB30S.north_res,dt,filt)';
% VB30S.east_subtidal=pl66tn(VB30S.east_res,dt,filt)';
% VB30S.vert_subtidal=pl66tn(VB30S.Wsig,dt,filt)';
tempo=VB30S.north_res-pl66tn(VB30S.north_res,dt,filt1)';
VB30S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB30S.east_res-pl66tn(VB30S.east_res,dt,filt1)';
VB30S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB30S.Wsig-pl66tn(VB30S.Wsig,dt,filt1)';
VB30S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB30S velocity filtering');

OC25NA.north_res=OC25NA.Vsig-repmat(nanmean(OC25NA.Vsig),length(OC25NA.sig),1);
OC25NA.east_res=OC25NA.Usig-repmat(nanmean(OC25NA.Usig),length(OC25NA.sig),1);
% OC25NA.north_subtidal=pl66tn(OC25NA.north_res,dt,filt)';
% OC25NA.east_subtidal=pl66tn(OC25NA.east_res,dt,filt)';
% OC25NA.vert_subtidal=pl66tn(OC25NA.Wsig,dt,filt)';
tempo=OC25NA.north_res-pl66tn(OC25NA.north_res,dt,filt1)';
OC25NA.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25NA.east_res-pl66tn(OC25NA.east_res,dt,filt1)';
OC25NA.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25NA.Wsig-pl66tn(OC25NA.Wsig,dt,filt1)';
OC25NA.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC25NA velocity filtering');

OC25NB.north_res=OC25NB.Vsig-repmat(nanmean(OC25NB.Vsig),length(OC25NB.sig),1);
OC25NB.east_res=OC25NB.Usig-repmat(nanmean(OC25NB.Usig),length(OC25NB.sig),1);
% OC25NB.north_subtidal=pl66tn(OC25NB.north_res,dt,filt)';
% OC25NB.east_subtidal=pl66tn(OC25NB.east_res,dt,filt)';
% OC25NB.vert_subtidal=pl66tn(OC25NB.Wsig,dt,filt)';
tempo=OC25NB.north_res-pl66tn(OC25NB.north_res,dt,filt1)';
OC25NB.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25NB.east_res-pl66tn(OC25NB.east_res,dt,filt1)';
OC25NB.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25NB.Wsig-pl66tn(OC25NB.Wsig,dt,filt1)';
OC25NB.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC25NB velocity filtering');

OC25M.north_res=OC25M.Vsig-repmat(nanmean(OC25M.Vsig),length(OC25M.sig),1);
OC25M.east_res=OC25M.Usig-repmat(nanmean(OC25M.Usig),length(OC25M.sig),1);
% OC25M.north_subtidal=pl66tn(OC25M.north_res,dt,filt)';
% OC25M.east_subtidal=pl66tn(OC25M.east_res,dt,filt)';
% OC25M.vert_subtidal=pl66tn(OC25M.Wsig,dt,filt)';
tempo=OC25M.north_res-pl66tn(OC25M.north_res,dt,filt1)';
OC25M.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25M.east_res-pl66tn(OC25M.east_res,dt,filt1)';
OC25M.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25M.Wsig-pl66tn(OC25M.Wsig,dt,filt1)';
OC25M.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC25M velocity filtering');

OC25SB.north_res=OC25SB.Vsig-repmat(nanmean(OC25SB.Vsig),length(OC25SB.sig),1);
OC25SB.east_res=OC25SB.Usig-repmat(nanmean(OC25SB.Usig),length(OC25SB.sig),1);
% OC25SB.north_subtidal=pl66tn(OC25SB.north_res,dt,filt)';
% OC25SB.east_subtidal=pl66tn(OC25SB.east_res,dt,filt)';
% OC25SB.vert_subtidal=pl66tn(OC25SB.Wsig,dt,filt)';
tempo=OC25SB.north_res-pl66tn(OC25SB.north_res,dt,filt1)';
OC25SB.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25SB.east_res-pl66tn(OC25SB.east_res,dt,filt1)';
OC25SB.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25SB.Wsig-pl66tn(OC25SB.Wsig,dt,filt1)';
OC25SB.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC25SB velocity filtering');

OC25SA.north_res=OC25SA.Vsig-repmat(nanmean(OC25SA.Vsig),length(OC25SA.sig),1);
OC25SA.east_res=OC25SA.Usig-repmat(nanmean(OC25SA.Usig),length(OC25SA.sig),1);
% OC25SA.north_subtidal=pl66tn(OC25SA.north_res,dt,filt)';
% OC25SA.east_subtidal=pl66tn(OC25SA.east_res,dt,filt)';
% OC25SA.vert_subtidal=pl66tn(OC25SA.Wsig,dt,filt)';
tempo=OC25SA.north_res-pl66tn(OC25SA.north_res,dt,filt1)';
OC25SA.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25SA.east_res-pl66tn(OC25SA.east_res,dt,filt1)';
OC25SA.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC25SA.Wsig-pl66tn(OC25SA.Wsig,dt,filt1)';
OC25SA.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC25SA velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
VB25N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
VB25N.north_res=VB25N.Vsig-repmat(nanmean(VB25N.Vsig),length(VB25N.sig),1);
VB25N.east_res=VB25N.Usig-repmat(nanmean(VB25N.Usig),length(VB25N.sig),1);
% VB25N.north_subtidal=pl66tn(VB25N.north_res,dt,filt)';
% VB25N.east_subtidal=pl66tn(VB25N.east_res,dt,filt)';
% VB25N.vert_subtidal=pl66tn(VB25N.Wsig,dt,filt)';
tempo=VB25N.north_res-pl66tn(VB25N.north_res,dt,filt1)';
VB25N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB25N.east_res-pl66tn(VB25N.east_res,dt,filt1)';
VB25N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=VB25N.Wsig-pl66tn(VB25N.Wsig,dt,filt1)';
VB25N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished VB25N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL20N.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
NRL20N.north_res=NRL20N.Vsig-repmat(nanmean(NRL20N.Vsig),length(NRL20N.sig),1);
NRL20N.east_res=NRL20N.Usig-repmat(nanmean(NRL20N.Usig),length(NRL20N.sig),1);
% NRL20N.north_subtidal=pl66tn(NRL20N.north_res,dt,filt)';
% NRL20N.east_subtidal=pl66tn(NRL20N.east_res,dt,filt)';
% NRL20N.vert_subtidal=pl66tn(NRL20N.Wsig,dt,filt)';
tempo=NRL20N.north_res-pl66tn(NRL20N.north_res,dt,filt1)';
NRL20N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL20N.east_res-pl66tn(NRL20N.east_res,dt,filt1)';
NRL20N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
%no vert vel
disp('Finished NRL20N velocity filtering');

Nz = 25 ;  % number of vertical levels 
dsig = 1/Nz ; % delta sigma
NRL20S.sig = (0:dsig:1)' ; % sigma 
clear Nz dsig 
NRL20S.north_res=NRL20S.Vsig-repmat(nanmean(NRL20S.Vsig),length(NRL20S.sig),1);
NRL20S.east_res=NRL20S.Usig-repmat(nanmean(NRL20S.Usig),length(NRL20S.sig),1);
% NRL20S.north_subtidal=pl66tn(NRL20S.north_res,dt,filt)';
% NRL20S.east_subtidal=pl66tn(NRL20S.east_res,dt,filt)';
% NRL20S.vert_subtidal=pl66tn(NRL20S.Wsig,dt,filt)';
tempo=NRL20S.north_res-pl66tn(NRL20S.north_res,dt,filt1)';
NRL20S.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL20S.east_res-pl66tn(NRL20S.east_res,dt,filt1)';
NRL20S.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NRL20S.Wsig-pl66tn(NRL20S.Wsig,dt,filt1)';
NRL20S.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished NRL20S velocity filtering');

OC10N.north_res=OC10N.Vsig-repmat(nanmean(OC10N.Vsig),length(OC10N.sig),1);
OC10N.east_res=OC10N.Usig-repmat(nanmean(OC10N.Usig),length(OC10N.sig),1);
% OC10N.north_subtidal=pl66tn(OC10N.north_res,dt,filt)';
% OC10N.east_subtidal=pl66tn(OC10N.east_res,dt,filt)';
% OC10N.vert_subtidal=pl66tn(OC10N.Wsig,dt,filt)';
tempo=OC10N.north_res-pl66tn(OC10N.north_res,dt,filt1)';
OC10N.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC10N.east_res-pl66tn(OC10N.east_res,dt,filt1)';
OC10N.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=OC10N.Wsig-pl66tn(OC10N.Wsig,dt,filt1)';
OC10N.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished OC10N velocity filtering');

STR3B.north_res=STR3B.Vsig-repmat(nanmean(STR3B.Vsig),length(STR3B.sig),1);
STR3B.east_res=STR3B.Usig-repmat(nanmean(STR3B.Usig),length(STR3B.sig),1);
% STR3B.north_subtidal=pl66tn(STR3B.north_res,dt,filt)';
% STR3B.east_subtidal=pl66tn(STR3B.east_res,dt,filt)';
% STR3B.vert_subtidal=pl66tn(STR3B.Wsig,dt,filt)';
tempo=STR3B.north_res-pl66tn(STR3B.north_res,dt,filt1)';
STR3B.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=STR3B.east_res-pl66tn(STR3B.east_res,dt,filt1)';
STR3B.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=STR3B.Wsig-pl66tn(STR3B.Wsig,dt,filt1)';
STR3B.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished STR3B velocity filtering');

% % % NO vertical velocity 
NW2.north_res=NW2.Vsig-repmat(nanmean(NW2.Vsig),length(NW2.sig),1);
NW2.east_res=NW2.Usig-repmat(nanmean(NW2.Usig),length(NW2.sig),1);
% NW2.north_subtidal=pl66tn(NW2.north_res,dt,filt)';
% NW2.east_subtidal=pl66tn(NW2.east_res,dt,filt)';
% NW2.vert_subtidal=pl66tn(NW2.Wsig,dt,filt)';
tempo=NW2.north_res-pl66tn(NW2.north_res,dt,filt1)';
NW2.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=NW2.east_res-pl66tn(NW2.east_res,dt,filt1)';
NW2.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
% tempo=NW2.Wsig-pl66tn(NW2.Wsig,dt,filt1)';
% NW2.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished NW2 velocity filtering');

% % % NO vertical velocity 
RW3.north_res=RW3.Vsig-repmat(nanmean(RW3.Vsig),length(RW3.sig),1);
RW3.east_res=RW3.Usig-repmat(nanmean(RW3.Usig),length(RW3.sig),1);
% RW3.north_subtidal=pl66tn(RW3.north_res,dt,filt)';
% RW3.east_subtidal=pl66tn(RW3.east_res,dt,filt)';
% RW3.vert_subtidal=pl66tn(RW3.Wsig,dt,filt)';
tempo=RW3.north_res-pl66tn(RW3.north_res,dt,filt1)';
RW3.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=RW3.east_res-pl66tn(RW3.east_res,dt,filt1)';
RW3.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
% tempo=RW3.Wsig-pl66tn(RW3.Wsig,dt,filt1)';
% RW3.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished RW3 velocity filtering');

% % % NO vertical velocity 
BW1.north_res=BW1.Vsig-repmat(nanmean(BW1.Vsig),length(BW1.sig),1);
BW1.east_res=BW1.Usig-repmat(nanmean(BW1.Usig),length(BW1.sig),1);
% BW1.north_subtidal=pl66tn(BW1.north_res,dt,filt)';
% BW1.east_subtidal=pl66tn(BW1.east_res,dt,filt)';
% BW1.vert_subtidal=pl66tn(BW1.Wsig,dt,filt)';
tempo=BW1.north_res-pl66tn(BW1.north_res,dt,filt1)';
BW1.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=BW1.east_res-pl66tn(BW1.east_res,dt,filt1)';
BW1.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
% tempo=BW1.Wsig-pl66tn(BW1.Wsig,dt,filt1)';
% BW1.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished BW1 velocity filtering');

% % % NO vertical velocity 
PW5.north_res=PW5.Vsig-repmat(nanmean(PW5.Vsig),length(PW5.sig),1);
PW5.east_res=PW5.Usig-repmat(nanmean(PW5.Usig),length(PW5.sig),1);
% PW5.north_subtidal=pl66tn(PW5.north_res,dt,filt)';
% PW5.east_subtidal=pl66tn(PW5.east_res,dt,filt)';
% PW5.vert_subtidal=pl66tn(PW5.Wsig,dt,filt)';
tempo=PW5.north_res-pl66tn(PW5.north_res,dt,filt1)';
PW5.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=PW5.east_res-pl66tn(PW5.east_res,dt,filt1)';
PW5.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
% tempo=PW5.Wsig-pl66tn(PW5.Wsig,dt,filt1)';
% PW5.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished PW5 velocity filtering');


STR5B.north_res=STR5B.Vsig-repmat(nanmean(STR5B.Vsig),length(STR5B.sig),1);
STR5B.east_res=STR5B.Usig-repmat(nanmean(STR5B.Usig),length(STR5B.sig),1);
% STR5B.north_subtidal=pl66tn(STR5B.north_res,dt,filt)';
% STR5B.east_subtidal=pl66tn(STR5B.east_res,dt,filt)';
% STR5B.vert_subtidal=pl66tn(STR5B.Wsig,dt,filt)';
tempo=STR5B.north_res-pl66tn(STR5B.north_res,dt,filt1)';
STR5B.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=STR5B.east_res-pl66tn(STR5B.east_res,dt,filt1)';
STR5B.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=STR5B.Wsig-pl66tn(STR5B.Wsig,dt,filt1)';
STR5B.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished STR5B velocity filtering');


STR6B.north_res=STR6B.Vsig-repmat(nanmean(STR6B.Vsig),length(STR6B.sig),1);
STR6B.east_res=STR6B.Usig-repmat(nanmean(STR6B.Usig),length(STR6B.sig),1);
% STR6B.north_subtidal=pl66tn(STR6B.north_res,dt,filt)';
% STR6B.east_subtidal=pl66tn(STR6B.east_res,dt,filt)';
% STR6B.vert_subtidal=pl66tn(STR6B.Wsig,dt,filt)';
tempo=STR6B.north_res-pl66tn(STR6B.north_res,dt,filt1)';
STR6B.north_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=STR6B.east_res-pl66tn(STR6B.east_res,dt,filt1)';
STR6B.east_sd=pl66tn(tempo,dt,filt3)'; clear tempo
tempo=STR6B.Wsig-pl66tn(STR6B.Wsig,dt,filt1)';
STR6B.vert_sd=pl66tn(tempo,dt,filt3)'; clear tempo
disp('Finished STR6B velocity filtering');

%% Plot the filtered velocities at OC50 
if 0
    close all
    mapb=flipud(cbrewer('div','PuOr',40));
    map=flipud(cbrewer('div','RdBu',40));
    
    figure('position',[118          52        1266         753]);
    subplot(4,2,1);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.north_res); colormap(gca,mapb); caxis([-.5 0.5]); colorbar;
    title(['North' char(10) 'Residual (total-depthmean)']);
    subplot(4,2,3);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.north_subtidal);colormap(gca,mapb); caxis([-.5 0.5]);colorbar;
    title('> 33hrs')
    subplot(4,2,5);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.north_semidiurnal);colormap(gca,mapb); caxis([-.5 0.5]);colorbar;
    title('1 hr -16hrs')
    subplot(4,2,7);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.north_hf); colormap(gca,mapb); caxis([-.5 0.5]);colorbar;
    title('3 min - 1 hr')
    
    subplot(4,2,2);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.east_res); colormap(gca,map); caxis([-.5 0.5]);colorbar;
    title(['East' char(10) 'Residual (total-depthmean)']);
    subplot(4,2,4);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.east_subtidal);colormap(gca,map); caxis([-.5 0.5]);colorbar;
    title('> 33hrs')
    subplot(4,2,6);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.east_semidiurnal);colormap(gca,map); caxis([-.5 0.5]);colorbar;
    title('1 hr -16hrs')
    subplot(4,2,8);
    pcolorjw(OC50.dn_sig,OC50.H,OC50.east_hf); colormap(gca,map); caxis([-.5 0.5]);colorbar;
    title('3 min - 1 hr')
    
    if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir 'KineticEnergy/OC50Vel']);
    end
end

%% Calculate Kinetic Energy Denisty
MS100.rho=nanmean(nanmean(MS100.Rho_sig));
ke_sd.MS100=nanmean(0.5*MS100.rho*(MS100.east_sd.^2+MS100.north_sd.^2+MS100.vert_sd.^2));

OC50.rho=nanmean(nanmean(OC50.Rho_sig));
ke_sd.OC50=nanmean(0.5*OC50.rho*(OC50.east_sd.^2+OC50.north_sd.^2+OC50.vert_sd.^2));

sal=33.4188991428255;
NRL50N.Psig = ones(size(sal))*nanmax(NRL50N.H) - NRL50N.H ; % pressure 
NRL50N.Rho_sig=sw_dens(0*NRL50N.Tsig+sal,NRL50N.Tsig,NRL50N.Psig) ;
NRL50N.rho=nanmean(nanmean(NRL50N.Rho_sig));
ke_sd.NRL50N=nanmean(0.5*NRL50N.rho*(NRL50N.east_sd.^2+NRL50N.north_sd.^2+NRL50N.vert_sd.^2));

sal=33.4188991428255;
PS50.Psig = ones(size(sal))*nanmax(PS50.H) - PS50.H ; % pressure 
PS50.Rho_sig=sw_dens(0*PS50.Tsig+sal,PS50.Tsig,PS50.Psig) ;
PS50.rho=nanmean(nanmean(PS50.Rho_sig));
ke_sd.PS50=nanmean(0.5*PS50.rho*(PS50.east_sd.^2+PS50.north_sd.^2+PS50.vert_sd.^2));

sal=33.4188991428255;
VB50N.Psig = ones(size(sal))*nanmax(VB50N.H) - VB50N.H ; % pressure 
VB50N.Rho_sig=sw_dens(0*VB50N.Tsig+sal,VB50N.Tsig,VB50N.Psig) ;
VB50N.rho=nanmean(nanmean(VB50N.Rho_sig));
ke_sd.VB50N=nanmean(0.5*VB50N.rho*(VB50N.east_sd.^2+VB50N.north_sd.^2+VB50N.vert_sd.^2));

sal=33.4188991428255;
VB50S.Psig = ones(size(sal))*nanmax(VB50S.H) - VB50S.H ; % pressure 
VB50S.Rho_sig=sw_dens(0*VB50S.Tsig+sal,VB50S.Tsig,VB50S.Psig) ;
VB50S.rho=nanmean(nanmean(VB50S.Rho_sig));
ke_sd.VB50S=nanmean(0.5*VB50S.rho*(VB50S.east_sd.^2+VB50S.north_sd.^2+VB50S.vert_sd.^2));

sal=33.4188991428255;
NRL50S.Psig = ones(size(sal))*nanmax(NRL50S.H) - NRL50S.H ; % pressure 
NRL50S.Rho_sig=sw_dens(0*NRL50S.Tsig+sal,NRL50S.Tsig,NRL50S.Psig) ;
NRL50S.rho=nanmean(nanmean(NRL50S.Rho_sig));
ke_sd.NRL50S=nanmean(0.5*NRL50S.rho*(NRL50S.east_sd.^2+NRL50S.north_sd.^2+NRL50S.vert_sd.^2));

OC40N.rho=nanmean(nanmean(OC40N.Rho_sig));
ke_sd.OC40N=nanmean(0.5*OC40N.rho*(OC40N.east_sd.^2+OC40N.north_sd.^2+OC40N.vert_sd.^2));

OC40S.rho=nanmean(nanmean(OC40S.Rho_sig));
ke_sd.OC40S=nanmean(0.5*OC40S.rho*(OC40S.east_sd.^2+OC40S.north_sd.^2+OC40S.vert_sd.^2));

sal=33.4188991428255;
PS40N.Psig = ones(size(sal))*nanmax(PS40N.H) - PS40N.H ; % pressure 
PS40N.Rho_sig=sw_dens(0*PS40N.Tsig+sal,PS40N.Tsig,PS40N.Psig) ;
PS40N.rho=nanmean(nanmean(PS40N.Rho_sig));
ke_sd.PS40N=nanmean(0.5*PS40N.rho*(PS40N.east_sd.^2+PS40N.north_sd.^2+PS40N.vert_sd.^2));

sal=33.4188991428255;
PS40M.Psig = ones(size(sal))*nanmax(PS40M.H) - PS40M.H ; % pressure 
PS40M.Rho_sig=sw_dens(0*PS40M.Tsig+sal,PS40M.Tsig,PS40M.Psig) ;
PS40M.rho=nanmean(nanmean(PS40M.Rho_sig));
ke_sd.PS40M=nanmean(0.5*PS40M.rho*(PS40M.east_sd.^2+PS40M.north_sd.^2+PS40M.vert_sd.^2));

sal=33.4188991428255;
PS40S.Psig = ones(size(sal))*nanmax(PS40S.H) - PS40S.H ; % pressure 
PS40S.Rho_sig=sw_dens(0*PS40S.Tsig+sal,PS40S.Tsig,PS40S.Psig) ;
PS40S.rho=nanmean(nanmean(PS40S.Rho_sig));
ke_sd.PS40S=nanmean(0.5*PS40S.rho*(PS40S.east_sd.^2+PS40S.north_sd.^2+PS40S.vert_sd.^2));

sal=33.4188991428255;
PS30M.Psig = ones(size(sal))*nanmax(PS30M.H) - PS30M.H ; % pressure 
PS30M.Rho_sig=sw_dens(0*PS30M.Tsig+sal,PS30M.Tsig,PS30M.Psig) ;
PS30M.rho=nanmean(nanmean(PS30M.Rho_sig));
ke_sd.PS30M=nanmean(0.5*PS30M.rho*(PS30M.east_sd.^2+PS30M.north_sd.^2+PS30M.vert_sd.^2));

sal=33.4188991428255;
PS30S.Psig = ones(size(sal))*nanmax(PS30S.H) - PS30S.H ; % pressure 
PS30S.Rho_sig=sw_dens(0*PS30S.Tsig+sal,PS30S.Tsig,PS30S.Psig) ;
PS30S.rho=nanmean(nanmean(PS30S.Rho_sig));
ke_sd.PS30S=nanmean(0.5*PS30S.rho*(PS30S.east_sd.^2+PS30S.north_sd.^2+PS30S.vert_sd.^2));

sal=33.4188991428255;
VB30N.Psig = ones(size(sal))*nanmax(VB30N.H) - VB30N.H ; % pressure 
VB30N.Rho_sig=sw_dens(0*VB30N.Tsig+sal,VB30N.Tsig,VB30N.Psig) ;
VB30N.rho=nanmean(nanmean(VB30N.Rho_sig));
ke_sd.VB30N=nanmean(0.5*VB30N.rho*(VB30N.east_sd.^2+VB30N.north_sd.^2+VB30N.vert_sd.^2));

sal=33.4188991428255;
VB30S.Psig = ones(size(sal))*nanmax(VB30S.H) - VB30S.H ; % pressure 
VB30S.Rho_sig=sw_dens(0*VB30S.Tsig+sal,VB30S.Tsig,VB30S.Psig) ;
VB30S.rho=nanmean(nanmean(VB30S.Rho_sig));
ke_sd.VB30S=nanmean(0.5*VB30S.rho*(VB30S.east_sd.^2+VB30S.north_sd.^2+VB30S.vert_sd.^2));

OC25NA.rho=nanmean(nanmean(OC25NA.Rho_sig));
ke_sd.OC25NA=nanmean(0.5*OC25NA.rho*(OC25NA.east_sd.^2+OC25NA.north_sd.^2+OC25NA.vert_sd.^2));

OC25NB.rho=nanmean(nanmean(OC25NB.Rho_sig));
ke_sd.OC25NB=nanmean(0.5*OC25NB.rho*(OC25NB.east_sd.^2+OC25NB.north_sd.^2+OC25NB.vert_sd.^2));

OC25M.rho=nanmean(nanmean(OC25M.Rho_sig));
ke_sd.OC25M=nanmean(0.5*OC25M.rho*(OC25M.east_sd.^2+OC25M.north_sd.^2+OC25M.vert_sd.^2));

OC25SB.rho=nanmean(nanmean(OC25SB.Rho_sig));
ke_sd.OC25SB=nanmean(0.5*OC25SB.rho*(OC25SB.east_sd.^2+OC25SB.north_sd.^2+OC25SB.vert_sd.^2));

OC25SA.rho=nanmean(nanmean(OC25SA.Rho_sig));
ke_sd.OC25SA=nanmean(0.5*OC25SA.rho*(OC25SA.east_sd.^2+OC25SA.north_sd.^2+OC25SA.vert_sd.^2));

sal=33.4188991428255;
VB25N.Psig = ones(size(sal))*nanmax(VB25N.H) - VB25N.H ; % pressure 
VB25N.Rho_sig=sw_dens(0*VB25N.Tsig+sal,VB25N.Tsig,VB25N.Psig) ;
VB25N.rho=nanmean(nanmean(VB25N.Rho_sig));
ke_sd.VB25N=nanmean(0.5*VB25N.rho*(VB25N.east_sd.^2+VB25N.north_sd.^2+VB25N.vert_sd.^2));

sal=33.4188991428255;
NRL20N.Psig = ones(size(sal))*nanmax(NRL20N.H) - NRL20N.H ; % pressure 
NRL20N.Rho_sig=sw_dens(0*NRL20N.Tsig+sal,NRL20N.Tsig,NRL20N.Psig) ;
NRL20N.rho=nanmean(nanmean(NRL20N.Rho_sig));
ke_sd.NRL20N=nanmean(0.5*NRL20N.rho*(NRL20N.east_sd.^2+NRL20N.north_sd.^2));

sal=33.4188991428255;
NRL20S.Psig = ones(size(sal))*nanmax(NRL20S.H) - NRL20S.H ; % pressure 
NRL20S.Rho_sig=sw_dens(0*NRL20S.Tsig+sal,NRL20S.Tsig,NRL20S.Psig) ;
NRL20S.rho=nanmean(nanmean(NRL20S.Rho_sig));
ke_sd.NRL20S=nanmean(0.5*NRL20S.rho*(NRL20S.east_sd.^2+NRL20S.north_sd.^2+NRL20S.vert_sd.^2));

OC10N.rho=nanmean(nanmean(OC10N.Rho_sig));
ke_sd.OC10N=nanmean(0.5*OC10N.rho*(OC10N.east_sd.^2+OC10N.north_sd.^2+OC10N.vert_sd.^2));

STR3B.rho=nanmean(nanmean(STR3B.Rho_sig));
ke_sd.STR3B=nanmean(0.5*STR3B.rho*(STR3B.east_sd.^2+STR3B.north_sd.^2+STR3B.vert_sd.^2));

% no vert 
NW2.rho=nanmean(nanmean(NW2.Rho_sig));
ke_sd.NW2=nanmean(0.5*NW2.rho*(NW2.east_sd.^2+NW2.north_sd.^2));

% no vert 
RW3.rho=nanmean(nanmean(RW3.Rho_sig));
ke_sd.RW3=nanmean(0.5*RW3.rho*(RW3.east_sd.^2+RW3.north_sd.^2));

% no vert 
BW1.rho=nanmean(nanmean(BW1.Rho_sig));
ke_sd.BW1=nanmean(0.5*BW1.rho*(BW1.east_sd.^2+BW1.north_sd.^2));

% no vert 
PW5.rho=nanmean(nanmean(PW5.Rho_sig));
ke_sd.PW5=nanmean(0.5*PW5.rho*(PW5.east_sd.^2+PW5.north_sd.^2));


STR5B.rho=nanmean(nanmean(STR5B.Rho_sig));
ke_sd.STR5B=nanmean(0.5*STR5B.rho*(STR5B.east_sd.^2+STR5B.north_sd.^2+STR5B.vert_sd.^2));

STR6B.rho=nanmean(nanmean(STR6B.Rho_sig));
ke_sd.STR6B=nanmean(0.5*STR6B.rho*(STR6B.east_sd.^2+STR6B.north_sd.^2+STR6B.vert_sd.^2));

%% 
dn=datenum(2017,9,6,0,0,0):1/24:datenum(2017,11,3,0,0,0);
KEsd.ms100=interp1(MS100.dn,ke_sd.MS100,dn);
KEsd.oc50=interp1(OC50.dn,ke_sd.OC50,dn);
KEsd.nrl50n=interp1(NRL50N.dn,ke_sd.NRL50N,dn);
KEsd.ps50=interp1(PS50.dn,ke_sd.PS50,dn);
KEsd.vb50n=interp1(VB50N.dn,ke_sd.VB50N,dn);
KEsd.vb50s=interp1(VB50S.dn,ke_sd.VB50S,dn);
KEsd.nrl50s=interp1(NRL50S.dn,ke_sd.NRL50S,dn);
KEsd.oc40n=interp1(OC40N.dn,ke_sd.OC40N,dn);
KEsd.oc40s=interp1(OC40S.dn,ke_sd.OC40S,dn);
KEsd.ps40n=interp1(PS40N.dn,ke_sd.PS40N,dn);
KEsd.ps40m=interp1(PS40M.dn,ke_sd.PS40M,dn);
KEsd.ps40s=interp1(PS40S.dn,ke_sd.PS40S,dn);
KEsd.ps30m=interp1(PS30M.dn,ke_sd.PS30M,dn);
KEsd.ps30s=interp1(PS30S.dn,ke_sd.PS30S,dn);
KEsd.vb30n=interp1(VB30N.dn,ke_sd.VB30N,dn);
KEsd.vb30s=interp1(VB30S.dn,ke_sd.VB30S,dn);
KEsd.oc25na=interp1(OC25NA.dn,ke_sd.OC25NA,dn);
KEsd.oc25nb=interp1(OC25NB.dn,ke_sd.OC25NB,dn);
KEsd.oc25m=interp1(OC25M.dn,ke_sd.OC25M,dn);
KEsd.oc25sb=interp1(OC25SB.dn,ke_sd.OC25SB,dn);
KEsd.oc25sa=interp1(OC25SA.dn,ke_sd.OC25SA,dn);
KEsd.vb25n=interp1(VB25N.dn,ke_sd.VB25N,dn);
KEsd.nrl20n=interp1(NRL20N.dn,ke_sd.NRL20N,dn);
KEsd.nrl20s=interp1(NRL20S.dn,ke_sd.NRL20S,dn);
KEsd.oc10n=interp1(OC10N.dn,ke_sd.OC10N,dn);
KEsd.str3b=interp1(STR3B.dn,ke_sd.STR3B,dn);
KEsd.nw2=interp1(NW2.dn,ke_sd.NW2,dn);
KEsd.rw3=interp1(RW3.dn,ke_sd.RW3,dn);
KEsd.bw1=interp1(BW1.dn,ke_sd.BW1,dn);
KEsd.pw5=interp1(PW5.dn,ke_sd.PW5,dn);
KEsd.str5b=interp1(STR5B.dn,ke_sd.STR5B,dn);
KEsd.str6b=interp1(STR6B.dn,ke_sd.STR6B,dn);

ii=find(~isnan(KEsd.oc50) & ~isnan(KEsd.nrl50n) & ~isnan(KEsd.ps50) & ...
    ~isnan(KEsd.vb50n) & ~isnan(KEsd.vb50s) & ~isnan(KEsd.nrl50s));

if 0
    close all;
    figure;
    plot(1:6,[nanmean(KEsd.oc50(ii)) nanmean(KEsd.nrl50n(ii)) nanmean(KEsd.ps50(ii))...
        nanmean(KEsd.vb50n(ii)) nanmean(KEsd.vb50s(ii)) nanmean(KEsd.nrl50s(ii))],...
        '.','markersize',20);
    set(gca,'xtick',1:1:6,'xticklabels',[{'OC50'};{'NRL50N'};{'PS50'};{'VB50N'};{'VB50S'};{'NRL50S'}]);
    ylabel('time averaged, depth avg KE(sd)');
end

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'KineticEnergy/TimeAvgKE_sd']);
end
%% 

load('/Volumes/InnerShelf1/OC1709A/bathy/PtSalcoast.mat');
bathy=load('/Volumes/InnerShelf1/Moorings/gathered_grids.mat');
load('/Volumes/InnerShelf1/Moorings/MoorLocs_AlongshoreVar.mat');
[coast.x, coast.y]=lltoxy_PtSal(coast.lat,coast.lon);
coast.patchx=[10*1000;  coast.x;  10*1000]; 
coast.patchy=[-31259.1251413617; coast.y; 38034.4551997054]; 
[bathy.x, bathy.y]=lltoxy_PtSal(bathy.G.g200m.lat,bathy.G.g200m.lon);
[x, y]=lltoxy_PtSal(lat,lon);

%% Average over time period all the moorings are in the water 

ii=find(~isnan(KEsd.ms100) & ~isnan(KEsd.oc50) & ~isnan(KEsd.nrl50n) & ~isnan(KEsd.ps50) & ...
    ~isnan(KEsd.vb50n) & ~isnan(KEsd.vb50s) & ~isnan(KEsd.nrl50s) & ~isnan(KEsd.oc40n)...
    & ~isnan(KEsd.oc40s) & ~isnan(KEsd.ps40n) & ~isnan(KEsd.ps40m)...
    & ~isnan(KEsd.ps40s)  & ~isnan(KEsd.ps30m) & ~isnan(KEsd.ps30s) ...
    & ~isnan(KEsd.vb30n) & ~isnan(KEsd.vb30s) & ~isnan(KEsd.oc25na)...
    & ~isnan(KEsd.oc25nb) & ~isnan(KEsd.oc25m) & ~isnan(KEsd.oc25sa)...
    & ~isnan(KEsd.vb25n) & ~isnan(KEsd.nrl20n) & ~isnan(KEsd.nrl20s)...
    & ~isnan(KEsd.oc10n) & ~isnan(KEsd.str3b) & ~isnan(KEsd.nw2)...
    & ~isnan(KEsd.rw3) & ~isnan(KEsd.bw1) & ~isnan(KEsd.pw5)...
    & ~isnan(KEsd.str5b) & ~isnan(KEsd.str6b));

KEsd_timemean=nan*x;
KEsd_timemean(1)=nanmean(KEsd.ms100(ii));
KEsd_timemean(2)=nanmean(KEsd.oc50(ii));
KEsd_timemean(3)=nanmean(KEsd.nrl50n(ii));
KEsd_timemean(4)=nanmean(KEsd.ps50(ii));
KEsd_timemean(5)=nanmean(KEsd.vb50n(ii));
KEsd_timemean(6)=nanmean(KEsd.vb50s(ii));
KEsd_timemean(7)=nanmean(KEsd.nrl50s(ii));
KEsd_timemean(8)=nanmean(KEsd.oc40n(ii));
KEsd_timemean(9)=nanmean(KEsd.oc40s(ii));
KEsd_timemean(10)=nanmean(KEsd.ps40n(ii));
KEsd_timemean(11)=nanmean(KEsd.ps40m(ii));
KEsd_timemean(12)=nanmean(KEsd.ps40s(ii));
KEsd_timemean(19)=nanmean(KEsd.ps30m(ii));
KEsd_timemean(20)=nanmean(KEsd.ps30s(ii));
KEsd_timemean(21)=nanmean(KEsd.vb30n(ii));
KEsd_timemean(22)=nanmean(KEsd.vb30s(ii));
KEsd_timemean(23)=nanmean(KEsd.oc25na(ii));
KEsd_timemean(24)=nanmean(KEsd.oc25nb(ii));
KEsd_timemean(25)=nanmean(KEsd.oc25m(ii));
KEsd_timemean(26)=nanmean(KEsd.oc25sb(ii));
KEsd_timemean(27)=nanmean(KEsd.oc25sa(ii));
KEsd_timemean(28)=nanmean(KEsd.vb25n(ii));
KEsd_timemean(30)=nanmean(KEsd.nrl20n(ii));
KEsd_timemean(31)=nanmean(KEsd.nrl20s(ii));
KEsd_timemean(34)=nanmean(KEsd.oc10n(ii));
KEsd_timemean(35)=nanmean(KEsd.str3b(ii));
KEsd_timemean(36)=nanmean(KEsd.nw2(ii));
KEsd_timemean(37)=nanmean(KEsd.rw3(ii));
KEsd_timemean(38)=nanmean(KEsd.bw1(ii));
KEsd_timemean(39)=nanmean(KEsd.pw5(ii));
KEsd_timemean(41)=nanmean(KEsd.str5b(ii));
KEsd_timemean(42)=nanmean(KEsd.str6b(ii));

ind=ii;
xx=ind(end);


%% check the time peridod we are averaging over 
close all 
figure('position',[64    54   946   751]);
plot(dn(ind), KEsd.ms100(ind));
hold on 
plot(dn(ind), KEsd.oc50(ind)+10);
plot(dn(ind), KEsd.nrl50n(ind)+10*2);
plot(dn(ind), KEsd.ps50(ind)+10*3);
plot(dn(ind), KEsd.vb50n(ind)+10*4);
plot(dn(ind), KEsd.vb50s(ind)+10*5);
plot(dn(ind), KEsd.nrl50s(ind)+10*6);
plot(dn(ind), KEsd.oc40n(ind)+10*7);
plot(dn(ind), KEsd.oc40s(ind)+10*8);
plot(dn(ind), KEsd.ps40n(ind)+10*9);
plot(dn(ind), KEsd.ps40m(ind)+10*10);
plot(dn(ind), KEsd.ps40s(ind)+10*11);
plot(dn(ind), KEsd.ps30m(ind)+10*12);
plot(dn(ind), KEsd.ps30s(ind)+10*13);
plot(dn(ind), KEsd.vb30n(ind)+10*14);
plot(dn(ind), KEsd.vb30s(ind)+10*15);
plot(dn(ind), KEsd.oc25na(ind)+10*16);
plot(dn(ind), KEsd.oc25nb(ind)+10*17);
plot(dn(ind), KEsd.oc25m(ind)+10*18);
plot(dn(ind), KEsd.oc25sb(ind)+10*19); % no data in this time frame 
plot(dn(ind), KEsd.oc25sa(ind)+10*20);
plot(dn(ind), KEsd.vb25n(ind)+10*21);
plot(dn(ind), KEsd.nrl20n(ind)+10*22);
plot(dn(ind), KEsd.nrl20s(ind)+10*23);
plot(dn(ind), KEsd.oc10n(ind)+10*24);
plot(dn(ind), KEsd.str3b(ind)+10*25);
plot(dn(ind), KEsd.nw2(ind)+10*26);
plot(dn(ind), KEsd.rw3(ind)+10*27);
plot(dn(ind), KEsd.bw1(ind)+10*28);
plot(dn(ind), KEsd.pw5(ind)+10*29);
plot(dn(ind), KEsd.str5b(ind)+10*30);
plot(dn(ind), KEsd.str6b(ind)+10*31);
datetick('x','keepticks');
set(gca,'fontweight','bold','fontsize',12);
ylabel('Semidiurnal Kinetic Energy (J m^-^3 )');

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'Timeseries_Allmoorings']);
end


%% Map of time averaged KE denisty 
close all

figure('position',[440   100   479   698]);
hold on
patch(coast.patchx/1000,coast.patchy/1000,[1 1 1].*.8,'HandleVisibility','on');
plot(coast.x/1000,coast.y/1000,'-','color',[1 1 1].*0.5,'linewidth',2);
[c,h] = contour(bathy.x/1000,bathy.y/1000,-bathy.G.g200m.h,[-120 -100:10:10],...
    'color',[1 1 1].*0.8,'HandleVisibility','off','labelspacing',375);

h1=scatter(x/1000,y/1000,2000,KEsd_timemean,'.');
c1=colorbar; colormap(flipud(get(gca,'colormap'))); caxis([0.5 5.5]);
set(c1,'position',[0.73121085594989         0.563037249283667      ...
    0.0579331941544898         0.2753030073009]);
c1.Label.String='(J m^-^3)';

yaxis(-15,15);xaxis(-15,10);
aspect_jack(gca,lat(1));
title(['Depth- and Time- Averaged Semidiurnal Kinetic Energy (J m^-^3)' char(10) ...
    datestr(dn(ind(1)),'dd') ' - ' datestr(dn(xx),'dd') ' September 2017'])
box on
set(gca,'fontsize',12,'fontweight','bold','layer','top','TickDir', 'both',...
    'TickLength',[0.0075 0.0075]);
ylabel('Northing (km)');
xlabel('Easting (km)');

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'TimeAvgKE_sd_MAP_91919']);
end
    
%% 
dn=datenum(2017,9,6,0,0,0):1/24:datenum(2017,11,3,0,0,0);

tempavg.oc50=interp1(OC50.dn,nanmean(OC50.Tsig),dn);
temp_s.oc50=interp1(OC50.dn,OC50.Tsig(end,:),dn);
temp_b.oc50=interp1(OC50.dn,OC50.Tsig(1,:),dn);
tempavg.nrl50n=interp1(NRL50N.dn,nanmean(NRL50N.Tsig),dn);
temp_s.nrl50n=interp1(NRL50N.dn,NRL50N.Tsig(end,:),dn);
temp_b.nrl50n=interp1(NRL50N.dn,NRL50N.Tsig(1,:),dn);
tempavg.ps50=interp1(PS50.dn,nanmean(PS50.Tsig),dn);
temp_s.ps50=interp1(PS50.dn,PS50.Tsig(end,:),dn);
temp_b.ps50=interp1(PS50.dn,PS50.Tsig(1,:),dn);
tempavg.vb50n=interp1(VB50N.dn,nanmean(VB50N.Tsig),dn);
temp_s.vb50n=interp1(VB50N.dn,VB50N.Tsig(end,:),dn);
temp_b.vb50n=interp1(VB50N.dn,VB50N.Tsig(1,:),dn);
tempavg.vb50s=interp1(VB50S.dn,nanmean(VB50S.Tsig),dn);
temp_s.vb50s=interp1(VB50S.dn,VB50S.Tsig(end,:),dn);
temp_b.vb50s=interp1(VB50S.dn,VB50S.Tsig(1,:),dn);
tempavg.nrl50s=interp1(NRL50S.dn,nanmean(NRL50S.Tsig),dn);
temp_s.nrl50s=interp1(NRL50S.dn,NRL50S.Tsig(end,:),dn);
temp_b.nrl50s=interp1(NRL50S.dn,NRL50S.Tsig(1,:),dn);


%% Conditions for JPO paper
dn1=datenum(2017,9,6);
dn2=datenum(2017,11,2);
dt=7;


close all;
figure('position',[ 185    52   694   753]);
hax=tight_subplot(6,1,[0.02 0.1],[0.04 0.02],[0.13 0.02]);

axes(hax(1));
% plot(dn_atm,wspd,'-k','linewidth',2); hold on 
% plot(dn_atm,wdir,'.k','linewidth',2); hold on 
ylim=20;
patch('XData',[datenum(2017,9,10,21,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,10,21,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',[0.94 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
hold on;
patch('XData',[datenum(2017,9,28,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,9,28,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,15,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,15,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,25,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,25,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));quiver(dn_atm2,dn_atm2.*0,weast,wnorth,'k'); %,'autoscale','off');%,'maxheadsize',0)
% quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,24),pl33tn(wind_n,1,24),'k')
% quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,33),pl33tn(wind_n,1,33),'k','autoscale','off','maxheadsize',0)
xaxis(dn1,dn2);
yaxis(-2,2);
ylabel(['Wind speed' char(10) '(m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',-2:1:2,'yticklabels',[-2:1:2]*7.5,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 

axes(hax(2));
ylim=1.25;
patch('XData',[datenum(2017,9,10,21,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,10,21,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',[0.94 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
hold on     
patch('XData',[datenum(2017,9,28,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,9,28,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,15,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,15,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,25,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,25,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
plot(tide.dn,fillmissing(tide.wl,'linear'),'-','color',[0.65 0.65 0.65],'linewidth',1); hold on 
xaxis(dn1,dn2);
yaxis(-1.25,1.25);
ylabel('BT tide (m)');
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',-1:1:1,'yticklabels',-1:1:1,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 


axes(hax(3));
ylim=14;
patch('XData',[datenum(2017,9,10,21,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,10,21,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',[0.94 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
hold on  
patch('XData',[datenum(2017,9,28,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,9,28,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,15,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,15,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,25,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,25,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));% filt=33;
filt=24*4;
plot(OC50.dn,pl33tn(ke_sd.OC50,1/60,filt),'k','linewidth',2);
plot(NRL50N.dn,pl33tn(ke_sd.NRL50N,1/60,filt),'r','linewidth',2);
plot(PS50.dn,pl33tn(ke_sd.PS50,1/60,filt),'color',hex2rgb('#1682E0'),'linewidth',2);
plot(VB50N.dn,pl33tn(ke_sd.VB50N,1/60,filt),'color',hex2rgb('#45C162'),'linewidth',2);
plot(VB50S.dn,pl33tn(ke_sd.VB50S,1/60,filt),'color',[0.7 0.7 0.7],'linewidth',2);
plot(NRL50S.dn,pl33tn(ke_sd.NRL50S,1/60,filt),'color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2); yaxis(0,9);
ylabels = 2:2:14;
ylabels(2:2:end,:) = nan; % remove every other one
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',2:2:14,'yticklabel',ylabels,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel(['Semidiurnal KE' char(10) '(J m^-^3 )'])
% ylabel(['Semidiurnal $\overline{KE}$' char(10) '(kg/m s2)'],'interpreter','latex','fontsize',14)

axes(hax(4));
ylim=24;
patch('XData',[datenum(2017,9,10,21,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,10,21,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',[0.94 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
hold on  
patch('XData',[datenum(2017,9,28,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,9,28,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,15,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,15,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,25,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,25,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
filt=33;
plot(dn,pl33tn(temp_s.oc50-temp_b.oc50,1,filt),'k','linewidth',2);
plot(dn,pl33tn(temp_s.nrl50n-temp_b.nrl50n,1,filt),'r','linewidth',2);
plot(dn,pl33tn(temp_s.ps50-temp_b.ps50,1,filt),'color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,pl33tn(temp_s.vb50n-temp_b.vb50n,1,filt),'color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,pl33tn(temp_s.vb50s-temp_b.vb50s,1,filt),'color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,pl33tn(temp_s.nrl50s-temp_b.nrl50s,1,filt),'color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2); yaxis(0, 8);
ylabels = (0:1:8);
ylabels(2:2:end) = nan; % remove every other one
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',2:2:6,'yticklabel',2:2:6,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel(['Stratification' char(10) '\DeltaT'])

axes(hax(5));
patch('XData',[datenum(2017,9,10,21,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,10,21,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',[0.94 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
hold on   
patch('XData',[datenum(2017,9,28,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,9,28,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,15,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,15,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,25,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,25,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
h1=plot(OC50.DNrho,OC50.Cw,'-k','linewidth',2);  
h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.Cw,'-r','linewidth',2);  
h2b=plot(PS50.subtidal.DNrho,PS50.subtidal.Cw,'-b','color',hex2rgb('#1682E0'),'linewidth',2);  
h3=plot(VB50N.subtidal.DNrho,VB50N.subtidal.Cw,'-','color',hex2rgb('#45C162'),'linewidth',2); 
h4=plot(VB50S.subtidal.DNrho,VB50S.subtidal.Cw,'-','color',[0.7 0.7 0.7],'linewidth',2); 
h5=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.Cw,'-','color',hex2rgb('#FF9456'),'linewidth',2); 
ylabel (['c (m s^-^1)' ]);
grid on; hh=gca; hh.XMinorGrid='on'; 

% h1=plot(OC50.DNrho,OC50.Cw./OC50.alpha,'-k','linewidth',2);  
% h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.Cw./NRL50N.subtidal.alpha,'-r','linewidth',2);  
% h2b=plot(PS50.subtidal.DNrho,PS50.subtidal.Cw./PS50.subtidal.alpha,'-b','color',hex2rgb('#1682E0'),'linewidth',2);  
% h3=plot(VB50N.subtidal.DNrho,VB50N.subtidal.Cw./VB50N.subtidal.alpha,'-','color',hex2rgb('#45C162'),'linewidth',2); 
% h4=plot(VB50S.subtidal.DNrho,VB50S.subtidal.Cw./VB50S.subtidal.alpha,'-','color',[0.7 0.7 0.7],'linewidth',2); 
% h5=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.Cw./NRL50S.subtidal.alpha,'-','color',hex2rgb('#FF9456'),'linewidth',2); 
% ylabel (['c/\alpha (m)' ]);
xaxis(dn1,dn2);
yaxis(0.13,.33);
datetick('x','keeplimits','keepticks');
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',.15:.05:.45,'yticklabel',.15:.05:.45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 

axes(hax(6));
ylim=0.07;
patch('XData',[datenum(2017,9,10,21,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,22,6,0,0) datenum(2017,9,10,21,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',[0.94 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
hold on     
patch('XData',[datenum(2017,9,28,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,10,3,0,0,0) datenum(2017,9,28,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,15,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,19,0,0,0) datenum(2017,10,15,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));
patch('XData',[datenum(2017,10,25,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,31,0,0,0) datenum(2017,10,25,0,0,0)], ...
    'YData',[-ylim; -ylim; ylim; ylim], ...
    'FaceColor',hex2rgb('#F1FFD6'),'edgecolor',hex2rgb('#F1FFD6'));plot([datenum(2017,9,5) datenum(2017,11,3)],[0 0],'--k','linewidth',1) ;hold on 
h1=plot(OC50.DNrho,OC50.alpha,'-k','linewidth',2); 
h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.alpha,'-r','linewidth',2); 
h2b=plot(PS50.subtidal.DNrho,PS50.subtidal.alpha,'b-','color',hex2rgb('#1682E0'),'linewidth',2);  
h3=plot(VB50N.subtidal.DNrho,VB50N.subtidal.alpha,'-g','color',hex2rgb('#45C162'),'linewidth',2); 
h4=plot(VB50S.subtidal.DNrho,VB50S.subtidal.alpha,'-','color',[0.7 0.7 0.7],'linewidth',2); 
h5=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.alpha,'-','color',hex2rgb('#FF9456'),'linewidth',2); 
ylabel ('\alpha (s^-^1)');
yaxis(-0.01,0.015);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',-0.02:0.005:0.02,'yticklabels',-0.02:0.005:0.02,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);

% h1=plot(OC50.DNrho,OC50.alpha./OC50.Cw,'-k','linewidth',2); 
% h2=plot(NRL50N.subtidal.DNrho,NRL50N.subtidal.alpha./NRL50N.subtidal.Cw,'-r','linewidth',2); 
% h2b=plot(PS50.subtidal.DNrho,PS50.subtidal.alpha./PS50.subtidal.Cw,'b-','color',hex2rgb('#1682E0'),'linewidth',2);  
% h3=plot(VB50N.subtidal.DNrho,VB50N.subtidal.alpha./VB50N.subtidal.Cw,'-g','color',hex2rgb('#45C162'),'linewidth',2); 
% h4=plot(VB50S.subtidal.DNrho,VB50S.subtidal.alpha./VB50S.subtidal.Cw,'-','color',[0.7 0.7 0.7],'linewidth',2); 
% h5=plot(NRL50S.subtidal.DNrho,NRL50S.subtidal.alpha./NRL50S.subtidal.Cw,'-','color',hex2rgb('#FF9456'),'linewidth',2); 
% ylabel ('\alpha/c_o (m^-^1)');
% % yaxis(-0.04,0.06);
% set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
%     'box','on','layer','top','ytick',-0.04:0.02:0.06,'yticklabels',-0.04:0.02:0.06'TickDir', 'both','TickLength',[0.0075 0.0075]);
xaxis(dn1,dn2);
leg=legend([h1 h2 h2b h3 h4 h5],'OC50','NRL50N','PS50','VB50N','VB50S','NRL50S');
set(leg,'Position',[0.561597305475493         0.437959197375965         0.115994236311239         0.108666089678692],...
    'edgecolor','k','fontsize',10);

grid on; hh=gca; hh.XMinorGrid='on'; 


if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'Conditions']);
end 

% Joule is kg m^2 s^-2

% units here are kg m^-3 m^2 s^2 

%% 
figure; 
quiver(dn_atm2,dn_atm2.*0,weast,wnorth,'k'); %,'autoscale','off');%,'maxheadsize',0)
% quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,33),pl33tn(wind_n,1,33),'k','autoscale','off','maxheadsize',0)
xaxis(dn1,dn2);
yaxis(-2,2);
hold on
plot(dn_atm,-pl33tn(wspd,1,33)/7.5,'r','linewidth',2);
ylabel(['Wind speed' char(10) '(m/s)']);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',-2:1:2,'yticklabels',[-2:1:2]*7.5);
grid on; hh=gca; hh.XMinorGrid='on'; 
ax1=gca;
set(ax1,'color','none');

% axes;
% ax2=gca;
% set(ax2,'position',get(ax1,'position'));
% set(ax2,'color','none');
% plot(dn_atm,pl33tn(wspd,1,33),'r','linewidth',2);
% xaxis(dn1,dn2);


%% 

% dn1=datenum(2017,9,6);
% dn2=datenum(2017,11,2);
% dt=7;

 close all 
 
% dn1=datenum(2017,9,11);
% dn2=datenum(2017,9,18);
dn1=dn(ind(1));
dn2=dn(ind(end));
dt=0.5;

labels=datestr(dn1:dt:dn2,'mm/dd');
labels(2:2:end,:) = nan; % remove every other one

figure('position',[440   100   923   698]);
hax=tight_subplot(6,1,[0.04 0.1],[0.04 0.02],[0.05 0.02]);

axes(hax(1))
plot(OC50.dn,ke_sd.OC50,'-','color',0.7*[1 1 1],'linewidth',1.5); 
hold on
plot(OC50.dn,pl33tn(ke_sd.OC50,1/60,1),'k','linewidth',3); 
yaxis(0,50);xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',15:15:45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on

axes(hax(2))
plot(NRL50N.dn,ke_sd.NRL50N,'-','color',0.7*[1 1 1],'linewidth',1.5);
hold on
plot(NRL50N.dn,pl33tn(ke_sd.NRL50N,1/60,1),'r','linewidth',3);
yaxis(0,50);xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',15:15:45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on

axes(hax(3))
plot(PS50.dn,ke_sd.PS50,'-','color',0.7*[1 1 1],'linewidth',1.5);
hold on 
plot(PS50.dn,pl33tn(ke_sd.PS50,1/60,1),'color',hex2rgb('#1682E0'),'linewidth',3);
yaxis(0,50);xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',15:15:45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on

axes(hax(4))
plot(VB50N.dn,ke_sd.VB50N,'-','color',0.7*[1 1 1],'linewidth',1.5);
hold on; 
plot(VB50N.dn,pl33tn(ke_sd.VB50N,1/60,1),'color',hex2rgb('#45C162'),'linewidth',3);
yaxis(0,50);xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',15:15:45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on

axes(hax(5))
plot(VB50S.dn,ke_sd.VB50S,'-','color',0.7*[1 1 1],'linewidth',1.5);
hold on 
plot(VB50S.dn,pl33tn(ke_sd.VB50S,1/60,1),'color',0.5*[ 1 1 1],'linewidth',3);
yaxis(0,50);xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',15:15:45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on

axes(hax(6))
plot(NRL50S.dn,ke_sd.NRL50S,'-','color',0.7*[1 1 1],'linewidth',1.5);
hold on 
plot(NRL50S.dn,pl33tn(ke_sd.NRL50S,1/60,1),'color',hex2rgb('#FF9456'),'linewidth',3);
yaxis(0,50);xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',labels,...
    'box','on','layer','top','ytick',15:15:45,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir 'KEsd_timeseries10719']);
end 

%% 

dn=datenum(2017,9,6,0,0,0):1/24/60:datenum(2017,11,3,0,0,0);

ii=find(~isnan(MS100.DNrho));
Alpha.ms100=interp1(MS100.DNrho(ii),MS100.alpha(ii),dn);
speed.ms100=interp1(MS100.DNrho(ii),MS100.Cw(ii),dn);
ii=find(isnan(MS100.DNrho));
jj=findnearest_JACK(dn,MS100.DNrho(107)+1/24):findnearest_JACK(dn,MS100.DNrho(117)-1/24);
Alpha.ms100(jj)=nan;
speed.ms100(jj)=nan;
tempavg.ms100=interp1(MS100.dn,nanmean(MS100.Tsig),dn);
temp_s.ms100=interp1(MS100.dn,MS100.Tsig(end,:),dn);
temp_b.ms100=interp1(MS100.dn,MS100.Tsig(1,:),dn);


ii=find(~isnan(OC50.DNrho));
Alpha.oc50=interp1(OC50.DNrho(ii),OC50.alpha(ii),dn);
speed.oc50=interp1(OC50.DNrho(ii),OC50.Cw(ii),dn);
ii=find(isnan(OC50.DNrho));
jj=findnearest_JACK(dn,OC50.DNrho(108)+1/24):findnearest_JACK(dn,OC50.DNrho(118)-1/24);
Alpha.oc50(jj)=nan;
speed.oc50(jj)=nan;
tempavg.oc50=interp1(OC50.dn,nanmean(OC50.Tsig),dn);
temp_s.oc50=interp1(OC50.dn,OC50.Tsig(end,:),dn);
temp_b.oc50=interp1(OC50.dn,OC50.Tsig(1,:),dn);


Alpha.nrl50n=interp1(NRL50N.subtidal.DNrho,NRL50N.subtidal.alpha,dn);
Alpha.ps50=interp1(PS50.subtidal.DNrho,PS50.subtidal.alpha,dn);
Alpha.vb50n=interp1(VB50N.subtidal.DNrho,VB50N.subtidal.alpha,dn);
Alpha.vb50s=interp1(VB50S.subtidal.DNrho,VB50S.subtidal.alpha,dn);
Alpha.nrl50s=interp1(NRL50S.subtidal.DNrho,NRL50S.subtidal.alpha,dn);

speed.nrl50n=interp1(NRL50N.subtidal.DNrho,NRL50N.subtidal.Cw,dn);
speed.ps50=interp1(PS50.subtidal.DNrho,PS50.subtidal.Cw,dn);
speed.vb50n=interp1(VB50N.subtidal.DNrho,VB50N.subtidal.Cw,dn);
speed.vb50s=interp1(VB50S.subtidal.DNrho,VB50S.subtidal.Cw,dn);
speed.nrl50s=interp1(NRL50S.subtidal.DNrho,NRL50S.subtidal.Cw,dn);

tempavg.nrl50n=interp1(NRL50N.dn,nanmean(NRL50N.Tsig),dn);
temp_s.nrl50n=interp1(NRL50N.dn,NRL50N.Tsig(end,:),dn);
temp_b.nrl50n=interp1(NRL50N.dn,NRL50N.Tsig(1,:),dn);
tempavg.ps50=interp1(PS50.dn,nanmean(PS50.Tsig),dn);
temp_s.ps50=interp1(PS50.dn,PS50.Tsig(end,:),dn);
temp_b.ps50=interp1(PS50.dn,PS50.Tsig(1,:),dn);
tempavg.vb50n=interp1(VB50N.dn,nanmean(VB50N.Tsig),dn);
temp_s.vb50n=interp1(VB50N.dn,VB50N.Tsig(end,:),dn);
temp_b.vb50n=interp1(VB50N.dn,VB50N.Tsig(1,:),dn);
tempavg.vb50s=interp1(VB50S.dn,nanmean(VB50S.Tsig),dn);
temp_s.vb50s=interp1(VB50S.dn,VB50S.Tsig(end,:),dn);
temp_b.vb50s=interp1(VB50S.dn,VB50S.Tsig(1,:),dn);
tempavg.nrl50s=interp1(NRL50S.dn,nanmean(NRL50S.Tsig),dn);
temp_s.nrl50s=interp1(NRL50S.dn,NRL50S.Tsig(end,:),dn);
temp_b.nrl50s=interp1(NRL50S.dn,NRL50S.Tsig(1,:),dn);

ii=find(~isnan(OC40N.DNrho));
Alpha.oc40n=interp1(OC40N.DNrho(ii),OC40N.alpha(ii),dn);
speed.oc40n=interp1(OC40N.DNrho(ii),OC40N.Cw(ii),dn);
ii=find(isnan(OC40N.DNrho));
jj=findnearest_JACK(dn,OC40N.DNrho(125)+1/24):findnearest_JACK(dn,OC40N.DNrho(136)-1/24);
Alpha.oc40n(jj)=nan;
speed.oc40n(jj)=nan;
tempavg.oc40n=interp1(OC40N.dn,nanmean(OC40N.Tsig),dn);
temp_s.oc40n=interp1(OC40N.dn,OC40N.Tsig(end,:),dn);
temp_b.oc40n=interp1(OC40N.dn,OC40N.Tsig(1,:),dn);

ii=find(~isnan(OC40S.DNrho));
Alpha.oc40s=interp1(OC40S.DNrho(ii),OC40S.alpha(ii),dn);
speed.oc40s=interp1(OC40S.DNrho(ii),OC40S.Cw(ii),dn);
ii=find(isnan(OC40S.DNrho));
jj=findnearest_JACK(dn,OC40S.DNrho(113)+1/24):findnearest_JACK(dn,OC40S.DNrho(121)-1/24);
Alpha.oc40s(jj)=nan;
speed.oc40s(jj)=nan;
tempavg.oc40s=interp1(OC40S.dn,nanmean(OC40S.Tsig),dn);
temp_s.oc40s=interp1(OC40S.dn,OC40S.Tsig(end,:),dn);
temp_b.oc40s=interp1(OC40S.dn,OC40S.Tsig(1,:),dn);


Alpha.ps40n=interp1(PS40N.subtidal.DNrho,PS40N.subtidal.alpha,dn);
Alpha.ps40m=interp1(PS40M.subtidal.DNrho,PS40M.subtidal.alpha,dn);
Alpha.ps40s=interp1(PS40S.subtidal.DNrho,PS40S.subtidal.alpha,dn);
Alpha.nrl35n=interp1(NRL35N.subtidal.DNrho,NRL35N.subtidal.alpha,dn);
Alpha.ps35m=interp1(PS35M.subtidal.DNrho,PS35M.subtidal.alpha,dn);
Alpha.nrl35s=interp1(NRL35S.subtidal.DNrho,NRL35S.subtidal.alpha,dn);

speed.ps40n=interp1(PS40N.subtidal.DNrho,PS40N.subtidal.Cw,dn);
speed.ps40m=interp1(PS40M.subtidal.DNrho,PS40M.subtidal.Cw,dn);
speed.ps40s=interp1(PS40S.subtidal.DNrho,PS40S.subtidal.Cw,dn);
speed.nrl35n=interp1(NRL35N.subtidal.DNrho,NRL35N.subtidal.Cw,dn);
speed.ps35m=interp1(PS35M.subtidal.DNrho,PS35M.subtidal.Cw,dn);
speed.nrl35s=interp1(NRL35S.subtidal.DNrho,NRL35S.subtidal.Cw,dn);

tempavg.ps40n=interp1(PS40N.dn,nanmean(PS40N.Tsig),dn);
temp_s.ps40n=interp1(PS40N.dn,PS40N.Tsig(end,:),dn);
temp_b.ps40n=interp1(PS40N.dn,PS40N.Tsig(1,:),dn);
tempavg.ps40m=interp1(PS40M.dn,nanmean(PS40M.Tsig),dn);
temp_s.ps40m=interp1(PS40M.dn,PS40M.Tsig(end,:),dn);
temp_b.ps40m=interp1(PS40M.dn,PS40M.Tsig(1,:),dn);
tempavg.ps40s=interp1(PS40S.dn,nanmean(PS40S.Tsig),dn);
temp_s.ps40s=interp1(PS40S.dn,PS40S.Tsig(end,:),dn);
temp_b.ps40s=interp1(PS40S.dn,PS40S.Tsig(1,:),dn);
tempavg.nrl35n=interp1(NRL35N.dn,nanmean(NRL35N.Tsig),dn);
temp_s.nrl35n=interp1(NRL35N.dn,NRL35N.Tsig(end,:),dn);
temp_b.nrl35n=interp1(NRL35N.dn,NRL35N.Tsig(1,:),dn);
tempavg.ps35m=interp1(PS35M.dn,nanmean(PS35M.Tsig),dn);
temp_s.ps35m=interp1(PS35M.dn,PS35M.Tsig(end,:),dn);
temp_b.ps35m=interp1(PS35M.dn,PS35M.Tsig(1,:),dn);
tempavg.nrl35s=interp1(NRL35S.dn,nanmean(NRL35S.Tsig),dn);
temp_s.nrl35s=interp1(NRL35S.dn,NRL35S.Tsig(end,:),dn);
temp_b.nrl35s=interp1(NRL35S.dn,NRL35S.Tsig(1,:),dn);

ii=find(~isnan(OC32N.DNrho));
Alpha.oc32n=interp1(OC32N.DNrho(ii),OC32N.alpha(ii),dn);
speed.oc32n=interp1(OC32N.DNrho(ii),OC32N.Cw(ii),dn);
tempavg.oc32n=interp1(OC32N.dn,nanmean(OC32N.Tsig),dn);
temp_s.oc32n=interp1(OC32N.dn,OC32N.Tsig(end,:),dn);
temp_b.oc32n=interp1(OC32N.dn,OC32N.Tsig(1,:),dn);

ii=find(~isnan(OC32S.DNrho));
Alpha.oc32s=interp1(OC32S.DNrho(ii),OC32S.alpha(ii),dn);
speed.oc32s=interp1(OC32S.DNrho(ii),OC32S.Cw(ii),dn);
tempavg.oc32s=interp1(OC32S.dn,nanmean(OC32S.Tsig),dn);
temp_s.oc32s=interp1(OC32S.dn,OC32S.Tsig(end,:),dn);
temp_b.oc32s=interp1(OC32S.dn,OC32S.Tsig(1,:),dn);

Alpha.ps30n=interp1(PS30N.subtidal.DNrho,PS30N.subtidal.alpha,dn);
Alpha.ps30m=interp1(PS30M.subtidal.DNrho,PS30M.subtidal.alpha,dn);
Alpha.ps30s=interp1(PS30S.subtidal.DNrho,PS30S.subtidal.alpha,dn);
Alpha.vb30n=interp1(VB30N.subtidal.DNrho,VB30N.subtidal.alpha,dn);
Alpha.vb30s=interp1(VB30S.subtidal.DNrho,VB30S.subtidal.alpha,dn);

speed.ps30n=interp1(PS30N.subtidal.DNrho,PS30N.subtidal.Cw,dn);
speed.ps30m=interp1(PS30M.subtidal.DNrho,PS30M.subtidal.Cw,dn);
speed.ps30s=interp1(PS30S.subtidal.DNrho,PS30S.subtidal.Cw,dn);
speed.vb30n=interp1(VB30N.subtidal.DNrho,VB30N.subtidal.Cw,dn);
speed.vb30s=interp1(VB30S.subtidal.DNrho,VB30S.subtidal.Cw,dn);

tempavg.ps30n=interp1(PS30N.dn,nanmean(PS30N.Tsig),dn);
temp_s.ps30n=interp1(PS30N.dn,PS30N.Tsig(end,:),dn);
temp_b.ps30n=interp1(PS30N.dn,PS30N.Tsig(1,:),dn);
tempavg.ps30m=interp1(PS30M.dn,nanmean(PS30M.Tsig),dn);
temp_s.ps30m=interp1(PS30M.dn,PS30M.Tsig(end,:),dn);
temp_b.ps30m=interp1(PS30M.dn,PS30M.Tsig(1,:),dn);
tempavg.ps30s=interp1(PS30S.dn,nanmean(PS30S.Tsig),dn);
temp_s.ps30s=interp1(PS30S.dn,PS30S.Tsig(end,:),dn);
temp_b.ps30s=interp1(PS30S.dn,PS30S.Tsig(1,:),dn);
tempavg.vb30n=interp1(VB30N.dn,nanmean(VB30N.Tsig),dn);
temp_s.vb30n=interp1(VB30N.dn,VB30N.Tsig(end,:),dn);
temp_b.vb30n=interp1(VB30N.dn,VB30N.Tsig(1,:),dn);
tempavg.vb30s=interp1(VB30S.dn,nanmean(VB30S.Tsig),dn);
temp_s.vb30s=interp1(VB30S.dn,VB30S.Tsig(end,:),dn);
temp_b.vb30s=interp1(VB30S.dn,VB30S.Tsig(1,:),dn);

ii=find(~isnan(OC25NA.DNrho));
Alpha.oc25na=interp1(OC25NA.DNrho(ii),OC25NA.alpha(ii),dn);
speed.oc25na=interp1(OC25NA.DNrho(ii),OC25NA.Cw(ii),dn);
tempavg.oc25na=interp1(OC25NA.dn,nanmean(OC25NA.Tsig),dn);
temp_s.oc25na=interp1(OC25NA.dn,OC25NA.Tsig(end,:),dn);
temp_b.oc25na=interp1(OC25NA.dn,OC25NA.Tsig(1,:),dn);

ii=find(~isnan(OC25NB.DNrho));
Alpha.oc25nb=interp1(OC25NB.DNrho(ii),OC25NB.alpha(ii),dn);
speed.oc25nb=interp1(OC25NB.DNrho(ii),OC25NB.Cw(ii),dn);
tempavg.oc25nb=interp1(OC25NB.dn,nanmean(OC25NB.Tsig),dn);
temp_s.oc25nb=interp1(OC25NB.dn,OC25NB.Tsig(end,:),dn);
temp_b.oc25nb=interp1(OC25NB.dn,OC25NB.Tsig(1,:),dn);

ii=find(~isnan(OC25M.DNrho));
Alpha.oc25m=interp1(OC25M.DNrho(ii),OC25M.alpha(ii),dn);
speed.oc25m=interp1(OC25M.DNrho(ii),OC25M.Cw(ii),dn);
tempavg.oc25m=interp1(OC25M.dn,nanmean(OC25M.Tsig),dn);
temp_s.oc25m=interp1(OC25M.dn,OC25M.Tsig(end,:),dn);
temp_b.oc25m=interp1(OC25M.dn,OC25M.Tsig(1,:),dn);

ii=find(~isnan(OC25SB.DNrho));
Alpha.oc25sb=interp1(OC25SB.DNrho(ii),OC25SB.alpha(ii),dn);
speed.oc25sb=interp1(OC25SB.DNrho(ii),OC25SB.Cw(ii),dn);
tempavg.oc25sb=interp1(OC25SB.dn,nanmean(OC25SB.Tsig),dn);
temp_s.oc25sb=interp1(OC25SB.dn,OC25SB.Tsig(end,:),dn);
temp_b.oc25sb=interp1(OC25SB.dn,OC25SB.Tsig(1,:),dn);

ii=find(~isnan(OC25SA.DNrho));
Alpha.oc25sa=interp1(OC25SA.DNrho(ii),OC25SA.alpha(ii),dn);
speed.oc25sa=interp1(OC25SA.DNrho(ii),OC25SA.Cw(ii),dn);
tempavg.oc25sa=interp1(OC25SA.dn,nanmean(OC25SA.Tsig),dn);
temp_s.oc25sa=interp1(OC25SA.dn,OC25SA.Tsig(end,:),dn);
temp_b.oc25sa=interp1(OC25SA.dn,OC25SA.Tsig(1,:),dn);

Alpha.vb25n=interp1(VB25N.subtidal.DNrho,VB25N.subtidal.alpha,dn);
Alpha.vb25s=interp1(VB25S.subtidal.DNrho,VB25S.subtidal.alpha,dn);
Alpha.nrl20n=interp1(NRL20N.subtidal.DNrho,NRL20N.subtidal.alpha,dn);
Alpha.nrl20s=interp1(NRL20S.subtidal.DNrho,NRL20S.subtidal.alpha,dn);

speed.vb25n=interp1(VB25N.subtidal.DNrho,VB25N.subtidal.Cw,dn);
speed.vb25s=interp1(VB25S.subtidal.DNrho,VB25S.subtidal.Cw,dn);
speed.nrl20n=interp1(NRL20N.subtidal.DNrho,NRL20N.subtidal.Cw,dn);
speed.nrl20s=interp1(NRL20S.subtidal.DNrho,NRL20S.subtidal.Cw,dn);

tempavg.vb25n=interp1(VB25N.dn,nanmean(VB25N.Tsig),dn);
temp_s.vb25n=interp1(VB25N.dn,VB25N.Tsig(end,:),dn);
temp_b.vb25n=interp1(VB25N.dn,VB25N.Tsig(1,:),dn);
tempavg.vb25s=interp1(VB25S.dn,nanmean(VB25S.Tsig),dn);
temp_s.vb25s=interp1(VB25S.dn,VB25S.Tsig(end,:),dn);
temp_b.vb25s=interp1(VB25S.dn,VB25S.Tsig(1,:),dn);
tempavg.nrl20n=interp1(NRL20N.dn,nanmean(NRL20N.Tsig),dn);
temp_s.nrl20n=interp1(NRL20N.dn,NRL20N.Tsig(end,:),dn);
temp_b.nrl20n=interp1(NRL20N.dn,NRL20N.Tsig(1,:),dn);
tempavg.nrl20s=interp1(NRL20S.dn,nanmean(NRL20S.Tsig),dn);
temp_s.nrl20s=interp1(NRL20S.dn,NRL20S.Tsig(end,:),dn);
temp_b.nrl20s=interp1(NRL20S.dn,NRL20S.Tsig(1,:),dn);

ii=find(~isnan(OC17N.DNrho));
Alpha.oc17n=interp1(OC17N.DNrho(ii),OC17N.alpha(ii),dn);
speed.oc17n=interp1(OC17N.DNrho(ii),OC17N.Cw(ii),dn);
tempavg.oc17n=interp1(OC17N.dn,nanmean(OC17N.Tsig),dn);
temp_s.oc17n=interp1(OC17N.dn,OC17N.Tsig(end,:),dn);
temp_b.oc17n=interp1(OC17N.dn,OC17N.Tsig(1,:),dn);

ii=find(~isnan(OC17S.DNrho));
Alpha.oc17s=interp1(OC17S.DNrho(ii),OC17S.alpha(ii),dn);
speed.oc17s=interp1(OC17S.DNrho(ii),OC17S.Cw(ii),dn);
tempavg.oc17s=interp1(OC17S.dn,nanmean(OC17S.Tsig),dn);
temp_s.oc17s=interp1(OC17S.dn,OC17S.Tsig(end,:),dn);
temp_b.oc17s=interp1(OC17S.dn,OC17S.Tsig(1,:),dn);

ii=find(~isnan(OC10N.DNrho));
Alpha.oc10n=interp1(OC10N.DNrho(ii),OC10N.alpha(ii),dn);
speed.oc10n=interp1(OC10N.DNrho(ii),OC10N.Cw(ii),dn);
tempavg.oc10n=interp1(OC10N.dn,nanmean(OC10N.Tsig),dn);
temp_s.oc10n=interp1(OC10N.dn,OC10N.Tsig(end,:),dn);
temp_b.oc10n=interp1(OC10N.dn,OC10N.Tsig(1,:),dn);

ii=find(~isnan(STR3B.DNrho));
Alpha.str3b=interp1(STR3B.DNrho(ii),STR3B.alpha(ii),dn);
speed.str3b=interp1(STR3B.DNrho(ii),STR3B.Cw(ii),dn);
tempavg.str3b=interp1(STR3B.dn,nanmean(STR3B.Tsig),dn);
temp_s.str3b=interp1(STR3B.dn,STR3B.Tsig(end,:),dn);
temp_b.str3b=interp1(STR3B.dn,STR3B.Tsig(1,:),dn);

ii=find(~isnan(NW2.subtidal.DNrho));
Alpha.nw2=interp1(NW2.subtidal.DNrho(ii),NW2.subtidal.alpha(ii),dn);
speed.nw2=interp1(NW2.subtidal.DNrho(ii),NW2.subtidal.Cw(ii),dn);
tempavg.nw2=interp1(NW2.dn,nanmean(NW2.Tsig),dn);
temp_s.nw2=interp1(NW2.dn,NW2.Tsig(end,:),dn);
temp_b.nw2=interp1(NW2.dn,NW2.Tsig(1,:),dn);

ii=find(~isnan(RW3.subtidal.DNrho));
Alpha.rw3=interp1(RW3.subtidal.DNrho(ii),RW3.subtidal.alpha(ii),dn);
speed.rw3=interp1(RW3.subtidal.DNrho(ii),RW3.subtidal.Cw(ii),dn);
tempavg.rw3=interp1(RW3.dn,nanmean(RW3.Tsig),dn);
temp_s.rw3=interp1(RW3.dn,RW3.Tsig(end,:),dn);
temp_b.rw3=interp1(RW3.dn,RW3.Tsig(1,:),dn);

ii=find(~isnan(BW1.subtidal.DNrho));
Alpha.bw1=interp1(BW1.subtidal.DNrho(ii),BW1.subtidal.alpha(ii),dn);
speed.bw1=interp1(BW1.subtidal.DNrho(ii),BW1.subtidal.Cw(ii),dn);
tempavg.bw1=interp1(BW1.dn,nanmean(BW1.Tsig),dn);
temp_s.bw1=interp1(BW1.dn,BW1.Tsig(end,:),dn);
temp_b.bw1=interp1(BW1.dn,BW1.Tsig(1,:),dn);

ii=find(~isnan(PW5.subtidal.DNrho));
Alpha.pw5=interp1(PW5.subtidal.DNrho(ii),PW5.subtidal.alpha(ii),dn);
speed.pw5=interp1(PW5.subtidal.DNrho(ii),PW5.subtidal.Cw(ii),dn);
tempavg.pw5=interp1(PW5.dn,nanmean(PW5.Tsig),dn);
temp_s.pw5=interp1(PW5.dn,PW5.Tsig(end,:),dn);
temp_b.pw5=interp1(PW5.dn,PW5.Tsig(1,:),dn);

ii=find(~isnan(STR4F.subtidal.DNrho));
Alpha.str4f=interp1(STR4F.subtidal.DNrho(ii),STR4F.subtidal.alpha(ii),dn);
speed.str4f=interp1(STR4F.subtidal.DNrho(ii),STR4F.subtidal.Cw(ii),dn);
tempavg.str4f=interp1(STR4F.dn,nanmean(STR4F.Tsig),dn);
temp_s.str4f=interp1(STR4F.dn,STR4F.Tsig(end,:),dn);
temp_b.str4f=interp1(STR4F.dn,STR4F.Tsig(1,:),dn);

ii=find(~isnan(STR5B.subtidal.DNrho));
Alpha.str5b=interp1(STR5B.subtidal.DNrho(ii),STR5B.subtidal.alpha(ii),dn);
speed.str5b=interp1(STR5B.subtidal.DNrho(ii),STR5B.subtidal.Cw(ii),dn);
tempavg.str5b=interp1(STR5B.dn,nanmean(STR5B.Tsig),dn);
temp_s.str5b=interp1(STR5B.dn,STR5B.Tsig(end,:),dn);
temp_b.str5b=interp1(STR5B.dn,STR5B.Tsig(1,:),dn);

ii=find(~isnan(STR6B.subtidal.DNrho));
Alpha.str6b=interp1(STR6B.subtidal.DNrho(ii),STR6B.subtidal.alpha(ii),dn);
speed.str6b=interp1(STR6B.subtidal.DNrho(ii),STR6B.subtidal.Cw(ii),dn);
tempavg.str6b=interp1(STR6B.dn,nanmean(STR6B.Tsig),dn);
temp_s.str6b=interp1(STR6B.dn,STR6B.Tsig(end,:),dn);
temp_b.str6b=interp1(STR6B.dn,STR6B.Tsig(1,:),dn);


%% Plot Temp Map

    dn1=datenum(2017,10,19,0,0,0); 
    alpha=lon.*nan;
    Cw=lon.*nan;
    temps=lon.*nan;
    
    
    ii=findnearest_JACK(dn,dn1); 
    
    jj=find(ismember(names,'MS100'));
    alpha(jj)=Alpha.ms100(ii); Cw(jj)=speed.ms100(ii); temps(jj)=temp_s.ms100(ii); 
    tempb(jj)=temp_b.ms100(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'OC50'));
    alpha(jj)=Alpha.oc50(ii);  Cw(jj)=speed.oc50(ii);temps(jj)=temp_s.oc50(ii);
    tempb(jj)=temp_b.oc50(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC40N'));
    alpha(jj)=Alpha.oc40n(ii);   Cw(jj)=speed.oc40n(ii);temps(jj)=temp_s.oc40n(ii);   
    tempb(jj)=temp_b.oc40n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC40S'));
    alpha(jj)=Alpha.oc40s(ii);   Cw(jj)=speed.oc40s(ii);temps(jj)=temp_s.oc40s(ii);    
    tempb(jj)=temp_b.oc40n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC32N'));
    alpha(jj)=Alpha.oc32n(ii);   Cw(jj)=speed.oc32n(ii);temps(jj)=temp_s.oc32n(ii);
    tempb(jj)=temp_b.oc32n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC32S'));
    alpha(jj)=Alpha.oc32s(ii);   Cw(jj)=speed.oc32s(ii);temps(jj)=temp_s.oc32s(ii);
    tempb(jj)=temp_b.oc32s(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC25NA'));
    alpha(jj)=Alpha.oc25na(ii);  Cw(jj)=speed.oc25na(ii);temps(jj)=temp_s.oc25na(ii);
    tempb(jj)=temp_b.oc25na(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC25NB'));
    alpha(jj)=Alpha.oc25nb(ii);  Cw(jj)=speed.oc25nb(ii);temps(jj)=temp_s.oc25nb(ii);
    tempb(jj)=temp_b.oc25nb(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC25M'));
    alpha(jj)=Alpha.oc25m(ii);   Cw(jj)=speed.oc25m(ii);temps(jj)=temp_s.oc25m(ii);
    tempb(jj)=temp_b.oc25m(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC25SB'));
    alpha(jj)=Alpha.oc25sb(ii);  Cw(jj)=speed.oc25sb(ii);temps(jj)=temp_s.oc25sb(ii);
    tempb(jj)=temp_b.oc25sb(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC25SA'));
    alpha(jj)=Alpha.oc25sa(ii);  Cw(jj)=speed.oc25sa(ii);temps(jj)=temp_s.oc25sa(ii);
    tempb(jj)=temp_b.oc25sa(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC17N'));
    alpha(jj)=Alpha.oc17n(ii);   Cw(jj)=speed.oc17n(ii);temps(jj)=temp_s.oc17n(ii);
    tempb(jj)=temp_b.oc17n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC17S'));
    alpha(jj)=Alpha.oc17s(ii);   Cw(jj)=speed.oc17s(ii);temps(jj)=temp_s.oc17s(ii);
    tempb(jj)=temp_b.oc17s(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'OC10N'));
    alpha(jj)=Alpha.oc10n(ii);   Cw(jj)=speed.oc10n(ii);temps(jj)=temp_s.oc10n(ii);
    tempb(jj)=temp_b.oc10n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'STR3B'));
    alpha(jj)=Alpha.str3b(ii);   Cw(jj)=speed.str3b(ii);temps(jj)=temp_s.str3b(ii);
    tempb(jj)=temp_b.str3b(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'NRL50N'));
    alpha(jj)=Alpha.nrl50n(ii);  Cw(jj)=speed.nrl50n(ii);temps(jj)=temp_s.nrl50n(ii);
    tempb(jj)=temp_b.nrl50n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'NRL50S'));
    alpha(jj)=Alpha.nrl50s(ii);  Cw(jj)=speed.nrl50s(ii);temps(jj)=temp_s.nrl50s(ii);
    tempb(jj)=temp_b.nrl50s(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PS50'));
    alpha(jj)=Alpha.ps50(ii);    Cw(jj)=speed.ps50(ii);temps(jj)=temp_s.ps50(ii);
    tempb(jj)=temp_b.ps50(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'VB50N'));
    alpha(jj)=Alpha.vb50n(ii);    Cw(jj)=speed.vb50n(ii);temps(jj)=temp_s.vb50n(ii);
    tempb(jj)=temp_b.vb50n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'VB50S'));
    alpha(jj)=Alpha.vb50s(ii);    Cw(jj)=speed.vb50s(ii);temps(jj)=temp_s.vb50s(ii);
    tempb(jj)=temp_b.vb50s(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PS40N'));
    alpha(jj)=Alpha.ps40n(ii);    Cw(jj)=speed.ps40n(ii);temps(jj)=temp_s.ps40n(ii);
    tempb(jj)=temp_b.ps40n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PS40M'));
    alpha(jj)=Alpha.ps40m(ii);    Cw(jj)=speed.ps40m(ii);temps(jj)=temp_s.ps40m(ii);
    tempb(jj)=temp_b.ps40m(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PS40S'));
    alpha(jj)=Alpha.ps40s(ii);    Cw(jj)=speed.ps40s(ii);temps(jj)=temp_s.ps40s(ii);
    tempb(jj)=temp_b.ps40s(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'NRL35N'));
    alpha(jj)=Alpha.nrl35n(ii);   Cw(jj)=speed.nrl35n(ii);  temps(jj)=temp_s.nrl35n(ii);
    tempb(jj)=temp_b.nrl35n(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PS35M'));
    alpha(jj)=Alpha.ps35m(ii);    Cw(jj)=speed.ps35m(ii);temps(jj)=temp_s.ps35m(ii);
    tempb(jj)=temp_b.ps35m(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'NRL35S'));
    alpha(jj)=Alpha.nrl35s(ii);   Cw(jj)=speed.nrl35s(ii);temps(jj)=temp_s.nrl35s(ii);
    tempb(jj)=temp_b.nrl35s(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PS30N'));
    alpha(jj)=Alpha.ps30n(ii);   Cw(jj)=speed.ps30n(ii);temps(jj)=temp_s.ps30n(ii);
    tempb(jj)=temp_b.ps30n(ii); strat(jj)=temps(jj)-tempb(jj);

   
    jj=find(ismember(names,'PS30M'));
    alpha(jj)=Alpha.ps30m(ii);   Cw(jj)=speed.ps30m(ii);temps(jj)=temp_s.ps30m(ii);
    tempb(jj)=temp_b.ps30m(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'PS30S'));
    alpha(jj)=Alpha.ps30s(ii);   Cw(jj)=speed.ps30s(ii);temps(jj)=temp_s.ps30s(ii);
    tempb(jj)=temp_b.ps30s(ii); strat(jj)=temps(jj)-tempb(jj);
        
    jj=find(ismember(names,'VB30N'));
    alpha(jj)=Alpha.vb30n(ii);   Cw(jj)=speed.vb30n(ii);temps(jj)=temp_s.vb30n(ii);
    tempb(jj)=temp_b.vb30n(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'VB30S'));
    alpha(jj)=Alpha.vb30s(ii);   Cw(jj)=speed.vb30s(ii);temps(jj)=temp_s.vb30s(ii);
    tempb(jj)=temp_b.vb30s(ii); strat(jj)=temps(jj)-tempb(jj);
       
    jj=find(ismember(names,'VB25N'));
    alpha(jj)=Alpha.vb25n(ii);   Cw(jj)=speed.vb25n(ii);temps(jj)=temp_s.vb25n(ii);
    tempb(jj)=temp_b.vb25n(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'VB25S'));
    alpha(jj)=Alpha.vb25s(ii);   Cw(jj)=speed.vb25s(ii);temps(jj)=temp_s.vb25s(ii);
    tempb(jj)=temp_b.vb25s(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'NRL20N'));
    alpha(jj)=Alpha.nrl20n(ii);   Cw(jj)=speed.nrl20n(ii);temps(jj)=temp_s.nrl20n(ii);
    tempb(jj)=temp_b.nrl20n(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'NRL20S'));
    alpha(jj)=Alpha.nrl20s(ii);   Cw(jj)=speed.nrl20s(ii);temps(jj)=temp_s.nrl20s(ii);
    tempb(jj)=temp_b.nrl20s(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'NW2'));
    alpha(jj)=Alpha.nw2(ii);   Cw(jj)=speed.nw2(ii);temps(jj)=temp_s.nw2(ii);
    tempb(jj)=temp_b.nw2(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'RW3'));
    alpha(jj)=Alpha.rw3(ii);   Cw(jj)=speed.rw3(ii);temps(jj)=temp_s.rw3(ii);
    tempb(jj)=temp_b.rw3(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'BW1'));
    alpha(jj)=Alpha.bw1(ii);   Cw(jj)=speed.bw1(ii);temps(jj)=temp_s.bw1(ii);
    tempb(jj)=temp_b.bw1(ii); strat(jj)=temps(jj)-tempb(jj);

    jj=find(ismember(names,'PW5'));
    alpha(jj)=Alpha.pw5(ii);   Cw(jj)=speed.pw5(ii);temps(jj)=temp_s.pw5(ii);
    tempb(jj)=temp_b.pw5(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'STR4F'));
    alpha(jj)=Alpha.str4f(ii);   Cw(jj)=speed.str4f(ii);temps(jj)=temp_s.str4f(ii);
    tempb(jj)=temp_b.str4f(ii); strat(jj)=temps(jj)-tempb(jj);
    
    jj=find(ismember(names,'STR5B'));
    alpha(jj)=Alpha.str5b(ii);   Cw(jj)=speed.str5b(ii);temps(jj)=temp_s.str5b(ii);
    tempb(jj)=temp_b.str5b(ii); strat(jj)=temps(jj)-tempb(jj);
    
     jj=find(ismember(names,'STR6B'));
    alpha(jj)=Alpha.str6b(ii);   Cw(jj)=speed.str6b(ii); temps(jj)=temp_s.str6b(ii);
    tempb(jj)=temp_b.str6b(ii); strat(jj)=temps(jj)-tempb(jj);
    
    
    
    close all;
    load('/Volumes/InnerShelf1/JackAnalysis/AlongShelfVariability/AlphaColormap.mat');
    map=flipud(map);
%     map=flipud(cbrewer('div','PRGn',40));
%     jet_JACK
    
    figure('position',[88   239   948   564]);
    hax=tight_subplot(1,3,[0.01 0.001],[0.04 0.1],[0.03 0.02]);
    
    axes(hax(1))
    hold on
    patch(coast.patchx/1000,coast.patchy/1000,[1 1 1].*.8,'HandleVisibility','on');
    plot(coast.x/1000,coast.y/1000,'-','color',[1 1 1].*0.5,'linewidth',2);
    [c,h] = contour(bathy.x/1000,bathy.y/1000,-bathy.G.g200m.h,[-120 -100:10:10],...
        'color',[1 1 1].*0.8,'HandleVisibility','off','labelspacing',375);
    h1=scatter(x/1000,y/1000,2000,temps,'.');
%     set(h1,'markeredgecolor','k')        
    cc=colorbar; colormap(flipud(map)); caxis([10 20]);
    set(cc,'position',[0.272151898734177         0.528724420329614        0.0245522933678615         0.285105366904428]);
%     cl=xlabel(cc,'Surface Temp');
%     set(cl,'position',[-0.782666714986165      -0.00142130072346799                         0],...
%         'fontweight','bold','fontsize',18);
    yaxis(-20,15);xaxis(-15,10);
    aspect_jack(gca,lat(1));
    title(['Surface Temp' char(10) datestr(dn1,'mmm dd HH:MM')]);
    set(gca,'fontsize',14,'fontweight','bold','box','on','layer','top');
    
    axes(hax(2))
    hold on
    patch(coast.patchx/1000,coast.patchy/1000,[1 1 1].*.8,'HandleVisibility','on');
    plot(coast.x/1000,coast.y/1000,'-','color',[1 1 1].*0.5,'linewidth',2);
    [c,h] = contour(bathy.x/1000,bathy.y/1000,-bathy.G.g200m.h,[-120 -100:10:10],...
        'color',[1 1 1].*0.8,'HandleVisibility','off','labelspacing',375);
    h1=scatter(x/1000,y/1000,2000,tempb,'.');
%     set(h1,'markeredgecolor','k')        
    cc=colorbar; colormap(flipud(map)); caxis([10 20]);
    set(cc,'position',[0.59     0.528724420329614        0.0245522933678615         0.285105366904428]);
%     cl=xlabel(cc,'Surface Temp');
%     set(cl,'position',[-0.782666714986165      -0.00142130072346799                         0],...
%         'fontweight','bold','fontsize',18);
    yaxis(-20,15);xaxis(-15,10);
    aspect_jack(gca,lat(1));
    title(['Bottom Temp' char(10) datestr(dn1,'mmm dd HH:MM')]);
    set(gca,'fontsize',14,'fontweight','bold','box','on','layer','top');
    
    axes(hax(3))
    hold on
    patch(coast.patchx/1000,coast.patchy/1000,[1 1 1].*.8,'HandleVisibility','on');
    plot(coast.x/1000,coast.y/1000,'-','color',[1 1 1].*0.5,'linewidth',2);
    [c,h] = contour(bathy.x/1000,bathy.y/1000,-bathy.G.g200m.h,[-120 -100:10:10],...
        'color',[1 1 1].*0.8,'HandleVisibility','off','labelspacing',375);
    h1=scatter(x/1000,y/1000,2000,strat,'.');
    cc=colorbar; colormap(flipud(map)); caxis([0 5]);
    set(cc,'position',[0.91     0.528724420329614        0.0245522933678615         0.285105366904428]);
    yaxis(-20,15);xaxis(-15,10);
    aspect_jack(gca,lat(1));
    title(['\DeltaT' char(10) datestr(dn1,'mmm dd HH:MM')]);
    set(gca,'fontsize',14,'fontweight','bold','box','on','layer','top');
    
    
    
    if 0
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r300',[anadir  datestr(dn1,'mmmdd_HHMM')]);
    end
    
%% 

dn1=datenum(2017,10,15,0,0,0); 
dn2=datenum(2017,10,19,0,0,0); 

dt=0.25;
labels=datestr(dn1:dt:dn2,'mmmdd HH:MM');
labels(2:2:end,:) = nan; % remove every other one

close all 
figure('position',[ 133         185        1118         600]); 
subplot(3,1,1);
plot(dn,temp_s.oc50,'-k','linewidth',2);
hold on 
% plot(dn,pl33tn(temp_s.oc50,1/60,16),'-k','linewidth',3);
plot(dn,temp_s.nrl50n,'-r','linewidth',2);
plot(dn,temp_s.ps50,'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,temp_s.vb50n,'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,temp_s.vb50s,'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,temp_s.nrl50s,'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels);
title('surface temp')

subplot(3,1,2);
hold on 
plot(dn,temp_b.oc50,'-k','linewidth',2);
plot(dn,temp_b.nrl50n,'-r','linewidth',2);
plot(dn,temp_b.ps50,'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,temp_b.vb50n,'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,temp_b.vb50s,'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,temp_b.nrl50s,'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels);
title('bottom temp')


subplot(3,1,3);
plot(dn,temp_s.oc50-temp_b.oc50,'-k','linewidth',2);
hold on 
plot(dn,temp_s.nrl50n-temp_b.nrl50n,'-r','linewidth',2);
plot(dn,temp_s.ps50-temp_b.ps50,'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,temp_s.vb50n-temp_b.vb50n,'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,temp_s.vb50s-temp_b.vb50s,'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,temp_s.nrl50s-temp_b.nrl50s,'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels);
title('\DeltaT')

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir  'Case2']);
end


%%
dn1=datenum(2017,10,15,0,0,0); 
dn2=datenum(2017,10,19,0,0,0); 

dt=0.25;
labels=datestr(dn1:dt:dn2,'mmmdd HH:MM');
labels(2:2:end,:) = nan; % remove every other one

filt=33;

close all 
figure('position',[ 133         185        1118         600]); 
subplot(3,1,1);
plot(dn,pl33tn(temp_s.oc50,1/60,filt),'-k','linewidth',2);
hold on 
plot(dn,pl33tn(temp_s.nrl50n,1/60,filt),'-r','linewidth',2);
plot(dn,pl33tn(temp_s.ps50,1/60,filt),'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,pl33tn(temp_s.vb50n,1/60,filt),'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,pl33tn(temp_s.vb50s,1/60,filt),'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,pl33tn(temp_s.nrl50s,1/60,filt),'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels);
title('surface temp')
grid on; 

subplot(3,1,2);
hold on 
plot(dn,pl33tn(temp_b.oc50,1/60,filt),'-k','linewidth',2);
plot(dn,pl33tn(temp_b.nrl50n,1/60,filt),'-r','linewidth',2);
plot(dn,pl33tn(temp_b.ps50,1/60,filt),'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,pl33tn(temp_b.vb50n,1/60,filt),'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,pl33tn(temp_b.vb50s,1/60,filt),'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,pl33tn(temp_b.nrl50s,1/60,filt),'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels);
title('bottom temp')
grid on; 


subplot(3,1,3);
plot(dn,pl33tn((temp_s.oc50-temp_b.oc50),1/60,filt),'-k','linewidth',2);
hold on 
plot(dn,pl33tn((temp_s.nrl50n-temp_b.nrl50n),1/60,filt),'-r','linewidth',2);
plot(dn,pl33tn((temp_s.ps50-temp_b.ps50),1/60,filt),'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,pl33tn((temp_s.vb50n-temp_b.vb50n),1/60,filt),'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,pl33tn((temp_s.vb50s-temp_b.vb50s),1/60,filt),'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,pl33tn((temp_s.nrl50s-temp_b.nrl50s),1/60,filt),'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'xtick',dn1:dt:dn2,'xticklabels',labels);
title('\DeltaT')
grid on; 

if 0
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300',[anadir  'Case2_filt' num2str(filt)]);
end


%% Stratification during a feature 
period= 4; % 1, 2, 3, 4 (4 is the grey period - where bores are trackable

switch period
    case 1
        dn1=datenum(2017,9,27);
        dn2=datenum(2017,10,4);
    case 2
        dn1=datenum(2017,10,14);
        dn2=datenum(2017,10,20);    
    case 3
        dn1=datenum(2017,10,25);
        dn2=datenum(2017,10,31); 
    case 4
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);     
end
dt=1;

map=flipud(cbrewer('div','RdYlBu',100));


switch period
    case 1
        cmin=1024;
        cmax=1025.75;
        isobold=1025.25;
    case 2
        cmin=1024.5;
        cmax=1025.75;
        isobold=1025.25;
    case 3
        cmin=1024.5;
        cmax=1025.75;
        isobold=1025.25;
    case 4
        cmin=1023.75;
        cmax=1025.75;
        isobold=1025.25;
end

close all;
figure('position',[185    51   684   754]);
hax=tight_subplot(7,1,[0.02 0.1],[0.04 0.02],[0.1 0.12]);


axes(hax(1))
% % % plot(dn_atm,wspd,'-k','linewidth',2); hold on 
% % % plot(dn_atm,wdir,'.k','linewidth',2); hold on 
% % ylim=20;
% % hold on   
% % patch('XData',[datenum(2017,9,11,9,0,0) datenum(2017,9,12,12,0,0) datenum(2017,9,12,12,0,0) datenum(2017,9,11,9,0,0)], ...
% %     'YData',[-ylim; -ylim; ylim; ylim], ...
% %     'FaceColor',[0.94 0.74 0.94],'edgecolor',[0.94 0.94 0.94]);
% % % patch('XData',[datenum(2017,9,27,0,0,0) datenum(2017,9,28,0,0,0) datenum(2017,9,28,0,0,0) datenum(2017,9,27,0,0,0)], ...
% % %     'YData',[-ylim; -ylim; ylim; ylim], ...
% % %     'FaceColor',[0.74 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
% % 
% % % patch('XData',[datenum(2017,9,29,0,0,0) datenum(2017,9,30,0,0,0) datenum(2017,9,30,0,0,0) datenum(2017,9,29,0,0,0)], ...
% % %     'YData',[-ylim; -ylim; ylim; ylim], ...
% % %     'FaceColor',[0.74 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
% % patch('XData',[datenum(2017,10,1,0,0,0) datenum(2017,10,2,0,0,0) datenum(2017,10,2,0,0,0) datenum(2017,10,1,0,0,0)], ...
% %     'YData',[-ylim; -ylim; ylim; ylim], ...
% %     'FaceColor',[0.74 0.94 0.94],'edgecolor',[0.94 0.94 0.94]);
% % quiver(dn_atm2,dn_atm2.*0,weast,wnorth,'k'); %,'autoscale','off');%,'maxheadsize',0)
% % % quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,24),pl33tn(wind_n,1,24),'k')
% % % quiver(dn_atm,dn_atm.*0,pl33tn(wind_e,1,33),pl33tn(wind_n,1,33),'k','autoscale','off','maxheadsize',0)
% % xaxis(dn1,dn2);
% % yaxis(-2,2);
% % ylabel(['Wind speed' char(10) '(m/s)']);
% % set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
% %     'box','on','layer','top','ytick',-2:1:2,'yticklabels',[-2:1:2]*7.5,...
% %     'TickDir', 'both','TickLength',[0.0075 0.0075]);
% % grid on; hh=gca; hh.XMinorGrid='on'; 


plot(dn,pl33tn((temp_s.oc50-temp_b.oc50),1/60,filt),'-k','linewidth',2);
hold on 
plot(dn,pl33tn((temp_s.nrl50n-temp_b.nrl50n),1/60,filt),'-r','linewidth',2);
plot(dn,pl33tn((temp_s.ps50-temp_b.ps50),1/60,filt),'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,pl33tn((temp_s.vb50n-temp_b.vb50n),1/60,filt),'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,pl33tn((temp_s.vb50s-temp_b.vb50s),1/60,filt),'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,pl33tn((temp_s.nrl50s-temp_b.nrl50s),1/60,filt),'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top',...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);title('\DeltaT')
grid on; 

axes(hax(2)); hold on;
set(gca,'position',[  0.1         0.712451686244789                      0.78         0.117142857142857]);
pcolorjw(OC50.dn_sig,OC50.H,OC50.Rho_sig);
colormap(flipud(map)); 
cc=colorbar;
caxis([cmin cmax]);
contour(OC50.DNrho(2:108),OC50.Z,OC50.Rho(:,2:108),[1021:0.25:1030],'k');
contour(OC50.DNrho(118:end),OC50.Z,OC50.Rho(:,118:end),[1021:0.25:1030],'k');
contour(OC50.DNrho(2:108),OC50.Z,OC50.Rho(:,2:108),[isobold isobold],'k','linewidth',1.5);
contour(OC50.DNrho(118:end),OC50.Z,OC50.Rho(:,118:end),[isobold isobold],'k','linewidth',1.5);
set(cc,'position',[0.892768292922024        0.0411140583554377        0.0355942801773912         0.789124668435013]);
xaxis(dn1,dn2); yaxis(0,50);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('OC50');
title('Density (kg/m^3)');

axes(hax(3));hold on;
set(gca,'position',[  0.1         0.579287608942781                      0.78         0.117142857142857]);
pcolorjw(NRL50N.dn_sig,NRL50N.H,NRL50N.Rho_sig);
colormap(flipud(map));
caxis([cmin cmax]);
contour(NRL50N.subtidal.DNrho,NRL50N.subtidal.Z,NRL50N.subtidal.Rho,[1021:0.25:1030],'k');
contour(NRL50N.subtidal.DNrho,NRL50N.subtidal.Z,NRL50N.subtidal.Rho,[isobold isobold],'k','linewidth',1.5);
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('NRL50N');

axes(hax(4));hold on; 
set(gca,'position',[  0.1         0.443471011746874                      0.78         0.117142857142857]);
pcolorjw(PS50.dn_sig,PS50.H,PS50.Rho_sig);
colormap(flipud(map));
caxis([cmin cmax]);
contour(PS50.subtidal.DNrho,PS50.subtidal.Z,PS50.subtidal.Rho,[1021:0.25:1030],'k');
contour(PS50.subtidal.DNrho,PS50.subtidal.Z,PS50.subtidal.Rho,[isobold isobold],'k','linewidth',1.5);
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('PS50');

axes(hax(5));hold on; 
set(gca,'position',[  0.1         0.308980674497916                      0.78         0.117142857142857]);
pcolorjw(VB50N.dn_sig,VB50N.H,VB50N.Rho_sig);
colormap(flipud(map));
caxis([cmin cmax]);
contour(VB50N.subtidal.DNrho,VB50N.subtidal.Z,VB50N.subtidal.Rho,[1021:0.25:1030],'k');
contour(VB50N.subtidal.DNrho,VB50N.subtidal.Z,VB50N.subtidal.Rho,[isobold isobold],'k','linewidth',1.5);
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('VB50N');

axes(hax(6));hold on; 
set(gca,'position',[  0.1         0.174490337248958                      0.78         0.117142857142857]);
pcolorjw(VB50S.dn_sig,VB50S.H,VB50S.Rho_sig);
colormap(flipud(map));
caxis([cmin cmax]);
contour(VB50S.subtidal.DNrho,VB50S.subtidal.Z,VB50S.subtidal.Rho,[1021:0.25:1030],'k');
contour(VB50S.subtidal.DNrho,VB50S.subtidal.Z,VB50S.subtidal.Rho,[isobold isobold],'k','linewidth',1.5);
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('VB50S');


axes(hax(7));hold on; 
set(gca,'position',[  0.1         0.0400000000000001                      0.78         0.117142857142857]);
pcolorjw(NRL50S.dn_sig,NRL50S.H,NRL50S.Rho_sig);
colormap(flipud(map));
caxis([cmin cmax]);
contour(NRL50S.subtidal.DNrho,NRL50S.subtidal.Z,NRL50S.subtidal.Rho,[1021:0.25:1030],'k');
contour(NRL50S.subtidal.DNrho,NRL50S.subtidal.Z,NRL50S.subtidal.Rho,[isobold isobold],'k','linewidth',1.5);
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('NRL50S');

if 0
    set(gcf,'PaperPositionMode','auto');
    switch period 
        case 1
            print(gcf,'-dpng','-r300',[anadir 'Case1_strat']);
        case 2
            print(gcf,'-dpng','-r300',[anadir 'Case2_strat']);
        case 3
            print(gcf,'-dpng','-r300',[anadir 'Case3_strat']);
        case 4
            print(gcf,'-dpng','-r300',[anadir 'GreyPeriod_strat']);
    end
end 

%% Stratification during a feature 
period= 3; % 1, 2, 3, 4 (4 is the grey period - where bores are trackable

switch period
    case 1
        dn1=datenum(2017,9,27);
        dn2=datenum(2017,10,4);
    case 2
        dn1=datenum(2017,10,14);
        dn2=datenum(2017,10,20);    
    case 3
        dn1=datenum(2017,10,25);
        dn2=datenum(2017,10,31); 
    case 4
        dn1=datenum(2017,9,11);
        dn2=datenum(2017,9,22);     
end
dt=1;

map=(cbrewer('div','RdBu',100));

cmin=-0.3;
cmax=0.3;
        
close all;
figure('position',[185    51   684   754]);
hax=tight_subplot(7,1,[0.02 0.1],[0.04 0.02],[0.1 0.12]);

filt =33; 

axes(hax(1))
plot(dn,pl33tn((temp_s.oc50-temp_b.oc50),1/60,filt),'-k','linewidth',2);
hold on 
plot(dn,pl33tn((temp_s.nrl50n-temp_b.nrl50n),1/60,filt),'-r','linewidth',2);
plot(dn,pl33tn((temp_s.ps50-temp_b.ps50),1/60,filt),'-','color',hex2rgb('#1682E0'),'linewidth',2);
plot(dn,pl33tn((temp_s.vb50n-temp_b.vb50n),1/60,filt),'-','color',hex2rgb('#45C162'),'linewidth',2);
plot(dn,pl33tn((temp_s.vb50s-temp_b.vb50s),1/60,filt),'-','color',[0.7 0.7 0.7],'linewidth',2);
plot(dn,pl33tn((temp_s.nrl50s-temp_b.nrl50s),1/60,filt),'-','color',hex2rgb('#FF9456'),'linewidth',2);
xaxis(dn1,dn2);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top',...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);title('\DeltaT')
grid on; 

axes(hax(2)); hold on;
set(gca,'position',[  0.1         0.712451686244789                      0.78         0.117142857142857]);
pcolorjw(OC50.dn_sig,OC50.H,OC50.east_sd);
colormap(flipud(map)); 
cc=colorbar;
caxis([cmin cmax]);
contour(OC50.dn_sig,OC50.H,OC50.Rho_sig,[1021:0.25:1030],'k');
set(cc,'position',[0.892768292922024        0.0411140583554377        0.0355942801773912         0.789124668435013]);
xaxis(dn1,dn2); yaxis(0,50);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('OC50');
title('Density (kg/m^3)');

axes(hax(3));hold on;
set(gca,'position',[  0.1         0.579287608942781                      0.78         0.117142857142857]);
pcolorjw(NRL50N.dn_sig,NRL50N.H,NRL50N.east_sd);
colormap(flipud(map));
caxis([cmin cmax]);
contour(NRL50N.dn_sig,NRL50N.H,NRL50N.Rho_sig,[1021:0.25:1030],'k');
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('NRL50N');

axes(hax(4));hold on; 
set(gca,'position',[  0.1         0.443471011746874                      0.78         0.117142857142857]);
pcolorjw(PS50.dn_sig,PS50.H,PS50.east_sd);
colormap(flipud(map));
caxis([cmin cmax]);
contour(PS50.dn_sig,PS50.H,PS50.Rho_sig,[1021:0.25:1030],'k');
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('PS50');

axes(hax(5));hold on; 
set(gca,'position',[  0.1         0.308980674497916                      0.78         0.117142857142857]);
pcolorjw(VB50N.dn_sig,VB50N.H,VB50N.east_sd);
colormap(flipud(map));
caxis([cmin cmax]);
contour(VB50N.dn_sig,VB50N.H,VB50N.Rho_sig,[1021:0.25:1030],'k');
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('VB50N');

axes(hax(6));hold on; 
set(gca,'position',[  0.1         0.174490337248958                      0.78         0.117142857142857]);
pcolorjw(VB50S.dn_sig,VB50S.H,VB50S.east_sd);
colormap(flipud(map));
caxis([cmin cmax]);
contour(VB50S.dn_sig,VB50S.H,VB50S.Rho_sig,[1021:0.25:1030],'k');
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',[],...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('VB50S');


axes(hax(7));hold on; 
set(gca,'position',[  0.1         0.0400000000000001                      0.78         0.117142857142857]);
pcolorjw(NRL50S.dn_sig,NRL50S.H,NRL50S.east_sd);
colormap(flipud(map));
caxis([cmin cmax]);
contour(NRL50S.dn_sig,NRL50S.H,NRL50S.Rho_sig,[1021:0.25:1030],'k');
xaxis(dn1,dn2); yaxis(0,51);
set(gca,'fontweight','bold','fontsize',12,'xtick',dn1:dt:dn2,'xticklabels',datestr(dn1:dt:dn2,'mm/dd'),...
    'box','on','layer','top','ytick',10:10:40,'yticklabels',10:10:40,...
    'TickDir', 'both','TickLength',[0.0075 0.0075]);
grid on; hh=gca; hh.XMinorGrid='on'; 
ylabel('NRL50S');

if 1
    set(gcf,'PaperPositionMode','auto');
    switch period 
        case 1
            print(gcf,'-dpng','-r300',[anadir 'Case1_bores']);
        case 2
            print(gcf,'-dpng','-r300',[anadir 'Case2_bores']);
        case 3
            print(gcf,'-dpng','-r300',[anadir 'Case3_bores']);
        case 4
            print(gcf,'-dpng','-r300',[anadir 'GreyPeriod_bores']);
    end
end 
