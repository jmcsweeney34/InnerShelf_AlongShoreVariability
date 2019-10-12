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
