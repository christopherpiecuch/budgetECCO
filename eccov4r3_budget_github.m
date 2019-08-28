



%%%%%
% hello
clear all, close all, clc
tic % keep track of time

%%%%%
% cd to base directory
% (enter your base directory)
baseDir='/Users/christopherpiecuch/Desktop/budgetECCO/';
cd(baseDir)

%%%%%
% instantiate gcmfaces
% 1. cd to gcmfaces directory
% the gcmfaces package has been provided with the budgetECCO directory
% gcmfaces codes were produced by gael forget (mit) and downloaded 
% from https://github.com/gaelforget/gcmfaces
% make sure to unzip ``gcmfaces.zip'' directory included with this code
cd([baseDir,'gcmfaces/'])

% 2. define where grid files are
dirGrid=[baseDir,'gcmfaces/grid/'];

% 3. call global and load grid
% these lines below will probably generate the warning message
% ``Warning: mygrid has not yet been loaded to memory 
% > In gcmfaces_global (line 115)''
% but that's okay; don't worry
global mygrid; mygrid=[];
gcmfaces_global
nFaces=5;
grid_load(dirGrid,nFaces,'compact');

%%%%%
% define some useful space and time quantities
dzMatF = mk3D(mygrid.DRF, mygrid.hFacC);
dzMat = dzMatF.*mygrid.hFacC;
RACMat = mk3D(mygrid.RAC, mygrid.hFacC);
VVV = mygrid.mskC.*mygrid.hFacC.*mk3D(mygrid.RAC,mygrid.mskC).*mk3D(mygrid.DRF,mygrid.mskC);
nLevels = numel(mygrid.RC);
secPerHour = 3600;
hFacC = mygrid.hFacC;
dt=diff([1 (datenum(1992,2:288,1)-datenum(1992,1,1,12,0,0))*24 (datenum(1992,288,31,12,0,0)-datenum(1992,1,1,12,0,0))*24]);
% ecco version 4 release 3 has 288 time monthly time steps
% from january 1992 to december 2015
% here i choose 4 arbitrary months of model output 
% (january to april 2010) just to demonstrate the code
% and speed things up, alleviating longer calculation
% times and memory issues due to loading too much 
% model output into the workspace, e.g., if you were 
% to try to compute the budgets for the full 1992-2015 period
t0=217; % jan 2010
tf=220; % dec 2010
nTS=numel(t0:tf);

%%%%%
% load ocean-bottom geothermal flux
% this file can be downloaded from https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3/input_init/
fid=fopen([baseDir,'geothermalFlux.bin'],'r','b');
geoflx2d=convert2gcmfaces(reshape(fread(fid,'float32'),90,1170));
fclose(fid);
% take 2d geothermal flux field and apply at bottom of 3d ocean
mskc=mygrid.mskC;
mskc(isnan(mskc))=0;
mskcp1=mskc;
mskcp1(:,:,51)=0;
mskcp1(:,:,1)=[];
mskb=mskc-mskcp1;
geoflx3d=mk3D(geoflx2d,mskc).*mskb.*mygrid.mskC;
clear mskc mskcp1 mskb geoflx2d

%%%%%
% start loading fields

%%%%%
% define directories where outputs are
% this code assumes you have filled these two 
% subdirectories with the relevant "nctiles" model
% output files, as described in README and downloadable
% from https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3/nctiles_monthly/
% and https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3/nctiles_monthly_snapshots/
dirIn1=[baseDir,'nctiles_monthly/'];
dirIn2=[baseDir,'nctiles_monthly_snapshots/'];

%%%%%
% load two-dimensional monthly averaged fields
ETAN=read_nctiles([dirIn1,'ETAN/ETAN'],'ETAN',t0:tf);
SFLUX=read_nctiles([dirIn1,'SFLUX/SFLUX'],'SFLUX',t0:tf);
TFLUX=read_nctiles([dirIn1,'TFLUX/TFLUX'],'TFLUX',t0:tf);
oceQsw=read_nctiles([dirIn1,'oceQsw/oceQsw'],'oceQsw',t0:tf);
oceFWflx=read_nctiles([dirIn1,'oceFWflx/oceFWflx'],'oceFWflx',t0:tf);
oceSPflx=read_nctiles([dirIn1,'oceSPflx/oceSPflx'],'oceSPflx',t0:tf);

%%%%%
% load three-dimensional monthly averaged fields
THETA=read_nctiles([dirIn1,'THETA/THETA'],'THETA',t0:tf);
SALT=read_nctiles([dirIn1,'SALT/SALT'],'SALT',t0:tf);
ADVr_TH=read_nctiles([dirIn1,'ADVr_TH/ADVr_TH'],'ADVr_TH',t0:tf);
ADVx_TH=read_nctiles([dirIn1,'ADVx_TH/ADVx_TH'],'ADVx_TH',t0:tf);
ADVy_TH=read_nctiles([dirIn1,'ADVy_TH/ADVy_TH'],'ADVy_TH',t0:tf);
DFrI_TH=read_nctiles([dirIn1,'DFrI_TH/DFrI_TH'],'DFrI_TH',t0:tf);
DFrE_TH=read_nctiles([dirIn1,'DFrE_TH/DFrE_TH'],'DFrE_TH',t0:tf);
DFxE_TH=read_nctiles([dirIn1,'DFxE_TH/DFxE_TH'],'DFxE_TH',t0:tf);
DFyE_TH=read_nctiles([dirIn1,'DFyE_TH/DFyE_TH'],'DFyE_TH',t0:tf);
ADVr_SLT=read_nctiles([dirIn1,'ADVr_SLT/ADVr_SLT'],'ADVr_SLT',t0:tf);
ADVx_SLT=read_nctiles([dirIn1,'ADVx_SLT/ADVx_SLT'],'ADVx_SLT',t0:tf);
ADVy_SLT=read_nctiles([dirIn1,'ADVy_SLT/ADVy_SLT'],'ADVy_SLT',t0:tf);
DFrI_SLT=read_nctiles([dirIn1,'DFrI_SLT/DFrI_SLT'],'DFrI_SLT',t0:tf);
DFrE_SLT=read_nctiles([dirIn1,'DFrE_SLT/DFrE_SLT'],'DFrE_SLT',t0:tf);
DFxE_SLT=read_nctiles([dirIn1,'DFxE_SLT/DFxE_SLT'],'DFxE_SLT',t0:tf);
DFyE_SLT=read_nctiles([dirIn1,'DFyE_SLT/DFyE_SLT'],'DFyE_SLT',t0:tf);
oceSPtnd=read_nctiles([dirIn1,'oceSPtnd/oceSPtnd'],'oceSPtnd',t0:tf);
UVELMASS=read_nctiles([dirIn1,'UVELMASS/UVELMASS'],'UVELMASS',t0:tf);
VVELMASS=read_nctiles([dirIn1,'VVELMASS/VVELMASS'],'VVELMASS',t0:tf);
WVELMASS=read_nctiles([dirIn1,'WVELMASS/WVELMASS'],'WVELMASS',t0:tf);

%%%%%
% load two- and three-dimensional monthly snapshots
THETA_SNAP=read_nctiles([dirIn2,'THETA/THETA'],'THETA',(t0-1):tf);
SALT_SNAP=read_nctiles([dirIn2,'SALT/SALT'],'SALT',(t0-1):tf);
ETAN_SNAP=read_nctiles([dirIn2,'ETAN/ETAN'],'ETAN',(t0-1):tf+1);

%%%%%
% START COMPUTATIONS
%%%%%
% for details see Piecuch (2017) PDF document 
% ``A Note on Practical Evaluation of Budgets in 
% ECCO Version 4 Release 3'' provided with code package

rhoconst=1029; % constant reference ocean density kg/m^3
heatcap=3994; % specific heat capacity of seawater J/kg/K
rcp=rhoconst*heatcap;

%%%%%
% MASS
%%%%%
tendM=0*THETA;
hConvM=0*THETA;
vConvM=0*THETA;
forcM=0*THETA;
for nt=1:nTS, disp(num2str(nt))
% total tendency
 tendM(:,:,:,nt)=(1./mk3D(mygrid.Depth,mygrid.mskC)).*mk3D((ETAN_SNAP(:,:,nt+1)-ETAN_SNAP(:,:,nt))/(secPerHour*dt(nt)),mygrid.mskC);
% horizontal convergence
 hConvM(:,:,:,nt)=mygrid.mskC.*calc_UV_conv(UVELMASS(:,:,:,nt),VVELMASS(:,:,:,nt),{'dh'})./(RACMat.*hFacC);
% vertical divergence
 for nz=1:nLevels, disp(num2str(nz))
  nzp1=min([nz+1,nLevels]);
  vConvM(:,:,nz,nt)=squeeze(WVELMASS(:,:,nzp1,nt)*double(nz~=nLevels)-WVELMASS(:,:,nz,nt)*double(nz~=1))./(dzMat(:,:,nz));
 end
% forcing
 forcM(:,:,:,nt)=mygrid.mskC.*mk3D(oceFWflx(:,:,nt),mygrid.mskC)./(dzMat*rhoconst); forcM(:,:,2:50,nt)=0*mygrid.mskC(:,:,2:50);
end
TOT_V=(tendM);
FRC_V=(forcM);
CON_V=(hConvM+vConvM);
clear tendM forcM *ConvM
%%%%%

%%%%%
% HEAT  
%%%%% 
% total tendency
HC_snap1=0*THETA;
HC_snap2=0*THETA;
tendHC=0*THETA;
for nt=1:nTS
 HC_snap1(:,:,:,nt)=(THETA_SNAP(:,:,:,nt).*(1+mk3D(ETAN_SNAP(:,:,nt)./mygrid.Depth,dzMat)));
 HC_snap2(:,:,:,nt)=(THETA_SNAP(:,:,:,nt+1).*(1+mk3D(ETAN_SNAP(:,:,nt+1)./mygrid.Depth,dzMat)));
 tendHC(:,:,:,nt)=(HC_snap2(:,:,:,nt)-HC_snap1(:,:,:,nt))/(secPerHour*(dt(nt)));
end

% horizontal divergences
adv_hConvHT=0*THETA;
dif_hConvHT=0*THETA;
for nt=1:nTS
 adv_hConvHT(:,:,:,nt)=calc_UV_conv(ADVx_TH(:,:,:,nt),ADVy_TH(:,:,:,nt))./VVV;
 dif_hConvHT(:,:,:,nt)=calc_UV_conv(DFxE_TH(:,:,:,nt),DFyE_TH(:,:,:,nt))./VVV;
end

% vertical divergences
adv_vConvHT=0*ADVx_TH;
dif_vConvHT=0*ADVx_TH;
for nt=1:nTS, disp(num2str(nt))
for nz=1:nLevels, disp(['... ',num2str(nz)])
 nzp1=min([nz+1,nLevels]);
 adv_vConvHT(:,:,nz,nt)=squeeze(ADVr_TH(:,:,nzp1,nt)*double(nz~=nLevels)-ADVr_TH(:,:,nz,nt));
 dif_vConvHT(:,:,nz,nt)=squeeze(DFrI_TH(:,:,nzp1,nt)*double(nz~=nLevels)-DFrI_TH(:,:,nz,nt)+DFrE_TH(:,:,nzp1,nt)*double(nz~=nLevels)-DFrE_TH(:,:,nz,nt));
end
adv_vConvHT(:,:,:,nt)=adv_vConvHT(:,:,:,nt)./VVV;
dif_vConvHT(:,:,:,nt)=dif_vConvHT(:,:,:,nt)./VVV;

end
% surface heat flux
R=0.62;
zeta1=0.6;
zeta2=20;
q1=R*exp(1/zeta1*mygrid.RF(1:nLevels))+(1-R)*exp(1/zeta2*mygrid.RF(1:nLevels));
q2=R*exp(1/zeta1*mygrid.RF(2:(nLevels+1)))+(1-R)*exp(1/zeta2*mygrid.RF(2:(nLevels+1)));
% correction for the 200m cutoff
zCut = find(mygrid.RC<-200,1);
q1(zCut:nLevels)=0;
q2((zCut-1):nLevels)=0;
HF=0*ADVx_TH;
msk=mygrid.mskC; msk(isnan(msk))=0;
for nt=1:nTS, disp(num2str(nt))
for nz=1:nLevels, disp(['... ',num2str(nz)])
 if nz==1
  HF(:,:,nz,nt)=TFLUX(:,:,nt)-(1-(q1(nz)-q2(nz)))*oceQsw(:,:,nt);
 else
  nzp1=min([nz+1,nLevels]);
  HF(:,:,nz,nt)=HF(:,:,nz,nt)+( (mygrid.mskC(:,:,nz)==1).*q1(nz)-(mygrid.mskC(:,:,nzp1)==1).*q2(nz)).*oceQsw(:,:,nt);
 end
end
% add geothermal
HF(:,:,:,nt) = HF(:,:,:,nt) + geoflx3d;
HF(:,:,:,nt)=HF(:,:,:,nt)./(rcp*dzMat);
end

TOT_T=(tendHC);
ADV_T=(adv_hConvHT+adv_vConvHT);
DIF_T=(dif_hConvHT+dif_vConvHT);
FRC_T=(HF);
clear tendHC adv*T dif*T HF

%%%%%
% SALT (NOT salinity ... yet)
%%%%% 
% total tendency
SL_snap1=0*SALT;
SL_snap2=0*SALT;
tendSL=0*SALT;
for nt=1:nTS
 SL_snap1(:,:,:,nt)=(SALT_SNAP(:,:,:,nt).*(1+mk3D(ETAN_SNAP(:,:,nt)./mygrid.Depth,dzMat)));
 SL_snap2(:,:,:,nt)=(SALT_SNAP(:,:,:,nt+1).*(1+mk3D(ETAN_SNAP(:,:,nt+1)./mygrid.Depth,dzMat)));
 tendSL(:,:,:,nt)=(SL_snap2(:,:,:,nt)-SL_snap1(:,:,:,nt))/(secPerHour*(dt(nt)));
end

% horizontal divergences
adv_hConvSL=0*SALT;
dif_hConvSL=0*SALT;
for nt=1:nTS
 adv_hConvSL(:,:,:,nt)=calc_UV_conv(ADVx_SLT(:,:,:,nt),ADVy_SLT(:,:,:,nt))./VVV;
 dif_hConvSL(:,:,:,nt)=calc_UV_conv(DFxE_SLT(:,:,:,nt),DFyE_SLT(:,:,:,nt))./VVV;
end

% vertical divergences
adv_vConvSL=0*ADVx_SLT;
dif_vConvSL=0*ADVx_SLT;
for nt=1:nTS, disp(num2str(nt))
for nz=1:nLevels, disp(['... ',num2str(nz)])
 nzp1=min([nz+1,nLevels]);
 adv_vConvSL(:,:,nz,nt)=squeeze(ADVr_SLT(:,:,nzp1,nt)*double(nz~=nLevels)-ADVr_SLT(:,:,nz,nt));
 dif_vConvSL(:,:,nz,nt)=squeeze(DFrI_SLT(:,:,nzp1,nt)*double(nz~=nLevels)-DFrI_SLT(:,:,nz,nt)+DFrE_SLT(:,:,nzp1,nt)*double(nz~=nLevels)-DFrE_SLT(:,:,nz,nt));
end
adv_vConvSL(:,:,:,nt)=adv_vConvSL(:,:,:,nt)./VVV;
dif_vConvSL(:,:,:,nt)=dif_vConvSL(:,:,:,nt)./VVV;

end

% surface salt flux
SF=0*oceSPtnd;
for nt=1:nTS
for nz=1:nLevels
 if nz==1
  SF(:,:,nz,nt)=SFLUX(:,:,nt)/rhoconst;
 end
 SF(:,:,nz,nt)=SF(:,:,nz,nt)+oceSPtnd(:,:,nz,nt)/rhoconst;
end
SF(:,:,:,nt)=SF(:,:,:,nt)./(dzMat);
end

TOT_S=(tendSL);
ADV_S=(adv_hConvSL+adv_vConvSL);
DIF_S=(dif_hConvSL+dif_vConvSL);
FRC_S=(SF);
clear tendSL adv*SL dif*SL SF

%%%%%
% NOW DO SALINITY
%%%%%
SSTAR=0*SALT;
TOT=0*SALT;
for nt=1:nTS
 rstarfac=(mygrid.Depth+ETAN(:,:,nt))./mygrid.Depth;
 SSTAR(:,:,:,nt) = mk3D(rstarfac,mygrid.mskC); 
 TOT(:,:,:,nt)=(SALT_SNAP(:,:,:,nt+1)-SALT_SNAP(:,:,:,nt))/(secPerHour*dt(nt));
end
FRC = -SALT.*FRC_V./SSTAR + FRC_S./SSTAR;
ADV = -SALT.*CON_V./SSTAR + ADV_S./SSTAR;
DIF = DIF_S./SSTAR;

toc % stop keeping track of time

clearvars -except mygrid *_T *_S *_V TOT ADV DIF FRC