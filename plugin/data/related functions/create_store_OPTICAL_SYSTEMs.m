%create_store_OPTICAL_SYSTEMs.m

%With this script is possible to create and stored a set of possible Optical Systems
%that can be uploaded when needed.

%% # D-GAUSS 
% %# 	
% %# Lens Design Foundammental, p.41"  {Its given also abcd matrix=	

OptStruct.name='Cemented double objective';
OptStruct.source='Lens Design Foundammental, p.41';

OptStruct.unit='mm';

OptStruct.Radius=[7.3895,-5.1784,-16.225];% surface radius of curvature
OptStruct.DistZ=[0,1.05,0.4]; 
OptStruct.N_type={'set'};
OptStruct.N=[1.517,1.649,1];
% OptStruct.N=repmat(N,length(wavelength),1);
OptStruct.type={'refr','refr','refr'};
for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.Diam=[4,4,4];% Aperture Diameter

OptStruct.features.efl=12; %mm
OptStruct.features.Fnum=3;%mm
OptStruct.features.ABCD=[0.9404865,0.9390067;-0.0833332,0.9800774]; %ABCD matrix coeff
OptStruct.features.note={'unit in mm', 'marginal ray at heigth 2 mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('DoubletObjective.mat')
    %Delete existing file
    delete('DoubletObjective.mat')
end
save('DoubletObjective.mat','OptStruct','-mat')
cd(cd_old)

% Clear work space 
clearvars 'OptStruct'

%% # D-GAUSS F/2 22deg HFOV	
% %# US patent 2,673,491 Tronnier"	
% %# Modern Lens Design, p.312"	
% SCALED VERSION TO 50 mm of focal length

OptStruct.name='D-GAUSS F/2 22deg HFOV, scaled version	';
OptStruct.source='Modern Lens Design, p.312';

OptStruct.unit='mm';

OptStruct.Radius=[29.475,84.83,19.275,40.77,12.75,0,-14.495,40.77,-20.385,437.065,-39.73];% surface radius of curvature
OptStruct.DistZ=[0,3.76,0.12,4.025,3.275,5.705,4.5,1.18,6.065,0.19,3.22]; 
OptStruct.N_type={'set'};
OptStruct.N=[1.67,1,1.67,1.699,1,1,1.603,1.658,1,1.717,1];
% OptStruct.abbeNumber=[47.1,89.3,47.1,30.1,89.3,89.3,38,57.3,89.3,48,89.3]; %Air abbeNumber 89.3 {http://refractiveindex.info/legacy/?group=GASES&material=Air}

OptStruct.type={'refr','refr','refr','refr','refr','diaphragm','refr','refr','refr','refr','refr'};
for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.Diam=[25.2,25.2,23,23,18,17.1,17,20,20,20,12];% Aperture Diameter

OptStruct.features.efl=50; %mm
OptStruct.features.Fnum=2;%mm
OptStruct.features.halfFoVdeg=[22]; %Field of View in deg
OptStruct.features.note={'unit in mm', 'scaled to have a focal length of 50mm '};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('D-GAUSS_F_2_22deg_HFOV_scaled.mat')
    %Delete existing file
    delete('D-GAUSS_F_2_22deg_HFOV_scaled.mat')
end
save('D-GAUSS_F_2_22deg_HFOV_scaled.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'

%% # D-GAUSS F/2 22deg HFOV	
% %# US patent 2,673,491 Tronnier"	
% %# Modern Lens Design, p.312"	
% ORIGINAL VERSION VERSION TO 50 mm of focal length

OptStruct.name='D-GAUSS F/2 22deg HFOV, ORIGINAL VERSION ';
OptStruct.source='Modern Lens Design, p.312';
OptStruct.source_link='http://graphics.stanford.edu/courses/cs348b-03/lectures/camera.pdf';

OptStruct.unit='mm';

OptStruct.Radius=[58.950,169.660,38.550,81.540,25.500,0,-28.990,81.540,-40.770,874.130,-79.460];% surface radius of curvature
OptStruct.DistZ=[0,7.520,0.240,8.050,6.550,11.410,9,2.360,12.130,0.380,6.440]; 
OptStruct.N_type={'set'};
OptStruct.N=[1.67,1,1.67,1.699,1,1,1.603,1.658,1,1.717,1];
% OptStruct.abbeNumber=[47.1,89.3,47.1,30.1,89.3,89.3,38,57.3,89.3,48,89.3]; %Air abbeNumber 89.3 {http://refractiveindex.info/legacy/?group=GASES&material=Air}

OptStruct.type={'refr','refr','refr','refr','refr','diaphragm','refr','refr','refr','refr','refr'};
for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.Diam=[50.4,50.4,46.0,46.0,36.0,34.2,34.0,40,40,40,40];% Aperture Diameter

%Set Film position
OptStruc.film_z=OptStruct.PosZ(end)+72.228;

OptStruct.features.efl=10; %mm
OptStruct.features.Fnum=2;%mm
OptStruct.features.halfFoVdeg=[22]; %Field of View in deg
OptStruct.features.note={'unit in mm', 'ORIGINAL VESTION '};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('D-GAUSS_F_2_22deg_HFOV_original.mat')
    %Delete existing file
    delete('D-GAUSS_F_2_22deg_HFOV_original.mat')
end
save('D-GAUSS_F_2_22deg_HFOV_original.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'

%% # F/2.0, 15 HFOV	split-front crown triplet
%# EP#237, 212–1925.	
%# Moden Optical Engineering , W Smith 4th edition, pag 485"	

OptStruct.name='F/2.0, 150 HFOV	split-front crown triplet ';
OptStruct.source=' Moden Optical Engineering , W Smith, pag 485';

OptStruct.unit='mm';

OptStruct.Radius=[51,-441,35.3,47.8,-254.8,28.30,0,107.8,-60.3];% surface radius of curvature
OptStruct.DistZ=[0,8.80,0.03,7.8,8.4,2,10,19.4,4.9]; 


OptStruct.Diam=2*[25.1,25.1,22,20,18,16,15.69,16,16];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr','refr','refr','refr','diaphragm','refr','refr'};

OptStruct.N_type={'SK11','air0','SK11','air0','SF2','air0','air0','SK11','air0'};

% for j=1:length(OptStruct.N_type)
%     OptStruct.N(:,j)=RefractiveIndexDispersion(wavelength,unit,OptStruct.N_type{j});
% end
%FILM
OptStruct.film_z=OptStruct.PosZ(end)+56.89;
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=99.79; %mm
OptStruct.features.bfl=56.89; %mm
OptStruct.features.NA=0.251; %mm
OptStruct.features.Fnum=2;%mm
OptStruct.features.halfFoVdeg=[15]; %Field of View in deg

OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('split_front_crown_triplet_F_2_15deg_HFOV_original.mat')
    %Delete existing file
    delete('split_front_crown_triplet_F_2_15deg_HFOV_original.mat')
end
save('split_front_crown_triplet_F_2_15deg_HFOV_original.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'

%% # F/8 90 HFOV Fisheye Miyamoto Josa 1964
%# 	
%# Moden Optical Engineering , W Smith, pag 547"	

OptStruct.name='F/8 90 HFOV Fisheye Miyamoto Josa 1964 ';
OptStruct.source=' Moden Optical Engineering , W Smith, pag 547';

OptStruct.unit='mm';

OptStruct.Radius=[599.383,235.825,605.513,111.094,-452.384,127.733,462.892,inf,inf,inf,38507.649,95.081,-162.638,1376.167,177.275,-400.339,-337.536];% surface radius of curvature
OptStruct.DistZ=[0,35.030,190.161,30.025,120.102,10.008,45.038,25.021,15.013,36.281,13.762,10.008,110.093,130.110,20.017,150.127,18.766]; 


OptStruct.Diam=2*[448.4,234,251.8,110.1,93.5,93.5,93.5,65.4,65.5,15.8,84.1,84.1,84.1,139.0,139,139,139];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr','refr','refr','refr','refr','flat','flat','diaphragm','refr','refr','refr','refr','refr','refr','refr'};

% OptStruct.N_type={'BK7','air0','FK5','air0','FK5','SF56','air0','K3','air0','air0','SF56','LAF2','air0','SF56','BSF52','BASF6','air0'};  % BSF52 is not available sub with BASF7 (similar index and Abbe number
OptStruct.N_type={'BK7','air0','FK5','air0','FK5','SF56','air0','K3','air0','air0','SF56','LAF2','air0','SF56','BASF7','BASF6','air0'}; 
% OptStruct.N=[1.517,1,1.487,1,1.487,1.785,1,1.518,1,1,1.785,1.744,1,1.785,1.702,1.668,1];
% OptStruct.abbeNumber=[64.2,89.3,70.4,89.3,70.4,26.1,89.3,59,89.3,89.3,26.1,44.7,89.3,26.1,41.0,41.9,89.3];%Air abbeNumber 89.3 {http://refractiveindex.info/legacy/?group=GASES&material=Air}

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+150.119;
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=100; %mm
OptStruct.features.bfl=150.1; %mm
OptStruct.features.NA=-0.0626; %mm
OptStruct.features.Fnum=8;%mm
OptStruct.features.FoVdeg=[180]; %Field of View in deg

OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('Fisheye_Miyamoto_Josa_1964.mat')
    %Delete existing file
    delete('Fisheye_Miyamoto_Josa_1964.mat')
end
save('Fisheye_Miyamoto_Josa_1964.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'


%% # F/4.0 Telephoto Camera Lens
%# 	
%# Moden Optical Engineering , W Smith, pag 549"	

OptStruct.name='F/4.0 Telephoto Camera Lens ';
OptStruct.source=' Moden Optical Engineering , W Smith, pag 549';

OptStruct.unit='mm';

OptStruct.Radius=[27.03,-176.93,30.66,76.46,-212.41,36.22,506.55,-67.74,-20.57,-78.32];% surface radius of curvature
OptStruct.DistZ=[0,4.5,0.1,3,1.4,2,30.84,2.5,1.66,1.5]; 


OptStruct.Diam=2*[12.6,12.6,12.6,12.6,12.6,11,9.3,9.3,9.3,9.3];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr','refr','refr','refr','refr','refr','refr','refr'};

OptStruct.N_type={'BK10','air0','BK10','air0','SF5','air0','SF57','air0','LaK8','air0'};

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+32.46;
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=99.9; %mm
OptStruct.features.bfl=32.46; %mm
OptStruct.features.NA=0.125; %mm
OptStruct.features.Fnum=4;%mm


OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('F_4_Telephoto Camera Lens.mat')
    %Delete existing file
    delete('F_4_Telephoto Camera Lens.mat')
end
save('F_4_Telephoto Camera Lens.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'


%% # F/3.5, 11.9° HFOV	Cooke Triplet
%# EP#237, 212–1925.	
%# Moden Optical Engineering , W Smith, pag 118"	

OptStruct.name='F/3.5, 11.9° HFOV	Cooke Triplet';
OptStruct.source=' Moden Optical Engineering , W Smith, pag 118';

OptStruct.unit='mm';

OptStruct.Radius=[37.4,-341.48,-42.65,36.4,0,204.52,-37.05];%% surface radius of curvature
OptStruct.DistZ=[0,5.9,12.93,2.5,2,9.85,5.9]; 


OptStruct.Diam=2*[14.7,14.7,10.8,10.8,10.3,11.6,11.6];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr','refr','diaphragm','refr','refr'};

OptStruct.N_type={'SK4','air0','SF2','air0','air0','SK4','air0'};

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+[77.405];
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=101.181; %mm
OptStruct.features.bfl=77.405; %mm
OptStruct.features.NA=0.144; %mm
OptStruct.features.Fnum=3.47;%mm
OptStruct.features.FoVdeg=11.9; %[

OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('F3_5_Cooke_Triplet.mat')
    %Delete existing file
    delete('F3_5_Cooke_Triplet.mat')
end
save('F3_5_Cooke_Triplet.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'

%% # Cooke Triplet SAMPLE for ABERRATION
%# L3_OPTI517_Aberration .	
%# Prof. Sasian {Arizona University} ;http://fp.optics.arizona.edu/Sasian/2013opti517/L3_OPTI517_Aberrations.pdf	
% and Lens Design Fundamental- Chapter 10-pag 143 
OptStruct.name='Cooke Triplet sample for Aberration';
OptStruct.source1=' http://fp.optics.arizona.edu/Sasian/2013opti517/L3_OPTI517_Aberrations.pdf';
OptStruct.source2=' Chapter 10_ Abberation Coeffs .pdf';

OptStruct.unit='mm';

OptStruct.Radius=[22.05,371.58,-30.1,20.01,64.47,-23.48];%% surface radius of curvature
OptStruct.DistZ=[0,4.83,5.86,0.98,4.82,5]; 


OptStruct.Diam=2*[10,10,9,9,10,10];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr','refr','refr','refr'};

OptStruct.N_type={'LAK9','air0','SF5','air0','LAK9','air0'};

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+[77.405];
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=50; %mm
% OptStruct.features.bfl=77.405; %mm
% OptStruct.features.NA=0.144; %mm
OptStruct.features.Fnum=4;%mm
OptStruct.features.FoVdeg=40; 
OptStruct.features.waveCoeff587nm=[2.66;0.02;-11.13;10.77;6.19]; %wavefront aberration coeff at 587 nm

OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('Cooke_Triplet_Sample_Aberr.mat')
    %Delete existing file
    delete('Cooke_Triplet_Sample_Aberr.mat')
end
save('Cooke_Triplet_Sample_Aberr.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'

%% # Triplet object SAMPLE for ABERRATION
%#  .	
%# Kidget MJ  Fundamental Optical Design 2002- Chapter 6	
% {Also in "Lens Design Fundamental (2nd ed) pag 420- Kingslake}

OptStruct.name='Triplet Object';
OptStruct.source= 'Kidget MJ  Fundamental Optical Design 2002- Chapter 6';

OptStruct.unit='mm';

OptStruct.Radius=[42.98790,-248.07740,-38.21035,43.95894,656.66349,-33.50754];%% surface radius of curvature
OptStruct.DistZ=[0,4,10.51018,2.5,9.86946,4.5,86.48643]; 


OptStruct.Diam=2*[11.5,11.5,9.852,8.8885,11,11];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr','refr','refr','refr'};

OptStruct.N_type={'SK16','air0','S-F4','air0','SK16','air0'};

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+[86.48643];
OptStruct.film_Diam=2*[37.166];
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=100; %mm
OptStruct.features.Fnum=4.5;%mm
OptStruct.features.FoVdeg=40; %deg
OptStruct.features.EnP_Diam=22.22; %mm


OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('Triplet_Object.mat')
    %Delete existing file
    delete('Triplet_Object.mat')
end
save('Triplet_Object.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'


%% # Crown-first Fraunhofer doublet
%#  .	
%# Kidget MJ  Fundamental Optical Design 2002- Chapter 8-pag 172	
% 
OptStruct.name='Crown-first Fraunhofer doublet';
OptStruct.source= 'Kidget MJ  Fundamental Optical Design 2002- Chapter 6';

OptStruct.unit='mm';

OptStruct.Radius=[64.10,-43.249,-183];%% surface radius of curvature
OptStruct.DistZ=[0,3.5,1.5]; 


OptStruct.Diam=2*[10,9.926,9.86];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'refr','refr','refr'};

OptStruct.N_type={'SK11','F19','air0'};

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+[97.44950];
OptStruct.film_Diam=2*[1.763];
OptStruct.augParam_Film=[0;0];

%Optical system features
OptStruct.features.efl=100; %mm
% OptStruct.features.Fnum=4.5;%mm
% OptStruct.features.FoVdeg=40; %deg
OptStruct.features.LagInv=-0.1746; %mm
% Parax Ray Tracing
OptStruct.features.paraxRayTracing.marginalRay.y=[10.00000,9.80313,9.74498 ];%height
OptStruct.features.paraxRayTracing.marginalRay.u=[0.00000,-0.05625,-0.03877,-0.100000 ];%angle
OptStruct.features.paraxRayTracing.marginalRay.Delta=[-0.03597,0.01271,-0.07674];% Delta(u/n)
OptStruct.features.paraxRayTracing.marginalRay.A=[0.15601,-0.44243,-0.15305]; % Refractive Invariant
OptStruct.features.paraxRayTracing.chiefRay.y=[0.00000,0.03907,0.05486 ];%height
OptStruct.features.paraxRayTracing.chiefRay.u=[0.01746,0.01116,0.01053,0.017349 ];%angle
OptStruct.features.paraxRayTracing.chiefRay.A=[0.01746,0.01604,0.01705]; % Refractive Invariant
%Seidel Coeffs
OptStruct.features.SeidelCoeffs.SI=0.001889;OptStruct.features.SeidelCoeffs.SII=-0.000088; 
OptStruct.features.SeidelCoeffs.SIII=0.000295;OptStruct.features.SeidelCoeffs.SIV=0.000210; 
OptStruct.features.SeidelCoeffs.SV=0.000002; 
%Surface contribution
OptStruct.features.SeidelCoeffs.log.SI=[0.008754,-0.024383,0.017518];
OptStruct.features.SeidelCoeffs.log.SII=[0.000979,0.000884,-0.001952];
OptStruct.features.SeidelCoeffs.log.SIII=[0.000110,-0.000032,0.000217];
OptStruct.features.SeidelCoeffs.log.SIV=[0.000171,-0.000028,0.000066];
OptStruct.features.SeidelCoeffs.log.SV=[0.000031 ,0.000002,-0.000032];

OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('Crown_First_Fraunhofer_Doublet.mat')
    %Delete existing file
    delete('Crown_First_Fraunhofer_Doublet.mat')
end
save('Crown_First_Fraunhofer_Doublet.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'


%% # Singlet plus Diaphragm
%# L3_OPTI517_Aberration .	
%# Prof. Sasian {Arizona University} ;http://fp.optics.arizona.edu/Sasian/2013opti517/L3_OPTI517_Aberrations.pdf	
% and Lens Design Fundamental- Chapter 2-pag 143 
OptStruct.name='Singlet plus Diaphragm';
OptStruct.source1=' Chapter2_Basice concepts in geometrical optics.pdf';


OptStruct.unit='mm';

OptStruct.Radius=[Inf,Inf,-51.680];%% surface radius of curvature
OptStruct.DistZ=[0,30.775,5]; 


OptStruct.Diam=[12.5,30,30];% Aperture Diameter

for i=1:length(OptStruct.DistZ)
    OptStruct.PosZ(i)=sum(OptStruct.DistZ(1:i));%Element position along optical axis
end
OptStruct.type={'diaphragm','refr','refr'};

OptStruct.N_type={'air0','BK7','air0'};

%FILM
OptStruct.film_z=OptStruct.PosZ(end)+[100];
OptStruct.augParam_Film=[0;0];

%Optical system features
% OptStruct.features.efl=50; %mm
% OptStruct.features.bfl=77.405; %mm
% OptStruct.features.NA=0.144; %mm
% OptStruct.features.Fnum=4;%mm
OptStruct.features.FoVdeg=30; %[

OptStruct.features.note={'unit in mm'};%mm
cd_old=cd('E:\Michael\PSF3D\data');

if exist('Singlet_plus_Diaphragm.mat')
    %Delete existing file
    delete('Singlet_plus_Diaphragm.mat')
end
save('Singlet_plus_Diaphragm.mat','OptStruct','-mat')
cd(cd_old)

clearvars 'OptStruct'
