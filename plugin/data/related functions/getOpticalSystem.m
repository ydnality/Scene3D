% This function allows to select whihc optical system use
% lecelting amoung a predefinited set , created earlier


function [OptStruct]=getOpticalSystem(wave,file_name)

%% INPUT
% wave: sampling wavelegth for the refractive indices in [unit]
% file_name: name of the optical system to load

%% OUTPUT
% OptStruct: structure containing several feature of the optical system



%% IMPORT DATA through selection

if exist('file_name')
    switch file_name
        case {'D-GAUSS_F_2_22deg_HFOV_original';'D-GAUSS_F_2_22deg_HFOV_original.mat';'D-Gauss original' ;'Doublet Guass Original'}
            load('D-GAUSS_F_2_22deg_HFOV_original.mat');
        case {'D-GAUSS_F_2_22deg_HFOV_scaled';'D-GAUSS_F_2_22deg_HFOV_scaled.mat';'D-Gauss scaled';'Doublet Guass Scaled'}
            load('D-GAUSS_F_2_22deg_HFOV_scaled.mat');
        case {'DoubletObjective';'DoubletObjective.mat';'Doublet Objective';'Doublet Obj'}
            load('DoubletObjective.mat');
        case {'F_4_Telephoto Camera Lens';'Telephoto Camera Lens.mat';'TelephotoCameraLens';'TelCamLens'}
            load('F_4_Telephoto Camera Lens.mat');
        case {'Fisheye_Miyamoto_Josa_1964';'Fisheye_Miyamoto_Josa_1964.mat';'Fisheye';'Fisheye1964';'Miyamoto1964'}
            load('Fisheye_Miyamoto_Josa_1964.mat');
        case {'split_front_crown_triplet_F_2_15deg_HFOV_original';'split_front_crown_triplet_F_2_15deg_HFOV_original.mat';'split_front_crown_triplet';'Split Front Crown Triplet'}
            load('split_front_crown_triplet.mat');
        case {'Cooke_Triplet_Sample_Aberr';'Cooke_Triplet_Sample_Aberr.mat';'Cooke Triplet Aberration'}
            load('Cooke_Triplet_Sample_Aberr.mat');
        case {'Crown_First_Fraunhofer_Doublet';'Fraunhofer Doublet'; 'Crown First'}
            load('Crown_First_Fraunhofer_Doublet.mat');
         case {'Triplet_Object.mat';'Triplet_Object';'Triplet Object'}
            load('Triplet_Object.mat'); 
        case {'F3_5_Cooke_Triplet.mat';'F3_5_Cooke_Triplet';'F3_5 Cooke Triplet'}
            load('F3_5_Cooke_Triplet.mat'); 
            
        otherwise
            cd_old=cd('E:\Michael\PSF3D\data');
            load(uigetfile); %all Optical Systems have been save as OptStruct
            cd(cd_old);
            
    end
            
else
    cd_old=cd('E:\Michael\PSF3D\data');
    load(uigetfile); %all Optical Systems have been save as OptStruct
    cd(cd_old);
end

unit=OptStruct.unit;

%% Check Fields

if (isfield(OptStruct,'N')) 
    if size(OptStruct.N,1)==1
        %N is constant 
        OptStruct.N=repmat(OptStruct.N,size(wave,1),1); %repmat
    else
        warning('The given refr. indices do not match with sampling wavelength! Refr. indices are re-computed!')
        for j=1:length(OptStruct.N_type)
            [OptStruct.N(:,j),OptStruct.abbeNumber(:,j)]=RefractiveIndexDispersion(wave,unit,OptStruct.N_type{j});
         end
    end
else %if the refr. indices don't exist 
    for j=1:length(OptStruct.N_type)
        [OptStruct.N(:,j),OptStruct.abbeNumber(:,j)]=RefractiveIndexDispersion(wave,unit,OptStruct.N_type{j});
    end
end

% Film position
if (isfield(OptStruct,'film_z')) %Check if the film structure exist
    warning ('Film position is given!!')
else
    if exist('film_z')
        OptStruct.film_z=film_z;
    else
        prompt = ['Insert film position along z_axis [',unit,'] \n last vertex at ', num2str(OptStruct.PosZ(end)), unit,' :'];
        OptStruct.film_z = input(prompt);
    end
end

%Film imperfection from rotational symmetry [augment parameter]
if not(isfield(OptStruct,'augParam_Film')) %Check if the film structure exist
    OptStruct.augParam_Film=[0;0];
    warning ('augmented parameter for the Film set to zeros')
end
