% Generate an imaging system according to the DATABASE 
%with specific feature about object (source) and acquisition surface (film)

function [ImagSyst]=paraxGenerateImagSystem(file_name,Aper_Diam,Obj,filmParam,n_ob,n_im,wave,unit)

%INPUT
%file_name: optical structure available among the DB
%Aper_Diam: Aperture diameter
%Obj: struct about the object  .z= position along the optical axis
%                               .y= eccentricity position
%filmParam: struct about the image sensor  .mode : default or in-focus
                                           %.waveRef: reference wave if in focus
                                           %.dim : [wdim,hdim] in unit
                                           %.pitch: [w_pixel_distance,                                           %h_pixel_distance]
%n_ob: refractive index of the object space
%n_im: refractive index of the image space
%wave: column vector of sampling wavelength
%unit: for all distance, included wavelengt [mm]


%OUTPUT
%ImagSyst: struct about the imaging System including the optical element,
%the film and the image sensor

%% LOAD THE OPTICAL SYSTEM FROM THE DATABASE
OptStruct=getOpticalSystem(wave,file_name);

% Built  the surfaces

for k=1:length(OptStruct.Radius)
    R(k)=OptStruct.Radius(k);
    if (abs(R(k))==Inf) && not(strcmp(OptStruct.type{k},'diaphragm')||strcmp(OptStruct.type{k},'diaph'))
        OptStruct.type{k}='flat';
    end
    z_pos(k)=OptStruct.PosZ(k);
    n(:,k)=OptStruct.N(:,k);
    diam(k)=OptStruct.Diam(k);
    abbeNumber(k)=OptStruct.abbeNumber(k);
    switch OptStruct.type{k}
        case {'refr','flat'}          
            surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,OptStruct.type{k},R(k),n(:,k),abbeNumber(k));
            diam(k)=OptStruct.Diam(k);
        case {'diaphragm','diaph'}
            surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,OptStruct.type{k});
            if isempty (Aper_Diam)
                diam(k)=OptStruct.Diam(k);
            else
                diam(k)=Aper_Diam; %set desidered diaphragm aperture
            end
    end
    
end
%% ASSEMPLE THE SURFACE AND CREATE  THE OPTICAL SYSTEM
[OptSyst]=paraxCreateOptSyst(surf,n_ob,n_im,unit,wave);

%% CREATE an IMAGING SYSTEM : Optical System + Film
profile='flat';
[film]=paraxCreateFilm(0,profile,filmParam.dim,filmParam.pitch,unit);
switch filmParam.mode
    case{'default'}        
        film_zpos=OptStruct.film_z;
    case {'in-focus'}
        % FAKE 
        film_zpos1=OptStruct.film_z;
        [ImagSyst1]=paraxCreateImagSyst(OptSyst,film,film_zpos1,OptStruct.augParam_Film);
        profile1='point';
        [source1]=paraxCreateObject(Obj.z,Obj.y,profile1,unit);
        %Add to the Imaging System
        [ImagSyst1]=paraxAddObject2ImagSyst(ImagSyst1,source1);
%         inWref=find(wave==filmParam.waveRef);
%         film_zpos=OptSyst.cardPoints.dFi(wave==filmParam.waveRef)+OptSyst.cardPoints.lastVertex;       
        film_zpos=ImagSyst1.object{end}.ConjGauss.z_im(wave==filmParam.waveRef);
end

[ImagSyst]=paraxCreateImagSyst(OptSyst,film,film_zpos,OptStruct.augParam_Film);




if not(isempty(Obj))
    %% Create a Object and add to the Imaging system
    profile1='point';
    [source1]=paraxCreateObject(Obj.z,Obj.y,profile1,unit);
    %Add to the Imaging System
    [ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,source1);
end
