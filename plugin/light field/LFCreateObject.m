% This function is aimed at build a object decribing a source/reflectance of light
% field that can be differentiated for -type,-position,etc.

function [LFobject]=LFCreateObject(type,unit, wave,varargin)

%INPUT
%Type:  of light field source/reflectance surface
%        -point source  varargin{1}:=z_pos (positio along the optical axis) ....
           %;varargin{2}:=y_ecc (eccentricity from the optical axis);
           %varargin{3}:= type of radiance
%        -plane wave  varargin{1}:=angle of eccentricity; varargin{2}=type of radiance
%        -flat source  varargin{1}:=z_pos; varargin{2}:= y eccentricity [vector of sampled eccentricity];
%                       varargin{3}:= type of radiance

%unit: unit witch referes to eache distance e.g. 'mm' ,'m'
%wavelength: (column vector)  wavelength for sampling the radiance spectrum

%OUTPUT
%LFobject: light field object

error('LFCreateObject is deprecated');

end

%Append fields to the OUTPUT
LFobject.wave=wave;
LFobject.unit=unit;


switch type
    case {'point source';'point';'point reflectance'}
        LFobject.type='point';
        LFobject.z_pos=varargin{1}; % position along the optical axis [unit]
        LFobject.y_ecc=varargin{2}; % eccentricity from the optical axis [unit]
        switch varargin{3}
            case {'Lambertian';'Lamb';'isotropic';'isotr';'iso'}
                LFobject.rad_type='isotropic';
            case {'modulate','mask','customized'}
                LFobject.rad_type='customized';
            otherwise
                warning('Not valid the type of radiance emitted or reflected by the object')
                LFobject.rad_type='Not valid';
        end
           
    case {'plane wave';'plane'}
        LFobject.type={'plane'};
        LFobject.z_pos=Inf; % position along the optical axis [unit]
        LFobject.angle_ecc=varargin{1}; % eccentricity from the optical axis [unit]
        switch varargin{2}
            case {'Lambertian';'Lamb';'isotropic';'isotr';'iso'}
                LFobject.rad_type='isotropic';
            otherwise
                warning('Not valid the type of radiance emitted or reflected by the object')
        end
        
    case {'flat source','flat'}
        LFobject.type='flat';
        LFobject.z_pos=varargin{1}; % position along the optical axis [unit]
        LFobject.y_ecc=varargin{2}; % vector describing the  eccentricity [unit]
        switch varargin{3}
            case {'Lambertian';'Lamb';'isotropic';'isotr';'iso'}
                LFobject.rad_type='isotropic';
            case {'modulate','mask','customized'}
                LFobject.rad_type='customized';
            otherwise
                warning('Not valid the type of radiance emitted or reflected by the object')
                LFobject.rad_type='Not valid';
        end
        
    otherwise
        warning ('The type of object for Light Field is DEFINED appropriately!')
end
