function [n]=SFSL5index(w,unit)

% Index of refraction for ____ (): S-FSL5 {http://refractiveindex.info/legacy/?group=OHARA&material=S-FSL5}

%Info at 587.6 nm
%n_ref= 1.48749; 
%chromatic dispersion (dn/dw)=-0.0361 um^(-1)
%Abbe numbers: Vd=70.24; Ve=70.04     {http://refractiveindex.info/legacy/include/formulae-abbe.html}


%% INPUT
% w: (1xP) p different wavelength in [mm]


%OUTPUT
%n: (Px1) refr. index at p different wavelength [adim] NB set in column


%Adapt wavelengths to coeff. unit

switch unit
        case 'mm'
            w_um=w*1e3;
        case 'm'
            w_um=w*1e6;
        
end

%Sellmeier equation parameters for unit= um
n_ord=3; %order number
B=[1.17447043,1.40056154*1e-2,1.19272435]; %numerator coeff.s [ref. to um]
C=[8.41855181*1e-3,-5.81790767*1e-2,1.29599726*1e2]; %den coeff.s [ref. to um]

[n]=sellmeier_func(w_um,n_ord,B,C);

if size(n,1)~=length(w) 
    n=n'; %set in column
end