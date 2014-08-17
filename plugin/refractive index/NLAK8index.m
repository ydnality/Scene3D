function [n]=NLAK8index(w,unit)

% Index of refraction for Lanthanum crown (LaK): N-LAK8  {http://refractiveindex.info/legacy/?group=SCHOTT&material=N-LAK8}

%Info at 587.6 nm
%n_ref= 1.713; 
%chromatic dispersion (dn/dw)=-0.0684 um^(-1)
%Abbe numbers: Vd=53.83; Ve=53.61     {http://refractiveindex.info/legacy/include/formulae-abbe.html}


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
B=[1.33183167, 0.546623206, 1.19084015]; %numerator coeff.s [ref. to um]
C=[0.00620023871,0.0216465439, 82.5827736]; %den coeff.s [ref. to um]

[n]=sellmeier_func(w_um,n_ord,B,C);

if size(n,1)~=length(w) 
    n=n'; %set in column
end