function [n]=STIH13index(w,unit)

% Index of refraction for ____ (): S-TIH13 {http://refractiveindex.info/legacy/?group=OHARA&material=S-TIH13}

%Info at 587.6 nm
%n_ref= 1.74076; 
%chromatic dispersion (dn/dw)=-0.134 um^(-1)
%Abbe numbers: Vd=27.79; Ve=2    7.56 {http://refractiveindex.info/legacy/include/formulae-abbe.html}


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
B=[1.62224674,2.93844586*1e-1,1.99225164]; %numerator coeff.s [ref. to um]
C=[1.18368386*1e-2,5.90208025*1e-2,1.71959976*1e2]; %den coeff.s [ref. to um]

[n]=sellmeier_func(w_um,n_ord,B,C);

if size(n,1)~=length(w) 
    n=n'; %set in column
end