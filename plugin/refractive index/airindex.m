function [n]=airindex(w,unit)

% Index of refraction for standard air @15°C   {http://refractiveindex.info/?shelf=other&book=air&page=Ciddor}
%Info at 587.6 nm
%n_ref= 1.0002772; 
%chromatic dispersion (dn/dw)=-0.000015908 um^(-1)
%Abbe numbers: Vd=89.30    {http://refractiveindex.info/legacy/include/formulae-abbe.html}


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
        case {'um'}
            w_um=w;
        
end

% equation parameters for unit= um

B=[5792105*1e-8, 167917*1e-8]; %numerator coeff.s [ref. to um]
C=[238.0185,57.362]; %den coeff.s [ref. to um]

n=(1+(B(1))./(C(1)-w_um.^(-2))+(B(2))./(C(2)-w_um.^(-2)));


if size(n,1)~=length(w) 
    n=n'; %set in column
end