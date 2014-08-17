function [n]=NSF57index(w,unit)

% Index of refraction for flit materials N-SF57   {http://refractiveindex.info/legacy/?group=SCHOTT&material=N-SF57}
%Info at 587.6 nm
%n_ref= 1.84665; 
%chromatic dispersion (dn/dw)=-0.178 um^(-1)
%Abbe numbers: Vd=23.78  Ve=23.59    {http://refractiveindex.info/legacy/include/formulae-abbe.html}


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
B=[1.87543831,0.37375749,2.30001797]; %numerator coeff.s [ref. to um]
C=[0.0141749518,0.0640509927,177.389795]; %den coeff.s [ref. to um]

[n]=sellmeier_func(w_um,n_ord,B,C);

if size(n,1)~=length(w) 
    n=n'; %set in column
end