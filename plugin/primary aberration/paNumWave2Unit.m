%Convert  Wavefront Error coeffs from number of wavelength to unit [mm or um or nm]: 
%Smith, George. The eye and visual optical instruments. Cambridge University Press, 1997." -Chapter 33- 33.1.2.


function  [Coeff_unit]=paNumWave2Unit(Coeff)

%INPUT
%Coeff: structure with wavefront coeff normalized for wavelength [.W040,.W131,.W222,.W220,.W.311,.W400]



%OUTPUT
%Coeff_unit: structure with 4th order wavefront coeff
%[.W040,.W131,.W222,.W220,.W.311] ...
... or Seidel coeff [.SI,.SII,.SIII,.SIV,.SV] ...
... or Peak Value coeff [.W40,.W31,.W22,.W20,.W11]



%% GET WAVELENGTH
wave=getfield(Coeff,'wave');
unit=getfield(Coeff,'unit');


%% GET FIELD names of the INPUT
fnames=fieldnames(Coeff);

%% COMPUTE WAVEFRONT COEFFs
Coeff_unit=struct;
for si=1:length(fnames)
    switch fnames{si}
    case {'W040';'W131';'W222';'W220';'W311';'SI';'SII';'SIII';'SIV';'SV';'W40';'W31';'W22';'W20';'W11'}
        wX=getfield(Coeff,fnames{si});
        for ii=1:size(wX,2)
            wN(:,ii)=wX(:,ii).*wave;
        end
    Coeff_unit=setfield(Coeff_unit,fnames{si},wN);
    case {'unit'}
        Coeff_unit=setfield(Coeff_unit,'unit',unit);
    case {'wave'}
        Coeff_unit=setfield(Coeff_unit,'wave',wave);
    otherwise 
    end
end