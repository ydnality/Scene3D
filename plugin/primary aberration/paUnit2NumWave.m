%Convert  Wavefront Error coeffs  to the wavlenght [change of unit]: 
%Smith, George. The eye and visual optical instruments. Cambridge University Press, 1997." -Chapter 33- 33.1.2.


function  [Coeff_norm]=paUnit2NumWave(Coeff)

%INPUT
%Coeff: structure with 4th order wavefront coeff
%[.W040,.W131,.W222,.W220,.W.311] ...
... or Seidel coeff [.SI,.SII,.SIII,.SIV,.SV] ...
... or Peak Value coeff [.W40,.W31,.W22,.W20,.W11]


%OUTPUT
%waveCoeff_normalized: structure with wavefront coeff normalized [.W040,.W131,.W222,.W220,.W.311,.W400]



%% GET WAVELENGTH
wave=getfield(Coeff,'wave');


%% GET FIELD names of the INPUT
fnames=fieldnames(Coeff);

%% COMPUTE WAVEFRONT COEFFs
Coeff_norm=struct;
for si=1:length(fnames)
    switch fnames{si}
    case {'W040';'W131';'W222';'W220';'W311';'SI';'SII';'SIII';'SIV';'SV';'W40';'W31';'W22';'W20';'W11'}
        wX=getfield(Coeff,fnames{si});
        for ii=1:size(wX,2)
            wN(:,ii)=wX(:,ii)./wave;
        end
    Coeff_norm=setfield(Coeff_norm,fnames{si},wN);
    case {'unit'}
        Coeff_norm=setfield(Coeff_norm,'unit','#wave');
    otherwise 
    end
end