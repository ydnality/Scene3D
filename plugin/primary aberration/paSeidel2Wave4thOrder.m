%Convert the Seidel coeffs to the Wavefront Error coeffs according to the article: 
%Sasián, José. "Introduction to Aberration in Optical Imaging System" -Chapter 10- Table 10.1.




function  [waveCoeff]=paSeidel2Wave4thOrder(SeidelCoeff)

%INPUT
%SeidelCoeff: structure with seidel coeffs [.SI,.SII,.SIII,.SIV,.SV,.SVI]


%OUTPUT
%peakCoeff: structure with the corresponding wavefront coeffs [.W040,.W131,.W222,.W220,.W.311,.W400]


%NB: Conversion is independent on the unit

%% GET FIELD names of the INPUT
fnames=fieldnames(SeidelCoeff);

%% COMPUTE WAVEFRONT COEFFs

for si=1:length(fnames)
    switch fnames{si}
        case{'SI';'S1'}
            SI=getfield(SeidelCoeff,fnames{si});
            waveCoeff.W040=SI/8;
        case{'SII';'S2'}
            SII=getfield(SeidelCoeff,fnames{si});
            waveCoeff.W131=SII/2;
        case{'SIII';'S3'}
            SIII=getfield(SeidelCoeff,fnames{si});
            waveCoeff.W222=SIII/2;
        case{'SIV';'S4'}
            SIV=getfield(SeidelCoeff,fnames{si});
            if isfield(SeidelCoeff,'SIII')
                SIII=getfield(SeidelCoeff,'SIII');
                waveCoeff.W220=(SIV+SIII)/4;
            elseif isfield(SeidelCoeff,'S3')
                SIII=getfield(SeidelCoeff,'S3');
                waveCoeff.W220=(SIV+SIII)/4;
            else
                warning('For coeff W220 (Field of Curvature) is mandatory also the SIII!!')
                waveCoeff.W220=[];
            end          
            
        case{'SV';'S5'}
            SV=getfield(SeidelCoeff,fnames{si});
            waveCoeff.W311=SV/2;
        case {'unit'}
            waveCoeff.unit=getfield(SeidelCoeff,fnames{si});
        case {'wave'}
            waveCoeff.wave=getfield(SeidelCoeff,fnames{si});
         otherwise
            warning (['The field ',fnames{si}, 'is unknown!'])
    end
end