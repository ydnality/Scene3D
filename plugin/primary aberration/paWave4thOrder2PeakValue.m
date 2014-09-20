% Adjust wavefront coeff according to the real field height since they are
% compute for the max field height at the given depth

function [peakCoeff]=paWave4thOrder2PeakValue(waveCoeff,y,y_max)


%INPUT:
%waveCoeff: [.W040,.W131,.W222,.W220,.W.311,.W400] "forth order coeffs"
%      are accepted also Seidel Coeff [.SI,.SII,.SIII,.SIV,.SV]
%y: real field height
%y_max: max field height used to compute Seidel and converted Wavefront Coeffs


%OUTPUT:
% peakCoeff

%% CHECK is are available Wavefront 4th Order Coeffs or Seidel
fnames=fieldnames(waveCoeff);
type='w4th'; % Ground Hypothesis
for ci=1:length(fnames)
    switch fnames{ci}
    case {'SI';'SII';'SIII';'SIV';'SV'}
        type='seidel';
    case {'W040';'W131';'W222';'W220';'W311'}
    case {'wave','unit'}
    otherwise
    end
end

if strcmp (type,'seidel')
    waveCoeff=paSeidel2Wave4thOrder(waveCoeff); %Compute 4th order wavefront coeffs
end
%% SET wave and unit field in output
peakCoeff.wave=waveCoeff.wave;
peakCoeff.unit=waveCoeff.unit;

%% Compute ray field ratio
ry=y./y_max;

%Flag
PDT_flag=0;
ind_name=[];
ind=1;



for ci=1:length(fnames)
    
    switch fnames{ci} %From "Smith, George. The eye and visual optical instruments. Cambridge University Press, 1997."
                       % Chapter 33- eq. 33.60                                 
                       
        % Fourth order coeff
        case {'W040'}
           W40=getfield(waveCoeff,fnames{ci});
           peakCoeff.W40=W40;
        case {'W131'}
           W31=getfield(waveCoeff,fnames{ci});
           peakCoeff.W31=W31.*ry;
        case {'W222'}
            W22=getfield(waveCoeff,fnames{ci});
            peakCoeff.W22=W22.*ry.^2;
        case {'W220'}
            W20=getfield(waveCoeff,fnames{ci});
            peakCoeff.W20=W20.*ry.^2;
        case {'W311'}
            W11=getfield(waveCoeff,fnames{ci});
            peakCoeff.W11=W11.*ry.^3;
        otherwise
          PDT_flag=1; 
          ind_name(ind)=ci;
          ind=ind+1;
    end
    
%     SI=getfield(SeidelCoeff,fnames{si});
end


if PDT_flag==1
    for ii=1:length(ind_name)
        switch fnames{ind_name(ii)}
%% If Piston-Defocus-Tilt add to realtive value
        case {'W000'} %piston
           piston=getfield(waveCoeff,fnames{ind_name(ii)});
           peakCoeff.W000=piston;
        case {'W020'} %Defocus
           defocus=getfield(waveCoeff,fnames{ind_name(ii)});
           peakCoeff.W20=peakCoeff.W20+defocus;
       case {'W011'} %Tilt
           tilt=getfield(waveCoeff,fnames{ind_name(ii)});
           peakCoeff.W11=peakCoeff.W11+tilt;
        end
    end
end
