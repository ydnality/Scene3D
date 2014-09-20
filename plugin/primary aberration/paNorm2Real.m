% Convert aberration coeffs from normalized polar coordinate to cartesian
% polar coordinate

function  [C]=paNorm2Real(Cnorm,ExP_diam)

%INPUT
%Cnorm: struct with normalized coeffs for aberration [ for now Seidel
%coefs (SI,SII,SIII,SIV,SV), peak value (W40,....) and 4th Order Wave Aberr
%(W040.....)
% ExP_diam: exit pupil diameter

%OUTPUT
% C: struct with the same kind of coeffs (Seidel/PeakValue/4thOrdWaveAberr)

% 
%NOTE:
% THE ZERNIKE COEFFs are just draft!! This section has to be checked!


%% Check and get unit
unit=Cnorm.unit;
wave=Cnorm.wave;
%number of wavelength sample
nw=size(wave,1);


    


%% GET  field name
fnames=paGetFieldNames(Cnorm);


% FIND TYPE OF COEFFs and the index for the radius order
switch fnames{1}(1)
    case {'W','w'}
        switch length(fnames{1})
            case 3
                C_type='peak'; %W40, W31,....
                Rind=2;
            case 4
                C_type='4thwave'; %W040, W131,....
                Rind=3;
            otherwise
                error (['This coeffs are not 4thOr WAVE\Peak Value, sample ', fnames{1}])
        end
    case {'S'}
        C_type='seidel'; %SI,SII,...
        %convert Seidel Romanic indices to Numeric  SI-> S1; SIII->S3
        Cnorm0=struct;
        for ci=1:length(fnames)
            value=getfield(Cnorm,fnames{ci});
            switch fnames{ci}
                case {'SI'}
                    Cnorm0=setfield(Cnorm0,'S1',value);
                case {'SII'}
                    Cnorm0=setfield(Cnorm0,'S2',value);                    
                case {'SIII'}
                    Cnorm0=setfield(Cnorm0,'S3',value);                    
                case {'SIV'}                    
                    Cnorm0=setfield(Cnorm0,'S4',value);                    
                case {'SV'}
                    Cnorm0=setfield(Cnorm0,'S5',value);                    
                otherwise
                    error ([fnames{ci}, 'is not accepted as Seidel Coeffs'])                    
            end
        end
        Cnorm=Cnorm0; %substitute original input   S1,S2,S3,....
        Rind=2;
    case{'C'}
        C_type='zernike';  %C00,C11,C20,...
        Rind=2; 
    otherwise
        error (['Not accepted this type of coeffs, sample ', fnames{1}])
end

%% CONVERT the Coeff according to the change of coordiantes
switch C_type;
    case {'seidel'}
        fnamesNEW=paGetFieldNames(Cnorm);
    case {'zernike';'peak';'4thwave'}
        fnamesNEW=fnames;
    otherwise
        error(['Not available ', C_type, 'as type of coeffs' ])
end


C=struct;
ExP_radius=ExP_diam/2; %get the exit pupil diameter
for ci=1:length(fnamesNEW)
    val=getfield(Cnorm,fnamesNEW{ci});
    R_ord=str2num(fnamesNEW{ci}(Rind)); %Radius coordinate order 
    valNew=val./(ExP_radius.^R_ord);  % W°= W/(R^ord)
    C=setfield(C,fnames{ci},valNew); %notice here used fname and not fnamesNEW because the field names ... of the output struct
    ... have to be the same of the input [NOTICE in case of SEIDEL they have been change]

end

%% Then append to the OUTPUT extra field (wave and unit)
%Unit has to change because its coeffs NOW is not measured as distance but
%in distance over the distance^ (order of radius )

C.wave=wave;
C.unit=[unit,'\',unit,'^radius order'];
