% Get the field names for aberration of different structure [Seidel
% coefficients, Peak Values,  4th order wavefront coeffs, zernike coeffs]
% avoid 'wave' and 'unit' field


function [fnames]=paGetFieldNames(Coeff)

%INPUT
%Coeff: struct contain different type of aberration coeff 

%OUTPUT
%fnames: names of the Coeff struct which really contain coeff




%number of Primary aberration coeffs
fnames=fieldnames(Coeff);
nfname=length(fnames); 
% eliminate field of 'wave' and 'unit'

nfs=struct;
for ai=1:nfname
    switch fnames{ai}
        case{'unit';'wave';'note';'type'}
        otherwise
            value=getfield(Coeff,fnames{ai});
            nfs=setfield(nfs,fnames{ai},value);
    end
end

%% SET OUTPUT

fnames=fieldnames(nfs);
