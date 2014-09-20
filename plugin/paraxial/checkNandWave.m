function [newN]=checkNandWave(N,wavelength)

%The function checks the size of N and wavelength and corrects if possible

%INPUT
%N: scalar or column vector
%wavelength: scalar o column vector

%OUTPUT
%newN: corrected refractive index vector


if size(N,1)==size(wavelength,1)
    newN=N;
else
    if (size(wavelength,1)>1) && (size(N,1)==1) && (size(N,2)==1) %case of no dispersive material
        newN=repmat(N,size(wavelength,1),1);
    else
        warning('Refractive index and sampling wavelength do not match!! The surface struc misses this information')
        newN=[];
    end
end