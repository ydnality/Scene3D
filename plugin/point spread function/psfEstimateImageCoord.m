% Estimate image coordinate for the given pupil coordinate (that can be
% polar, cartesian and real or normalized


function [x_im,y_im]=psfEstimateImageCoord(coordType,P1,P2,wave,NA,vararging)

% INPUT
% coordType: specify the type of input pupil coordinate {'normalized cartesian','normalized polar',''real cartesian','real polar'}
% P1: first coordinate of pupil function  {ro or x}
% P2: second coordinate of pupil function  {thetha or y}
% wave: wavelenght vector
% NA: wavelength dependent Numerical Aperture ( column vector)
% vararging : {1}: Exit Pupil Diameter (wavelenght dependent) [needed for convert real to normalized pupil coordinate]

% OUTPUT
% x_im: wavelength dependent x-coordinate in image space 
% y_im: wavelength dependent y-coordinate in image space


%% Get some relevant parameter

nW=size(NA,1); %number of wavelength

if not(size(NA,1)==nW)
    error('Numerical Aperture number of elements is different to wavelenght!!!')
end

switch coordType
    case{'normalized cartesian';'normalizedcartesian';'norm cartesian';'norm cart';'normcart'}
        x_p=P1;y_p=P2; 
        %Estimate 'real' image coordinate
         for li=1:nW       
            if (size(P1,1)==1) || (size(P1,2)==1) %case of INPUT as vector
                [x_im(li,:),y_im(li,:)]=psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            else % case of INPUT as matrix
                [x_im(:,:,li),y_im(:,:,li)]=psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            end
        end
    case{'normalized polar';'normalizedpolar';'norm pol';'norm pol';'normpol'}
        ro=P1;theta=P2; 
        % Convert polar to cartesian coordinate
        x_p=ro.*cos(theta);y_p=ro.*sin(theta);
        %Estimate 'real' image coordinate
        for li=1:nW       
            if (size(P1,1)==1) || (size(P1,2)==1) %case of INPUT as vector
                [x_im(li,:),y_im(li,:)]=psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            else % case of INPUT as matrix
                [x_im(:,:,li),y_im(:,:,li)]=psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            end
        end

    case{'real cartesian';'realcartesian';'real cartesian';'real cart';'realcart'}
        
    case{'real polar';'realpolar';'real pol';'real pol';'realpol'}
        
    otherwise
        error ([coordType,' not valid as pupil coordinate type !'])
end






