% Estimate the Apodization Function of the pupil function for the specified
% type 

function [ApodW]=psfGetApodFun(typeApod,ro,theta,debugF,varargin)

% INPUT
% typeApod: type of aberration {'uniform','defocus'}
% ro: vector or matrix for radial coordinate
% theta: vector of matrix for angle coordinate
% varargin  {1}: specify if plot the result for DEBUG

% OUTPUT
% W: Apodization function



switch typeApod
    case {'uniform';'unif';'default'}          
        ApodW=circMask(ro);
    case {'Stiles Crawford';'StilesCrawford';'Stiles-Crawford'}
        % Human eye-Chapter 36- Handbook of optical system
        % {http://www.wiley-vch.de/books/sample/3527403809_c01.pdf}
        % I(r)=10^(-p*r^2); p: parameter (0.066 is default)
         if nargin>4
             p=varargin{1}; % get new parameter             
         else
             p=0.066; % default parameter [mm^-1]
         end
         if nargin>5
             Diam=varargin{2}; % aperture diameter
         else
%              Diam=8; %mm Assumed aperture of 8 mm
             Diam=6; %mm Assumed aperture of 6 mm
         end
         
         % Compute effective parameter
         p_eff=p*(Diam/2).^2;
         
         Mask=circMask(ro);
%          ApodW=10.^(-p.*ro.^2).*Mask; %Gaussian aposization
         ApodW=10.^(-p_eff.*ro.^2).*Mask; %Gaussian aposization
         
         % REFERENCE: D. A. Atchison, A. Joblin and G. Smith,
                       % Influence of the Stiles–Crawford effect apodization on spatial visual performance, JOSA A 15, 2545 (1998).
         
    case {'customized'}
        error ('This section has to be filled')
    otherwise
        error(['Specify a valid Apodization type']) 
end

%% DEBUG
if not(isempty(debugF)) && (debugF==1)
    x_p=ro.*cos(theta); % x-axis normalized pupil coordinate
    y_p= ro.*sin(theta); % y-axis normalized pupil coordinate
     surf(x_p,y_p,ApodW)
%     contour(x_p,y_p,PhaseW(:,:, nW0))
%     contourf(x_p,y_p,PhaseW(:,:, nW0))
    xlabel('normalized x-axis '),ylabel('normalized y-axis'),zlabel('Relative intensity')
    title(['Apodization ', typeApod, ' function'])
    
end