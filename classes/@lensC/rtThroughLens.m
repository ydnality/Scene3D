function obj = rtThroughLens(obj, rays, nLines, rtType)
% Gateway routine for tracing ray through the lens
%
%   lens.rtThroughLens(rays,nLines,rtType)
%
% Inputs the rays at the entrance aperture and calculates the rays
% at the exit aperture.
%
% Inputs:  
%    rays -   Rays at the entrance aperture
%    nLines - Number of lines to draw (0 means don't draw).
%             This can be a number, or it can be a
%             struct with .spacing and .numLines
%    rtType - 'ideal' or 'realistic'.  that switches between
%             rtIdealThroughLens and rtRealisticThroughLens
%
% Is this right (BW)? 
% On return, the input variable rays, which starts out representing the
% rays at the entrance aperture, is changed to be the position and
% direction of the rays at the exit aperture.
%
% Simplify the code.
% Why is there similar code in ppsfCamera?

% The order is from furthest from film to film, which is how the rays pass
% through the optics.

if (ieNotDefined('nLines')), nLines = false; end
if (ieNotDefined('rtType')), rtType = 'realistic'; end

rtType = ieParamFormat(rtType);
switch rtType
    case 'ideal'
        obj = obj.rtIdealThroughLens(rays, nLines);
    case 'realistic'
        obj = obj.rtRealisticThroughLens(rays, nLines);
    case 'linear'
        error ('not implemented yet');
    otherwise
        error ('unknown ray trace type');
end

end

