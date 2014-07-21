function obj = rtThroughLens(obj, rays, nLines, rtType)
% Ray trace rays through the lens
%
% lens.rtThroughLens(rays,nLines,rtType)
%
% Inputs:  lens object,rays at the entrance aperture,
% calculates the rays at the exit aperture.  The exiting rays
% are represented by the
%
% On return, the input variable rays, which starts out
% representing the rays at the entrance aperture, is changed to
% be the position and direction of the rays at the exit
% aperture.
%
% Simplify the code.
% Also, why is there similar code in ppsfCamera?

% The order is from furthest from film to film, which is also
% how the rays pass through the optics.

if (ieNotDefined('nLines')), nLines = false; end
if (ieNotDefined('rtType')), rtType = 'realistic'; end
rtType = ieParamFormat(rtType);
switch rtType
    case 'ideal'
        obj = obj.rtIdealThroughLens(rays, nLines);
    case 'realistic'
        obj = obj.rtRealisticThroughLens( rays, nLines);
    case 'linear'
        error ('not implemented yet');
    otherwise
        error ('unknown ray trace type');
end

end

