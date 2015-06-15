function obj =  draw(obj, fHdl)
% Draw the the multi-element lens surfaces in a graph window
%
%   lens.draw
%
% fHdl:  A figure handle to use.  If not passed, then vcNewGraphWin is
%        called and that figure handle is set in the lens object field fHdl
%
% See also:  psfCamera.draw
%            That calls this lens draw and also draws the rays for the
%            point spread. 
%
% AL/BW Vistasoft Team, Copyright 2014

%% Create the figure and set the parameters

if ieNotDefined('figureHandle'), fHdl = vcNewGraphWin; axis equal;
else                             figure(fHdl);
end
obj.fHdl = fHdl;

lWidth = 2; lColor = 'k';  % Drawing parameters
xlabel('mm'); ylabel('mm');

% Make sure that 
minx = 0; maxx = 0;
miny = 0; maxy = 0;

%% We draw one surface/aperture at a time
nSurfaces = obj.get('n surfaces');
for lensEl = 1:nSurfaces
    
    % Get the current surface element
    curEl = obj.surfaceArray(lensEl);
    apertureD = curEl.apertureD;
    c = curEl.sCenter(3);
    r = curEl.sRadius;
    if (r ~= 0)
         
        % One edge of the lens sits at this z intercept
        zIntercept = curEl.get('z intercept');
        minx = min(minx,zIntercept);
        maxx = max(maxx,zIntercept);

        % Solve for the points on the curve for this lens
        % surface element.  We are drawing in the z-y plane
        % because the z-axis is the horizontal axis, and the
        % y-axis is the vertical. The center of the sphere is
        % at (0,0,z), so the formula is
        %
        %     r^2 = (x)^2 + (y)^2 + (z - c)^2
        %
        % But since we are in the x=0 plane this simplifies to
        %
        %   r^2 = (y)^2 + (z - c)^2
        %
        % We solve for y in terms of z.
        
        % Solve for the Z-range we will need
        %  r^2 = (apertureD)^2 + (Zedge - c)^2
        %  zEdge = sqrt(r^2 - (apertureD)^2 ) + c
        %
        % There will be a positive and negative solution
        zEdgeN = -sqrt(r^2 - (apertureD/2)^2 ) + c;
        zEdgeP =  sqrt(r^2 - (apertureD/2)^2 ) + c;
        
        % Choose the zEdge that is closer to the zIntercept.
        zEdge = zEdgeN;
        if abs(zEdgeN - zIntercept) > abs(zEdgeP - zIntercept)
            zEdge = zEdgeP;
        end
        
        % This is the range of z values we will consider.
        zPlot = linspace(zIntercept, zEdge, 100);
        
        % We get the positive and negative y-values
        yPlot  =  sqrt(r^2 - (zPlot - c) .^2);
        yPlotN = -sqrt(r^2 - (zPlot - c) .^2);
        
        
        % The positive solutions
        line('xData',zPlot, 'yData',yPlot,...
            'color',lColor,'linewidth',lWidth);
        
        % The negative solutions
        line('xData',zPlot, 'yData',yPlotN,...
            'color',lColor,'linewidth',lWidth);
        
        miny = min(miny, min(yPlotN(:)));
        maxy = max(maxy, max(yPlot(:)));
        
    else
        %Draw the aperture opening if radius = 0
        
        %TODO: draw the difference between specified aperture
        %from file and specified aperture from object instance
        
        %right now: take the minimum value
        curAperture = min(curEl.apertureD/2, obj.apertureMiddleD/2);
        zIntercept = curEl.sCenter(3);
        
        l = line(zIntercept * ones(2,1), -1*[curEl.apertureD/2 curAperture]);
        set(l,'linewidth',lWidth,'color',lColor);
        l = line(zIntercept * ones(2,1), [curAperture curEl.apertureD/2]);
        set(l,'linewidth',lWidth,'color',lColor);
        
        minx = min(minx,zIntercept);
        maxx = max(maxx,zIntercept);
        miny = min(miny,-1*curEl.apertureD/2);
        maxy = max(maxy,   curEl.apertureD/2);
    end
    
end

% Make sure the surfaces are all shown within the range
% We believe that maxx and maxy are always positive, and minx and miny are
% always negative.  But maybe we should deal with the sign issue here for
% generality in the future.  If you have a bug, that might be.
set(gca,'xlim',1.1*[minx,maxx],'ylim',1.1*[miny,maxy])

end
