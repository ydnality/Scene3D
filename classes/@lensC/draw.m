function obj =  draw(obj)
% Draw the the multi-element lens in a graph window
%
%   lens.draw
%
% See also:  psfCamera.draw
%  That calls this lens draw and also draws the rays for the point spread,
%  I think.
%
% AL/BW Vistasoft Team, Copyright 2014

% Create the figure and set the parameters
vcNewGraphWin;
axis equal;
lWidth = 2; lColor = 'k';  % Drawing parameters

% We draw one surface/aperture at a time
nSurfaces = obj.get('n surfaces');
for lensEl = 1:nSurfaces
    
    % Get the current surface element
    curEl = obj.surfaceArray(lensEl);
    
    if (curEl.sRadius ~= 0)
        %% Simpler code
        testme = false;
        if testme
            [z,y] = circlePoints([curEl.sCenter,0],curEl.sRadius);
            if curEl.sRadius > 0
                l = (abs(y) < curEl.apertureD/2) & (z < curEl.sCenter(3));
            else
                l = (abs(y) < curEl.apertureD/2) & (z > curEl.sCenter(3));
            end
            p = plot(z(l),y(l),'k.'); set(p,'markersize',1); hold on;
        else
            % Draw a spherical element
            
            % Get previous and next elements, within bounds
            nextEl = obj.surfaceArray(min(lensEl+1, end));
            prevEl = obj.surfaceArray(max(lensEl-1, 1));
            
            % Lens elements do NOT always end when the neighboring
            % element begins.  this allows for a fudge factor.  This
            % won't matter too much because the aperture radius will
            % be the limiting factor. (I don't understand this
            % comment, other than this is a fudge factor somehwere.
            % BW).
            delta = 10;
            
            % Determine the minimum/maximum z-positions for the
            % curve. Which side has the position depends on the
            % sign of the curvature
            if (curEl.sRadius > 0 )
                % The fudge factor is used here
                leftBoundary  = curEl.get('zIntercept');
                rightBoundary = nextEl.get('zIntercept') + delta;
            else
                % Here is fudge again
                leftBoundary  = prevEl.get('zIntercept') - delta;
                rightBoundary = curEl.get('zIntercept');
            end
            
            % This is the range of z values we will consider.
            zPlot = linspace(leftBoundary, rightBoundary, 100);
            
            
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
            % We get the positive and negative y-values
            yPlot  =  sqrt(curEl.sRadius^2 - (zPlot - curEl.sCenter(3)) .^2);
            yPlotN = -sqrt(curEl.sRadius^2 - (zPlot - curEl.sCenter(3)) .^2);
            
            % NOTE: We may have problems with a concave lens
            % because of how we are choosing the z range - but
            % these are rare.
            %
            % We will plot for the range of y values that are less
            % than the radius of the spherical surface.
            withinRange = (yPlot < curEl.apertureD/2);
            
            % The positive solutions
            l = line(zPlot(withinRange), yPlot(withinRange));
            set(l,'linewidth',lWidth,'color',lColor);
            
            % The negative solutions
            l = line(zPlot(withinRange), yPlotN(withinRange));
            set(l,'linewidth',lWidth,'color',lColor);
        end
        
    else
        %Draw the aperture opening if radius = 0
        
        %TODO: draw the difference between specified aperture
        %from file and specified aperture from object instance
        
        %right now: take the minimum value
        curAperture = min(curEl.apertureD/2, obj.apertureMiddleD/2);
        
        l = line(curEl.sCenter(3) * ones(2,1), -1*[curEl.apertureD/2 curAperture]);
        set(l,'linewidth',lWidth,'color',lColor);
        l = line(curEl.sCenter(3) * ones(2,1), [curAperture curEl.apertureD/2]);
        set(l,'linewidth',lWidth,'color',lColor);
        
    end
    
end

end
