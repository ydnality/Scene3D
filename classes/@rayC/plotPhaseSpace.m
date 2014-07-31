function obj = plotPhaseSpace(obj)
% Plots the phase-space representation of the rays

    angles = toProjectedAngles(obj);
    
    %% plot the projection onto the x-z plane first
    position = obj.origin(:,1);
    theta_x = angles(:,1);
    
    vcNewGraphWin;
    plot(position, theta_x, '.');  
    xlabel('x');
    ylabel('theta_x');
    
    %% plot the projection onto the x-z plane first
    position = obj.origin(:,2);
    theta_y = angles(:,2);
    
    vcNewGraphWin;
    plot(position, theta_y, '.');  
    xlabel('y');
    ylabel('theta_y');    
    
end
