function [ x0, y0, z0, r ] = fitConcentricSpheres( xCoord, yCoord, zCoord, wght )


x = ones(4,1);

for i=1:1000
    grad = zeros(4,1);
    J = zeros(4);
   
    d = sqrt((x(1) - xCoord).^2 + (x(2) - yCoord).^2 + (x(3)-zCoord).^2);
    
    grad(1) = 2*sum((d-wght*x(4)).*(1./d).*(x(1)-xCoord));
    grad(2) = 2*sum((d-wght*x(4)).*(1./d).*(x(2)-yCoord));
    grad(3) = 2*sum((d-wght*x(4)).*(1./d).*(x(3)-zCoord));
    grad(4) = -2*sum(wght.*(d - wght.*x(4)));
    
    J(1,1) = 2*sum(1 - (wght.*x(4).*d - (wght*x(4).*((x(1) - xCoord).^2))./d)./(d.^2));
    J(1,2) = 2*sum((wght*x(4).*(x(1)-xCoord).*(x(2)-yCoord))./(d.^(3)));
    J(1,3) = 2*sum((wght*x(4).*(x(1)-xCoord).*(x(3)-zCoord))./(d.^(3)));
    J(1,4) = -2*sum(wght.*(x(1)-xCoord)./d);
        
    J(2,1) = 2*sum((wght*x(4).*(x(1)-xCoord).*(x(2)-yCoord))./(d.^(3)));
    J(2,2) = 2*sum(1 - (wght.*x(4).*d - (wght*x(4).*((x(2) - yCoord).^2))./d)./(d.^2));
    J(2,3) = 2*sum((wght*x(4).*(x(2)-yCoord).*(x(3)-zCoord))./(d.^(3)));
    J(2,4) = -2*sum(wght.*(x(2)-yCoord)./d);
    
    J(3,1) = 2*sum((wght*x(4).*(x(1)-xCoord).*(x(3)-zCoord))./(d.^(3)));
    J(3,2) = 2*sum((wght*x(4).*(x(3)-zCoord).*(x(2)-yCoord))./(d.^(3)));
    J(3,3) = 2*sum(1 - (wght.*x(4).*d - (wght*x(4).*((x(3) - zCoord).^2))./d)./(d.^2));
    J(3,4) = -2*sum(wght.*(x(3)-zCoord)./d);

    J(4,1) = -2*sum(wght.*(x(1)-xCoord)./d);
    J(4,2) = -2*sum(wght.*(x(2)-yCoord)./d);
    J(4,3) = -2*sum(wght.*(x(3)-zCoord)./d);
    J(4,4) = 2*sum(wght.^2);
    
    x = x - (J\grad);
    % x = x - gamma*grad;
    
    fprintf('Iter %i (%f,%f,%f,%f), gradient %f\n',i,x(1),x(2),x(3),x(4),norm(grad));
    if norm(grad) < eps, break; end
    
end

x0 = x(1);
y0 = x(2);
z0 = x(3);
r = x(4);

end

