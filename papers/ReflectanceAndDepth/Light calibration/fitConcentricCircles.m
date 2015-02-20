function [ x0, y0, r ] = fitConcentricCircles( xCoord, yCoord, label )

nCircles = length(unique(label));

x = zeros(2+nCircles,1);
gamma = 0.0001;

for i=1:10000
    J = [];
    for n=1:nCircles

        m = length(label == n);
        d = sqrt((x(1) - xCoord(label==n)).^2 + (x(2) - yCoord(label == n)).^2);

        rCurr = x(2+n);
        xCurr = x(1);
        yCurr = x(2);
        
        da = rCurr*(- xCurr + xCoord(label == n))./d - mean(xCoord(label == n)) + xCurr;
        db = rCurr*(- yCurr + yCoord(label == n))./d - mean(yCoord(label == n)) + yCurr;
        
        Jcircle = [da db zeros(m,n-1) -d+rCurr zeros(m,nCircles-n)];
        
        J = [J; Jcircle];
    end

    grad = sum(J)';
    x = x - gamma*grad;

    if norm(grad) < eps, break; end
    
end

norm(grad)

x0 = x(1);
y0 = x(2);
r = x(3:end);

end

