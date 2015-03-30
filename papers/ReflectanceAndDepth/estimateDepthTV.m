function [ est, initEst, hist ] = estimateDepthTV( A, B, C, alpha, varargin )

% [ est, hist ] = estimateDepthTV( A, B, C, alpha, varargin )
%
% This function accepts matrices A,B,C that are used to compute an initial
% estimate of the depth map. The pixel depth d is obtained by solving
% 
%   a*d^2 + b*d + c =0
%
% The initial estimate is then refined by solving with ADMM the following TV
% penalized optimization problem
%
% minimize |d-dref| + alha*TV(d')
% subject to d>=0
%
% At this point the TV constriant is applied to the anisotropic gradient of
% d in the x,y spatial dimensions.
%
% Inputs:
%   A,B,C - matrices with 2nd order polynominal coefficients. They have to
%       have the same size, which defines the spatial arrangement of
%       pixels.
%   alpha - the TV penalty strength
%
% Outputs:
%   est - a matrix with the smoothed depth map estimate.
%   initEst - a matrix with the initial depth map estimate.
%   hist - parameters and convergence data of the ADMM solver.
%
% Copyright, Henryk Blasinski 2015

p = inputParser;

p.addParamValue('maxIter',1000,@isscalar);  % Maximum number of iterations
p.addParamValue('tol',1e-6,@isscalar);      % Per pixel sub-optimiality 
p.addParamValue('verbose',false);           % Display messages

% These are some ADMM tuning parameters, not really necessary to change
% these.
p.addParamValue('mu',10,@isscalar);
p.addParamValue('tauIncr',2,@isscalar);
p.addParamValue('tauDecr',2,@isscalar);
p.addParamValue('rhoInit',0.1,@isscalar);
p.parse(varargin{:});
inputs = p.Results;


% First compute the depth by solving the quadratic
delta = sqrt(B.^2 - 4.*A.*C);
d1 = (-B - delta)./(2*A);
d2 = (-B + delta)./(2*A);
dEst = max(d1,d2);

initEst = 1./dEst;
initEst(isinf(initEst)) = 0;

% And now the ADMM

[h, w ] = size(A);

% Initialize temporary and dual variables
Y1 = zeros(h,w-1);
Y2 = zeros(h-1,w);
Y3 = zeros(h,w);

Y1minus = Y1;
Y2minus = Y2;
Y3minus = Y3;

U1 = zeros(size(Y1));
U2 = zeros(size(Y2));
U3 = zeros(size(Y3));

est = zeros(h,w);

% Initialize some book keeping variables
rho = inputs.rhoInit*ones(inputs.maxIter+1,1);
prRes = zeros(inputs.maxIter+1,1);
dualRes = zeros(inputs.maxIter+1,1);
pcgIter = zeros(inputs.maxIter+1,1);

for i =1:inputs.maxIter
    
    % Solve a least-squares problem using conjugate gradient descent.
    Atb = applyAt(initEst,Y1-U1,Y2-U2,Y3-U3,h,w,rho(i)/2);
    AtAhndl = @(x) applyAtA(x,h,w,rho(i)/2);
    
    tic
    [est, ~, relres, pcgIter(i)] = pcg(AtAhndl,Atb,1e-6,10000,[],[],est(:));
    t1 = toc;
    
    
    % Y variable update
    est = reshape(est,[h,w]);
    
    tmpX = imfilter(est,[1 -1 0]);    
    tmpY = imfilter(est,[1 -1 0]');
   
    tmpX = tmpX(:,2:w,:);
    tmpY = tmpY(2:h,:,:);

    
    % Soft thresholding
    Y1 = max(U1+tmpX - alpha/rho(i),0) - max(-(U1+tmpX) - alpha/rho(i),0);
    Y2 = max(U2+tmpY - alpha/rho(i),0) - max(-(U2+tmpY) - alpha/rho(i),0);
    % Nonegativity thresholding
    Y3 = est; Y3(Y3<0) = 0;
    
    % Dual variable update
    U1 = U1 + tmpX - Y1;
    U2 = U2 + tmpY - Y2;
    U3 = U3 + est - Y3;
    
    % Compute residuals
    prRes(i) = sqrt(norm(tmpX-Y1,'fro')^2 + norm(tmpY-Y2,'fro')^2 + norm(est-Y3,'fro')^2); 
    dualRes(i) = sqrt(norm(Y1-Y1minus,'fro')^2 + norm(Y2-Y2minus,'fro')^2 + norm(Y3-Y3minus,'fro')^2); 
    if max(prRes(i)/(h*w),dualRes(i)/(h*w)) < inputs.tol, break; end;
    
    if inputs.verbose
        fprintf('Iter %i, pr res %f, dual res %f\n',i,prRes(i),dualRes(i));
        fprintf('     -> PCG (%f), err %f, nIter %i\n',t1,relres,pcgIter(i));
    end
    
    % Re-scale the parameter rho
    if prRes(i) > inputs.mu*dualRes(i),
        rho(i+1) = rho(i)*inputs.tauIncr;
    end;
    if dualRes(i) > inputs.mu*prRes(i),
        rho(i+1) = rho(i)/inputs.tauDecr;
    end;

    % Now we need to re-scale the scaled dual variable U as well as the
    % dual residuals
    U1 = U1*rho(i)/rho(i+1);
    U2 = U2*rho(i)/rho(i+1);
    U3 = U3*rho(i)/rho(i+1);

    Y1minus = Y1;
    Y2minus = Y2;
    Y3minus = Y3;
    
end

hist.rho = rho(1:i);
hist.prRes = prRes(1:i);
hist.dualRes = dualRes(1:i);
hist.pcgIter = pcgIter(1:i);

end



function res = applyAtA( meas, h , w, beta)

% meas will be passed as a vector
meas = reshape(meas,h,w);

% Derivatives in the x direction

% Fast
firstLastColumn = imfilter(meas(:,[1 2 w-1 w],:),[1 -1 0],'same');
interior = imfilter(meas,[-1 2 -1],'same');
interior(:,1) = firstLastColumn(:,2);
interior(:,w) = -firstLastColumn(:,4);

dx2ImgMat = interior;

% derivatives in the y direction

% Fast
firstLastRow = imfilter(meas([1 2 h-1 h],:,:), [1 -1 0]','same');
interior = imfilter(meas,[-1 2 -1]','same');
interior(1,:) = firstLastRow(2,:);
interior(h,:) = -firstLastRow(4,:);

dy2ImgMat = interior;

res = (1+beta)*meas + beta*(dx2ImgMat + dy2ImgMat);
res = res(:);

end

function res = applyAt( meas, smX, smY, nonneg, h, w, beta)


% Derivatives in the x direction
dx2Img = imfilter(smX,[-1 1 0]);
dx2Img = dx2Img(:,2:w-1);
dx2Img = cat(2,smX(:,1),dx2Img,-smX(:,w-1));

% derivatives in the y direction
dy2Img = imfilter(smY,[-1 1 0]');
dy2Img = dy2Img(2:h-1,:);
dy2Img = cat(1,smY(1,:),dy2Img,-smY(h-1,:));


s34 = beta*(dy2Img + dx2Img);

res = meas + beta*nonneg + s34;
res = res(:);

end
