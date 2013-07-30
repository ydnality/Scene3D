function [bases,sValues,mn] = lmComputeBases(samples,removeMean,nBases)
%
%   [bases, sValues] = lmComputeBases(samples,removeMean,nBases)
%
% Author: FX, BW
%  Compute the basis functions of the samples using Principal Components
%  Analysis (PCA) .  By default, we remove the mean as in conventional
%  principal components calculations.
% 
%  samples --         Each column is a data sample
%  removeMean:        Subtract out the mean  (default is to remove it)
%  nBases:            How many bases to compute
%
%  OUTPUT 
%  bases   -- bases functions
%  sValues -- singular values
%  mn      -- mean data vector
%
%  Maybe this should be called pca?  Maybe we should have a linear models
%  toolbox with routines like lmPCA and lmOneMode and so forth.

if ~exist('removeMean') | isempty(removeMean)
    removeMean = 1;
end

if removeMean
    mn = mean(samples);
%     samples = samples - repmat(mn,[1,size(samples,2)]);
    samples = samples - repmat(mn,[size(samples,1), 1]);
else
    mn = 0;
end

% By calculating the covariance first, we can shrink the size of the matrix used in the svd, 
% allowing some conditions to run that would otherwise be too large.

[bases,sValues,V] = svd(samples*samples',0);
sValues = diag(sValues);

if exist('nBases','var')
    sValues = sValues(1:nBases);
    bases = bases(:,1:nBases);
end

return;
