function d = ptDistances(origin,ePoints)
% Compute the distance from the origin to the endpoints
%
% AL/BW Vistasoft Team, 2015

% The origin is forced to a row vector, and we assume that the end
% points are in the rows of the matrix ePoints

% We could check that length(origin) == size(ePoints,2)
%
delta = bsxfun(@minus,ePoints,origin(:)');

% The two-norm of each column
d = sqrt(sum(abs(delta').^2,1))'; 

end