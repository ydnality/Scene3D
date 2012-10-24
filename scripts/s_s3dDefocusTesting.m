% Subbarao claims that the alpha in the first term is an error.  Without
% it, however, the OTF doesn't scale to (near) 1.  Even with it, the value
% is a little short of 1.  So ...

% We need s (reduced spatial frequency) and alpha (which is related to w20
% and the defocus parameter)
% This are computed in opticsDefocusCore
ii = 1;
a = alpha(ii,:);
nf = abs(abs(s(ii,:)))/2;   % Normalized spatial frequency from reduced SF
beta = acos(nf); %

H1 = beta .* besselj(1,a) + ...
    sin(2*beta)/2 .* (besselj(1,a) - besselj(3,a)) - ...
    sin(4*beta)/4 .* (besselj(3,a) - besselj(5,a));

H2 = sin(beta).*(besselj(0,a) - besselj(2,a)) + ...
    sin(3*beta)/3 .* (besselj(2,a) - besselj(4,a)) - ...
    sin(5*beta)/5 .* (besselj(4,a) - besselj(6,a));

OTF = (4 ./ (pi*a)) .* ((cos(a .* nf)) .* H1) - (sin(a.*nf) .* H2);

figure; plot(OTF); hold on


%%
%
tl = 0;
OTF = ...
    (4./(pi*a)).*cos(a.*nf).* ...
    (beta.*besselj(1,a)+ ...
    1/2*sin(2*beta).* (besselj(1,a)-besselj(3,a))...
    -1/4*sin(4*beta).*(besselj(3,a)-besselj(5,a))+tl)...
    -(4./(pi*a)).*sin(a.*s(ii)/2).*...
    (sin(beta).*(besselj(0,a)-besselj(2,a))...
    +1/3*sin(3*beta).*(besselj(2,a)-besselj(4,a))...
    -1/5*sin(5*beta).*(besselj(4,a)-besselj(6,a))-tl);

plot(OTF)





