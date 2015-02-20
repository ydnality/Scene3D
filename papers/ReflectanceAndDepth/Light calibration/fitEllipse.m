function [ Cmat, aEst, bEst, xcEst, ycEst, tauEst ] = fitEllipse( x,y )

D = [x.*x x.*y y.*y x y ones(length(x),1)];
S = D'*D;
C = zeros(6,6);
C(1,3) = -2;
C(2,2) = 1;
C(3,1) = -2;
[gevec, geval] = eig(S,C);
[~, negC] = find(geval < 0 & ~isinf(geval));
Cmat = gevec(:,negC);

A = Cmat(1);
B = Cmat(2);
C = Cmat(3);
D = Cmat(4);
E = Cmat(5);
F = Cmat(6);

M0 = [F D/2 E/2; D/2 A B/2; E/2 B/2 C];
M = [A B/2; B/2 C];
eVal = sort(eig(M));

aEst = sqrt(-det(M0)/(det(M)*eVal(1)));
bEst = sqrt(-det(M0)/(det(M)*eVal(2)));
xcEst = (B*E - 2*C*D)/(4*A*C - B^2);
ycEst = (B*D - 2*A*E)/(4*A*C - B^2);
tauEst = acotd((A-C)/B)/2;



end

