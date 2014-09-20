function [n]=cauchy_func(w,Num,n_ord)

% Sellmeier equation {http://en.wikipedia.org/wiki/Cauchy's_equation}

%INPUT
%w: (1xP) P different wavelength [nm]
%n_ord: order of the Cauchy. Eq terms 
%Num: (1 x n_ord) numerator coefficients [ for wavelength in nm]

%OUTPUT
%n: (1xP) refractive index at the P different wavelength [adim]

%NOTE: Coeff.s has expressed for lambda in [nm] 
if nargin<3
    n_ord=length(Num); %use all terms for compute the refractive index 
elseif length(Num)==length(n_ord);
    Num=Num(1:n_ord); %conside only the first n_ord    
end




for ni=1:n_ord    
    if ni==1
        n=Num(ni);% first term
    else
        exp_ord=2*(ni-1); %order of exponent
        w_exp_ord=w.^exp_ord; %lambda .^(2(ni-1))
        n=n+Num(ni)./w_exp_ord; %add to the previuos term
    end
end


