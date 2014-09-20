function [n]=sellmeier_func(w,n_ord,B,C)

% Sellmeier equation {http://en.wikipedia.org/wiki/Sellmeier_equation}

%INPUT
%w: (1xP) P different wavelength [unit]
%n_ord: order of the Sellm. Eq (typical num ord=3)
%B: (1 x n_ord) numerator coefficients 
%C: (1 x n_ord) denominator coefficients

%OUTPUT
%n: (1xP) refractive index at the P different wavelength [adim]

%NOTE: Coeff.s have to be consistent with the wavelength ( typically in um)

if (n_ord==length(B)) && (length(C)==length(B))
    n_2=1;
    for i1=1:n_ord
        n_2=n_2+(B(i1).*w.^2./(w.^2-C(i1))); %Iterative approach
    end
    n=sqrt(n_2); %index
else
    n=[];
end

