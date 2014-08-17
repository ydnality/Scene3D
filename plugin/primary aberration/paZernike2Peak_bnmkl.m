% Compute the weight [ b_nmkl ] for the specific Zernike coeff [c_nnm]  and the 
% relate peak coeff [a_kl]


function [b_nmkl]=paZernike2Peak_bnmkl(n,m,k,l,varargin)

%INPUT
%n: Zernike coeff- radius order
%m: Zernike coeff- cosine ordeer
%k: Peak Value coeff- radius order
%l: Peak coeff- cosine order
%varargin: specify the method implemented

%OUTPUT
%d_nmkl: weight     

% See formula (13.39)-Data from "Malacara, Daniel, ed. Optical shop testing. Vol. 59. John Wiley & Sons, 2007." 
% Alternative: 
% Tyson, Robert K. "Conversion of Zernike aberration coefficients to Seidel and higher-order power-series aberration coefficients." optics letters 7.6 (1982): 262-264.

%% BETWEEN  Malacara (Mahajan) and Conforti there is a difference in the definition of the Zernike orthonormal basis.
% different term is [2(n+1)(1+deltaK(m=0))]^(1/2)

%% Select method

if nargin>4
    switch varargin{1}
        case {'Mahajan'}
            method='Mahajan';
        case {'Tyson';'Tyson1982'}
            method='Tyson1982';
        otherwise
            error (['Not valid ',varargin{1},' as method'])            
    end
else
   method='Mahajan'; %by default
end



%% IMPLEMENT METHOD

% Delta Kronacher
        if m==0
            dK=1;
        else
            dK=0;
        end

switch method
    case {'Mahajan'}
        %% CHECK EXISTENCE
        c1= k-l;
        cond1=(c1>=0)& not(mod(c1,2));
        c2= n-m;
        cond2=(c2>=0)& not(mod(c2,2));
        c3= n-k;
        cond3=(c3>=0)& not(mod(c3,2));
        c4= k-m;
        cond4=(c4>=0)& not(mod(c4,2));
        c5= m-l;
        cond5=(c5>=0)& not(mod(c5,2));

        if (cond1 && cond2 && cond3 && cond4 && cond5)
            if (m==0)&&(n==0)&&(k==0)&&(l==0)
                b_nmkl=1;
            elseif m==l
                % Terms 
                T1=sqrt(2*(n+1)*(1+dK));
                T2num=(-1)^((n-k)/2)*factorial((n+k)/2)*2^(m-1);
                T2den=factorial((n-k)/2)*factorial((n+k)/2)*factorial ((k-m)/2);
                % output
                b_nmkl=T1*T2num/T2den;
            else
                % Terms 
                T1=sqrt(2*(n+1)/(1+dK));
                T2num=(-1)^((n-k)/2)*factorial((n+k)/2)*m*2^l*(-1)^((m-l)/2)*factorial((m+l)/2-1);
                T2den=factorial((n-k)/2)*factorial((k+m)/2)*factorial((k-m)/2)*factorial(l)*(m-l)*factorial((m-l)/2-1);
                % output
                b_nmkl=T1*T2num/T2den;
            end
        else
            b_nmkl=0;
        end
    case {'Tyson1982'}
            %% CHECK EXISTENCE
        c1= k-l;
        cond1=(c1>=0)& not(mod(c1,2));
        c2= n-m;
        cond2=(c2>=0)& not(mod(c2,2));
        c3= n-k;
        cond3=(c3>=0)& not(mod(c3,2));
        c4= k-m;
        cond4=(c4>=0)& not(mod(c4,2));
        c5= m-l;
        cond5=(c5>=0)& not(mod(c5,2));
        c6=[l,m,k,n];
        cond6=all(diff(c6)>=0);
        

        if (cond1 && cond2 && cond3 && cond4 && cond5 && cond6)
            if (m==0)&&(n==0)&&(k==0)&&(l==0)
                b_nmkl=1;
            elseif m==l
                % Terms 
                T1=(1+dK)*2^(l-1);
                T2num=(-1)^((n-k)/2)*factorial((n+k)/2);
                T2den=sqrt(1+dK)*factorial((n-k)/2)*factorial((k+m)/2)*factorial((k-m)/2);
                % output
                b_nmkl=T1*T2num/T2den;
            else
                % Terms 
                T1=1/(sqrt(1+dK));
                T2num=(-1)^((n-k)/2)*factorial((n+k)/2)*((-1)^((m-l)/2))*2^l*m*factorial((m+l)/2-1);
                T2den=factorial((n-k)/2)*factorial((k+m)/2)*factorial((k-m)/2)*factorial(m-l)*factorial((m-l)/2-1)*factorial(l);
                b_nmkl=T1*T2num/T2den;
            end
        else
            b_nmkl=0;
        end
        
end