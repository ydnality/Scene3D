% Compute the weight [ d_nmkl ] the specific peak coeff [a_kl] and the 
% relate Zernike coeff [c_nnm]


function [d_nmkl]=paPeak2Zernike_dnmkl(n,m,k,l,varargin)

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
% Conforti, G. "Zernike aberration coefficients from Seidel and higher-order power-series coefficients." Optics letters 8.7 (1983): 407-408.


%% BETWEEN  Malacara (Mahajan) and Conforti there is a difference in the definition of the Zernike orthonormal basis.
% different term is [2(n+1)(1+deltaK(m=0))]^(1/2)

%% Select method

if nargin>4
    switch varargin{1}
        case {'Mahajan'}
            method='Mahajan';
        case {'Conforti';'Conforti1983'}
            method='Conforti1983';
        case {'Conforti_mod';'Conforti1983 mod';'Conforti1983_mod';'OSA'}
            method='Conforti1983 mod';
        otherwise
            error (['Not valid ',varargin{1},' as method'])            
    end
else
   method='Conforti1983 mod'; %by default
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
        c1= n-m;
        cond1=(c1>=0)& not(mod(c1,2));
        c2= k-l;
        cond2=(c2>=0)& not(mod(c2,2));
        cond2=1;
        c3= l-m;
        cond3=(c3>=0)& not(mod(c3,2));


        if (cond1 && cond2 && cond3)
            % Terms not in the sum
%             if m==l
%                 T1=(2*(n+1)*(1+dK))^1/2;
%             else
                T1=(2*(n+1)*(1+dK))^1/2;
%             end
            T2=(factorial(l))/((2^l)*factorial((l-m)/2)*factorial((l+m)/2));
            T3=T1*T2;

            % Sum term
            lim=(n-m)/2;
            for s=0:lim
                num=((-1)^s)*factorial((n-s));
                den=factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s)*(n-2*s+k+2);
                TS(s+1)=num/den;
            end

            %% OUTPUT
            d_nmkl=T3*sum(TS);
        else
            d_nmkl=0;
        end
    case {'Conforti1983'}
            %% CHECK EXISTENCE
            c1= n-m;
            cond1=(c1>=0)& not(mod(c1,2));
            % c2= k-l;
            % cond2=(c2>=0)& not(mod(c2,2));
%             cond2=1;
            c3= l-m;
            cond3=(c3>=0)& not(mod(c3,2));
%             if (cond1 && cond2 && cond3)
            if (cond1 && cond3)
                %% Terms not in the sum
                % Conforti, G. "Zernike aberration coefficients from Seidel and higher-order power-series coefficients." Optics letters 8.7 (1983): 407-408.
                T3_1=(n+1)*factorial(l)*2^(2-l);
                T3_2=(1+dK)*factorial((l+m)/2)*factorial((l-m)/2);
                T3=T3_1/T3_2;


                %% Sum term
                lim=(n-m)/2;
                for s=0:lim
                    num=(-1)^s*factorial(n-s);
                    den=factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s)*(n-2*s+k+2);
                    TS(s+1)=num/den;
                end

                %% OUTPUT
                d_nmkl=T3*sum(TS);
            else
                d_nmkl=0;
            end
        case {'Conforti1983 mod'}
            % Here is implemented the method described by Conforti but , to
            % have consistent result with Mahajan the weight are normalize
            % to coeff. used by Mahajan (and in standard) for each base
            % [2(n+1)/(1+dK)]^(1/2)
            %% CHECK EXISTENCE
            c1= n-m;
            cond1=(c1>=0)& not(mod(c1,2));
            c2= k-l;
            cond2=(c2>=0)& not(mod(c2,2));            
            c3= l-m;
            cond3=(c3>=0)& not(mod(c3,2));
            if (cond1 && cond2 && cond3)
%             if (cond1 && cond3)
                %% Terms not in the sum
                % Conforti, G. "Zernike aberration coefficients from Seidel and higher-order power-series coefficients." Optics letters 8.7 (1983): 407-408.
                T3_1=(n+1)*factorial(l)*2^(2-l);
                T3_2=(1+dK)*factorial((l+m)/2)*factorial((l-m)/2);
                T3=T3_1/T3_2;


                %% Sum term
                lim=(n-m)/2;
                for s=0:lim
                    num=(-1)^s*factorial(n-s);
                    den=factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s)*(n-2*s+k+2);
                    TS(s+1)=num/den;
                end
                 %% OUTPUT
                d_nmkl=T3*sum(TS);
                %% Modification from "Conforti1983"
                koeff=(2*(n+1)/(1+dK))^(1/2); %[2(n+1)/(1+dK)]^(1/2)
                d_nmkl=d_nmkl/koeff;

               
            else
                d_nmkl=0;
            end
        
end




