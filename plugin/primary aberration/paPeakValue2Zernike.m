% Convert Peak Wavefront Coeffs to Zernike Coeffs 
%                       (with rotational  symmetry)
%Data from "Malacara, Daniel, ed. Optical shop testing. Vol. 59. John Wiley & Sons, 2007." 
% Are defined as peak value:
% .W40;.W31;.W22;.W20;.W11 [Hopkins and Sasian notation]
% or .a40;.a31;.a22;.a20;.a11 [Mahajan notation]


function [ZernikeCoeff]=paPeakValue2Zernike(peakCoeff)

%INPUT
%PeaCoeff: this coeffs has been adapted to real field height [.W40;.W31;.W22;.W20;.W11]

%OUTPUT
%ZernikeCoeff: struct with several term





%% INITIALIZE WEIGHT



ZernikeList=struct('C00',[],'C11',[],'C20',[],'C22',[],'C31',[],...
    'C33',[],'C40',[],'C42',[],'C44',[],'C51',[],'C53',[],'C55',[],...
    'C60',[],'C62',[],'C64',[],'C66',[],'C71',[],'C73',[],'C75',[],...
    'C77',[],'C80',[]);

%% COMPUTE ZERNIKE coeffs

%number of Primary aberration coeffs
% fnames=fieldnames(peakCoeff);
% nfname=length(fnames); 
% % eliminate field of 'wave' and 'unit'
% inF=1;
% nfs=struct;
% for ai=1:nfname
%     switch fnames{ai}
%         case{'unit';'wave'}
%         otherwise
%             value=getfield(nfs,fnames{ai});
%             nfs=setfield(nfs,fnames{ai},value);
%     end
% end
[fnames]=paGetFieldNames(peakCoeff);
 nfname=length(fnames); 

%number of wavelengths
dummy=getfield(peakCoeff,fnames{1});
nw=size(dummy,1);

%num  of Zernike coeff evaluated
Znames=fieldnames(ZernikeList);
nZc=size(Znames,1);

for ai=1:nfname
    peakC=getfield(peakCoeff,fnames{ai}); %get value from wavefront coeff
    kc=str2num(fnames{ai}(2)); %k index
    lc=str2num(fnames{ai}(3)); %l index
    for ci=1:nZc
        nc=str2num(Znames{ci}(2)); %k index
        mc=str2num(Znames{ci}(3)); %l index
        %get specific weight
        dnmkl=paPeak2Zernike_dnmkl(nc,mc,kc,lc,'Conforti1983 mod');
        ZC(ci,ai,:)=dnmkl.*peakC;
    end    
end

%Get output sum all all element for different Zernike 
ZernC(:,:)=sum(ZC,2);


%% SET OUTPUT
ZernikeCoeff=struct;
for ci=1:nZc
    value=ZernC(ci,:)';
   ZernikeCoeff=setfield(ZernikeCoeff,Znames{ci},value); %Set the specific Zernike Coeff
end
ZernikeCoeff.wave=peakCoeff.wave;
ZernikeCoeff.unit=peakCoeff.unit;


% for li=1:nw %for each wavelength
%     for ci=1:nfname  %for each wavef. coeffs   
%          wCoeff=getfield(peakCoeff,fnames{ci}); %get value from wavefront coeff
%          k_coeff=str2num(fnames{ci}(2)); %k index
%         l_coeff=str2num(fnames{ci}(3)); %l index
%         
%         for wi=1:ncoeff
%             n_coeff=Z_nm(1,wi);m_coeff=Z_nm(2,wi);
%             weight=paPeak2Zernike_dnmkl(n_coeff,m_coeff,k_coeff,l_coeff,'Conforti1983 mod');
%             Zc(:,wi)=weight.*wCoeff(li,1);
%         end
% %         
% %         for wi=1:ncoeff %Find weight for zernike coeffs for the given wavef. coeff
% %             if (k_coeff==a_kl(1,wi))&&(l_coeff==a_kl(2,wi))
% %                 weight=wM(:,wi);
% %                 weight=paPeak2Zernike_dnmkl(n,m,k,l,'Conforti1983 mod')
% %                 Zc(:,wi)=weight.*wCoeff(li,1);
% %             end
% %         end        
%     end  
%     Zcoeff(:,li)=sum(Zc,2);
% end
% 
% 
% %% CREATE  OUTPUT and APPEND VALUEs to ZERNIKE STRUCTURE
% 
% 
% Zfnames=fieldnames(ZernikeCoeff);
% for zi=1:nZc    
%     ZernikeCoeff=setfield(ZernikeCoeff,Zfnames{zi},Zcoeff(zi,:)'); %Set the specific Zernike Coeff
% end
% 
% 
% 
