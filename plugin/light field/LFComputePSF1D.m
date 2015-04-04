% Computer 1D PSF  along the directional-axis of a 2D Light Field

function [PSF1D]=LFComputePSF1D(LF,method,vectY)


%INPUT
%LFobj: it an object describing a light light field 
%method: specify the method of computation for the PSF
%vectY: sampling vector for the PSF


%OUTPUT
%PSF1D: is the point spread funtion along the axis of the given LF object 
%       

error('Deprecated')

end


% %% CHECK whether the Y vector is sampled uniformally
dy=vectY(2:end)-vectY(1:end-1);

if not(all(dy==dy(1)))
    %Resampling
    warning (' Y vector for compute PSF profile is RESAMPLED with same sample number')
    num_elem=numel(vectY);
    deltaY=(max(vectY)-min(vectY))/(num_elem-1);
    vectY=min(vectY):deltaY:max(vectY);
%     dY=[vectY-deltaY/2,vectY(end)+deltaY/2]; % vector for histogram
    dY=[vectY-deltaY/2,vectY(end)+deltaY/2]; % vector for histogram
else
    dY=[vectY-dy./2,vectY(end)+dy(end)./2]; % so length(dY)=length(vectY)  
    deltaY=abs(dy(2)-dy(1));
end




%% SELECT COMPUTATIONAL METHOD

switch method
    
    case {'standard';'default';'no intensity'}
        %All rays with all possible angle of incidence are considered with the same efficiency
        PSF1D.method='standard';
        for li=1:size(LF.wave,1) 
            vectorY=reshape(LF.Y(:,:,li),1,[]); 
            [vectorY,isY]=sort(vectorY);
            %matrix of intensity
            vectorI=reshape(LF.intensity(:,:,li),1,[]);
            matrixI(li,:)=vectorI(isY);
            %count number of rays reaching the same position
            vettC(li,:)=histc(vectorY,dY);
            %normalize the number of rays for the total number
            Cnorm=vettC(li,:)./(deltaY);
            vettCnorm(li,:)=Cnorm./(sum(Cnorm));
            
        end
       %Histogram
       PSF1D.hist.Y=dY;PSF1D.hist.I=vettC; %not normalized
       PSF1D.hist.normI=vettCnorm;
       % Histogram derived
%        PSF1D.raw.Y=vectY;
%        PSF1D.raw.normI=vettCnorm(:,1:end-1); 
       
       %More reliable
       PSF1D.raw.Y=vectorY;
       PSF1D.raw.normI=matrixI;
       
    case {'cosine law','Lambertian cosine law'} 
        %All rays are weighted with the efficensy give by the  Lambert's cosine law , Rad_teta=Rad * cos(teta)
        PSF1D.method='cosine law';
        histI=zeros(size(LF.wave,1),size(dY,2));
        for li=1:size(LF.wave,1) 
            % Take value for specific wavelength
            Ulam=LF.U(:,:,li);
            Ylam=LF.Y(:,:,li);
            Ilam=LF.intensity(:,:,li);
            Wlam=cos(Ulam);% weight vector according to Lambert's cosine law
                        
            % Arrange in a vector
            vecUlab=reshape(Ulam,1,[]);
            vecYlam=reshape(Ylam,1,[]);
            vecIlam=reshape(Ilam,1,[]);
            vecWlam=reshape(Wlam,1,[]);
            %order along Y
            [vecYlam,indYlam]=sort(vecYlam);
            vecUlab=vecUlab(indYlam);
            vecIlam=vecIlam(indYlam);
            vecWlam=vecWlam(indYlam);
            %Overall intensity
            vecWI(li,:)=vecIlam.*vecWlam;
            %normalized
            vecWI_norm(li,:)=vecWI(li,:)/sum(vecWI(li,:));
            
            %DEBUG
%             figure
%             subplot (1,2,1)
%             scatter(vecYlam,vecUlab)
%             hold all
%             stem(dY,max(vecUlab)*ones(1,length(dY)),'--*r')
            
            
            %Compute histogram weighted by radiance intensity
            for di=2:length(dY)
                indY=find(vecYlam>=dY(di-1) & vecYlam<dY(di));
%                 histI(li,di-1)=sum(vecIlam(indY));
                histI(li,di-1)=sum(vecWlam(indY).*vecIlam(indY));
                if di==length(dY) %to be similar to histc
                    indY0=find(vecYlam==dY(di));
                    histI(li,di)=sum(vecWlam(indY0).*vecIlam(indY0));
                end
            end
            
%             %DEBUG
%             subplot(1,2,2)
%             bar(dY,histI)


        end
        
       %Histogram
       PSF1D.hist.Y=dY;PSF1D.hist.I=histI; %not normalized
%        PSF1D.hist.normI=vettCnorm;
       
       PSF1D.raw.Y=vecYlam;
       PSF1D.raw.normI= vecWI_norm; 
       
    otherwise
        warning('Not Valid Method of Computation for PSF1D!')
        PSF1D=[];
end


%% FIX the UN-MATCHING between Ray Angle and Film sampling

for li=1:size(LF.wave,1) 
    % ADD EXTRA ZEROS to force interpolation to zero out of known edge
    %Get known values
    k_Y=PSF1D.raw.Y;
    k_I=PSF1D.raw.normI(li,:);
    if (min(dY)<min(k_Y))||(max(dY)>max(k_Y))
        deltaY=abs(diff(dY)); %sampling vector for dY
        delta_k_Y=abs(diff(k_Y)); %sampling vector for k_Y
        dymin=min([min(deltaY),min(delta_k_Y)]); %select smaller sampling step
        if dymin<=0 %to avoid numeric 'error' that can give dymin=0 or negative
            
%             ind_k_Y=find(delta_k_Y<1e-5);
%             k_Y1=k_Y(ind_k_Y);
%             k_I1=k_I(ind_k_Y);
              [k_Y1,r1,c1]=unique(k_Y);
              k_I1=k_I(r1);
              k_Y=k_Y1;k_I=k_I1;
              
              % add 10 zero values (out of the known edge) sampled with dymin/10 to the known data
               dymin=mean(diff(k_Y));
        end
            
        if not((dymin>0) && (dymin<Inf))  %include the case of NaN
             dymin=LF.wave/100; %limit set to 100 times smaller that the current wavelent
        end
       % add 10 zero values (out of the known edge) sampled with dymin/10 to the known data
        Nsamples=10; % 10 samples are selected beacuse the intepolation algorithms need less samples
        dymin=dymin/10;
        vector_dymin=dymin:dymin:dymin*Nsamples;
        k_Y=[k_Y(1)-vector_dymin(end:-1:1),k_Y,k_Y(end)+vector_dymin];
        k_I=[zeros(1,Nsamples),k_I,zeros(1,Nsamples)];
    
        
    end
    
   %interpolation
    normI_res(li,:)=interp1(k_Y,k_I,dY);
    % check NaN value 
    I_nan=isnan(normI_res(li,:));
    ind_I_nan=find(I_nan);
    normI_res(li,ind_I_nan)=0; %subj NaN with Zero
    if all(I_nan)
        % the peak of the PSF has been missed, it reintroduced by the
        % mean value of the positive raw data
        ind_I_pos=find(k_I>0);
        if isempty(ind_I_pos)
            
        else
            [peak_value]=mean(k_I(ind_I_pos));
            [peak_pos]=mean(k_Y(ind_I_pos));
            %add to the exit
            dY(end+1)=peak_pos;
            normI_res(end+1)=peak_value;
            %sort exit
            [dY,ind_sortY]=sort(dY);
            normI_res(li,:)= normI_res(li,ind_sortY);
            warning('A value is added to the PSF to show the peak. NOT UNIFORM SAMPLING ANYMORE!!')
        end
        
        
        
    end
    %normalized
    normI_res_norm(li,:)=normI_res(li,:)/sum(normI_res(li,:));
end
PSF1D.profile.Y=dY;
PSF1D.profile.normI=normI_res_norm;

%DEBUG
% figure
% stem(PSF1D.raw.Y,PSF1D.raw.normI(1,:))
% hold all
% stem(PSF1D.profile.Y,PSF1D.profile.normI(1,:))



% for li=1:size(LF.wave,1) 
%     inds=find(PSF1D.raw.normI(li,:)>0); %indices different to zero
% %     profCnom=vettCnorm(li,1:end-1);
%     profCnom=PSF1D.raw.normI(li,:);
%     vettY_res=PSF1D.raw.Y(inds); %resampling vector
%     if length(inds)>1
%         vett_nI_res=interp1(vettY_res,profCnom(inds),[PSF1D.raw.Y([inds(1):inds(end)])]); %resampled profile
%         vett_normI_res(li,:)=[profCnom(1:inds(1)-1),vett_nI_res,profCnom(inds(end)+1:end)]; % complete profile   
%     else
%         vett_normI_res(li,:)=PSF1D.raw.normI(li,:); % NO CHANGE
%     end
%     % PSF can not have negative value, if present for interpolation move to
%     % zero
%     vett_normI_res(li,find(vett_normI_res(li,:)<0))=0;    
% end
% PSF1D.profile.Y=vectY;
% PSF1D.profile.normI=vett_normI_res;



