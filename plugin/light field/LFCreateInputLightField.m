% Create the INPUT light field due to an object point
% for a given Imaging system.
%There are several method to compute the light field defined as [Y;U] (eccentricity, direction)


function [LFinput,varargout]=LFCreateInputLightField(LFobject,ImagSyst,LFunit,type_sampl,v_sampl)

%INPUT
%LFobject: struct describing the source of light field
%           -point   v_sampl= [dU] or [#U]
%           -wave
%           -flat, if length(varargin {1})=1-> varargin{1}=[#U];if length(varargin {1})>1-> varargin{1}=[#Y,#U]
%ImagSyst: imaging systen
%LFunit: specify if the LF is expressed in paraxial or real unit,('paraxial','real')
%type_sampl: type of sampling of the Light Field input    
%       -> spatialuniform  (varargin{1}:=[dY;dU])
%       -> range uniform    (varargin{1}:=[#Y;#U])
%v_sampl: vector for sampling, can describe sampling step[dY,dU], or
%       desired number of samples in all the range [#Y,#U]
%       if length(v_sampl)==1; v_sampl=[dU] or [#U]

%OUTPUT
%LFinput: input light field for the imaging system



%% CREATE LIGHT FIELD INPUT SIGNAL

% Append type of unit for LIGHT FIELD
LFinput.LFunit=LFunit;
LFinput.wave=ImagSyst.wave;
switch LFobject.type
    case {'point'}
        
        
        %Create intermediate object
        [Obj]=paraxCreateObject(LFobject.z_pos,LFobject.y_ecc,'point',LFobject.unit);
        [ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,Obj);
        LFinput.ISobj{1}=ImagSyst.object{end}; %append the generated object in the imaging system
        LFinput.z_pos=ImagSyst.object{end}.z_pos; %set position along optical axis
        
        %FIND direction limits
        Yup=repmat(LFobject.y_ecc,size(ImagSyst.wave,1),1);
        Ylo=repmat(LFobject.y_ecc,size(ImagSyst.wave,1),1);
        switch LFunit
            case {'paraxial';'parax'}
                Uup=ImagSyst.object{end}.meridionalPlane.comaRay.upper.angleParax;
                Ulo=ImagSyst.object{end}.meridionalPlane.comaRay.lower.angleParax;                
            case {'real'}
                Uup=ImagSyst.object{end}.meridionalPlane.comaRay.upper.angle;
                Ulo=ImagSyst.object{end}.meridionalPlane.comaRay.lower.angle;
            otherwise 
                error ('Light Field unit is not valid! Specify between paraxial and real')
        end
        % SAMPLING and create the the INPUT
        switch type_sampl
            case {'spatialuniform';'spatunif';'spatial'}
                %uniform sampling for the given values
                if length(v_sampl)>1
                    dU=v_sampl(2); %delta U
                    warning('Selected 2nd element for sampling of U!')
                else
                    dU=v_sampl(1); %delta U
                end
                
                                                
            case {'range';'range unif';'range uniform'}
                %uniform sampling on the variable range
                if length(v_sampl)>1
                    nU=v_sampl(2); % #sample of U
                    warning('Selected 2nd element for sampling of U!')
                else
                     nU=v_sampl(1); %delta U
                end
                
                dU=(Uup-Ulo)./(nU-1); %delta U
                               
            otherwise
                error ('NOT valid type of LF input sampling')                
        
        end
        % CREATE VECTOR   and MESH
        
         LFinput.vettY=Ylo; %only a value        
        for li=1:size(ImagSyst.wave,1)
            LFinput.vettU(li,:)=[Ulo(li):dU(li):Uup(li)];
            [LFinput.U(:,:,li),LFinput.Y(:,:,li)]=meshgrid(LFinput.vettU(li,:),LFinput.vettY(li,:));
        end                          
             
       
    
   
        
    case {'wave'}
        
        [Obj]=paraxCreateObject(Inf,0,'point',LFobject.unit);
        [ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,Obj);
        LFobject.ISobj{1}=ImagSyst.object{end}; %append the generated object in the imaging system
        LFinput.z_pos=Inf; %set position along optical axis
        %FIND direction limits
        Uup=repmat(LFobject.angle_ecc,size(ImagSyst.wave,1),1);
        Ulo=repmat(LFobject.angle_ecc,size(ImagSyst.wave,1),1);
        switch LFunit
            case {'paraxial';'parax'}
                Yup=ImagSyst.object{end}.Radiance.EnP.diam(:,1);
                Ylo=ImagSyst.object{end}.Radiance.EnP.diam(:,2);                
            case {'real'}
                Yup=ImagSyst.object{end}.Radiance.EnP.diam(:,1);
                Ylo=ImagSyst.object{end}.Radiance.EnP.diam(:,2);
            otherwise 
                error ('Light Field unit is not valid! Specify between paraxial and real')
        end
        % SAMPLING and create the the INPUT
        switch type_sampl
            case {'spatialuniform';'spatunif';'spatial'}
                %uniform sampling for the given values
                
                dY=v_sampl(1); %delta Y
                                                
            case {'range';'range unif';'range uniform'}
                %uniform sampling on the variable range
                
                nY=v_sampl(1); % #sample of Y
                dY=(Yup-Ylo).\(nY-1); %delta Y
                               
            otherwise
                error ('NOT valid type of LF input sampling')                
        
        end
        % CREATE VECTOR        
        LFinput.vettU=Ulo; %only a value
        
        for li=1:size(ImagSyst.wave,1)
            LFinput.vettY(li,:)=[Ylo(li):dY:Yup(li)];
        end
       
    case {'flat source','flat'}
        %Append Z position
        LFinput.z_pos=LFobject.z_pos; %set position along optical axis 
        
        %Evaluate the need of resampling
        if length(v_sampl)>1
            % Resampling Y axis 
            Yup=max(LFobject.y_ecc); Ylo=min(LFobject.y_ecc);
            switch type_sampl
                case {'spatialuniform';'spatunif';'spatial'}
                    % Along Y
                    dY=v_sampl(1);%delta Y
                    vettY=[Ylo:dY:Yup]; % Y profile resampled 
                    % Along  U
                    dU=v_sampl(2); %delta U
                    
                case {'range';'range unif';'range uniform'}
                    %Along Y
                    nY=v_sampl(1); % #samples along Y
                    dY=(Yup-Ylo)./(nY-1); %delta Y
                    vettY=[Ylo:dY:Yup]; % Y profile resampled 
                    
                    %Along U
                    nU=v_sampl(2);% #samples along Y
                    
                otherwise
                    error ('NOT valid type of LF input sampling')                
            end
            
        else
            vettY=LFobject.y_ecc;% Y profile NOT RESAMPLED 
            switch type_sampl
                case {'spatialuniform';'spatunif';'spatial'}
                    % Along Y
                    % NOT Re-Sampling
                    % Along  U
                    dU=v_sampl(1); %delta U
                    
                case {'range';'range unif';'range uniform'}
                    %Along Y
                    % NOT Re-Sampling                    
                    %Along U
                    nU=v_sampl(1);% #samples along Y                    
                otherwise
                    error ('NOT valid type of LF input sampling')                
            end
        end
        
        
        %Create intermediate object at each sampled eccentricity
        for yi=1:length(vettY)
            [Obj]=paraxCreateObject(LFobject.z_pos,vettY(yi),'point',LFobject.unit);
            [ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,Obj);
            LFinput.ISobj{yi}=ImagSyst.object{end}; %append the generated object in the imaging system
            
            %FIND direction limits      
            switch LFunit
                case {'paraxial';'parax'}
                    Uup=ImagSyst.object{end}.meridionalPlane.comaRay.upper.angleParax;
                    Ulo=ImagSyst.object{end}.meridionalPlane.comaRay.lower.angleParax;                
                case {'real'}
                    Uup=ImagSyst.object{end}.meridionalPlane.comaRay.upper.angle;
                    Ulo=ImagSyst.object{end}.meridionalPlane.comaRay.lower.angle;
                otherwise 
                    error ('Light Field unit is not valid! Specify between paraxial and real')
            end
            % SAMPLING and create the the INPUT
            switch type_sampl
                case {'spatialuniform';'spatunif';'spatial'}
                    warning('For have uniform sampling of U and and sam sample at different wavelengths, used the MAX range for U')
                    vettU=[min(min(Ulo)):dU:max(max(Uup))];
                    vettU=repmat(vettU,1,size(ImagSyst.wave,1)); %repmat for spectrum
                    
                case {'range';'range unif';'range uniform'}
                    %uniform sampling on the variable range
                    dU=(Uup-Ulo)./(nU-1); %delta U
                    for li=1:size(ImagSyst.wave,1)
                        vettU(:,li)=[Ulo(li):dU(li):Uup(li)];
                    end 
                otherwise
                    error ('NOT valid type of LF input sampling')                

            end
            % CREATE MESH of THE LIGHT FIELD
            for li=1:size(ImagSyst.wave,1)
                mU(yi,:,li)=vettU(:,li); % spectral-dependent ray angles spread from specific eccentricity position                 
            end            
        end
        
        % CREATE VECTORs and MESH
        for li=1:size(ImagSyst.wave,1)
            [Y(:,:,li),U(:,:,li)]=meshgrid(mU(end,:,li),vettY); % Create MESH
        end
        
        %CREATE OUTPUT
        %vector
        LFinput.vettY=repmat(vettY,size(ImagSyst.wave,1),1);
        LFinput.vettU=mU; 
        %matrix for mesh
        LFinput.U=mU;
        LFinput.Y=Y;
        
        
        
    case{'customized'}
        error ('To be completed!!')
end



for li=1:size(ImagSyst.wave,1)
     % RADIANCE INTENSITY
    switch LFobject.rad_type
        case {'isotropic'}
            LFinput.intensity(:,:,li)=ones(size(LFinput.U(:,:,li),1),size(LFinput.U(:,:,li),2)); %all rays with same intensity
        case {'customized'}
            error ('*************TO BE COMPLETED********************')
        otherwise
            error ('NOT valid type of object source/reflectance ! Suggestion: include also this solution')
    end
end


%% COMPUTE DIGITAL WIGNER FUNCTION
for li=1:size(ImagSyst.wave,1)
    % ALONG U
    Umax=max(max(LFinput.U(:,:,li)));
    Umin=min(min(LFinput.U(:,:,li)));
    if length(dU)>1
        rangeU=[Umin:dU(li):Umax]; % range along U
    else
        rangeU=[Umin:dU:Umax]; % range along U
    end
    % ALONG Y
    Ymax=max(max(LFinput.Y(:,:,li)));
    Ymin=min(min(LFinput.Y(:,:,li)));
    
    if (Ymax==Ymin)
        rangeY(li,:)=Ymax;        
    else
        if length(dY)>1
            rangeY=[Ymin:dY(li):Ymax]; % range along U
        else
            rangeY=[Ymin:dY:Ymax]; % range along U
        end
    end
    
    
    % FILL the matrix
    M=ones(numel(rangeY),numel(rangeU)); %empty
    resU=reshape(LFinput.U(:,:,li),1,[]);
    resY=reshape(LFinput.Y(:,:,li),1,[]);
    resI=reshape(LFinput.intensity(:,:,li),1,[]);
    
%     for ni=1:numel(resU)
%         M(resY(ni),resU(ni))=resI(ni); %fill the spot with the intensity values
%     end
%     
%     for ri=1:size(LFinput.U(:,:,li),1)
%         for ci=1:size(LFinput.U(:,:,li),2)
%             
%         end
%     end
    
    %OUTPUT
    LFinput.WD(:,:,li).U=resU;
    LFinput.WD(:,:,li).Y=resY;    
    LFinput.WD(:,:,li).normI=resI;
    
end


