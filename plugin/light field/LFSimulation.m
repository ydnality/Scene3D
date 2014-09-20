% Compute the light field transformation occuring from a object point
% to a film for a given imaging system.
%There are several method to compute the light field defined as [Y;U] (eccentricity, direction)


function [LFout,varargout]=LFSimulation(LFinput,ImagSyst,film_index,comput_type,varargin)

%INPUT
%LFobject: struct describing the source of light field
%ImagSyst: imaging systen
%film_index: number of the film to image the LF (can be a vector to image
%the point to different image planes
%LFunit: specify if the LF is expressed in paraxial or real unit,('paraxial','real')
%compu_type : several type of computation ('ABCD','ABCD sequential', )
%type_sampl: type of sampling of the Light Field input    
%       -> spatialuniform  (varargin{1}:=[dY;dU])
%       -> range uniform    (varargin{1}:=[#Y;#U])

%OUTPUT
%LFout: output light field




%% CHECK

if any(film_index>length(ImagSyst.film))
    ind_over=find (film_index>length(ImagSyst.film))
    film_index(ind_over)=1;
    warning ('Not valid film number! Set to the first film of the Imaging System')
end

%Relevant parameter
wave=LFinput.wave;


%% COMPUTE LIGHT FIELD TRANSFORMATION

%Prepare output
for nfilm=1:length(film_index)
    LFout{nfilm}.LFunit=LFinput.LFunit;
    LFout{nfilm}.wave=LFinput.wave;
    LFout{nfilm}.ISfilm=ImagSyst.film{nfilm}; %attached film
    LFout{nfilm}.z_pos=ImagSyst.film{nfilm}.z_pos; %get the position
end


switch comput_type
    
    case {'ABCD';'standard';'standard ABCD';'abcd';'matrix'}
        log.FLlist={'source';'input OS';'output OS'; 'film' }; %list of position at which LF is stored
        log.FL{1}.Y=[LFinput.Y];log.FL{1}.U=[LFinput.U]; %input LF
        
        % Check if the sistem can be computed by 2x2 matrix or by 3x3
        % matrix (for quasi-rotationally symmetric)
        if all(all(ImagSyst.matrix.augParam==0))
            m2=1; %true
        else
            m2=0; %false
        end
        
        for li=1:size(wave,1)
            switch LFinput.LFunit
                case {'parax'}                    %NO MODULATION OF INTENSITY (4 now)!
                    for ri=1:size(LFinput.Y,1)
                         for ci=1:size(LFinput.U,2) 
                             
                            % STEP 1: From source of light field
                            LF1in=[LFinput.Y(ri,ci,li);LFinput.U(ri,ci,li)]; %create input for STEP 1
                            LF1out=LFinput.ISobj{ri}.matrix.abcd(:,:,li)*LF1in;
                            %Store in log file
                            log.FL{2}.Y(ri,ci,li)=LF1out(1); log.FL{2}.U(ri,ci,li)=LF1out(2);
                            % STEP 2: Through the Optical System
                            if m2
                                LF2out=ImagSyst.matrix.abcd(:,:,li)*LF1out;
                                %Store in log file
                                log.FL{3}.Y(ri,ci,li)=LF2out(1); log.FL{2}.U(ri,ci,li)=LF2out(2);
                            else
                                LF2out_aug=ImagSyst.matrix.abcdef(:,:,li)*[LF1out;1];
                                LF2out=LF2out_aug(1:2,1); % remove 3th value                                
                                log.FL{3}.Y(ri,ci,li)=LF2out(1); log.FL{2}.U(ri,ci,li)=LF2out(2);                                
                            end
                            % STEP 3: From the last surface of the optical
                            % system to the film/s
                            
                            for nfilm=1:length(film_index) %Imaging into the different FILM
                                 if all(all(ImagSyst.film{nfilm}.matrix.abcdef(1:2,3,:)==0)) %Check if the film is rotationally symmetric or quasi-rotationally symmetric
                                    LF3out=ImagSyst.film{nfilm}.matrix.abcd(:,:,li)*LF2out;
                                    %Store in log file
                                    log.FL{nfilm,4}.Y(ri,ci,li)=LF3out(1); log.FL{nfilm,4}.U(ri,ci,li)=LF3out(2);
                                else
                                    LF3out_aug=ImagSyst.matrix.abcdef(:,:,li)*[LF2out_aug];
                                    LF3out=LF2out_aug(1:2,1); % remove 3th value                                
                                    log.FL{3}.Y(ri,ci,li)=LF3out(1); log.FL{2}.U(ri,ci,li)=LF3out(2);                                
                                 end
                                 %SET OUTPUT
                                 LFout{nfilm}.Y(ri,ci,li)=LF3out(1);LFout{nfilm}.U(ri,ci,li)=LF3out(2);
                                 
                                 %NO INTENSITY MODULATION (4 now)!
                                 LFout{nfilm}.intensity(ri,ci,li)=LFinput.intensity(ri,ci,li);
                            end
                           
                            
                             
                         end
                    end
                case {'real'}
                    error ('Code has to be completed for REAL coordinates')
            end
                
            
            
        end
    case {'ABCD seq';'ABCD sequential';'matrix seq'}
    case {''}
    otherwise
        error ('Do not specify a VALID way to COMPUTE LIGHT FIELD!')
end




%% APPEND the log-file to the VARAGOUT
varargout{1}=log;

