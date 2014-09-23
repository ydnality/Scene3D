classdef LTC < clonableHandleObject
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wave;
        AInterp; 
        A1stInterp; 
        A2ndInterp; 
    end
    
    methods
        function obj = LTC(varargin)
        % Initialization of the Volume of Linear Transformations object    
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'wave'  %somehow syncrhonize this...
                        obj.wave = varargin{ii+1};
                    case 'ainterp'
                        obj.AInterp = varargin{ii+1};
                    case 'a1stinterp'
                        obj.A1stInterp = varargin{ii+1};
                    case 'a2ndinterp'
                        obj.A2ndInterp = varargin{ii+1};   
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
        end    
        
        function [outLFObject] = applyOnLF(obj, inputLFObject, apertureRadius)
            %applies the linear transform on an input light-field, and
            %returns the output light=field
            %
            %

            %withinAperture = middleXY(:,1).^2 + middleXY(:,2).^2 <= adjustedMiddleApertureRadius.^2;%apertureMiddleD/2;
            outputLF = [];
            outWaveIndex = [];
            inputLF = inputLFObject.get('LF');
            waveIndex = inputLFObject.get('waveIndex');
            %loops through all different waves
            for w = 1:length(obj.wave)
                inCurrentWaveBand = (waveIndex == w);   %todo: perhaps rethink this whole waveband thing...
                
                xCurrentWave = inputLF(:,inCurrentWaveBand);  %sifts out rays that only of the current wave band
                
                
                %middleXYCurrentWave = middleXY(inCurrentWaveBand,:);
                %withinAperture = middleXYCurrentWave(:,1).^2 + middleXYCurrentWave(:,2).^2 <= apertureRadius.^2;
                
                A1stInterpCurrentWave = obj.A1stInterp(:,:,w);  %use 2 dimensions of the matrices for multiplication
                A2ndInterpCurrentWave = obj.A2ndInterp(:,:,w);
                firstHalf = A1stInterpCurrentWave * xCurrentWave;
                
                withinAperture = firstHalf(1,:).^2 + firstHalf(2,:).^2 <= apertureRadius.^2;
                
                firstHalfBlock = firstHalf(:, withinAperture); %Apply aperture to rays from the first half
                bEstInterp = A2ndInterpCurrentWave * firstHalfBlock;
                
                %bOrigMaskedCurrentWave = bOrig(:, withinAperture); %ground
                %truth rays - for debug
                outputLF = cat(2, outputLF, bEstInterp);   %concatenate interpolatedB to final B list
                outWaveIndex = cat(1, outWaveIndex, ones(size(bEstInterp,2), 1)* w); %keep track of waveIndex
            end
            
            %create output light field based off the lf matrix and the
            %waveindex vector
            outLFObject = LFC('LF', outputLF, 'waveIndex', outWaveIndex, 'wave', obj.wave);
            
            
        end
        
    end
    
end

