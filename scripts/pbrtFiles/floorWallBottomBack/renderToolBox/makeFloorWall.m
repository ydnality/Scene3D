%%% RenderToolbox3 Copyright (c) 2012 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.
%
%% Render a shiny sphere sitting on a table.
clear;
clc;


%% Choose example files, make sure they're on the Matlab path.
AddWorkingPath(mfilename('fullpath'));
parentSceneFile = 'floorWall259.dae';
mappingsFile = 'floorWallMappings.txt';
conditionsFile = 'floorWallConditions.txt';

%% Choose batch renderer options.
hints.imageWidth = 320;
hints.imageHeight = 240;
hints.outputSubfolder = mfilename();

%% Render with Mitsuba and PBRT
toneMapFactor = 10;
isScale = true;
for renderer = {'PBRT'}
    hints.renderer = renderer{1};
    nativeSceneFiles = MakeSceneFiles(parentSceneFile, conditionsFile, mappingsFile, hints);
    radianceDataFiles = BatchRenderISET(nativeSceneFiles, hints);
%     montageName = sprintf('%s (%s)', 'TableSphere', hints.renderer);
%     montageFile = [montageName '.png'];
%     [SRGBMontage, XYZMontage] = ...
%         MakeMontage(radianceDataFiles, montageFile, toneMapFactor, isScale, hints);
%     ShowXYZAndSRGB([], SRGBMontage, montageName);
end

oiWindow;