%% s_pbrtExample
%
%
%
% Illustrate how to write a batch script to render a small scene multiple
% times with pbrt, using pbrtCreate, pbrtSet, pbrtWrite and so forth
%
%

%%
pbrt=pbrtCreate
pbrtWrite(pbrt, 'deleteMe.pbrt');

oi = s3dRenderScene('deleteMe.pbrt', 50, [dataPath '/tmp/'])
%% End
