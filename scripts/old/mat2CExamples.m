%% Example fread/fwrite for Matlab to C data conversion
%
% Defaults are for 32 bit format, we think.
%
% See Matlab fopen/fread/fwrite pages for manyh more details.


%% Create some photons
scene = sceneCreate;
oi = oiCreate;
oi = oiCompute(scene,oi);
p = oiGet(oi,'photons');

% Store their sizes
[r,c,w] = size(p);
sizeP = numel(p);


%% Dump them to a file
fID = fopen('scratch.dat','w','b');
fwrite(fID,p,'double');
fclose(fID);
size(p);

%% Read them back and check
fID = fopen('scratch.dat','r','b');
p2 = fread(fID,sizeP,'double');
fclose(fID);
p2 = reshape(p2,r,c,w);

vcNewGraphWin; imagesc(p2(:,:,10))
vcNewGraphWin; imagesc(p(:,:,10))

%% Write a mix in big-endian 32 bit format
fID = fopen('scratch.dat','w','b');
fprintf(fID,'%d %d %d\n',r,c,w);
fwrite(fID,p,'double');
fclose(fID);

%%
function oi = pbrt2oi(fname)
%Convert pbrt multispectral irradiance data into an ISET scene
%
%   oi = pbrt2oi(fname)
%
%

if ieNotDefined('fname'), error('File name required.'); end

%% Load the stored photons produced by AL's pbrt code
fID = fopen(fname,'r','l');
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 Inf]);
photons = fread(fID,prod(A),'double');

%A(2) = rows, A(1) = columns
photons = reshape(photons,A(2),A(1),A(3));
fclose(fID);

% Set the OI data
oi = oiCreate;
oi = initDefaultSpectrum(oi);
% Set the optics parameters from somewhere
% optics = oiGet(oi,'optics');

oi = oiSet(oi,'cphotons',photons);

return


end

function pbrtRead(fName)



end

function pbrtWrite(photons)
%
%  RGB Format of photons to be saved as a pbrt formatted spectral file.
%
% 

end


p2 = reshape(p2,A(1),A(2),A(3));

vcNewGraphWin; imagesc(p2(:,:,10))
vcNewGraphWin; imagesc(p(:,:,10))

%% Figure out the ordering of the write
q = 1:32;

% Down the row, then to columns, then to the 3rd dimension
q = reshape(q,4,4,2);

fID = fopen('scratch.dat','w','b');
fprintf(fID,'%d %d %d\n',4,4,2);
fwrite(fID,q,'double');
fclose(fID);

fID = fopen('scratch.dat','r','b');
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 Inf]);
p2 = fread(fID,prod(A),'double');
fclose(fID);

