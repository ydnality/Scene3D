function fileRead(obj, fullFileName)
% Reads PBRT lens file matrix of data
%
%   lens.fileRead(fullFileName)
%
% The file has focal length information added to the header
% This function converts the PBRT matrix of data into the format
% that we use for setting up the multielement lens.

% Open the lens file
%fid = fopen(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'));
fid = fopen(fullFileName);

% Read each of the lens and close the file
import = textscan(fid, '%s%s%s%s', 'delimiter' , '\t');
fclose(fid);

%first find the start of the lens line demarked by #   radius
firstColumn = import{1};

%Scan past the comment lines, looking for the data
% Check the first entry in each line for #
% When it is not, read the rest of the file as a
% d = fread(mumble,double)
% d= reshape(d,X,4);
%
% radius axpos N aperture
% The current test is just to see if the word radius is on the
% line.  This might be improved.
continu = true;
i = 1;
while(continu && i <= length(firstColumn))
    compare = regexp(firstColumn(i), 'radius');
    if(~(isempty(compare{1})))
        continu = false;
    end
    i = i+1;
end
% i is the row where the data begin

% The next lens are the matrix data
% put data into lens object
radius = str2double(import{1});
radius = radius(i:length(firstColumn));


% Change from pbrt Scene3D format to raytrace Scene3D format
% In PBRT, the row has the offset from the previous surface.  In
% PBRT the data are read from the bottom up.  The last row has no
% offset.
% In PBRT, we trace from the sensor to the scene.
% In Scene3D we trace from the scene to the sensor.
% So, the offsets are shifted down.  This means:
%
offset = str2double(import{2});
offset = offset(i:length(firstColumn));
offset = [0; offset(1:(end-1))];

% Index of refraction in the 3rd column
N = str2double(import{3});
N = N(i:length(firstColumn));

% Diameter of the aperture (or maybe radius.  Must determine).
aperture = str2double(import{4});
aperture = aperture(i:length(firstColumn));

%modify the object and reinitialize
obj.elementsSet(offset, radius, aperture, N);

end


