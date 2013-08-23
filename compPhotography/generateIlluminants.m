%%% This script generates illuminants and saves them to a text file, to be
%%% read by pbrt

%% generate illuminant files
d65 = illuminantCreate('d65', 400:10:700);
d65Photons = illuminantGet(d65, 'photons');
figure; plot(illuminantGet(d65, 'wave'), d65Photons');

tungsten = illuminantCreate('tungsten', 400:10:700);
tungstenPhotons = illuminantGet(tungsten, 'photons');
figure; plot(illuminantGet(tungsten, 'wave'), tungstenPhotons);

fluorescent = illuminantCreate('fluorescent', 400:10:700);
fluorescentPhotons = illuminantGet(fluorescent, 'photons');
figure; plot(illuminantGet(fluorescent, 'wave'), fluorescentPhotons);

% newIlluminant = illuminantCreate('vivitarFlash', 400:10:700);
% newIlluminantPhotons = illuminantGet(d65, 'photons');
% figure; plot(illuminantGet(newIlluminant, 'wave'), newIlluminantPhotons');


%% export to file in PBRT format

% open a file for writing (derived from:
% http://www.mathworks.com/help/matlab/import_export/writing-to-text-data-files-with-low-level-io.html)

illuminant = tungsten;
illuminantPhotons = tungstenPhotons;
outputName = 'tungsten.txt'

fid = fopen(outputName, 'w');
y = [illuminantGet(illuminant, 'wave'); illuminantPhotons' * 10^-16];
fprintf(fid, '%f  %f\n', y);
fclose(fid);


%% get reflectances from file

importedData = load([isetRootPath '/data/surfaces/reflectances/' 'DupontPaintChip_Vhrel.mat']);
figure; plot(importedData.data); title('DupontPaintChip Vhrel.mat reflectances');

importedData = load([isetRootPath '/data/surfaces/reflectances/' 'Nature_Vhrel.mat']);
figure; plot(importedData.data); title('Nature Vhrel.mat reflectances');

importedData = load([isetRootPath '/data/surfaces/reflectances/' 'Clothes_Vhrel.mat']);
figure; plot(importedData.data); title('Clothes Vhrel.mat reflectances');

importedData = load([isetRootPath '/data/surfaces/reflectances/' 'Food_Vhrel.mat']);
figure; plot(importedData.data); title('Food_Vhrel.mat reflectances');

importedData = load([isetRootPath '/data/surfaces/reflectances/' 'esserChart.mat']);
figure; plot(importedData.data); title('esserChart.mat reflectances');

importedData = load([isetRootPath '/data/surfaces/reflectances/' 'HyspexSkinReflectance.mat']);
figure; plot(importedData.data); title('HyspexSkinReflectance.mat reflectances');



%% export reflectance to a PBRT file

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 23);
figure; plot(reflectanceOut);
outputName = 'paintPurple.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%
wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 24);
figure; plot(reflectanceOut);
outputName = 'paintTorquoise.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);


%%%
wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 26);
figure; plot(reflectanceOut);
outputName = 'paintRed.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%
wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 48);
figure; plot(reflectanceOut);
outputName = 'paintGreen.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 4);
figure; plot(reflectanceOut);
outputName = 'paint4BlueGreen.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 80);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint80gray.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 50);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint50green.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 3);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint3red.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);


%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 100);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint100red.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 115);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint115Blue.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 97);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint97Yellow.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 77);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint77Orange.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

%%%

wave = importedData.wavelength';
reflectanceOut = importedData.data(:, 65);
xyz2srgb(reshape(ieXYZFromPhotons(reflectanceOut, 400:10:700), [1 1 3]))
figure; plot(reflectanceOut);
outputName = 'paint65Yellow.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/reflectances/' outputName], 'w');
y = [wave reflectanceOut]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);


%% export flash to PBRT format
importedData = load([isetRootPath '/data/lights/' 'VivitarFlash.mat']);
desiredWave = 400:10:700;
SPD = ieReadSpectra('data/lights/VivitarFlash', desiredWave);
figure; plot(desiredWave, SPD); title('VivitarFlash reflectances');
outputName = 'vivitarFlash.txt'

fid = fopen([s3dRootPath '/scripts/pbrtFiles/lights/' outputName], 'w');
y = [desiredWave' SPD * 5000]';  %the printing process transposes the matrix!
fprintf(fid, '%f  %f  \n', y);
fclose(fid);

