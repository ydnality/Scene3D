%must run this from the scripts directory for now!
%combines and displays resulting output from parallel job
%need to initially get the first OI to initialize tempPhotons
numProc = 24;
chdir('tempDir1');
firstOi = pbrt2oi('chromaticAberration_d.dat');
vcAddAndSelectObject(firstOi);
cd ..
%for some reason pbrt2oi does a cd.. already - should change this...
tempPhotons = oiGet(firstOi, 'photons');

for i = 2:24
          chdir(['tempDir' int2str(i)])
oi = pbrt2oi('chromaticAberration_d.dat');
tempPhotons = tempPhotons + oiGet(oi, 'photons');
    vcAddAndSelectObject(oi); 
    cd ..
end

tempPhotons = tempPhotons./numProc;
newOi = oi;
newOi = oiSet(newOi, 'photons', tempPhotons);
vcAddAndSelectObject(newOi); oiWindow;