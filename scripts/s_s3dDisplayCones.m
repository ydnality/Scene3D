

%% ISET will read the PBRT output
oi = pbrt2oi('cones_hq.dat');
vcAddAndSelectObject(oi);
oiWindow;

oi = pbrt2oi('cones_d_hq.dat');
vcAddAndSelectObject(oi);
oiWindow;

oi = pbrt2oi('cones_dca_hq.dat');
vcAddAndSelectObject(oi);
oiWindow;