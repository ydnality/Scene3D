

%% ISET will read the PBRT output
oi = pbrt2oi('benchScene_hq.dat');
vcAddAndSelectObject(oi);
oiWindow;

oi = pbrt2oi('benchScene_d_hq.dat');
vcAddAndSelectObject(oi);
oiWindow;

oi = pbrt2oi('benchScene_dca_hq.dat');
vcAddAndSelectObject(oi);
oiWindow;