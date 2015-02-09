Scene3D

This is the PBRT -> ISET interface. 

We use a modified (forked) form of PBRT that is in 

  https://github.com/ydnality/pbrt-v2-spectral/tree/spectral

which is derived from the master repository

  https://github.com/ydnality/pbrt-v2-spectral

The Scene3D repository contains 

	* Matlab files that call PBRT
        * Some PBRT data (examples)
	* Matlab files that read PBRT output and convert them to ISET objects

Andy added a lot of code about rays, meshes, forward calculations, curved sensors, diffraction, two flash, and so forth.  He also implemented the docker management for pbrt.

Michael Pieroni added some optics calculations in plugin.  This included the ABCD transformations and Black Box Model (BBM).  We need to integrate and maintain this in some way.

