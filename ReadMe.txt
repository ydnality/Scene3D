First version of renderToolBox to ISET interface. 

-makeDepthMap.m converts a binary depth map to a matlab format.  
-makeRadianceScene.m is an ISET function that takes in a given scene output from renderToolBox (picMat_orange_rad.mat) 
as well as a depth map (in matlab format), then converts this to an ISET scene.

-binaryRead.m is a helper function that allows makeDepthMap to do it's thing.
-depthmap.zbf is the original depth map (in binary format)
-picMac_orange_rad.mat is the spectral radiance, from renderToolBox.  MakeRadianceScene uses this and converts it into ISET format. 

These functions/scripts are specific to the data provided, for now.  
Also, the distance right now is very general for now (between 0 and 1, no units provided, 
where 0 is the minimum valid specified distance and 1 is the maximum valid specified distance). 
This was put into place so that the depth map would be easier to see as an image, 
but the absolute distance is important to be precise. 
In makeDepthMap.m, the minimum and maximum recognized distances can be changed.  


to run, run the following commands:

makeDepthMap
makeRadianceScene