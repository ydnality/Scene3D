
%Render cones scene
sceneName = 'cones autofocus';
oi = s3dRenderScene('../../autofocus/cones/autofocusTest.pbrt', .050);
oi = oiSet(oi, 'name', sceneName);
% vcAddAndSelectObject(oi);
% oiWindow;

rgbImage = oiGet(oi, 'rgb image');
figure; imshow(rgbImage);

% you can now interact with rgbImage, which is a 300 x 450 x 3 matrix

% useful functions:
% var(x) -> computes the variance of an array
% ex. var(rgbImage(:)) computes the variance of the complete image
% you should probably compute the variance of a PATCH of the image as a
% metric for contrast

% -----path to the .pbrt file-----
filepath = '/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/autofocusTest.pbrt';

% -----copy a temp .pbrt file-----
% -----create a old_temp.pbrt file
copyfile(filepath,'/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in original parameter
% -----find string in old_temp.pbrt-----
% -----write to new_temp.pbrt-----
fid = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt','r'); % open .pbrt file
nfid = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','w');
while ~feof(fid)
     tline = fgets(fid); % fgetl() will read row by row, and will not ignore enpty row
     if~ischar(tline),break,end % check if success
     %disp(tline);
     xresfind = strfind(tline,'xresolution'); % find keyword "xresolution"
     yresfind = strfind(tline,'yresolution'); % find keyword "yresolution"
     cwfind = strfind(tline,'cropwindow'); % find keyword "cropwindow"
     fdfind = strfind(tline,'filmdistance'); % find keyword "filmdistance"
     % strfind() Found: returns starting index. Not found: '[]' empty array
     if xresfind ~= 0
         xres = sscanf(tline,'%*s %*s %d'); % get the xreselution (int)     
         fprintf(nfid,'%s',tline);
     elseif yresfind ~= 0
         yres = sscanf(tline,'%*s %*s %d'); % get the xreselution (int)
         fprintf(nfid,'%s',tline);     
     elseif cwfind ~= 0
         orig_cw = sscanf(tline,'%*s %*s %*s %d %d %d %d'); % get cropwindow parameter [xmin xmax ymin ymax]
         fprintf(nfid,'%s',tline);     
     elseif fdfind ~= 0
         fd = sscanf(tline,'%*s %*s %f'); % get the xreselution (float)
         fprintf(nfid,'%s',tline);
     else
     fprintf(nfid,'%s',tline);
     end
end
fclose(nfid);
fclose(fid); % close .pbrt file

% -----replace old_temp.pbrt with new_temp.pbrt
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
copyfile('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt');
% -----END-----



[x,y] = ginput(1);
x = sscanf(num2str(x),'%d'); y = sscanf(num2str(y),'%d');
imagexmin = x-10; imagexmax = x+10;
imageymin = sscanf(num2str(y-10),'%d'); imageymax = sscanf(num2str(y+10),'%d');
%fprintf('x:%d , y:%d\n',x,y);
%fprintf('x:%d , y:%d\n',imageymax,imageymin);

if x<0 || y<0 || x>xres || y>yres
    fprintf('ERROR: Out of image!\n');
    break;
elseif 0<=x && x<10 || 0<=y && y<10 || xres-10<x && x<=xres || yres-10<y && y<=yres
    fprintf('ERROR: Can not get portion image!\n');
    break;
else
    smallrgbImage = rgbImage( imagexmin:imagexmax, imageymin:imageymax, :);
    %imshow(smallImage);
    tempVar = 0;%var(smallrgbImage(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace x,y resolution and cropwindow parameter
fid = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt','r'); % open .pbrt file
nfid = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','w');
while ~feof(fid)
     tline = fgets(fid); % fgetl() will read row by row, and will not ignore enpty row
     if~ischar(tline),break,end % check if success
     %disp(tline);
     xresfind = strfind(tline,'xresolution'); % find keyword "xresolution"
     yresfind = strfind(tline,'yresolution'); % find keyword "yresolution"
     cwfind = strfind(tline,'cropwindow'); % find keyword "cropwindow"
     fdfind = strfind(tline,'filmdistance'); % find keyword "filmdistance"
     % strfind() Found: returns starting index. Not found: '[]' empty array
     if xresfind ~= 0
         xres = sscanf(tline,'%*s %*s %d'); % get the xreselution (int)
         tline = regexprep(tline, num2str(xres), '60'); % replace original xresolution
         fprintf(nfid,'%s',tline);
     elseif yresfind ~= 0
         yres = sscanf(tline,'%*s %*s %d'); % get the xreselution (int)
         tline = regexprep(tline, num2str(yres),'40'); % replace original crowwindow
         fprintf(nfid,'%s',tline);     
     elseif cwfind ~= 0
         %orig_cw = sscanf(tline,'%*s %*s %*s %d %d %d %d'); % get cropwindow parameter [xmin xmax ymin ymax]
         %tline = regexprep(tline, num2str(orig_cw(1)), 'a');
         %tline = regexprep(tline, num2str(orig_cw(2)), 'b');
         new_cw = [(imagexmin)/xres; (imagexmax)/xres; (imageymin)/yres; (imageymax)/yres]; % image portion
         tline = regexprep(tline, '0 1 0 1', num2str(new_cw')); % replace original crowwindow
         fprintf(nfid,'%s',tline);     
%      elseif fdfind ~= 0
%          fd = sscanf(tline,'%*s %*s %f'); % get the xreselution (float)
%          %fprintf('filmdistance:%.2f\n',fd);
%          %tline = strrep(tline, num2str(fd), '36'); % replace 35 to 36
%          fprintf(nfid,'%s',tline);
     else
     fprintf(nfid,'%s',tline);
     end
end
fclose(nfid);
fclose(fid); % close .pbrt file

% -----replace old_temp.pbrt with new_temp.pbrt
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
copyfile('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt');
% -----END-----
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdmin = input('Please input min filmdistance:[1]','s');
fdmax = input('Please input max filmdistance:[100]','s');
fdmin = sscanf(fdmin,'%d'); fdmax = sscanf(fdmax,'%d');

optimalVar = tempVar; % Optimal var
opFilmdis = fd; % Optimal film distance

for n=fdmin:.1:fdmax
   
    [fid, msg] = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt','r'); % open .pbrt file
    [nfid, nmsg] = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','w');
    
    if msg == -1
        error(msg)
    end
    
    if nmsg == -1
        error(msg)
    end
    
    while ~feof(fid)
        %disp(ffid);
        tline = fgets(fid); % fgetl() will read row by row, and will not ignore enpty row
        if~ischar(tline),break,end % check if success
        fdfind = strfind(tline,'filmdistance'); % find keyword "filmdistance"
        % strfind() Found: returns starting index. Not found: '[]' empty array
        if fdfind ~= 0
            fd = sscanf(tline,'%*s %*s %f'); % get the xreselution (float)
            %fprintf('filmdistance:%.2f\n',fd);
            tline = strrep(tline, num2str(fd), num2str(n)); % replace fd to n
            fprintf(nfid,'%s',tline);
        else
            fprintf(nfid,'%s',tline);
        end
    end
    fclose(nfid);
    fclose(fid); % close .pbrt file
    
    % -----replace old_temp.pbrt with new_temp.pbrt
    delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
    %delete(filepath);
    copyfile('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
    %copyfile('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt',filepath');
    delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt');
    % -----END-----
    
    sceneName = 'cones autofocus';
    oi = s3dRenderScene('../../autofocus/cones/old_temp.pbrt', .050);
    oi = oiSet(oi, 'name', sceneName);
    % vcAddAndSelectObject(oi);
    % oiWindow;
    tempImage = oiGet(oi, 'rgb image');
    %smallTempimage = tempImage(imagexmin:imagexmax, imageymin:imageymax,:);
    %tempVar = var(smallTempimage(:));
    tempVar = var(tempImage(:));
    if tempVar > optimalVar
        optimalVar = tempVar;
        opFilmdis = n;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recover the file 
fid = fopen(filepath,'r'); % open .pbrt file
nfid = fopen('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','w');
while ~feof(fid)
     tline = fgets(fid); % fgetl() will read row by row, and will not ignore enpty row
     if~ischar(tline),break,end % check if success
     %disp(tline);
     xresfind = strfind(tline,'xresolution'); % find keyword "xresolution"
     yresfind = strfind(tline,'yresolution'); % find keyword "yresolution"
     cwfind = strfind(tline,'cropwindow'); % find keyword "cropwindow"
     fdfind = strfind(tline,'filmdistance'); % find keyword "filmdistance"
     % strfind() Found: returns starting index. Not found: '[]' empty array
     if xresfind ~= 0
         xres = sscanf(tline,'%*s %*s %d'); % get the xreselution (int)
         tline = regexprep(tline, '60', num2str(xres)); % replace original xresolution
         fprintf(nfid,'%s',tline);
     elseif yresfind ~= 0
         yres = sscanf(tline,'%*s %*s %d'); % get the xreselution (int)
         tline = regexprep(tline, '40',num2str(yres)); % replace original crowwindow
         fprintf(nfid,'%s',tline);     
     elseif cwfind ~= 0
         orig_cw = sscanf(tline,'%*s %*s %*s %d %d %d %d'); % get cropwindow parameter [xmin xmax ymin ymax]
         tline = regexprep(tline, num2str(orig_cw(1)), '0');
         tline = regexprep(tline, num2str(orig_cw(2)), '1');
         tline = regexprep(tline, num2str(orig_cw(3)), '0');
         tline = regexprep(tline, num2str(orig_cw(4)), '1'); % replace original crowwindow
         fprintf(nfid,'%s',tline);     
     elseif fdfind ~= 0
         fd = sscanf(tline,'%*s %*s %f'); % get the xreselution (float)
         %fprintf('filmdistance:%.2f\n',fd);
         tline = strrep(tline, num2str(fd), num2str(opFilmdis)); % replace 
         fprintf(nfid,'%s',tline);
     else
     fprintf(nfid,'%s',tline);
     end
end
fclose(nfid);
fclose(fid); % close .pbrt file

% -----replace old_temp.pbrt with new_temp.pbrt
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
copyfile('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt','/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt');
% -----END-----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% -----replace old_temp.pbrt with new_temp.pbrt
delete(filepath);
copyfile('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt',filepath);
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/old_temp.pbrt');
delete('/home/andy/Dropbox/Scene3D/Scene3DMatlab/autofocus/cones/new_temp.pbrt');
% -----END-----


sceneName = 'cones autofocus';
oi = s3dRenderScene('../../autofocus/cones/autofocusTest.pbrt', .050);
oi = oiSet(oi, 'name', sceneName);
% vcAddAndSelectObject(oi);
% oiWindow;
showImage = oiGet(oi, 'rgb image');
figure; imshow(showImage);
toc;




% user input
% str = input(userinput);