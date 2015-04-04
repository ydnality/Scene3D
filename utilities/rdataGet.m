function rdataGet(rurl,fname)
% Get a data file from the remote vista data set
%
%   rdataGet(rurl,fname)
%
% The data at the remote URL (rurl) are written to the fname inside the
% remote data directory stored in getpref

% Example
%
%  rurl = 'http://scarlet.stanford.edu/validation/SCIEN/LIGHTFIELD/scenes/benchLF.mat'
%  rdataGet(rurl,'benchLF.mat');
%
% getpref('scene3d','rdatadir')

% d = fullfile(s3dRootPath,'rdata');

% setpref('scene3d','rdatadir',
if ieNotDefined(rurl)
    web('http://scarlet.stanford.edu/validation/SCIEN')
    error('You must specify a remote url'); 
end

if ispref('scene3d','rdatadir'), 
    rdatadir = getpref('scene3d','rdatadir'); 
else
    d = fullfile(s3dRootPath,'rdata');
    setpref('scene3d','rdatadir',d);
    if ~exist(d,'dir')
        curd = pwd; 
        mkdir(d); chdir(d)
        s = sprintf('Remote data directory from \n');
        fid = fopen('Readme.txt','w');
        fprintf(fid,'%s',s); fclose(fid);
        chdir(curd);
    end
end

oname = fullfile(rdatadir,fname);

[s,r] = urlwrite(rurl,oname);

fprintf('Remote data written to %s\n',oname);

end