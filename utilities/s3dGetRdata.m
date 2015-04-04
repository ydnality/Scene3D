function s3dGetRdata(rurl,fname)
% Write data file from a remote data set to fname
%
%   s3dGetRdata(rurl,fname)
% 
% There is a version of this file for each different project. This is the
% version for
%
% The data at the remote URL (rurl) are written to the fname inside the
% remote data directory.  The remote data directory is stored in getpref.
% For each application this directory is set separately.  For example, in
% Scene3d this is set by default to setpref('scene3d','rdatadir',val).
% In Scene3d the default directory is fullfile(s3dRootPath,'rdata');
%
% In ISETBIO this is set to setpref('isetbio','rdatadir') and the default
% is fullfile(isetbioRootPath,'rdata').
%
% To check the remote data directory for scene3d use: 
%    getpref('scene3d','rdatadir')
%
% Example
%
%  rurl = 'http://scarlet.stanford.edu/validation/SCIEN/LIGHTFIELD/scenes/benchLF.mat'
%  s3dGetRdata(rurl,'benchLF.mat');
% 
% BW, Vistasoft Team, 2015

%%
if ieNotDefined('rurl')
    web('http://scarlet.stanford.edu/validation/SCIEN')
    error('You must specify a remote url'); 
end

% Check for a remote data directory.  If not set, use the default
if ispref('scene3d','rdatadir'), 
    rdatadir = getpref('scene3d','rdatadir'); 
else
    rdatadir = fullfile(s3dRootPath,'rdata');
    setpref('scene3d','rdatadir',rdatadir);
    if ~exist(rdatadir,'dir')
        curd = pwd; 
        mkdir(rdatadir); chdir(rdatadir)
        s = sprintf('Remote data directory from \n');
        fid = fopen('Readme.txt','w');
        fprintf(fid,'%s',s); fclose(fid);
        chdir(curd);
    end
end

% Make the call
oname = fullfile(rdatadir,fname);
[s,r] = urlwrite(rurl,oname);

if r
    fprintf('Remote data written to %s\n',oname);
else
    fprintf('Problem writing file\n');
end


end