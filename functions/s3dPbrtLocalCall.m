function s3dPbrtLocalCall(fullfname, generatedDir, outfile)

    validateattributes(fullfname, {'char'}, {'nonempty'});
    validateattributes(generatedDir, {'char'}, {'nonempty'});
    validateattributes(outfile, {'char'}, {'nonempty'});

    pbrtExe = 'pbrt';
   %if s, error('PBRT not found'); end
    
    % [p,n,ext] = fileparts(fullfname);
    %cmd = sprintf('%s %s --outfile %s\n',pbrtExe,fullfname,outfile);
    [~,n,e] = fileparts(fullfname); % Get name of pbrt input file
    tempInputFile = fullfile(generatedDir, [n e]);
    cmd = sprintf('%s %s --outfile %s\n',pbrtExe,tempInputFile,outfile);
    % chdir(p)
    unix(cmd)
end