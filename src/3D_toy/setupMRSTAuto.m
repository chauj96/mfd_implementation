function setupMRSTAuto(mrstRoot)
% Find or install MRST automatically using git.
% Case 1: If MRST is already on the MATLAB path, just activate needed modules.
% Case 2: Else if mrstRoot/startup.m exists, run it.
% Case 3: Else clone MRST from GitHub into mrstRoot, then run startup.m.

    if nargin < 1 || isempty(mrstRoot)
        mrstRoot = fullfile(fileparts(mfilename('fullpath')), 'external', 'mrst');
    end

    % Case 1: MRST already available on path
    if exist('mrstModule', 'file') == 2 && exist('computeGeometry', 'file') == 2
        fprintf('MRST already available on MATLAB path.\n');
        mrstModule add upr
        return;
    end

    startupFile = fullfile(mrstRoot, 'startup.m');

    % Case 2: Local MRST folder already exists
    if exist(startupFile, 'file') == 2
        fprintf('Found MRST at: %s\n', mrstRoot);
        run(startupFile);
        mrstModule add upr
        return;
    end

    % Case 3: Need to clone
    parentDir = fileparts(mrstRoot);
    if ~exist(parentDir, 'dir')
        mkdir(parentDir);
    end

    repoURL = 'https://github.com/SINTEF-AppliedCompSci/MRST.git';
    cmd = sprintf('git clone %s "%s"', repoURL, mrstRoot);

    fprintf('MRST not found. Cloning from GitHub...\n');
    [status, msg] = system(cmd);

    if status ~= 0
        error(['git clone failed.\n' ...
               'Command output:\n%s\n'], msg);
    end

    if exist(startupFile, 'file') ~= 2
        error('Clone completed, but startup.m was not found.');
    end

    run(startupFile);
    mrstModule add upr

    fprintf('MRST setup complete.\n');
end