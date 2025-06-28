function outName = GetNextFileName( prefix )
% Create another file name with given prefix, by adding an integer index

% Define directory and initialize index
dataDir = '../DATA/';
ind     = 1;
flag    = true;

% Get list of all .mat files
fileList = dir(fullfile(dataDir, '*.mat'));

% Extract only file names
fileNames = {fileList.name};

% Loop until an unused filename is found
while flag
    candidate = [prefix, num2str(ind), '.mat'];
    if any(strcmp(candidate, fileNames))
        ind = ind + 1;
    else
        flag = false;
    end
end

% Construct the output file name
outName = fullfile(dataDir, [prefix, num2str(ind), '.mat']);
fprintf('Next file: %s\n', outName);