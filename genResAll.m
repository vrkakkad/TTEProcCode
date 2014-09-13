function genResAll(DataDir)

if ~exist('DataDir','var')
    DataDir = pwd;
end

cd(DataDir)


nfiles = length(dir('arfi_par_*')); 

% Convert this to parfor or cluster eventually
for i=1:nfiles
    procTTE(DataDir,i)
end