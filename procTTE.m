function procTTE(DataDir,fidx,varargin)

if ~exist('DataDir','var')
    DataDir = pwd;
end

if ~exist('fidx','var')
    fidx = -1;
end

cd(DataDir)
%% Input Parameters
options.dataflow = struct(...
    'ARFI',1 ...
    ,'SWEI',1 ...
    ,'setID',fidx ...
    ,'saveRes',1 ...
    );

options.dispEst = struct(...
    'method','Loupas'...
    ,'ref_type','common' ...   % independent/common/progressive - indicates whether "no push" and "push" data sets will have independent references or a common reference, or use a moving reference
    ,'ref_idx',[] ...
    ,'noverlap',5 ...        % DO NOT CHANGE - number of time steps that are common between "no push" and "push" data sets (determined in sequenceParams file)
    ,'interpFactor',5 ...
    ,'kernelLength',4 ...
    ,'ccmode', 0 ...
    );

% options.motionFilter = struct(...
%     'enable',1 ...
%     ,'method','Polynomial' ...
%     ,'order',1 ...
%     ,'timeRange',[-inf -0.32 6.72 6.88] ...
%     ,'passBand',[20 2000] ...
%     );

%% Extract timeStamp
if ispc
    addpath C:\users\vrk4\Documents\GitHub\SC2000\arfiProcCode\
    addpath(genpath('C:\Users\vrk4\Tools\TransthoracicStudy\procCode')) % Change to GitHub path when it is ready
else
    addpath /emfd/vrk4/GitHub/SC2000/arfiProcCode
    % include GitHub path when it is ready
end

list = dir('arfi_par_*'); % get timeStamp based on existance of ARFI par files

if size(list,1)<options.dataflow.setID
    error('Data set index requested greater than number of data sets')
end

% Reading in timestamp for data set
if length(options.dataflow.setID)==14
    timeStamp = options.dataflow.setID;
elseif options.dataflow.setID == -1
    timeStamp = list(end).name(end-17:end-4);
else
    timeStamp = list(options.dataflow.setID).name(end-17:end-4);
end

fprintf('Loading data with timeStamp = %s\n', timeStamp);
%% Extract B-mode Data
fprintf(1,'Extracting B-mode Data...\n');
[bdata,bmodeSave] = extractBmode(timeStamp);

%% Extract ARFI/SWEI Data
fprintf(1,'Extracting ARFI/SWEI Data...\n');
if options.dataflow.ARFI
    [arfidata,arfiSave,options] = extractMmode(timeStamp,options,'ARFI');
    arfi_par = load(sprintf('arfi_par_%s.mat',timeStamp));
%     if options.motionFilter.enable
%         [arfidata] = motionFilter(arfidata,options);
%     end
end
if options.dataflow.SWEI
    [sweidata,sweiSave,options] = extractMmode(timeStamp,options,'SWEI');
    swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
%     if options.motionFilter.enable
%         [sweidata] = motionFilter(sweidata,options);
%     end
end

% Save time stamped results file
if options.dataflow.saveRes
    tic
    resfile = ['res_' timeStamp '.mat'];
    if (options.dataflow.ARFI && options.dataflow.SWEI)
        save(resfile,'bdata','arfidata','sweidata','options','-v7.3');
    elseif (options.dataflow.ARFI && ~options.dataflow.SWEI)
        save(resfile,'bdata','arfidata','options','-v7.3');
    elseif (~options.dataflow.ARFI && options.dataflow.SWEI)
        save(resfile,'bdata','sweidata','options','-v7.3');
    else
        save(resfile,'bdata','options','-v7.3');
    end
    fprintf(1,'Save Time = %2.2fs\n',toc)
end
