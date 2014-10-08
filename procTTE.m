function procTTE(DataDir,fidx,options)

if ~exist('DataDir','var')
    DataDir = pwd;
end
if ~exist('fidx','var')
    fidx = -1;
end
cd(DataDir)
% Input Parameters
if ~exist('options','var')
    options.dataflow = struct(...
        'display',1 ...
        ,'ARFI',1 ...
        ,'SWEI',1 ...
        ,'setID',fidx ...
        ,'saveRes',0 ...
        );
    options.dispEst = struct(...
        'method','Loupas'...
        ,'ref_type','progressive' ...   % independent/common/progressive - indicates whether "no push" and "push" data sets will have independent references or a common reference, or use a moving reference
        ,'ref_idx',[] ...
        ,'noverlap',5 ...        % DO NOT CHANGE - number of time steps that are common between "no push" and "push" data sets (determined in sequenceParams file)
        ,'nreverb',2 ...         % Not having the correct nreverb could mess up displacements when using progressive ref_type
        ,'interpFactor',5 ...
        ,'kernelLength',4 ...
        ,'ccmode', 0 ...
        );
end
% Extract timeStamp
if ispc
    addpath C:\users\vrk4\Documents\GitHub\SC2000\arfiProcCode\
    addpath(genpath('C:\users\vrk4\Documents\GitHub\TTEProcCode'))
elseif isunix
    addpath /emfd/vrk4/GitHub/SC2000/arfiProcCode
    addpath(genpath('/emfd/vrk4/GitHub/TTEProcCode'))
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
fprintf('Loading data with timeStamp = %s (Set # %d)\n', timeStamp,size(list,1));

% Extract ECG/Trigger Data
if exist(strcat('ECG_data_',timeStamp,'.mat'),'file')
    fprintf(1,'Loading ECG/Trigger Data...\n');
    ecgdata = extractECG(timeStamp,1);
else
    ecgdata = [];
end
% Extract B-mode Data
fprintf(1,'Extracting B-mode Data...\n');
[bdata,bmodeSave] = extractBmode(timeStamp);
% Extract ARFI/SWEI Data
fprintf(1,'Extracting ARFI/SWEI Data...\n');
if options.dataflow.ARFI
    [arfidata,arfiSave,options] = extractMmode(timeStamp,options,'ARFI');
    arfi_par = load(sprintf('arfi_par_%s.mat',timeStamp));
    arfidata = interpPushReverb(arfidata,options,arfi_par,'');
else
    arfidata = [];
    arfi_par = [];
end
if options.dataflow.SWEI
    [sweidata,sweiSave,options] = extractMmode(timeStamp,options,'SWEI');
    swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
    sweidata = interpPushReverb(sweidata,options,swei_par,'');
else
    sweidata = [];
    swei_par = [];
end

%%
if options.dataflow.display
dispTTE(ecgdata,bdata,arfidata,arfi_par,sweidata,swei_par,options);
end
%%
keyboard
%% Save time stamped results file
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
