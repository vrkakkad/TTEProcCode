function procTTE(DataDir,fidx,options)

if ~exist('DataDir','var')
    DataDir = pwd;
end
if ~exist('fidx','var')
    fidx = -1;
end

cudir = pwd;
cd(DataDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default Input Parameters
if ~exist('options','var')
    options.dataflow = struct(...
        'stream', stream_idx ... % [1] - realTime, [2] - review, [3] - cluster
        ,'display', 1 ...
        ,'ecg_test', 0 ...
        ,'ARFI', 1 ...
        ,'SWEI', 0 ...
        ,'oneSided', 1 ...
        ,'setID',fidx ...
        ,'saveRes', 0 ...
        );
    options.dispEst = struct(...
        'method','Loupas'...
        ,'ref_type','Progressive' ...   % anchored/progressive
        ,'ref_idx',[] ...
        ,'nreverb', 2 ...         % Not having the correct nreverb could mess up displacements when using progressive ref_type
        ,'interpFactor', 5 ...
        ,'kernelLength', 4 ... % 10 -> 2.5 mm
        ,'DOF_fraction', 2 ... % Fraction of Depth of Field (around focus) to compute displacements for. DOF = 9*lambda*F_num^2
        ,'ccmode', 1 ...
        );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Paths
if ispc
    addpath(genpath('C:\Users\vrk4\Documents\GiHub\SC2000\arfiProcCode\'))
    addpath(genpath('C:\Users\vrk4\Documents\GitHub\TTEProcCode'))
elseif isunix
    addpath(genpath('/emfd/vrk4/GitHub/SC2000/arfiProcCode'))
    addpath(genpath('/emfd/vrk4/GitHub/TTEProcCode'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract timeStamp
list = dir('arfi_par_*'); % get timeStamp based on existance of ARFI par files
if size(list,1)<options.dataflow.setID
    error('Data set index requested greater than number of data sets')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading in timestamp for data set
if length(options.dataflow.setID)==14
    timeStamp = options.dataflow.setID;
elseif options.dataflow.setID == -1
    timeStamp = list(end).name(end-17:end-4);
else
    timeStamp = list(options.dataflow.setID).name(end-17:end-4);
end

if fidx==-1
    fprintf('Loading data with timeStamp = %s (Set # %d of %d)\n', timeStamp,size(list,1),size(list,1));
else
    fprintf('Loading data with timeStamp = %s (Set # %d of %d)\n', timeStamp,fidx,size(list,1));
end

options.timeStamp = timeStamp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract B-mode Data
fprintf(1,'Extracting B-mode Data...\n');
[bdata,bmodeSave] = extractBmode(timeStamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract ARFI/SWEI Data
fprintf(1,'Extracting ARFI/SWEI Data...\n');
if options.dataflow.ARFI
    [arfidata,arfiSave,options] = extractMmode(timeStamp,options,'ARFI');
    arfi_par = load(sprintf('arfi_par_%s.mat',timeStamp));
    arfidata = interpPushReverb(arfidata,options,arfi_par,''); % Interpolates through push-reverb; so that LPF can to function
else
    arfidata = [];
    arfi_par = [];
end
if options.dataflow.SWEI
    [sweidata,sweiSave,options] = extractMmode(timeStamp,options,'SWEI');
    swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
    sweidata = interpPushReverb(sweidata,options,swei_par,''); % Interpolates through push-reverb; so that LPF can to function
else
    sweidata = [];
    swei_par = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract ECG/Trigger Data
if exist(strcat('ECG_data_',timeStamp,'.mat'),'file')
    fprintf(1,'Loading ECG/Trigger Data...\n');
    % Expected Length of Bmode, ARFI and SWEI acq respectively
    dt(1) = ceil(bdata.t(end));
    dt(2) = ceil(arfidata.acqTime(end)); 
    dt(3) = dt(2); % For now, this would need to be updated if ARFI/SWEI timing is independent
    
    ecgdata = extractECG(timeStamp,options.dataflow.ecg_test,dt);
else
    ecgdata = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save time stamped results file
if (options.dataflow.saveRes && ~options.dataflow.display)
    fprintf(1,'Saving Res files...\n');
    if options.dataflow.ARFI
        tic
        resfile = ['res_arfi_' timeStamp '.mat'];
        save(resfile,'bdata','ecgdata','arfidata','options','-v7.3');
        fprintf(1,'Save Time for ARFI = %2.2fs\n',toc)
    end
    if options.dataflow.SWEI
        tic
        resfile = ['res_swei_' timeStamp '.mat'];
        save(resfile,'bdata','ecgdata','sweidata','options','-v7.3');
        fprintf(1,'Save Time for SWEI = %2.2fs\n',toc)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display
if options.dataflow.display
    dispTTE(ecgdata,bdata,arfidata,arfi_par,sweidata,swei_par,options,timeStamp);
end

cd(cudir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
