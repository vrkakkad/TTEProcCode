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
        'stream', 'RT' ... % RT(realTime) || review || cluster
        ,'display', 1 ...
        ,'ecg_test', 0 ...
        ,'ARFI', 1 ...
        ,'DF_ARFI', 0 ...
        ,'SWEI', 0 ...
        ,'oneSided', [] ...
        ,'setID',fidx ...
        );
    
    options.dispEst = struct(...
        'method','Loupas'...
        ,'ref_type','Progressive' ...   % anchored || progressive
        ,'di', 1 ... Used only for progressive, 1 is A-B, B-C, 2 is A-C, B-D)
        ,'ref_idx', [1] ...
        ,'dims', [nan nan 2 nan nan] ... % [nref npush nreverb ntrack(total) ensemble] **Not having the correct nreverb could mess up displacements when using progressive ref_type**
        ,'interpFlag', 1 ... % Interpolate between push and reverb frames (1 = interpolate, 0 = nan out)
        ,'interpFactor', 5 ...
        ,'kernelLength', 5 ... % 10 -> 2.4 mm (Fundamental) 1.9 mm (Harmonic)
        ,'DOF_fraction', [] ... % Fraction of Depth of Field (around focus) to compute displacements for. DOF = 9*lambda*F_num^2
        ,'ccmode', 1 ...
        );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Paths
[~,temp] = fileparts(pwd);

if ispc
    addpath(genpath('C:\Users\vrk4\Documents\GiHub\SC2000\arfiProcCode\'))
    addpath(genpath('C:\Users\vrk4\Documents\GitHub\TTEProcCode'))
    resPath = strcat('C:\Users\vrk4\Tools\Transthoracic_Clinical\resDataArchive\',temp);
elseif isunix
    addpath(genpath('/emfd/vrk4/GitHub/SC2000/arfiProcCode'))
    addpath(genpath('/emfd/vrk4/GitHub/TTEProcCode'))
    resPath = [];
end

if ~exist(resPath,'dir'); mkdir(resPath); end
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
    fprintf('>>>>> Loading data with timeStamp = %s (Set # %d of %d)\n', timeStamp,size(list,1),size(list,1));
else
    fprintf('>>>>> Loading data with timeStamp = %s (Set # %d of %d)\n', timeStamp,fidx,size(list,1));
end

options.timeStamp = timeStamp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract B-mode Data
fprintf(1,'>>>>> Extracting B-mode Data...\n');
bdata = extractBmode(timeStamp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract ARFI/SWEI Data
fprintf(1,'>>>>> Extracting ARFI/SWEI Data...\n');

if options.dataflow.ARFI
    arfi_par = load(sprintf('arfi_par_%s.mat',timeStamp));
    options.dispEst.dims(1) = sum(arfi_par.nref);
    options.dispEst.dims(2) = arfi_par.npush;
    % nreverb is already defined
    options.dispEst.dims(4) = sum(arfi_par.ntrack);
    options.dispEst.dims(5) = arfi_par.ensemble;
    if isempty(options.dispEst.ref_idx)
        options.dispEst.ref_idx = arfi_par.nref;
    end
    [arfidata,options] = extractMmode(timeStamp,options,'ARFI');
    arfidata.df_flag = 0;
    if options.dispEst.interpFlag
        arfidata = interpPushReverb(arfidata,options,''); % Interpolate through push-reverb frames
    end
    % Compute Jitter
    arfidata = computeJitter(arfidata,options);
else
    arfidata = [];
    arfi_par = [];
end

if options.dataflow.DF_ARFI
    df_arfi_par = load(sprintf('swei_par_%s.mat',timeStamp));
    options.dispEst.dims(1) = sum(df_arfi_par.nref);
    options.dispEst.dims(2) = df_arfi_par.npush;
    % nreverb is already defined
    options.dispEst.dims(4) = sum(df_arfi_par.ntrack);
    options.dispEst.dims(5) = df_arfi_par.ensemble;
    if isempty(options.dispEst.ref_idx)
        options.dispEst.ref_idx = df_arfi_par.nref;
    end
    [df_arfidata,options] = extractMmode(timeStamp,options,'DF_ARFI');
    df_arfidata.df_flag = 1;
    if options.dispEst.interpFlag
        df_arfidata = interpPushReverb(df_arfidata,options,''); % Interpolate through push-reverb frames
    end
    % Compute Jitter
    df_arfidata = computeJitter(df_arfidata,options);
else
    df_arfidata = [];
    df_arfi_par = [];
end
 
if options.dataflow.SWEI
    %%%%%%%%% Check order of initializing options.dispEst.dims and loading par file
    resName = strcat(resPath,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_sweiRes.mat');
    optName = strcat(resPath,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_options.mat');
    swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
    options.dispEst.dims(1) = sum(swei_par.nref);
    options.dispEst.dims(2) = swei_par.npush;
    % nreverb is already defined
    options.dispEst.dims(4) = sum(swei_par.ntrack);
    options.dispEst.dims(5) = swei_par.ensemble;
    if isempty(options.dispEst.ref_idx)
        options.dispEst.ref_idx = swei_par.nref;
    end
    if (exist(resName,'file') && exist(optName,'file'))
        orig = load(optName);
        if isequal(options.dispEst,orig.options.dispEst)
            fprintf(1,'>>> Loading existing resFile...\n');
            load(resName);
            clear orig
        else
            fprintf(1,'>>> dispEst options do not match existing resFile. Computing displacements...\n');
            save(optName,'options','-v6');
            [sweidata,options] = extractMmode(timeStamp,options,'SWEI');
            swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
            options.dispEst.dims(1) = sum(swei_par.nref);
            options.dispEst.dims(2) = swei_par.npush;
            % nreverb is already defined
            options.dispEst.dims(4) = sum(swei_par.ntrack);
            options.dispEst.dims(5) = swei_par.ensemble;
            sweidata = interpPushReverb(sweidata,options,''); % Interpolates through push-reverb frames
            % Save sweidata and options to resData folder
            fprintf(1,'Saving sweidata and options file...\n');
            save(resName,'sweidata','-v6');
        end
    else
        fprintf(1,'>>> resFile does not exist. Computing displacements...\n');
        [sweidata,options] = extractMmode(timeStamp,options,'SWEI');
        swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
        sweidata = interpPushReverb(sweidata,options,''); % Interpolates through push-reverb frames
        % Save sweidata and options to resData folder
        fprintf(1,'Saving sweidata and options file...\n');
        save(resName,'sweidata','-v6');
        save(optName,'options','-v6');
        
    end
else
    sweidata = [];
    swei_par = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract ECG/Trigger Data
if exist(strcat('ECG_data_',timeStamp,'.mat'),'file')
    fprintf(1,'>>>>> Loading ECG/Trigger Data...\n');
    % Expected Length of Bmode, ARFI and SWEI acq respectively
    dt(1) = ceil(bdata.t(end));
    if options.dataflow.ARFI; dt(2) = ceil((arfi_par.numBeamGroups*arfi_par.numAcq-1)/arfi_par.pushPRF); dt(3) = dt(2); end
    if options.dataflow.DF_ARFI; dt(3) = ceil((df_arfi_par.numBeamGroups*df_arfi_par.numAcq-1)/df_arfi_par.pushPRF); dt(2) = dt(3); end
    if options.dataflow.SWEI; dt(3) = ceil((swei_par.numBeamGroups*swei_par.numAcq-1)/swei_par.pushPRF); dt(2) = dt(3); end
    
    [bECG,aECG,sECG,hr] = extractECG(timeStamp,options.dataflow.ecg_test,dt);
else
    bECG = []; aECG = []; sECG = []; hr = [];
end

bdata.ecg = bECG; bdata.hr = hr;
if options.dataflow.ARFI; arfidata.ecg = aECG; arfidata.hr = hr; end
if options.dataflow.DF_ARFI; df_arfidata.ecg = sECG; df_arfidata.hr = hr; end
if options.dataflow.SWEI; sweidata.ecg = sECG; sweidata.hr = hr; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defocused ARFI from SWEI Acquisition
% if options.dataflow.SWEI
%     if options.dataflow.oneSided
%         center_idx = 1;
%     else
%         center_idx = (swei_par.nBeams+1)/2;
%     end
%     df_arfidata = sweidata;
%     df_arfidata.disp = squeeze(df_arfidata.disp(:,center_idx,:,:));
%     df_arfidata.cc = squeeze(df_arfidata.cc(:,center_idx,:,:));
%     df_arfidata.df_flag = 1;
%     df_arfi_par = load(sprintf('swei_par_%s.mat',timeStamp));
%     %     df_arfidata = interpPushReverb(df_arfidata,options,df_arfi_par,'nan'); % NaN out push-reverb frames
% else
%     df_arfidata = [];
%     df_arfi_par = [];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display
if options.dataflow.display
    dispTTE(bdata,arfidata,arfi_par,sweidata,swei_par,df_arfidata,df_arfi_par,options);
end
cd(cudir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%