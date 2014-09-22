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
    ,'ref_type','independent' ...   % independent/common/progressive - indicates whether "no push" and "push" data sets will have independent references or a common reference, or use a moving reference
    ,'ref_idx',[] ...
    ,'noverlap',5 ...        % DO NOT CHANGE - number of time steps that are common between "no push" and "push" data sets (determined in sequenceParams file)
    ,'interpFactor',5 ...
    ,'kernelLength',4 ...
    ,'ccmode', 1 ...
    );

options.motionFilter = struct(...
    'enable',1 ...
    ,'method','BPF' ... % Polynomial/BPF/Both
    ... % Parameters for Polynomial filter
    ,'order',1 ... 
    ,'timeRange',[-inf -0.32] ... 
    ... % Parameters for Bandpass filter
    ,'passBand',[50 2000] ...
    );

%% Extract timeStamp
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

% Read in ARFI par file
arfi_par = load(sprintf('arfi_par_%s.mat',timeStamp));


fprintf('Loading data with timeStamp = %s\n', timeStamp);
%% Extract B-mode Data
fprintf(1,'Extracting B-mode Data...\n');
[bdata,bmodeSave] = extractBmode(timeStamp);

% dof = 7.22*1.540/arfi_par.pushFreq*(arfi_par.pushFnum)^2;
% edge = (arfi_par.pushFocalDepth + [-dof/2 dof/2])/10;
% 
% for i=1:81;
%     imagesc(bdata.blat,bdata.bax,bdata.bimg(:,:,i));
%     colormap(gray);axis image;
%     title(i);
%     hold on
%     plot(bdata.blat,edge(1)*ones(length(bdata.blat)),'b')
%     plot(bdata.blat,edge(2)*ones(length(bdata.blat)),'b')
%     plot(-0.5*ones(length(bdata.bax)),bdata.bax,'g')
%     plot(0.5*ones(length(bdata.bax)),bdata.bax,'g')
%     hold off
%     pause
% end
% 
% keyboard

%% Extract ARFI/SWEI Data
fprintf(1,'Extracting ARFI/SWEI Data...\n');
if options.dataflow.ARFI
    [arfidata,arfiSave,options] = extractMmode(timeStamp,options,'ARFI');
    arfi_par = load(sprintf('arfi_par_%s.mat',timeStamp));
    if options.motionFilter.enable
        [arfidata] = motionFilter(arfidata,options,arfi_par);
    end
end

if options.dataflow.SWEI
    [sweidata,sweiSave,options] = extractMmode(timeStamp,options,'SWEI');
    swei_par = load(sprintf('swei_par_%s.mat',timeStamp));
    if options.motionFilter.enable
        [sweidata] = motionFilter(sweidata,options,swei_par);
    end
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
