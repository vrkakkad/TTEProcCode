<<<<<<< HEAD
function runTTE_raw(DataDir,fidx,stream_idx)

close all

% TTE Wrapper for Raw Data

if ~exist('DataDir','var')
    DataDir = pwd;
end
if ~exist('fidx','var')
    fidx = -1;
end
if ~exist('stream_idx','var')
    stream_idx = 1;
end

% Sanity Checks
[status, hostname] = system('hostname');
hostname = deblank(hostname);
if (stream_idx ==3 && strcmpi(hostname,'patcuda1'))
    stream_idx = 2; 
    str = sprintf('Cluster mode [3] not supported on %s. Reverting to review mode [2]\n',hostname);
    warning(str)
end
clear status hostname
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to procTTE
options.dataflow = struct(...
    'stream', stream_idx ... % [1] - realTime, [2] - review, [3] - cluster 
    ,'display', 0 ...
    ,'ecg_test', 0 ...
    ,'ARFI', 1 ...
    ,'SWEI', 0 ...
    ,'oneSided', 1 ...
    ,'setID',fidx ...
    ,'saveRes', 1 ...
    );
options.dispEst = struct(...
    'method','Loupas'...
    ,'ref_type','Progressive' ...   % anchored/progressive
    ,'ref_idx',[] ...
    ,'nreverb', 2 ...         % Not having the correct nreverb could mess up displacements when using progressive ref_type
    ,'interpFactor', 5 ...
    ,'kernelLength', 10 ... % 10 -> 2.5 mm
    ,'DOF_fraction', 2.5 ... % Fraction of Depth of Field (around focus) to compute displacements for. DOF = 9*lambda*F_num^2
    ,'ccmode', 1 ...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to dispTTE
options.motionFilter = struct(...
    'enable', 1 ...
    ,'method','Both' ... % Polynomial/LPF/Both
    ... % Parameters for Polynomial filter
    ,'order', 2 ...
    ,'timeRange_push', [-1.5 -1 4.5 5] ... % [-1.5 -1 4.5 5]
    ,'pre_offset', -6.5 ... % [-6.5]
    ... % Parameters for Bandpass filter
    ,'LPF_Cutoff', 750 ...
    );
options.motionFilter.timeRange_pre = options.motionFilter.timeRange_push + options.motionFilter.pre_offset;

options.display = struct(...
    'theme', 'dark' ... % light/dark
    ,'IQrange',[-40 0] ...
    ,'gateWidth', 10 ...
    ,'gateOffset', 0 ...
    ,'n_pts', 5 ...
    ,'medfilt',[1 0.15] ... % median filter parameters - [axial (mm) acqTime (s)]
    ,'cc_filt', 1 ...
    ,'cc_thresh', 0.995 ...
    ... % ARFI Display Parameters
    ,'disprange',[-2 10] ...
    ,'normalize', 0 ...
    ,'t_disp_push', 0.75 ...
    ,'extras', 0 ... % [-1] - Suppress all extra plots; [0] - Asks for User Input on which plots to display
    ... % SWEI Display Parameters
    ,'velrange',[ ] ...
    ,'axial_scan',0 ...
    ,'sw_movie',0 ...
    ,'dvt_plots',0 ...
    ,'sw_display','disp' ... % Display displacements ('disp') or velocity ('vel') data
    );
options.display.t_disp_pre = options.motionFilter.timeRange_pre(1) + (options.display.t_disp_push - options.motionFilter.timeRange_push(1));

options.calcSWS = struct(...
    'enable',0 ...
    ,'method','LinReg' ... LinReg/ LatSum
    ,'metric', 'TTP' ...
    ,'r2_threshold',0.5 ...
    ,'SWSrange',[0 7] ...
    );

% Check input settings
switch options.dataflow.stream
    
    case 1 % realTime
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Real Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        options.dataflow.display = 1; fprintf(1,'\nDisplay: On');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dataflow.ARFI = 1; options.dataflow.SWEI = 0; fprintf(1,'\nARFI: On / SWEI: Off');
        options.dataflow.saveRes = 0; fprintf(1,'\nSaveRes: Off');
        options.display.extras = -1; fprintf(1,'\nExtras: Off [-1]');
        options.display.gateWidth = 20; fprintf(1,'\nGate Width = %d mm',options.display.gateWidth);
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
        
    case 2 % review
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Review %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        options.dataflow.display = 1; fprintf(1,'\nDisplay: On');
        options.dataflow.ecg_test = 1; fprintf(1,'\nECG Test: On');
        options.dataflow.saveRes = 0; fprintf(1,'\nSaveRes: Off');
        options.display.extras = 0; fprintf(1,'\nExtras: On [0]');
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
        
    case 3 % cluster
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        options.dataflow.display = 0; fprintf(1,'\nDisplay: Off');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dataflow.ARFI = 1; options.dataflow.SWEI = 1; fprintf(1,'\nARFI: On / SWEI: On');
        options.dataflow.saveRes = 1; fprintf(1,'\nSaveRes: On');
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Call
procTTE(DataDir,fidx,options)

clear
clc
=======
function runTTE_raw(DataDir,fidx,stream_idx)

close all

% TTE Wrapper for Raw Data

if ~exist('DataDir','var')
    DataDir = pwd;
end
if ~exist('fidx','var')
    fidx = -1;
end
if ~exist('stream_idx','var')
    stream_idx = 1;
end

% Sanity Checks
[status, hostname] = system('hostname');
hostname = deblank(hostname);
if (stream_idx ==3 && strcmpi(hostname,'patcuda1'))
    stream_idx = 2; 
    str = sprintf('Cluster mode [3] not supported on %s. Reverting to review mode [2]\n',hostname);
    warning(str)
end
clear status hostname
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to procTTE
options.dataflow = struct(...
    'stream', stream_idx ... % [1] - realTime, [2] - review, [3] - cluster 
    ,'display', 0 ...
    ,'ecg_test', 0 ...
    ,'ARFI', 1 ...
    ,'SWEI', 0 ...
    ,'oneSided', 1 ...
    ,'setID',fidx ...
    ,'saveRes', 1 ...
    );
options.dispEst = struct(...
    'method','Loupas'...
    ,'ref_type','Progressive' ...   % anchored/progressive
    ,'ref_idx',[] ...
    ,'nreverb', 2 ...         % Not having the correct nreverb could mess up displacements when using progressive ref_type
    ,'interpFactor', 5 ...
    ,'kernelLength', 10 ... % 10 -> 2.5 mm
    ,'DOF_fraction', 2.5 ... % Fraction of Depth of Field (around focus) to compute displacements for. DOF = 9*lambda*F_num^2
    ,'ccmode', 1 ...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to dispTTE
options.motionFilter = struct(...
    'enable', 1 ...
    ,'method','Both' ... % Polynomial/LPF/Both
    ... % Parameters for Polynomial filter
    ,'order', 2 ...
    ,'timeRange_push', [-1.5 -1 4.5 5] ... % [-1.5 -1 4.5 5]
    ,'pre_offset', -6.5 ... % [-6.5]
    ... % Parameters for Bandpass filter
    ,'LPF_Cutoff', 750 ...
    );
options.motionFilter.timeRange_pre = options.motionFilter.timeRange_push + options.motionFilter.pre_offset;

options.display = struct(...
    'theme', 'dark' ... % light/dark
    ,'IQrange',[-40 0] ...
    ,'gateWidth', 10 ...
    ,'gateOffset', 0 ...
    ,'n_pts', 5 ...
    ,'medfilt',[1 0.15] ... % median filter parameters - [axial (mm) acqTime (s)]
    ,'cc_filt', 1 ...
    ,'cc_thresh', 0.995 ...
    ... % ARFI Display Parameters
    ,'disprange',[-2 10] ...
    ,'normalize', 0 ...
    ,'t_disp_push', 0.75 ...
    ,'extras', 0 ... % [-1] - Suppress all extra plots; [0] - Asks for User Input on which plots to display
    ... % SWEI Display Parameters
    ,'velrange',[ ] ...
    ,'axial_scan',0 ...
    ,'sw_movie',0 ...
    ,'dvt_plots',0 ...
    ,'sw_display','disp' ... % Display displacements ('disp') or velocity ('vel') data
    );
options.display.t_disp_pre = options.motionFilter.timeRange_pre(1) + (options.display.t_disp_push - options.motionFilter.timeRange_push(1));

options.calcSWS = struct(...
    'enable',0 ...
    ,'method','LinReg' ... LinReg/ LatSum
    ,'metric', 'TTP' ...
    ,'r2_threshold',0.5 ...
    ,'SWSrange',[0 7] ...
    );

% Check input settings
switch options.dataflow.stream
    
    case 1 % realTime
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Real Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        options.dataflow.display = 1; fprintf(1,'\nDisplay: On');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dataflow.ARFI = 1; options.dataflow.SWEI = 0; fprintf(1,'\nARFI: On / SWEI: Off');
        options.dataflow.saveRes = 0; fprintf(1,'\nSaveRes: Off');
        options.display.extras = -1; fprintf(1,'\nExtras: Off [-1]');
        options.display.gateWidth = 20; fprintf(1,'\nGate Width = %d mm',options.display.gateWidth);
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
        
    case 2 % review
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Review %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        options.dataflow.display = 1; fprintf(1,'\nDisplay: On');
        options.dataflow.ecg_test = 1; fprintf(1,'\nECG Test: On');
        options.dataflow.saveRes = 0; fprintf(1,'\nSaveRes: Off');
        options.display.extras = 0; fprintf(1,'\nExtras: On [0]');
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
        
    case 3 % cluster
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        options.dataflow.display = 0; fprintf(1,'\nDisplay: Off');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dataflow.ARFI = 1; options.dataflow.SWEI = 1; fprintf(1,'\nARFI: On / SWEI: On');
        options.dataflow.saveRes = 1; fprintf(1,'\nSaveRes: On');
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Call
procTTE(DataDir,fidx,options)

clear
clc
>>>>>>> afa77557451dedcbdaaf013db29dc3d36822037b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%