function runTTE_raw(DataDir,stream,fidx)
% close all
% TTE Wrapper for Raw Data
if ~exist('DataDir','var')
    DataDir = pwd;
end
if ~exist('stream','var')
    stream = 'RT';
end
if ~exist('fidx','var')
    fidx = -1;
end

% Sanity Checks
[~, hostname] = system('hostname');
hostname = deblank(hostname);

if (strcmpi(stream,'cluster') && (strcmpi(hostname,'patcuda1') || strcmpi(hostname,'rizzo')))
    stream = 'review'; 
    warning('Cluster mode not supported on %sReverting to review mode\n',hostname)
end

clear hostname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to procTTE
options.dataflow = struct(...
    'stream', stream ... % RT (realTime) || review || cluster 
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
    ,'ref_type','progressive' ...   % anchored || progressive
    ,'di', 1 ... Used only for progressive, 1 is A-B, B-C, 2 is A-C, B-D)
    ,'ref_idx', [1] ...
    ,'dims', [nan nan 2 nan nan] ... % [nref npush nreverb ntrack(total) ensemble] **Not having the correct nreverb could mess up displacements when using progressive ref_type**
    ,'interpFlag', 1 ... % Interpolate between push and reverb frames (1 = interpolate, 0 = nan out)
    ,'interpFactor', 5 ...
    ,'kernelLength', 5 ... % 10 -> 2.4 mm (Fundamental) 1.9 mm (Harmonic)
    ,'DOF_fraction', [4] ... % Fraction of Depth of Field (around focus) to compute displacements for. DOF = 9*lambda*F_num^2
    ,'ccmode', 1 ...
    );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to dispTTE
options.motionFilter = struct(...
    'method','Poly' ... % Off | LPF | BPF | LPF+Poly | BPF+Poly | Poly
    ... % Parameters for Bandpass filter
    ,'Cutoff', [25 750] ... % Must be two entries: LPF uses second entry as cutoff/ BPF uses entries as passband cutoffs
    ... % Parameters for Polynomial filter
    ,'order', 1 ...
    ,'timeRange_push', [-1.5 -1 4.5 5] ... % [-1.5 -1 4.5 5] [-5 -1 7.5 8.5]
    ,'pre_offset', -7 ... % [-6.5]
    );
options.motionFilter.timeRange_pre = options.motionFilter.timeRange_push + options.motionFilter.pre_offset;
% options.motionFilter.timeRange_pre = options.motionFilter.timeRange_push;

options.display = struct(...
    'theme', 'dark' ... % light/dark
    ,'playCine', 1 ...
    ,'IQrange',[-30 0] ...
    ,'gateWidth', 10 ...
    ,'gateOffset', 0 ...
    ,'n_pts', 3 ...
    ,'medfilt',[1 0.15] ... % median filter parameters - [axial (mm) acqTime (s)]
    ,'cc_filt', 0 ...
    ,'cc_thresh', 0.99 ...
    ... % ARFI Display Parameters
    ,'dispRange',[-1 10] ...
    ,'showPre', 1 ...
    ,'autoRange', 0 ...
    ,'normalize', 0 ...
    ,'t_disp_push', [] ... % [] = first track after reverb - last reference
    ,'t_disp_pre', [] ... % [] = auto calculated based on time range of Polynomial Filter
    ,'sysRange', [17.5 27.5]/100 ...
    ,'diaRange',[60 85]/100 ...
    ,'extras', 0 ... % [-1] - Suppress all extra plots || [0] - Asks for User Input on which plots to display
    ... % SWEI Display Parameters
    ,'velRange',[-5 5] ...
    ,'axial_scan',0 ...
    ,'sw_movie',1 ...
    ,'dvt_plots',0 ...
    ,'sw_display','disp' ... % Display displacements ('disp') or velocity ('vel') data
    );
options.calcSWS = struct(...
    'enable',0 ...
    ,'method','LinReg' ... LinReg/ LatSum
    ,'metric', 'TTP' ...
    ,'r2_threshold',0.5 ...
    ,'SWSrange',[0 7] ...
    );

%% Check input settings
switch options.dataflow.stream
    case 'RT' % realTime
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Real Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        options.dataflow.display = 1; fprintf(1,'\nDisplay: On');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dataflow.ARFI = 1; options.dataflow.SWEI = 0; fprintf(1,'\nARFI: On / SWEI: Off');
        options.dispEst.DOF_fraction = 2; fprintf(1,'\nDOF Fraction  = %2.1f',options.dispEst.DOF_fraction);
        options.display.extras = -1; fprintf(1,'\nExtras: Off [-1]');
        options.display.gateWidth = 20; fprintf(1,'\nGate Width = %d mm',options.display.gateWidth);
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
    case 'review'
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Review %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        options.dataflow.display = 1; fprintf(1,'\nDisplay: On');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dispEst.DOF_fraction = 1.5; fprintf(1,'\nDOF Fraction  = %2.1f',options.dispEst.DOF_fraction);
        options.dsplay.extras = 0; fprintf(1,'\nExtras: On [0]');
        options.display.gateWidth = 10; fprintf(1,'\nGate Width = %d mm',options.display.gateWidth);
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
    case 'cluster'
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%% Cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        options.dataflow.display = 0; fprintf(1,'\nDisplay: Off');
        options.dataflow.ecg_test = 0; fprintf(1,'\nECG Test: Off');
        options.dataflow.ARFI = 1; options.dataflow.SWEI = 1; fprintf(1,'\nARFI: On / SWEI: On');
        options.dispEst.DOF_fraction = 10; fprintf(1,'\nDOF Fraction  = %2.1f',options.dispEst.DOF_fraction);
        options.display.extras = -1; fprintf(1,'\nExtras: Off [-1]');
        fprintf(1,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); fprintf(1,'\n\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Call
procTTE(DataDir,fidx,options)
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%