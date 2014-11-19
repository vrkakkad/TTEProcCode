function dispTTE(ecgdata,bdata,arfidata,arfi_par,sweidata,swei_par,options)

% Process Raw Displacements and Display

options.motionFilter = struct(...
    'enable',1 ...
    ,'method','Both' ... % Polynomial/LPF/Both
    ... % Parameters for Polynomial filter
    ,'order',1 ...
    ,'timeRange_pre',[-8 -7 -1 -0.5] ...
    ,'timeRange_push',[-1 -0.5 6 7] ...
    ... % Parameters for Bandpass filter
    ,'LPF_Cutoff',1000 ...
    );

options.display = struct(...
    'IQrange',[-50 0] ...
    ,'gateWidth',5 ...
    ,'gateOffset',0 ...
    ,'medfilt',[1 0.15] ... % median filter parameters - [axial (mm) acqTime (s)]
    ,'cc_filt',1 ...
    ,'cc_thresh', 0.999 ...
    ... % ARFI Display Parameters
    ,'disprange',[-5 20] ...
    ,'normalize',0 ...
    ,'t_disp_pre',-4 ...
    ,'t_disp_push',0.5 ...
    ,'IQtraces',0 ...
    ... % SWEI Display Parameters
    ,'velrange',[-2 10] ...
    ,'axial_scan',0 ...
    ,'sw_movie',0 ...
    ,'dvt_plots',0 ...
    ,'sw_display','disp' ... % Display displacements ('disp') or velocity ('vel') data
    );

options.calcSWS = struct(...
    'enable',1 ...
    ,'method','LinReg' ... LinReg/ LatSum
    ,'metric', 'TTP' ...
    ,'r2_threshold',0.9 ...
    ,'SWSrange',[0 6] ...
    );
 
if options.dataflow.ARFI
    if options.motionFilter.enable
        arfidata_mf_pre = motionFilter(arfidata,options,arfi_par,'pre');
        arfidata_mf_push = motionFilter(arfidata,options,arfi_par,'push');
    else
        arfidata_mf_pre = [];
        arfidata_mf_push = [];
    end
    dispARFI(ecgdata,bdata,arfidata,arfidata_mf_pre,arfidata_mf_push,options,arfi_par);
end

if options.dataflow.SWEI
    if options.motionFilter.enable
        sweidata_mf_pre = motionFilter(sweidata,options,swei_par,'pre');
        sweidata_mf_push = motionFilter(sweidata,options,swei_par,'push');
    else
        sweidata_mf_pre = [];
        sweidata_mf_push = [];
    end
    dispSWEI(ecgdata,bdata,sweidata,sweidata_mf_pre,sweidata_mf_push,options,swei_par);
end