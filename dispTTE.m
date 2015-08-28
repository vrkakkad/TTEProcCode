function dispTTE(bdata,arfidata,arfi_par,sweidata,swei_par,df_arfidata,df_arfi_par,options)
%% Default Input Parameters
if ~isfield(options,'motionFilter')
    options.motionFilter = struct(...
        'method','BPF' ... % Off | LPF | BPF | LPF+Poly | BPF+Poly | Poly
        ... % Parameters for Bandpass filter
        ,'Cutoff', [50 1000] ... % Must be two entries: LPF uses second entry as cutoff/ BPF uses entries as passband cutoffs
        ... % Parameters for Polynomial filter
        ,'order', 2 ...
        ,'timeRange_push', [-1.5 -1 4.5 5] ... % [-1.5 -1 4.5 5]
        ,'pre_offset', -6.5 ... % [-6.5]
        );
    options.motionFilter.timeRange_pre = options.motionFilter.timeRange_push + options.motionFilter.pre_offset;
end
if ~isfield(options,'display')
    options.display = struct(...
        'theme', 'dark' ... % light/dark
        ,'playCine', 0 ...
        ,'IQrange',[-40 0] ...
        ,'gateWidth', 2.5 ...
        ,'gateOffset', 0 ...
        ,'n_pts', 50 ...
        ,'medfilt',[2.5 0.25] ... % median filter parameters - [axial (mm) acqTime (s)]
        ,'cc_filt', 1 ...
        ,'cc_thresh', 0.975 ... % 30 db SNR (SNR = rho/(1-rho))
        ... % ARFI Display Parameters
        ,'disprange',[ ] ...
        ,'showPre', 1 ...
        ,'autoRange', 0 ...
        ,'normalize', 0 ...
        ,'t_disp_push', [] ... % [] = first track after reverb - last reference
        ,'t_disp_pre', [] ... % [] = auto calculated based on time range of Polynomial Filter
        ,'sysRange', [20 40]/100 ...
        ,'diaRange',[70 90]/100 ...
        ,'extras', 0 ...  % [-1] - Suppress all extra plots || [0] - Asks for User Input on which plots to display
        ... % SWEI Display Parameters
        ,'show_df_arfi',1 ...
        ,'show_swei',0 ...
        ,'velrange',[-5 15] ...
        ,'axial_scan',0 ...
        ,'sw_movie',0 ...
        ,'dvt_plots',0 ...
        ,'sw_display','disp' ... % Display displacements ('disp') or velocity ('vel') data
        );
end
if ~isfield(options,'calcSWS')
    options.calcSWS = struct(...
        'enable',0 ...
        ,'method','LinReg' ... LinReg/ LatSum
        ,'metric', 'TTP' ...
        ,'r2_threshold',0.5 ...
        ,'SWSrange',[0 7] ...
        );
end
%% Sequence specific updates
% try trackPRF = 1000/arfi_par.priusec(1); catch; trackPRF = 1000/swei_par.priusec(1); end % kHz
% if trackPRF<2.5
%     fprintf(1,'>>> Continuous track sequence: Changing Motion Filter Parameters...  \n');
%     options.motionFilter.timeRange_push = [-3 -1 20 22];
%     options.motionFilter.pre_offset = -23;
%     options.motionFilter.timeRange_pre = options.motionFilter.timeRange_push + options.motionFilter.pre_offset;
%     options.display.t_disp_push = 4;
%     options.display.t_disp_pre = options.motionFilter.timeRange_pre(1) + (options.display.t_disp_push - options.motionFilter.timeRange_push(1));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process and Display
if options.dataflow.ARFI
    if ~strcmpi(options.motionFilter.method,'Off')
        fprintf(1,'>>>>> Motion Filtering "pre" foc arfidata...\n');
        arfidata = motionFilter(arfidata,options,arfi_par,'pre');
        fprintf(1,'>>>>> Motion Filtering "push" foc arfidata...\n');
        arfidata = motionFilter(arfidata,options,arfi_par,'push');
        %         arfidata = interpPushReverb(arfidata,options,'nan'); % nan out push and reverb
    else
        arfidata.disp_mf_pre = [];
        arfidata.disp_mf_push = [];
    end
    dispARFI(bdata,arfidata,options,arfi_par);
end

if options.dataflow.DF_ARFI
    if ~strcmpi(options.motionFilter.method,'Off')
        fprintf(1,'>>>>> Motion Filtering "pre" defoc arfidata...\n');
        df_arfidata = motionFilter(df_arfidata,options,df_arfi_par,'pre');
        fprintf(1,'>>>>> Motion Filtering "push" defoc arfidata...\n');
        df_arfidata = motionFilter(df_arfidata,options,df_arfi_par,'push');
        %         df_arfidata = interpPushReverb(df_arfidata,options,'nan'); % nan out push and reverb
    else
        df_arfidata.disp_mf_pre = [];
        df_arfidata.disp_mf_push = [];
    end
    dispARFI(bdata,df_arfidata,options,df_arfi_par);
end

if options.dataflow.SWEI
    if ~strcmpi(options.motionFilter.method,'Off')
        
        fprintf(1,'>>>>> Motion Filtering "pre" sweidata...\n');
        sweidata = motionFilter(sweidata,options,swei_par,'pre');
        fprintf(1,'>>>>> Motion Filtering "push" sweidata...\n');
        sweidata = motionFilter(sweidata,options,swei_par,'push');
    else
        
        sweidata.disp_mf_pre = [];
        sweidata.disp_mf_push = [];
    end
    
    dispSWEI(bdata,sweidata,options,swei_par);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%