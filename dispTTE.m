function dispTTE(ecgdata,bdata,arfidata,arfi_par,sweidata,swei_par,options)

% Process Raw Displacements and Display

options.motionFilter = struct(...
    'enable',1 ...
    ,'method','Both' ... % Polynomial/BPF/Both
    ... % Parameters for Polynomial filter
    ,'order',1 ...
    ,'timeRange',[-8 -7 6 7] ...
    ... % Parameters for Bandpass filter
    ,'passBand',[10 500] ...
    );


options.display = struct(...
    'gateWidth',2.5 ...
    ,'gateOffset',0 ...
    ,'disprange',[-5 40] ...
    ,'t_disp',0.75 ...
    ,'velrange',[-10 10] ...
    ,'calcSWS',1 ...
    ,'sw_display','vel' ... % Display displacements ('disp') or velocity ('vel') data
    );

if options.dataflow.ARFI
    if options.motionFilter.enable
        arfidata_mf = motionFilter(arfidata,options,arfi_par);
    else
        arfidata_mf = [];
    end
    dispARFI(ecgdata,bdata,arfidata,arfidata_mf,options,arfi_par);
end

if options.dataflow.SWEI
    if options.motionFilter.enable
        sweidata_mf = motionFilter(sweidata,options,swei_par);
    else
        sweidata_mf = [];
    end
    dispSWEI(ecgdata,bdata,sweidata,sweidata_mf,options,swei_par);
end