function datastruct = motionFilter(datastruct,options,par,type)

% Reshape Displacements from 4D back to 3D in case of SWEI data
if (strcmpi(options.dispEst.ref_type,'anchored') || strcmpi(options.dispEst.ref_type,'progressive'))
    if size(size(datastruct.disp),2)==4
        [ax beam push tstep] = size(datastruct.disp);
        datastruct.disp = reshape(datastruct.disp,ax,beam*push,tstep);
        reshape_flag = 1;
    else
        reshape_flag = 0;
    end
    nref = par.nref;
    npush = par.npush;
else
    error('Reference Type not recognized or not supported')
end

% dims = size(datastruct.disp);
% test = 100*sin(2*pi*5000*datastruct.trackTime*1e-3);
% test = reshape(test,1,1,[]);
% test = repmat(test,[dims(1) dims(2) 1]);
% datastruct.disp = test;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(options.motionFilter.method,'LPF') || strcmpi(options.motionFilter.method,'Both'))
    
    fprintf(1,'Executing Lowpass filter: Cutoff set to %d Hz\n',options.motionFilter.LPF_Cutoff)
    
    [datastruct.trackTime, datastruct.disp] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime,datastruct.disp,nref,par,options.motionFilter.LPF_Cutoff);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(options.motionFilter.method,'Polynomial') || strcmpi(options.motionFilter.method,'Both'))
    
    if strcmpi(type,'pre')
        fprintf(1,'Executing Polynomial filter: TimeRange set to [%s] us\n',num2str(options.motionFilter.timeRange_pre,'%2.2f '))
        timeRange = options.motionFilter.timeRange_pre;
    elseif strcmpi(type,'push')
        fprintf(1,'Executing Polynomial filter: TimeRange set to [%s] us\n',num2str(options.motionFilter.timeRange_push,'%2.2f '))
        timeRange = options.motionFilter.timeRange_push;
    else
        error('Wrong type into motionFilter')
    end
    
    tmask = false(size(datastruct.trackTime));
    for i = 1:2:length(timeRange)
        tmask(datastruct.trackTime>timeRange(i) & datastruct.trackTime<timeRange(i+1)) = true;
    end
    tmask(nref+(1:npush)) = false;
    tmask(nref+(1:npush*length(par.pushFocalDepth))) = false;
    
    [datastruct.disp datastruct.motion] = linearmotionfilter(datastruct.disp,datastruct.trackTime,find(tmask),options.motionFilter.order);
%    [datastruct.disp datastruct.motion] = linearmotionfilter(datastruct.disp,datastruct.trackTime,find(tmask),options.motionFilter.order,datastruct.cc,options.display.cc_thresh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if reshape_flag
    if (strcmpi(options.dispEst.ref_type,'anchored') || strcmpi(options.dispEst.ref_type,'progressive'))
        datastruct.disp = reshape(datastruct.disp,ax,beam,push,tstep);
    else
        error('Reference Type not recognized or not supported')
    end
    
end