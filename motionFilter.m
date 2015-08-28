function datastruct = motionFilter(datastruct,options,par,type)

disp = datastruct.disp;
% Reshape Displacements from 4D back to 3D in case of SWEI data
if (strcmpi(options.dispEst.ref_type,'anchored') || strcmpi(options.dispEst.ref_type,'progressive'))
    if size(size(disp),2)==4
        [ax,beam,push,tstep] = size(disp);
        disp = reshape(disp,ax,beam*push,tstep);
        reshape_flag = 1;
    else
        reshape_flag = 0;
    end
    nref = options.dispEst.dims(1);
    npush = options.dispEst.dims(2);
    nreverb = options.dispEst.dims(3);
else
    error('Reference Type not recognized or not supported')
end
% dims = size(datastruct.disp);
% test = 100*sin(2*pi*5000*datastruct.trackTime*1e-3);
% test = reshape(test,1,1,[]);
% test = repmat(test,[dims(1) dims(2) 1]);
% datastruct.disp = test;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(options.motionFilter.method,'PF')
%     test = nan(size(disp));
    if strfind(options.motionFilter.method,'LPF')
        fprintf(1,'Executing Lowpass filter: Cutoff set to %d Hz\n',options.motionFilter.Cutoff(2));
        filt_type = 'LPF';
    end
    if strfind(options.motionFilter.method,'BPF')
        fprintf(1,'Executing Bandpass filter: Cutoff set to [%d %d] Hz\n',options.motionFilter.Cutoff(1),options.motionFilter.Cutoff(2));
        filt_type = 'BPF';
    end
    if reshape_flag
        % SWEI Data is interpolated through push-reverb frames
        [~, disp] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime,disp,nref,par,options.motionFilter.Cutoff,filt_type);
    else
        if options.dispEst.interpFlag
            [~, disp] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime,disp,nref,par,options.motionFilter.Cutoff,filt_type);
        else
            % ARFI data has NaN vales at push-reverb frames and needs to be filtered separately
            [~, disp(:,:,1:nref-1)] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime(1:nref-1),disp(:,:,1:nref-1),nref,par,options.motionFilter.Cutoff,filt_type);
            [~, disp(:,:,nref+npush+nreverb+1:end)] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime(nref+npush+nreverb+1:end),disp(:,:,nref+npush+nreverb+1:end),nref,par,options.motionFilter.Cutoff,filt_type);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(options.motionFilter.method,'Poly')
    if strcmpi(type,'pre')
        fprintf(1,'Executing Polynomial filter: TimeRange set to [%s] us\n',num2str(options.motionFilter.timeRange_pre,'%2.2f '));
        timeRange = options.motionFilter.timeRange_pre;
    elseif strcmpi(type,'push')
        fprintf(1,'Executing Polynomial filter: TimeRange set to [%s] us\n',num2str(options.motionFilter.timeRange_push,'%2.2f '));
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
[disp,motion] = linearmotionfilter(disp,datastruct.trackTime,find(tmask),options.motionFilter.order);
%    [datastruct.disp datastruct.motion] = linearmotionfilter(datastruct.disp,datastruct.trackTime,find(tmask),options.motionFilter.order,datastruct.cc,options.display.cc_thresh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if reshape_flag
    if (strcmpi(options.dispEst.ref_type,'anchored') || strcmpi(options.dispEst.ref_type,'progressive'))
        disp = reshape(disp,ax,beam,push,tstep);
        if exist('motion','var'); motion = reshape(motion,ax,beam,push,tstep); end
    else
        error('Reference Type not recognized or not supported')
    end
    
end
if strcmpi(type,'pre')
    datastruct.disp_mf_pre = disp;
    if exist('motion','var'); datastruct.motion_pre = motion; end
elseif strcmpi(type,'push')
    datastruct.disp_mf_push = disp;
    if exist('motion','var'); datastruct.motion_push = motion; end
end