function datastruct = motionFilter(datastruct,options)

% Check the math to see if order of filtering operations matters
% ie. Polynomial -> BPF vs. BPF -> Polynomial

% Reshape Displacements from 4D back to 3D in case of SWEI data
if size(size(datastruct.disp_off),2)==4
    [ax beam push tstep] = size(datastruct.disp_off);
    datastruct.disp_off = reshape(datastruct.disp_off,ax,beam*push,tstep);
    datastruct.disp_on = reshape(datastruct.disp_on,ax,beam*push,tstep);
    reshape_flag = 1;
else
    reshape_flag = 0;
end

nref = datastruct.tstep_type(1);
npush = datastruct.tstep_type(2);
ntrack = datastruct.tstep_type(3);
nensemble = datastruct.tstep_type(4);

if (strcmpi(options.motionFilter.method,'Polynomial') || strcmpi(options.motionFilter.method,'Both'))
    tmask = false(size(datastruct.trackTime));
    timeRange = options.motionFilter.timeRange;
    for i = 1:2:length(timeRange)
        tmask(datastruct.trackTime>timeRange(i) & datastruct.trackTime<timeRange(i+1)) = true;
    end
    tmask(nref+(1:npush)) = false;
%     tmask(par.nref+(1:par.npush*length(par.pushFocalDepth))) = false;
    [datastruct.mf_disp_off datastruct.motion_off] = linearmotionfilter(datastruct.disp_off,datastruct.trackTime,find(tmask),options.motionFilter.order);
    [datastruct.mf_disp_on datastruct.motion_on] = linearmotionfilter(datastruct.disp_on,datastruct.trackTime,find(tmask),options.motionFilter.order);
    
    % Cubic interpolate push and reverb
    if options.dispEst.ref_idx == -1
        tidx1 = [nref+[-1 0] nref+npush+[3:4]];
        tidx2 = [nref+[1:npush+2]];
    else
        tidx1 = [nref+[-1 0] nref+npush+[2:3]];
        tidx2 = [nref+[1:npush+1]];
    end
    [residtmp motion1] = linearmotionfilter(datastruct.disp_on,datastruct.trackTime,tidx1,3);
    datastruct.disp_on(:,:,tidx2) = motion1(:,:,tidx2);
    clear residtmp motion1;
    
    [residtmp motion1] = linearmotionfilter(datastruct.mf_disp_on,datastruct.trackTime,tidx1,3);
    datastruct.mf_disp_on(:,:,tidx2) = motion1(:,:,tidx2);
    clear residtmp motion1 tidx1 tidx2;
    
elseif (strcmpi(options.motionFilter.method,'BPF') || strcmpi(options.motionFilter.method,'Both'))
    error('Coming soon...')
else
    error('Motion Filter method not recognized or not supported')
end

if reshape_flag
    datastruct.disp_off = reshape(datastruct.disp_off,ax,beam,push,tstep);
    datastruct.mf_disp_off = reshape(datastruct.mf_disp_off,ax,beam,push,tstep);
    
    datastruct.disp_on = reshape(datastruct.disp_on,ax,beam,push,tstep);
    datastruct.mf_disp_on = reshape(datastruct.mf_disp_on,ax,beam,push,tstep);
end

end
