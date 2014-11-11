function datastruct = motionFilter(datastruct,options,par,type)

% Reshape Displacements from 4D back to 3D in case of SWEI data
if strcmpi(options.dispEst.ref_type,'independent')
    if size(size(datastruct.disp_off),2)==4
        [ax beam push tstep] = size(datastruct.disp_off);
        datastruct.disp_off = reshape(datastruct.disp_off,ax,beam*push,tstep);
        datastruct.disp_on = reshape(datastruct.disp_on,ax,beam*push,tstep);
        reshape_flag = 1;
    else
        reshape_flag = 0;
    end
else
    if size(size(datastruct.disp),2)==4
        [ax beam push tstep] = size(datastruct.disp);
        datastruct.disp = reshape(datastruct.disp,ax,beam*push,tstep);
        reshape_flag = 1;
    else
        reshape_flag = 0;
    end
end

if strcmpi(options.dispEst.ref_type,'independent')
    nref = options.dispEst.noverlap;
    npush = par.npush;
else
    nref = par.nref;
    npush = par.npush;
end


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
    
    if strcmpi(options.dispEst.ref_type,'independent')
        [datastruct.disp_off datastruct.motion_off] = linearmotionfilter(datastruct.disp_off,datastruct.trackTime,find(tmask),options.motionFilter.order);
        [datastruct.disp_on datastruct.motion_on] = linearmotionfilter(datastruct.disp_on,datastruct.trackTime,find(tmask),options.motionFilter.order);
    else
        [datastruct.disp datastruct.motion] = linearmotionfilter(datastruct.disp,datastruct.trackTime,find(tmask),options.motionFilter.order);
    end

%% interpPushReverb takes care of this part now; not required here

%     % Cubic interpolate push and reverb if only using Polynomial filter
%     % (BPF takes care of this as well)
%     if ~strcmpi(options.motionFilter.method,'Both')
%         if options.dispEst.ref_idx == -1
%             tidx1 = [nref+[-1 0] nref+npush+[3:4]];
%             tidx2 = [nref+[1:npush+2]];
%         else
%             tidx1 = [nref+[-1 0] nref+npush+[3:4]];
%             tidx2 = [nref+[1:npush+2]];
%         end
%         
%         if strcmpi(options.dispEst.ref_type,'independent')
%             [residtmp motion1] = linearmotionfilter(datastruct.disp_on,datastruct.trackTime,tidx1,3);
%             datastruct.disp_on(:,:,tidx2) = motion1(:,:,tidx2);
%         else
%             [residtmp motion1] = linearmotionfilter(datastruct.disp,datastruct.trackTime,tidx1,3);
%             datastruct.disp(:,:,tidx2) = motion1(:,:,tidx2);
%         end
%         clear residtmp motion1;
%     end
end


if (strcmpi(options.motionFilter.method,'LPF') || strcmpi(options.motionFilter.method,'Both'))
    
    fprintf(1,'Executing Lowpass filter: Cutoff set to %d Hz\n',options.motionFilter.LPF_Cutoff)
    
    if strcmpi(options.dispEst.ref_type,'independent')
        [datastruct.trackTime, datastruct.disp_off] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime,datastruct.disp_off,nref,par,options.motionFilter.LPF_Cutoff);
        [datastruct.trackTime, datastruct.disp_on] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime,datastruct.disp_on,nref,par,options.motionFilter.LPF_Cutoff);
    else
        [datastruct.trackTime, datastruct.disp] = filtArfiData_TTE(datastruct.axial,datastruct.trackTime,datastruct.disp,nref,par,options.motionFilter.LPF_Cutoff);
    end
end



if reshape_flag
    if strcmpi(options.dispEst.ref_type,'independent')
        datastruct.disp_off = reshape(datastruct.disp_off,ax,beam,push,tstep);
        datastruct.disp_on = reshape(datastruct.disp_on,ax,beam,push,tstep);
    else
        datastruct.disp = reshape(datastruct.disp,ax,beam,push,tstep);
    end
    
end
