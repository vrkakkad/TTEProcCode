function [datastruct] = interpPushReverb(datastruct,options,par,type)
% This function interpolates through the push and reveb frames of the
% displacement data.

if (strcmpi(options.dispEst.ref_type,'anchored') || strcmpi(options.dispEst.ref_type,'progressive'))
    if size(size(datastruct.disp),2)==4
        [ax beam push tstep] = size(datastruct.disp);
        datastruct.disp = reshape(datastruct.disp,ax,beam*push,tstep);
        if isfield(datastruct,'disp_mf_pre') datastruct.disp_mf_pre = reshape(datastruct.disp_mf_pre,ax,beam*push,tstep); end
        if isfield(datastruct,'disp_mf_push') datastruct.disp_mf_push = reshape(datastruct.disp_mf_push,ax,beam*push,tstep); end
        if isfield(datastruct,'vel') % May need to fix this later
            datastruct.vel = reshape(datastruct.vel,ax,beam*push,tstep-1);
        end
        reshape_flag = 1;
    else
        reshape_flag = 0;
    end
    
    nrev = options.dispEst.nreverb;
    
    fprintf(1,'Interpolating through Frames: %s (npush = %d, nrev = %d)\n',num2str(par.nref+1:par.nref+par.npush+nrev),par.npush,nrev);
    if strcmpi(type,'nan')
        datastruct.disp(:,:,par.nref+1:par.nref+par.npush+nrev) = nan(size(datastruct.disp,1),size(datastruct.disp,2),par.npush+nrev);
        if isfield(datastruct,'disp_mf_pre')
            datastruct.disp_mf_pre(:,:,par.nref+1:par.nref+par.npush+nrev) = nan(size(datastruct.disp_mf_pre,1),size(datastruct.disp_mf_pre,2),par.npush+nrev);
        end
        if isfield(datastruct,'disp_mf_push')
            datastruct.disp_mf_push(:,:,par.nref+1:par.nref+par.npush+nrev) = nan(size(datastruct.disp_mf_push,1),size(datastruct.disp_mf_push,2),par.npush+nrev);
        end
        if isfield(datastruct,'vel') % May need to fix this later
            datastruct.vel(:,:,par.nref+1:par.nref+par.npush+nrev) = nan(size(datastruct.vel,1),size(datastruct.vel,2),par.npush+nrev);
        end
    else
        tidx1 = [par.nref+[-1 0] par.nref+par.npush+nrev+[1:2]];
        tidx2 = [par.nref+[1:par.npush+nrev]];
        [residtmp motion1] = linearmotionfilter(datastruct.disp,datastruct.trackTime,tidx1,3);
        datastruct.disp(:,:,tidx2) = motion1(:,:,tidx2);
        if isfield(datastruct,'disp_mf_pre')
            [residtmp motion1] = linearmotionfilter(datastruct.disp_mf_pre,datastruct.trackTime,tidx1,3);
            datastruct.disp_mf_pre(:,:,tidx2) = motion1(:,:,tidx2);
        end
        if isfield(datastruct,'disp_mf_push')
            [residtmp motion1] = linearmotionfilter(datastruct.disp_mf_push,datastruct.trackTime,tidx1,3);
            datastruct.disp_mf_push(:,:,tidx2) = motion1(:,:,tidx2);
        end
        if isfield(datastruct,'vel') % May need to fix this later
            [residtmp motion1] = linearmotionfilter(datastruct.vel,datastruct.t_vel,tidx1,3);
            datastruct.vel(:,:,tidx2) = motion1(:,:,tidx2);
        end
    end
    
    if reshape_flag
        datastruct.disp = reshape(datastruct.disp,ax,beam,push,tstep);
        if isfield(datastruct,'disp_mf_pre') datastruct.disp_mf_pre = reshape(datastruct.disp_mf_pre,ax,beam,push,tstep); end
        if isfield(datastruct,'disp_mf_push') datastruct.disp_mf_push = reshape(datastruct.disp_mf_push,ax,beam,push,tstep); end
        if isfield(datastruct,'vel') % May need to fix this later
            datastruct.vel = reshape(datastruct.vel,ax,beam,push,tstep-1);
        end
    end
else
    error('Reference Type not recognized or not supported')
end

