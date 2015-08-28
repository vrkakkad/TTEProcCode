function [datastruct] = interpPushReverb(datastruct,options,type)
% This function interpolates through the push and reveb frames of the
% displacement data.

if (strcmpi(options.dispEst.ref_type,'anchored') || strcmpi(options.dispEst.ref_type,'progressive'))
    if size(size(datastruct.disp),2)==4
        [ax beam push tstep] = size(datastruct.disp);
        datastruct.disp = reshape(datastruct.disp,ax,beam*push,tstep);
        if (isfield(datastruct,'disp_mf_pre') && ~isempty(datastruct.disp_mf_pre)), datastruct.disp_mf_pre = reshape(datastruct.disp_mf_pre,ax,beam*push,tstep); end
        if (isfield(datastruct,'disp_mf_push') && ~isempty(datastruct.disp_mf_push)), datastruct.disp_mf_push = reshape(datastruct.disp_mf_push,ax,beam*push,tstep); end
        if (isfield(datastruct,'vel') && ~isempty(datastruct.vel)), datastruct.vel = reshape(datastruct.vel,ax,beam*push,tstep-1); end
        if (isfield(datastruct,'vel_mf_pre') && ~isempty(datastruct.vel_mf_pre)), datastruct.vel_mf_pre = reshape(datastruct.vel_mf_pre,ax,beam*push,tstep-1); end
        if (isfield(datastruct,'vel_mf_push') && ~isempty(datastruct.vel_mf_push)), datastruct.vel_mf_push = reshape(datastruct.vel_mf_push,ax,beam*push,tstep-1); end
        reshape_flag = 1;
    else
        reshape_flag = 0;
    end
    if strcmpi(type,'nan')
        fprintf(1,'>>>>> Removing Frames: %s (npush = %d, nrev = %d)\n',num2str(options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)),options.dispEst.dims(2),options.dispEst.dims(3));
        datastruct.disp(:,:,options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)) = nan(size(datastruct.disp,1),size(datastruct.disp,2),options.dispEst.dims(2)+options.dispEst.dims(3));
        if (isfield(datastruct,'disp_mf_pre') && ~isempty(datastruct.disp_mf_pre))
            datastruct.disp_mf_pre(:,:,options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)) = nan(size(datastruct.disp_mf_pre,1),size(datastruct.disp_mf_pre,2),options.dispEst.dims(2)+options.dispEst.dims(3));
        end
        if (isfield(datastruct,'disp_mf_push') && ~isempty(datastruct.disp_mf_push))
            datastruct.disp_mf_push(:,:,options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)) = nan(size(datastruct.disp_mf_push,1),size(datastruct.disp_mf_push,2),options.dispEst.dims(2)+options.dispEst.dims(3));
        end
        if (isfield(datastruct,'vel') && ~isempty(datastruct.vel))
            datastruct.vel(:,:,options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)) = nan(size(datastruct.vel,1),size(datastruct.vel,2),options.dispEst.dims(2)+options.dispEst.dims(3));
        end
        if (isfield(datastruct,'vel_mf_pre') && ~isempty(datastruct.vel_mf_pre))
            datastruct.vel_mf_pre(:,:,options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)) = nan(size(datastruct.vel_mf_pre,1),size(datastruct.vel_mf_pre,2),options.dispEst.dims(2)+options.dispEst.dims(3));
        end
        if (isfield(datastruct,'vel_mf_push') && ~isempty(datastruct.vel_mf_push))
            datastruct.vel_mf_push(:,:,options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)) = nan(size(datastruct.vel_mf_push,1),size(datastruct.vel_mf_push,2),options.dispEst.dims(2)+options.dispEst.dims(3));
        end
    else
        fprintf(1,'>>>>> Interpolating through Frames: %s (npush = %d, nrev = %d)\n',num2str(options.dispEst.dims(1)+1:options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)),options.dispEst.dims(2),options.dispEst.dims(3));
        tidx1 = [options.dispEst.dims(1)+[-1:0] options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)+[1:2]];
        tidx2 = [options.dispEst.dims(1)+[1:options.dispEst.dims(2)+options.dispEst.dims(3)]];
            [residtmp motion1] = linearmotionfilter(datastruct.disp,datastruct.trackTime,tidx1,3);
        datastruct.disp(:,:,tidx2) = motion1(:,:,tidx2);
        if (isfield(datastruct,'disp_mf_pre') && ~isempty(datastruct.disp_mf_pre))
            [residtmp motion1] = linearmotionfilter(datastruct.disp_mf_pre,datastruct.trackTime,tidx1,3);
            datastruct.disp_mf_pre(:,:,tidx2) = motion1(:,:,tidx2);
        end
        if (isfield(datastruct,'disp_mf_push') && ~isempty(datastruct.disp_mf_push))
            [residtmp motion1] = linearmotionfilter(datastruct.disp_mf_push,datastruct.trackTime,tidx1,3);
            datastruct.disp_mf_push(:,:,tidx2) = motion1(:,:,tidx2);
        end
%         if (isfield(datastruct,'vel') && ~isempty(datastruct.vel))
%             [residtmp motion1] = linearmotionfilter(datastruct.vel,datastruct.t_vel,tidx1,3);
%             datastruct.vel(:,:,tidx2) = motion1(:,:,tidx2);
%         end
%         if (isfield(datastruct,'vel_mf_pre') && ~isempty(datastruct.vel_mf_pre))
%             [residtmp motion1] = linearmotionfilter(datastruct.vel_mf_pre,datastruct.t_vel,tidx1,3);
%             datastruct.vel_mf_pre(:,:,tidx2) = motion1(:,:,tidx2);
%         end
%         if (isfield(datastruct,'vel_mf_push') && ~isempty(datastruct.vel_mf_push))
%             [residtmp motion1] = linearmotionfilter(datastruct.vel_mf_push,datastruct.t_vel,tidx1,3);
%             datastruct.vel_mf_push(:,:,tidx2) = motion1(:,:,tidx2);
%         end
    end
    
    if reshape_flag
        datastruct.disp = reshape(datastruct.disp,ax,beam,push,tstep);
        if (isfield(datastruct,'disp_mf_pre') && ~isempty(datastruct.disp_mf_pre)),datastruct.disp_mf_pre = reshape(datastruct.disp_mf_pre,ax,beam,push,tstep); end
        if (isfield(datastruct,'disp_mf_push') && ~isempty(datastruct.disp_mf_push)),datastruct.disp_mf_push = reshape(datastruct.disp_mf_push,ax,beam,push,tstep); end
        if (isfield(datastruct,'vel') && ~isempty(datastruct.vel)),datastruct.vel = reshape(datastruct.vel,ax,beam,push,tstep-1);end
        if (isfield(datastruct,'vel_mf_pre') && ~isempty(datastruct.vel_mf_pre)),datastruct.vel_mf_pre = reshape(datastruct.vel_mf_pre,ax,beam,push,tstep-1);end
        if (isfield(datastruct,'vel_mf_push') && ~isempty(datastruct.vel_mf_push)),datastruct.vel_mf_push = reshape(datastruct.vel_mf_push,ax,beam,push,tstep-1);end
    end
else
    error('Reference Type not recognized or not supported')
end

