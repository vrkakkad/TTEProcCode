function dispSWEI(ecgdata,bdata,sweidata,sweidata_mf_pre,sweidata_mf_push,options,par)

% Not currently compatible with ref_type = 'independent'

dof = 7.22*1.540/par.pushFreq*(par.pushFnum)^2;
edge = (par.pushFocalDepth + [-dof/2 dof/2]);

edge2 = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
if (edge2(1)<sweidata.axial(1) || edge2(2)>sweidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',edge2(1),edge2(2),sweidata.axial(1),sweidata.axial(end));
end
edge_idx = [find(sweidata.axial>edge2(1),1,'first') find(sweidata.axial<edge2(2),1,'last')];

% Display HQ Bmode Frames
figure
if isunix
    set(gcf,'Position',[1201 386 1920 1070])
elseif ispc
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end
for i=1:size(bdata.bimg,3);
    subplot('Position',[0.1 0.6 0.3 0.3])
    imagesc(bdata.blat,bdata.bax,bdata.bimg(:,:,i));
    colormap(gray);axis image; freezeColors;
    title(i);
    hold on
    rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2)
    rectangle('Position',[3 edge2(1) 3 edge2(2)-edge2(1)],'EdgeColor','r','Linewidth',2)
    hold off
    xlabel('Lateral (mm)','fontsize',10','fontweight','bold')
    ylabel('Axial (mm)','fontsize',10','fontweight','bold')
    title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)\n',i,bdata.t(i)),'fontsize',10','fontweight','bold')
    pause(0.025)
    
    if i==1
        % Display M-mode Data
        subplot('Position',[0.5 0.7 0.4 0.2])
        imagesc(sweidata.acqTime,sweidata.IQaxial,abs(db(sweidata.IQ(:,ceil(par.nBeams/2):par.nBeams:end,1))))
        hold on
        plot(sweidata.acqTime,edge2(1)*ones(length(sweidata.acqTime)),'r','Linewidth',2)
        plot(sweidata.acqTime,edge2(2)*ones(length(sweidata.acqTime)),'r','Linewidth',2)
        plot(sweidata.acqTime,sweidata.axial(1)*ones(length(sweidata.acqTime)),'b','Linewidth',2)
        plot(sweidata.acqTime,sweidata.axial(end)*ones(length(sweidata.acqTime)),'b','Linewidth',2)
        hold off
        ylabel('Axial (mm)','fontsize',10','fontweight','bold')
        title('M-Mode Frames','fontsize',10','fontweight','bold')
        colormap(gray); freezeColors;
        %         grid on
    end
end

% Compute velocity data from displacement data
[vel_temp, sweidata.t_vel] = differentiateDisplacements(permute(sweidata.disp(:,2:end,:,:),[1 2 4 3]),sweidata.trackTime,options.motionFilter.LPF_Cutoff);
sweidata.vel = permute(vel_temp,[1 2 4 3]);
sweidata = interpPushReverb(sweidata,options,par,'nan'); % NaN out push and reverb frames
clear vel_temp
if options.motionFilter.enable
    [vel_temp, sweidata_mf_pre.t_vel] = differentiateDisplacements(permute(sweidata_mf_pre.disp(:,2:end,:,:),[1 2 4 3]),sweidata_mf_pre.trackTime,options.motionFilter.LPF_Cutoff);
    sweidata_mf_pre.vel = permute(vel_temp,[1 2 4 3]);
    sweidata_mf_pre = interpPushReverb(sweidata_mf_pre,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
    [vel_temp, sweidata_mf_push.t_vel] = differentiateDisplacements(permute(sweidata_mf_push.disp(:,2:end,:,:),[1 2 4 3]),sweidata_mf_push.trackTime,options.motionFilter.LPF_Cutoff);
    sweidata_mf_push.vel = permute(vel_temp,[1 2 4 3]);
    sweidata_mf_push = interpPushReverb(sweidata_mf_push,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
end

%% Plot Shear Wave Movies
if options.display.sw_movie
    figure
    if isunix
        set(gcf,'Position',[1002 16 423 526])
    elseif ispc
        set(gcf,'Position',[1002 16 423 526])
    end
    for i=1:size(sweidata.vel,3)
        for j=par.nref:size(sweidata.vel,4)
            subplot(221)
            imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.disp(:,2:end,i,j),options.display.disprange);axis image; title(['Disp: Push #',num2str(i)])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            if options.motionFilter.enable
                subplot(222)
                imagesc(sweidata.lat(end,:),sweidata.axial,sweidata_mf_push.disp(:,2:end,i,j),options.display.disprange);axis image; title(['Disp MF: Push #',num2str(i)])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
            end
            
            subplot(223)
            imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.vel(:,:,i,j),options.display.velrange);axis image; title(['Vel: Push #',num2str(i)])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            if options.motionFilter.enable
                subplot(224)
                imagesc(sweidata.lat(end,:),sweidata.axial,sweidata_mf_push.vel(:,:,i,j),options.display.velrange);axis image; title(['Vel MF: Push #',num2str(i)])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
            end
            pause(0.001)
        end
        pause(0.001)
    end
    close(gcf)
%     return
end

%% Plot disp vs. time sw traces
if options.display.dvt_plots
    line_idx = 1+ [1 7 13];
    raw = squeeze(sweidata.disp(ceil(linspace(edge_idx(1),edge_idx(2),5)),line_idx,:,par.nref:end));
    if options.motionFilter.enable
        mf = squeeze(sweidata_mf_push.disp(ceil(linspace(edge_idx(1),edge_idx(2),5)),line_idx,:,par.nref:end));
    else
        mf = nan(size(raw));
    end
    figure
    for p = 1:size(raw,3)
        a1 = subplot(121);hold on;grid on;title(['Raw: ',num2str(p)])
        a2 = subplot(122);hold on;grid on;title(['MF: ',num2str(p)])
        for i=1:5
            plot(sweidata.trackTime(par.nref:end),-10*(i-1) + squeeze(raw(i,:,p,:))','Parent',a1)
            plot(sweidata.trackTime(par.nref:end),-10*(i-1) + squeeze(mf(i,:,p,:))','Parent',a2)
        end
        pause
        cla(a1);cla(a2);
    end
    close(gcf)
end


% Compute Axially Averaged Data
switch options.display.sw_display
    case 'disp'
        raw = squeeze(mean(sweidata.disp(edge_idx(1):edge_idx(2),2:end,:,:),1));
        raw = permute(raw, [1 3 2]);
        if options.motionFilter.enable
            mf = squeeze(mean(sweidata_mf_push.disp(edge_idx(1):edge_idx(2),2:end,:,:),1));
            mf = permute(mf, [1 3 2]);
        else
            mf = nan(size(raw));
        end
        rng = options.display.disprange;
    case 'vel'
        raw = squeeze(mean(sweidata.vel(edge_idx(1):edge_idx(2),:,:,:),1));
        raw = permute(raw, [1 3 2]);
        if options.motionFilter.enable
            mf = squeeze(mean(sweidata_mf_push.vel(edge_idx(1):edge_idx(2),:,:,:),1));
            mf = permute(mf, [1 3 2]);
        else
            mf = nan(size(raw));
        end
        rng = options.display.velrange;
end

% Calculate Shear Wave Speed
if options.calcSWS.enable
    switch options.calcSWS.SWSmethod
        case 'TTP'
            if options.motionFilter.enable
                [pk_pre, tpk_pre] = subsamplepeak(sweidata_mf_pre.trackTime(1:par.nref),sweidata_mf_pre.disp(:,2:end,:,1:par.nref),4);
                [pk_push, tpk_push] = subsamplepeak(sweidata_mf_push.trackTime(par.nref+par.npush+par.nreverb+1:end),sweidata_mf_push.disp(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
            else
                [pk_pre, tpk_pre] = subsamplepeak(sweidata.trackTime(1:par.nref),sweidata.disp(:,2:end,:,1:par.nref),4);
                [pk_push, tpk_push] = subsamplepeak(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),sweidata.disp(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
            end
        case 'TTPS'
            if options.motionFilter.enable
                [pk_pre, tpk_pre] = subsamplepeak(sweidata_mf_pre.t_vel(1:par.nref),sweidata_mf_pre.vel(:,:,:,1:par.nref),4);
                [pk_push, tpk_push] = subsamplepeak(sweidata_mf_push.t_vel(par.nref+par.npush+par.nreverb+1:end),sweidata_mf_push.vel(:,:,:,par.nref+par.npush+par.nreverb+1:end),4);
            else
                [pk_pre, tpk_pre] = subsamplepeak(sweidata.t_vel(1:par.nref),sweidata.vel(:,:,:,1:par.nref),4);
                [pk_push, tpk_push] = subsamplepeak(sweidata.t_vel(par.nref+par.npush+par.nreverb+1:end),sweidata.vel(:,:,:,par.nref+par.npush+par.nreverb+1:end),4);
            end
    end
    
    %     figure;set(gcf,'Position',[1070 50 361 345]);for i=1:size(tpk_push,3);imagesc(sweidata.lat(end,:),sweidata.axial,tpk_push(:,:,i),[0 7]);
    %     colorbar;xlabel('Lateral (mm)');ylabel('Axial (mm)');title(['TPK Map: Push #',num2str(i)]);axis image;pause;end
    
    klen = 13;
    dxdi= linreg(sweidata.lat,klen,2);
    dxdt = nan(size(tpk_push,1),size(tpk_push,3));
    for i=1:size(tpk_push,3)
        [dtdi_push, r2_push] = linreg(tpk_push(:,:,i),klen,2);
        temp =  dxdi./dtdi_push;
        dxdt_push(:,i) = temp(:,7).*(r2_push(:,7)>options.calcSWS.r2_threshold);
        [dtdi_pre, r2_pre] = linreg(tpk_pre(:,:,i),klen,2);
        temp =  dxdi./dtdi_pre;
        dxdt_pre(:,i) = temp(:,7).*(r2_pre(:,7)>options.calcSWS.r2_threshold);
    end
    
%     figure;plot(dxdt_push,'x');grid on;ylim([0 6])
    temp = dxdt_push;
    temp(temp==0) = nan;
    figure;hist(temp(:),100);title([nanmean(temp(:)), nanmedian(temp(:))])
    pause
    close(gcf)
%     keyboard
    % Indices corresponding to median filter parameters
    nax = double(ceil(options.display.medfilt(1)/(sweidata.axial(2) - sweidata.axial(1))));
    nt = double(ceil(options.display.medfilt(2)/(sweidata.acqTime(2) - sweidata.acqTime(1))));
    
    keyboard
    subplot('Position',[0.5 0.53 0.4 0.1]);
    imagesc(sweidata.acqTime,sweidata.axial,medfilt2(dxdt_pre,[nax nt]),options.calcSWS.SWSrange)
    ylabel('Axial (mm)','fontsize',10','fontweight','bold')
    title('SWS over DOF: Pre Push','fontsize',10','fontweight','bold')
    subplot('Position',[0.5 0.37 0.4 0.1]);
    imagesc(sweidata.acqTime,sweidata.axial,medfilt2(dxdt_push,[nax nt]),options.calcSWS.SWSrange)
    ylabel('Axial (mm)','fontsize',10','fontweight','bold')
    title('SWS over DOF: Post Push','fontsize',10','fontweight','bold')
    hcb = colorbar;
    set(hcb,'Position',[0.91 0.37 0.0187/2 0.26])
    ylabel(hcb,'SWS (m/s)','fontsize',10','fontweight','bold')
    colormap(hot)
    if isempty(ecgdata)
        xlabel('Acquisition Time (s)','fontsize',10','fontweight','bold')
    end
end

if ~isempty(ecgdata)
    samples = zeros(1,size(sweidata.disp,3));
    for i=1:size(sweidata.disp,3)
        samples(i) = ecgdata.swei(find(ecgdata.swei(:,1)>sweidata.acqTime(i),1,'first'),2);
    end
    ecgdata.swei(:,2) = ecgdata.swei(:,2)/max(ecgdata.swei(:,2));
    
    if options.calcSWS.enable
        h1 = subplot('Position',[0.5 0.1 0.4 0.2]);
    else
        h1 = subplot('Position',[0.5 0.25 0.4 0.2]);
    end
    plot(ecgdata.swei(:,1),ecgdata.swei(:,2),'Linewidth',2);
    hold on
    plot(sweidata.acqTime,samples,'kx','MarkerSize',8)
    pt = plot(sweidata.acqTime(1),samples(1),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    hold off
    grid on
    title('ECG Trace','fontsize',10','fontweight','bold')
    xlabel('Acquisition Time (s)','fontsize',10','fontweight','bold')
    axis tight
    hold(h1)
end

h2a = subplot('Position',[0.05 0.35 0.37 0.15]);
h2b = subplot('Position',[0.05 0.1 0.37 0.15]);

for rep =1:3
for i=1:size(sweidata.disp,3)
    if ~isempty(ecgdata)
        set(pt,'Visible','off')
    end
    colormap(hot)
    axes(h2a)
    imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(edge_idx)),:),raw(:,par.nref+1:end,i),rng/3)
    xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
    ylabel('Lateral (mm)','fontsize',10','fontweight','bold')
    grid on
    if strcmpi(options.display.sw_display,'disp')
        title(sprintf('SWEI Profiles averaged over %2.2f - %2.2f mm\n Raw Displacements: Push # %d (t = %2.2fs)',edge2(1),edge2(2),i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
    elseif strcmpi(options.display.sw_display,'vel')
        title(sprintf('SWEI Profiles averaged over %2.2f - %2.2f mm\n Raw Velocities: Push # %d (t = %2.2fs)',edge2(1),edge2(2),i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
    end
    axes(h2b)
    imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(edge_idx)),:),mf(:,par.nref+1:end,i),rng/3)
    xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
    ylabel('Lateral (mm)','fontsize',10','fontweight','bold')
    grid on
    
    hcb = colorbar;
    set(hcb,'Position',[0.43 0.1 0.0187/2 0.4])
    
    if strcmpi(options.display.sw_display,'disp')
        title(sprintf('Motion Filtered Displacements: Push # %d (t = %2.2f s)',i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
        ylabel(hcb,'Displacement (\mum)','fontsize',10','fontweight','bold')
    elseif strcmpi(options.display.sw_display,'vel')
        title(sprintf('Motion Filtered Velocities: Push # %d (t = %2.2f s)',i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
        ylabel(hcb,'Velocity (mm/s)','fontsize',10','fontweight','bold')
    end
    
    if ~isempty(ecgdata)
        pt = plot(sweidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    end
    
    %     frame = getframe(1);
    %     im = frame2im(frame);
    %     [imind,cm] = rgb2ind(im,256);
    %     if i == 1;
    %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    %     else
    %         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    %
    %     end
    pause(0.05)
end
end
