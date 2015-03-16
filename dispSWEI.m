<<<<<<< HEAD
function [traced_gate,trace_flag] = dispSWEI(ecgdata,bdata,sweidata,sweidata_mf_pre,sweidata_mf_push,options,par)

%% Check for existance of traced gate
if isfield(sweidata,'traced_gate') && ~isempty(sweidata.traced_gate)
    prior_trace = 1;
    traced_gate = sweidata.traced_gate;
else
    prior_trace = 0;
    traced_gate = [];
end

dims = size(sweidata.disp);
ndepth = dims(1); nacqT = dims(2); ntrackT = dims(3);

edge = [sweidata.axial(1) sweidata.axial(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Gate Parameters
if prior_trace
    gate = repmat(sweidata.traced_gate,1,2) + repmat(options.display.gateOffset,length(sweidata.traced_gate),2) + repmat([-options.display.gateWidth/2 options.display.gateWidth/2],length(sweidata.traced_gate),1);
elseif ~prior_trace
    gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
    gate = repmat(gate,[nacqT 1]);
    traced_gate = [];
    if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

if options.display.axial_scan
    ngates = round(2*(edge(2)-edge(1))/options.display.gateWidth);
    offsets = linspace(edge(1)+options.display.gateWidth/2,edge(2)-options.display.gateWidth/2,ngates) - par.pushFocalDepth;
    options.display.gateOffset = offsets(1);
end

% Set Display Parameters
figure(2)
if isunix
    set(gcf,'Position',[1201 386 1920 1070])
    dispPar.fsize = 16;
elseif ispc
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    dispPar.fsize = 8;
end

if strcmpi(options.display.theme,'light')
    dispPar.fig = [1 1 1]; dispPar.txt = 'k'; dispPar.ax = [0.5 0.5 0.5];
    for i=1:size(bdata.bimg,3)
        temp = bdata.bimg(:,:,i);
        temp(temp==bdata.bimg(1,1,i)) = 256;
        bdata.bimg(:,:,i) = temp;
    end
elseif strcmpi(options.display.theme,'dark')
    dispPar.fig = 'k'; dispPar.txt = 'w'; dispPar.ax = [0.25 0.25 0.25];
    for i=1:size(bdata.bimg,3)
        temp = bdata.bimg(:,:,i);
        temp(temp==bdata.bimg(1,1,i)) = 0;
        bdata.bimg(:,:,i) = temp;
    end
end

set(1,'Color',dispPar.fig)

                                        img_rng = options.display.disprange;
                                        plot_rng = options.display.disprange;

dispPar.cmap = colormap(parula);
dispPar.cmap(end,:) = dispPar.ax;

dispPar.corder = winter(options.display.n_pts);

    
for i=1:nacqT
    gate_idx(i,:) = [find(sweidata.axial>gate(i,1),1,'first') find(sweidata.axial<gate(i,2),1,'last')];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Bmode Cine
for i=1:size(bdata.bimg,3);
    p1 = subplot('Position',[0.1 0.6 0.3 0.3]);
    imagesc(bdata.blat,bdata.bax,bdata.bimg(:,:,i));
    colormap(gray);axis image; freezeColors;
    hold(p1,'on')
    r1 = rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2,'Parent',p1);
    if prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linestyle','--','Linewidth',2,'Parent',p1);
    elseif ~prior_trace
        r2 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:)))  options.display.gateWidth],'EdgeColor','g','Linewidth',2,'Parent',p1);
    end
    hold off
    xlabel('Lateral (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    title(sprintf('B-Mode Cine: Frame %d (t = %1.1f s)',i,bdata.t(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    set(gca,'xcolor',dispPar.fig,'ycolor',dispPar.fig,'fontweight','bold')
    
    %     xlim([-25 25]);ylim([max(edge(1)-7.5,sweidata.IQaxial(1)) min(edge(2)+7.5,sweidata.IQaxial(end))])
    
    pause(0.025)
    if i==size(bdata.bimg,3)
        % Display M-mode IQ
        p2 = axes('Position',[0.45 0.3 0.52 1]);
        frame = abs(sweidata.IQ(:,:,1)); % Display first frame only
        imagesc(linspace(0,range(sweidata.lat(:))*nacqT,size(sweidata.IQ,2)),sweidata.IQaxial,db(frame/max(frame(:))),options.display.IQrange)
        hold(p2,'on')
        l1 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(1)*ones(1,nacqT),'b','Linewidth',2,'Parent',p2)
        l2 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(end)*ones(1,nacqT),'b','Linewidth',2,'Parent',p2)
        l3 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2,'Parent',p2)
        l4 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2,'Parent',p2)
        hold(p2,'off')
        axis image
        ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        set(gca,'xTickLabel',[])
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d',par.isHarmonic),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylim([max(edge(1)-7.5,sweidata.IQaxial(1)) min(edge(2)+7.5,sweidata.IQaxial(end))])
        colormap(gray); freezeColors;
        set(gca,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trace out center of depth gate
if (isfield(sweidata,'traced_gate') && isempty(options.display.gateOffset))
    gate = sweidata.traced_gate;
    hold on
    l1 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2);
    l2 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2);
    r1 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) max(gate(:))-min(gate(:))],'EdgeColor','g','LineStyle','--','Linewidth',2,'Parent',p1);
    trace_flag = input('Previously traced gate exists, do you want to retrace the gate? (y/n) [n]: ','s');
    if ~strcmpi(trace_flag,'y')
        trace_flag = 'n';
    else
        delete(l1,l2,r1)
    end
elseif (~isfield(sweidata,'traced_gate') && isempty(options.display.gateOffset))
    trace_flag = 'y';
end
if ((isempty(options.display.gateOffset) && ~isfield(sweidata,'traced_gate')) || strcmpi(trace_flag,'y'))
    hold on
    fprintf(1,'Ready to trace gate (GateWidth = %2.1f mm)...\n',options.display.gateWidth)
    for i=1:nacqT
        [x,y] = ginput(1);
        mark(i) = plot(x,y,'yo','MarkerSize',10);
        gate(i,:) = (y + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        clear x y
    end
    delete(mark)
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2)
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2)
    rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) max(gate(:))-min(gate(:))],'EdgeColor','g','LineStyle','--','Linewidth',2,'Parent',p1)
    fprintf(1,'Gate Traced.\n')
elseif ((isempty(options.display.gateOffset) && ~isfield(sweidata,'traced_gate')) && strcmpi(trace,'n'))
    gate = repmat(sweidata.traced_gate',1,2) + repmat([-options.display.gateWidth/2 options.display.gateWidth/2],length(sweidata.traced_gate),1);
    hold on
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2)
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2)
    rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) max(gate(:))-min(gate(:))],'EdgeColor','g','LineStyle','--','Linewidth',2,'Parent',p1)
end

traced_gate = mean(gate');

if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
end

for i=1:nacqT
    gate_idx(i,:) = [find(sweidata.axial>gate(i,1),1,'first') find(sweidata.axial<gate(i,2),1,'last')];
end

if ~isempty(ecgdata)
    samples = zeros(1,size(sweidata.disp,3));
    for i=1:size(sweidata.disp,3)
        samples(i) = ecgdata.swei(find(ecgdata.swei(:,1)>sweidata.acqTime(i),1,'first'),2);
    end
    ecgdata.swei(:,2) = ecgdata.swei(:,2)/max(ecgdata.swei(:,2));
    
    if options.calcSWS.enable
        h1 = axes('Position',[0.5 0.1 0.4 0.2]);
    else
        h1 = axes('Position',[0.5 0.25 0.4 0.2]);
    end
    plot(ecgdata.swei(:,1),ecgdata.swei(:,2),'Linewidth',2);
    hold on
    plot(sweidata.acqTime,samples,'kx','MarkerSize',8)
    pt = plot(sweidata.acqTime(1),samples(1),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    hold off
    grid on
    title('ECG Trace','fontsize',fsize','fontweight','bold')
    xlabel('Acquisition Time (s)','fontsize',fsize','fontweight','bold')
    axis tight
    hold(h1)
end

% Compute velocity data from displacement data
[vel_temp, sweidata.t_vel] = differentiateDisplacements(permute(sweidata.disp,[1 2 4 3]),sweidata.trackTime,options.motionFilter.LPF_Cutoff);
sweidata.vel = 1*permute(vel_temp,[1 2 4 3]);
sweidata = interpPushReverb(sweidata,options,par,'nan'); % NaN out push and reverb frames
clear vel_temp
if options.motionFilter.enable
    [vel_temp, sweidata_mf_pre.t_vel] = differentiateDisplacements(permute(sweidata_mf_pre.disp,[1 2 4 3]),sweidata_mf_pre.trackTime,options.motionFilter.LPF_Cutoff);
    sweidata_mf_pre.vel = 1*permute(vel_temp,[1 2 4 3]);
    sweidata_mf_pre = interpPushReverb(sweidata_mf_pre,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
    [vel_temp, sweidata_mf_push.t_vel] = differentiateDisplacements(permute(sweidata_mf_push.disp,[1 2 4 3]),sweidata_mf_push.trackTime,options.motionFilter.LPF_Cutoff);
    sweidata_mf_push.vel = 1*permute(vel_temp,[1 2 4 3]);
    sweidata_mf_push = interpPushReverb(sweidata_mf_push,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
end

% Coorelation mask filter
if options.display.cc_filt
    mask = sweidata.cc>options.display.cc_thresh;
    sweidata.disp(mask==0) = nan;
    sweidata.vel(mask(:,:,:,1:end-1)==0) = nan;
    if options.motionFilter.enable
        sweidata_mf_pre.disp(mask==0) = nan;
        sweidata_mf_pre.vel(mask(:,:,:,1:end-1)==0) = nan;
        sweidata_mf_push.disp(mask==0) = nan;
        sweidata_mf_push.vel(mask(:,:,:,1:end-1)==0) = nan;
    end
end

%% Plot Axial Scan Data
if options.display.axial_scan
    for i=1:ngates
        set(0,'currentFigure',f1)
        options.display.gateOffset = offsets(i);
        gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        gate = repmat(gate,[nacqT 1]);
        if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
            warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
        end
        
        for j=1:nacqT
            gate_idx(j,:) = [find(sweidata.axial>gate(j,1),1,'first') find(sweidata.axial<gate(j,2),1,'last')];
        end
        
        delete(r1)
        r1 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) options.display.gateWidth ],'EdgeColor','g','Linewidth',2,'Parent',p1);
        title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)\n Gate Offset = %2.2f mm ',size(bdata.bimg,3),bdata.t(size(bdata.bimg,3)),offsets(i)),'fontsize',fsize','fontweight','bold','Parent',p1)
        p2 = axes('Position',[0.45 0.3 0.52 1]);
        frame = sweidata.IQ(:,:,1);
        imagesc(linspace(0,range(sweidata.lat(:))*nacqT,size(sweidata.IQ,2)),sweidata.IQaxial,db(frame/max(frame(:))),options.display.IQrange)
        hold on
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(1)*ones(1,nacqT),'b','Linewidth',2)
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(end)*ones(1,nacqT),'b','Linewidth',2)
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2)
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2)
        hold off
        axis image
        ylabel('Axial (mm)','fontsize',fsize','fontweight','bold')
        set(gca,'xTickLabel',[])
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d',par.isHarmonic),'fontsize',fsize,'fontweight','bold')
        ylim([max(edge(1)-15,sweidata.IQaxial(1)) min(edge(2)+15,sweidata.IQaxial(end))])
        
        % Compute Axially Averaged Data
        switch options.display.sw_display
            case 'disp'
                raw = squeeze(nanmean(sweidata.disp(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                raw = permute(raw, [1 3 2]);
                if options.motionFilter.enable
                    mf = squeeze(nanmean(sweidata_mf_push.disp(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                    mf = permute(mf, [1 3 2]);
                else
                    mf = nan(size(raw));
                end
                rng = options.display.disprange*0.5;
            case 'vel'
                raw = squeeze(nanmean(sweidata.vel(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                raw = permute(raw, [1 3 2]);
                if options.motionFilter.enable
                    mf = squeeze(nanmean(sweidata_mf_push.vel(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                    mf = permute(mf, [1 3 2]);
                else
                    mf = nan(size(raw));
                end
                rng = options.display.velrange;
        end
        
        figure(101)
        if isunix
            set(gcf,'Position',[1200 409 950 1050])
        elseif ispc
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
        end
        for i=1:nacqT
            subplot(nacqT/5,5,i)
            p = imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx(1,:))),:),raw(:,(par.nref+1:end),i),rng);
            set(p,'alphadata',~isnan(raw(:,(par.nref+1:end),i)))
            set(gca,'color',[0.4 0.4 0.4])
            title(sweidata.acqTime(i),'fontsize',fsize','fontweight','bold')
            xlim([0 7])
        end
        
        figure(102)
        if isunix
            set(gcf,'Position',[2150 409 950 1050])
        elseif ispc
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
        end
        for i=1:nacqT
            subplot(nacqT/5,5,i)
            p = imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx(1,:))),:),mf(:,(par.nref+1:end),i),rng);
            set(p,'alphadata',~isnan(raw(:,(par.nref+1:end),i)))
            set(gca,'color',[0.4 0.4 0.4])
            title(sweidata.acqTime(i),'fontsize',fsize','fontweight','bold')
            xlim([0 7])
        end
        pause
    end
    close 101 102
    clc
end

%% Plot Shear Wave Movies
if options.display.sw_movie
    figure(101)
    if isunix
        set(gcf,'Position',[69 435 414 906])
    elseif ispc
        set(gcf,'Position',[1002 16 423 526])
    end
    for i=1:nacqT
        for j=par.nref+1:ntrackT-1
            subplot(221)
            p1 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.disp(:,2:end,i,j),options.display.disprange);axis image; title(['Disp: Push #',num2str(i)]);
            set(p1,'alphadata',~isnan(sweidata.disp(:,2:end,i,j)))
            set(gca,'color',[0.4 0.4 0.4])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            if options.motionFilter.enable
                subplot(222)
                p2 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata_mf_push.disp(:,2:end,i,j),options.display.disprange);axis image; title(['Disp MF: Push #',num2str(i)]);
                set(p2,'alphadata',~isnan(sweidata_mf_push.disp(:,2:end,i,j)))
                set(gca,'color',[0.4 0.4 0.4])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
            end
            
            subplot(223)
            p3 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.vel(:,2:end,i,j),options.display.velrange);axis image; title(['Vel: Push #',num2str(i)]);
            set(p3,'alphadata',~isnan(sweidata.vel(:,2:end,i,j)))
            set(gca,'color',[0.4 0.4 0.4])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            if options.motionFilter.enable
                subplot(224)
                p4 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata_mf_push.vel(:,2:end,i,j),options.display.velrange);axis image; title(['Vel MF: Push #',num2str(i)]);
                set(p4,'alphadata',~isnan(sweidata_mf_push.vel(:,2:end,i,j)))
                set(gca,'color',[0.4 0.4 0.4])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
            end
            pause(0.01)
        end
        pause
    end
    close(gcf)
end

%% Plot disp vs. time sw traces
if options.display.dvt_plots
    % Calculate indices for disp. vs. time plots
    for i=1:nacqT
        idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts));
    end
    
    line_idx = 1+ [1  7  13];
    for i=1:nacqT
        raw(:,:,i,:) = squeeze(sweidata.disp(idx(i,:),line_idx,i,par.nref:end));
    end
    if options.motionFilter.enable
        for i=1:nacqT
            mf(:,:,i,:) = squeeze(sweidata_mf_push.disp(idx(i,:),line_idx,i,par.nref:end));
        end
    else
        mf = nan(size(raw));
    end
    figure(101)
    if isunix
        set(gcf,'Position',[3 284 1196 574])
    elseif ispc
        set(gcf,'Position',[1002 16 423 526])
    end
    for p = 1:nacqT
        a1 = subplot(121);hold on;grid on;title(['Raw: ',num2str(p)]);
        set(a1,'yTickLabel',[])
        a2 = subplot(122);hold on;grid on;title(['MF: ',num2str(p)]);
        set(a2,'yTickLabel',[])
        offset = options.display.disprange(2);
        for i=1:5
            plot(sweidata.trackTime(par.nref:end),-offset*(i-1) + squeeze(raw(i,:,p,:))','-x','Parent',a1)
            plot(sweidata.trackTime(par.nref:end),-offset*(i-1) + squeeze(mf(i,:,p,:))','-x','Parent',a2)
        end
        set(gcf,'currentAxes',a1)
        axis tight
        set(gcf,'currentAxes',a2)
        axis tight
        pause
        cla(a1);cla(a2);
    end
    close(gcf)
end

% Compute Axially Averaged SpatiotemporalData
switch options.display.sw_display
    case 'disp'
        for i=1:nacqT
            raw(:,i,:) = squeeze(nanmean(sweidata.disp(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
        end
        raw = permute(raw, [1 3 2]);
        if options.motionFilter.enable
            for i=1:nacqT
                mf(:,i,:) = squeeze(nanmean(sweidata_mf_push.disp(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
            end
            mf = permute(mf, [1 3 2]);
        else
            mf = nan(size(raw));
        end
        rng = options.display.disprange;
    case 'vel'
        for i=1:nacqT
            raw(:,i,:) = squeeze(nanmean(sweidata.vel(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
        end
        raw = permute(raw, [1 3 2]);
        if options.motionFilter.enable
            for i=1:nacqT
                mf(:,i,:) = squeeze(nanmean(sweidata_mf_push.vel(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
            end
            mf = permute(mf, [1 3 2]);
        else
            mf = nan(size(raw));
        end
        rng = options.display.velrange;
end

idx = ceil(mean(gate_idx'));

int_factor = 4;
min_lat = 1.5; 
max_lat = 6.5; 
t_int = interp(double(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end)),int_factor);

% Interpolate Spatiotemporal data
for i=1:nacqT
    raw_int(:,:,i) = interp1(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),raw(:,par.nref+par.npush+par.nreverb+1:end,i)',t_int,'spline');
    if options.motionFilter.enable
        mf_int(:,:,i) = interp1(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),mf(:,par.nref+par.npush+par.nreverb+1:end,i)',t_int,'spline');
    end
end
raw_int = permute(raw_int,[2 1 3]);
sws_raw = zeros(1,nacqT);
if options.motionFilter.enable
    mf_int = permute(mf_int,[2 1 3]);
    sws_mf = zeros(1,nacqT);
end

tic
for i=1:nacqT
    temp = CalcSWSfromLatsums(raw_int(:,:,i),sweidata.lat(idx(i),:),t_int,1000,min_lat,max_lat,300);
    sws_raw(1,i) = temp.speed;
    if options.motionFilter.enable
        temp = CalcSWSfromLatsums(mf_int(:,par.nref+par.npush+par.nreverb+1:end,i),sweidata.lat(idx(i),:),t_int,1000,min_lat,max_lat,400);
        sws_mf(1,i) = temp.speed;
    end
    toc
end

figure;plot(sws_raw,'.--');ylim([0 8]);

keyboard
end

% % Calculate Shear Wave Speed
% if options.calcSWS.enable
%     switch options.calcSWS.method
%         case 'LinReg'
%             switch options.calcSWS.metric
%                 case 'TTP'
%                     if options.motionFilter.enable
%                         [pk, tpk] = subsamplepeak(sweidata_mf_push.trackTime(par.nref+par.npush+par.nreverb+1:end),sweidata_mf_push.disp(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     else
%                         [pk, tpk] = subsamplepeak(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),sweidata.disp(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     end
%                 case 'TTPS'
%                     if options.motionFilter.enable
%                         [pk, tpk] = subsamplepeak(sweidata_mf_push.t_vel(par.nref+par.npush+par.nreverb+1:end),sweidata_mf_push.vel(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     else
%                         [pk, tpk] = subsamplepeak(sweidata.t_vel(par.nref+par.npush+par.nreverb+1:end),sweidata.vel(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     end
%             end
%             
%             
%             figure(101);for i=1:40;imagesc(linspace(min(sweidata.lat(:)),max(sweidata.lat(:)),14),sweidata.axial,tpk(:,:,i),[min(min((tpk(:,:,i)))) max(max((tpk(:,:,i))))]);axis image;colorbar;title(i);pause;end
%             %             figure(101);for i=1:nacqT;errorbar(linspace(min(sweidata.lat(:)),max(sweidata.lat(:)),14),mean(tpk(gate_idx(1):gate_idx(2),:,i)),std(tpk(gate_idx(1):gate_idx(2),:,i)));
%             %             hold all;title(i);grid on;xlabel('Lateral (mm)');ylabel('Time to Peak (ms)');pause;end;close(gcf)
%             
%             
%             keyboard
%             klen = 5;
%             dxdi= linreg(sweidata.lat,klen,2);
%             dxdt = nan(size(tpk,1),size(tpk,3));
%             for i=1:size(tpk,3)
%                 [dtdi, r2] = linreg(tpk(:,:,i),klen,2);
%                 temp =  dxdi./dtdi;
%                 figure(102);for j=1:14;plot(temp(:,j));ylim([0 5]);grid on;title(j);pause;end
%                 dxdt(:,i) = temp(:,9).*(r2(:,9)>options.calcSWS.r2_threshold);
%                 %                 figure(102);hold all;plot(temp(:,9));ylim([0 5]);grid on;title(i);pause
%             end
%             close(gcf)
%             clear temp
%             
%             temp = dxdt;
%             temp(temp==0) = nan;
%             figure;hist(temp(:),100);title([nanmean(temp(:)), nanmedian(temp(:))])
%             pause
%             close(gcf)
%             
%             % Indices corresponding to median filter parameters
%             nax = double(ceil(options.display.medfilt(1)/(sweidata.axial(2) - sweidata.axial(1))));
%             nt = double(ceil(options.display.medfilt(2)/(sweidata.acqTime(2) - sweidata.acqTime(1))));
%             
%             axes('Position',[0.5 0.5 0.4 0.2]);
%             imagesc(sweidata.acqTime,sweidata.axial,medfilt2(dxdt,[nax nt]),options.calcSWS.SWSrange)
%             ylabel('Axial (mm)','fontsize',10','fontweight','bold')
%             title('SWS over DOF: Post Push','fontsize',10','fontweight','bold')
%             hcb = colorbar;
%             set(hcb,'Position',[0.91 0.5 0.0187/2 0.20])
%             ylabel(hcb,'SWS (m/s)','fontsize',10','fontweight','bold')
%             colormap(hot)
%             if isempty(ecgdata)
%                 xlabel('Acquisition Time (s)','fontsize',10','fontweight','bold')
%             end
%             
%         case 'LatSum'
%             keyboard
%             
%     end
% end
% 
% h2a = subplot('Position',[0.05 0.35 0.37 0.15]);
% h2b = subplot('Position',[0.05 0.1 0.37 0.15]);
% 
% for rep =1:1
%     for i=1:size(sweidata.disp,3)
%         if ~isempty(ecgdata)
%             set(pt,'Visible','off')
%         end
%         colormap(hot)
%         axes(h2a)
%         imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx)),:),raw(:,par.nref+1:end,i),rng/3)
%         xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
%         ylabel('Lateral (mm)','fontsize',10','fontweight','bold')
%         grid on
%         if strcmpi(options.display.sw_display,'disp')
%             title(sprintf('SWEI Profiles averaged over %2.2f - %2.2f mm\n Raw Displacements: Push # %d (t = %2.2fs)',gate(1),gate(2),i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%         elseif strcmpi(options.display.sw_display,'vel')
%             title(sprintf('SWEI Profiles averaged over %2.2f - %2.2f mm\n Raw Velocities: Push # %d (t = %2.2fs)',gate(1),gate(2),i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%         end
%         axes(h2b)
%         imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx)),:),mf(:,par.nref+1:end,i),rng/3)
%         xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
%         ylabel('Lateral (mm)','fontsize',10','fontweight','bold')
%         grid on
%         
%         hcb = colorbar;
%         set(hcb,'Position',[0.43 0.1 0.0187/2 0.4])
%         
%         if strcmpi(options.display.sw_display,'disp')
%             title(sprintf('Motion Filtered Displacements: Push # %d (t = %2.2f s)',i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%             ylabel(hcb,'Displacement (\mum)','fontsize',10','fontweight','bold')
%         elseif strcmpi(options.display.sw_display,'vel')
%             title(sprintf('Motion Filtered Velocities: Push # %d (t = %2.2f s)',i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%             ylabel(hcb,'Velocity (mm/s)','fontsize',10','fontweight','bold')
%         end
%         
%         if ~isempty(ecgdata)
%             pt = plot(sweidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
%         end
%         
%         %     frame = getframe(1);
%         %     im = frame2im(frame);
%         %     [imind,cm] = rgb2ind(im,256);
%         %     if i == 1;
%         %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
%         %     else
%         %         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
%         %
%         %     end
%         pause(0.05)
%     end
% end
=======
function [traced_gate,trace_flag] = dispSWEI(ecgdata,bdata,sweidata,sweidata_mf_pre,sweidata_mf_push,options,par)

%% Check for existance of traced gate
if isfield(sweidata,'traced_gate') && ~isempty(sweidata.traced_gate)
    prior_trace = 1;
    traced_gate = sweidata.traced_gate;
else
    prior_trace = 0;
    traced_gate = [];
end

dims = size(sweidata.disp);
ndepth = dims(1); nacqT = dims(2); ntrackT = dims(3);

edge = [sweidata.axial(1) sweidata.axial(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Gate Parameters
if prior_trace
    gate = repmat(sweidata.traced_gate,1,2) + repmat(options.display.gateOffset,length(sweidata.traced_gate),2) + repmat([-options.display.gateWidth/2 options.display.gateWidth/2],length(sweidata.traced_gate),1);
elseif ~prior_trace
    gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
    gate = repmat(gate,[nacqT 1]);
    traced_gate = [];
    if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

if options.display.axial_scan
    ngates = round(2*(edge(2)-edge(1))/options.display.gateWidth);
    offsets = linspace(edge(1)+options.display.gateWidth/2,edge(2)-options.display.gateWidth/2,ngates) - par.pushFocalDepth;
    options.display.gateOffset = offsets(1);
end

% Set Display Parameters
figure(2)
if isunix
    set(gcf,'Position',[1201 386 1920 1070])
    dispPar.fsize = 16;
elseif ispc
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    dispPar.fsize = 8;
end

if strcmpi(options.display.theme,'light')
    dispPar.fig = [1 1 1]; dispPar.txt = 'k'; dispPar.ax = [0.5 0.5 0.5];
    for i=1:size(bdata.bimg,3)
        temp = bdata.bimg(:,:,i);
        temp(temp==bdata.bimg(1,1,i)) = 256;
        bdata.bimg(:,:,i) = temp;
    end
elseif strcmpi(options.display.theme,'dark')
    dispPar.fig = 'k'; dispPar.txt = 'w'; dispPar.ax = [0.25 0.25 0.25];
    for i=1:size(bdata.bimg,3)
        temp = bdata.bimg(:,:,i);
        temp(temp==bdata.bimg(1,1,i)) = 0;
        bdata.bimg(:,:,i) = temp;
    end
end

set(1,'Color',dispPar.fig)

                                        img_rng = options.display.disprange;
                                        plot_rng = options.display.disprange;

dispPar.cmap = colormap(parula);
dispPar.cmap(end,:) = dispPar.ax;

dispPar.corder = winter(options.display.n_pts);

    
for i=1:nacqT
    gate_idx(i,:) = [find(sweidata.axial>gate(i,1),1,'first') find(sweidata.axial<gate(i,2),1,'last')];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Bmode Cine
for i=1:size(bdata.bimg,3);
    p1 = subplot('Position',[0.1 0.6 0.3 0.3]);
    imagesc(bdata.blat,bdata.bax,bdata.bimg(:,:,i));
    colormap(gray);axis image; freezeColors;
    hold(p1,'on')
    r1 = rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2,'Parent',p1);
    if prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linestyle','--','Linewidth',2,'Parent',p1);
    elseif ~prior_trace
        r2 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:)))  options.display.gateWidth],'EdgeColor','g','Linewidth',2,'Parent',p1);
    end
    hold off
    xlabel('Lateral (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    title(sprintf('B-Mode Cine: Frame %d (t = %1.1f s)',i,bdata.t(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    set(gca,'xcolor',dispPar.fig,'ycolor',dispPar.fig,'fontweight','bold')
    
    %     xlim([-25 25]);ylim([max(edge(1)-7.5,sweidata.IQaxial(1)) min(edge(2)+7.5,sweidata.IQaxial(end))])
    
    pause(0.025)
    if i==size(bdata.bimg,3)
        % Display M-mode IQ
        p2 = axes('Position',[0.45 0.3 0.52 1]);
        frame = abs(sweidata.IQ(:,:,1)); % Display first frame only
        imagesc(linspace(0,range(sweidata.lat(:))*nacqT,size(sweidata.IQ,2)),sweidata.IQaxial,db(frame/max(frame(:))),options.display.IQrange)
        hold(p2,'on')
        l1 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(1)*ones(1,nacqT),'b','Linewidth',2,'Parent',p2)
        l2 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(end)*ones(1,nacqT),'b','Linewidth',2,'Parent',p2)
        l3 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2,'Parent',p2)
        l4 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2,'Parent',p2)
        hold(p2,'off')
        axis image
        ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        set(gca,'xTickLabel',[])
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d',par.isHarmonic),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylim([max(edge(1)-7.5,sweidata.IQaxial(1)) min(edge(2)+7.5,sweidata.IQaxial(end))])
        colormap(gray); freezeColors;
        set(gca,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trace out center of depth gate
if (isfield(sweidata,'traced_gate') && isempty(options.display.gateOffset))
    gate = sweidata.traced_gate;
    hold on
    l1 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2);
    l2 = plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2);
    r1 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) max(gate(:))-min(gate(:))],'EdgeColor','g','LineStyle','--','Linewidth',2,'Parent',p1);
    trace_flag = input('Previously traced gate exists, do you want to retrace the gate? (y/n) [n]: ','s');
    if ~strcmpi(trace_flag,'y')
        trace_flag = 'n';
    else
        delete(l1,l2,r1)
    end
elseif (~isfield(sweidata,'traced_gate') && isempty(options.display.gateOffset))
    trace_flag = 'y';
end
if ((isempty(options.display.gateOffset) && ~isfield(sweidata,'traced_gate')) || strcmpi(trace_flag,'y'))
    hold on
    fprintf(1,'Ready to trace gate (GateWidth = %2.1f mm)...\n',options.display.gateWidth)
    for i=1:nacqT
        [x,y] = ginput(1);
        mark(i) = plot(x,y,'yo','MarkerSize',10);
        gate(i,:) = (y + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        clear x y
    end
    delete(mark)
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2)
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2)
    rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) max(gate(:))-min(gate(:))],'EdgeColor','g','LineStyle','--','Linewidth',2,'Parent',p1)
    fprintf(1,'Gate Traced.\n')
elseif ((isempty(options.display.gateOffset) && ~isfield(sweidata,'traced_gate')) && strcmpi(trace,'n'))
    gate = repmat(sweidata.traced_gate',1,2) + repmat([-options.display.gateWidth/2 options.display.gateWidth/2],length(sweidata.traced_gate),1);
    hold on
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2)
    plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2)
    rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) max(gate(:))-min(gate(:))],'EdgeColor','g','LineStyle','--','Linewidth',2,'Parent',p1)
end

traced_gate = mean(gate');

if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
end

for i=1:nacqT
    gate_idx(i,:) = [find(sweidata.axial>gate(i,1),1,'first') find(sweidata.axial<gate(i,2),1,'last')];
end

if ~isempty(ecgdata)
    samples = zeros(1,size(sweidata.disp,3));
    for i=1:size(sweidata.disp,3)
        samples(i) = ecgdata.swei(find(ecgdata.swei(:,1)>sweidata.acqTime(i),1,'first'),2);
    end
    ecgdata.swei(:,2) = ecgdata.swei(:,2)/max(ecgdata.swei(:,2));
    
    if options.calcSWS.enable
        h1 = axes('Position',[0.5 0.1 0.4 0.2]);
    else
        h1 = axes('Position',[0.5 0.25 0.4 0.2]);
    end
    plot(ecgdata.swei(:,1),ecgdata.swei(:,2),'Linewidth',2);
    hold on
    plot(sweidata.acqTime,samples,'kx','MarkerSize',8)
    pt = plot(sweidata.acqTime(1),samples(1),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    hold off
    grid on
    title('ECG Trace','fontsize',fsize','fontweight','bold')
    xlabel('Acquisition Time (s)','fontsize',fsize','fontweight','bold')
    axis tight
    hold(h1)
end

% Compute velocity data from displacement data
[vel_temp, sweidata.t_vel] = differentiateDisplacements(permute(sweidata.disp,[1 2 4 3]),sweidata.trackTime,options.motionFilter.LPF_Cutoff);
sweidata.vel = 1*permute(vel_temp,[1 2 4 3]);
sweidata = interpPushReverb(sweidata,options,par,'nan'); % NaN out push and reverb frames
clear vel_temp
if options.motionFilter.enable
    [vel_temp, sweidata_mf_pre.t_vel] = differentiateDisplacements(permute(sweidata_mf_pre.disp,[1 2 4 3]),sweidata_mf_pre.trackTime,options.motionFilter.LPF_Cutoff);
    sweidata_mf_pre.vel = 1*permute(vel_temp,[1 2 4 3]);
    sweidata_mf_pre = interpPushReverb(sweidata_mf_pre,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
    [vel_temp, sweidata_mf_push.t_vel] = differentiateDisplacements(permute(sweidata_mf_push.disp,[1 2 4 3]),sweidata_mf_push.trackTime,options.motionFilter.LPF_Cutoff);
    sweidata_mf_push.vel = 1*permute(vel_temp,[1 2 4 3]);
    sweidata_mf_push = interpPushReverb(sweidata_mf_push,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
end

% Coorelation mask filter
if options.display.cc_filt
    mask = sweidata.cc>options.display.cc_thresh;
    sweidata.disp(mask==0) = nan;
    sweidata.vel(mask(:,:,:,1:end-1)==0) = nan;
    if options.motionFilter.enable
        sweidata_mf_pre.disp(mask==0) = nan;
        sweidata_mf_pre.vel(mask(:,:,:,1:end-1)==0) = nan;
        sweidata_mf_push.disp(mask==0) = nan;
        sweidata_mf_push.vel(mask(:,:,:,1:end-1)==0) = nan;
    end
end

%% Plot Axial Scan Data
if options.display.axial_scan
    for i=1:ngates
        set(0,'currentFigure',f1)
        options.display.gateOffset = offsets(i);
        gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        gate = repmat(gate,[nacqT 1]);
        if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
            warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
        end
        
        for j=1:nacqT
            gate_idx(j,:) = [find(sweidata.axial>gate(j,1),1,'first') find(sweidata.axial<gate(j,2),1,'last')];
        end
        
        delete(r1)
        r1 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) options.display.gateWidth ],'EdgeColor','g','Linewidth',2,'Parent',p1);
        title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)\n Gate Offset = %2.2f mm ',size(bdata.bimg,3),bdata.t(size(bdata.bimg,3)),offsets(i)),'fontsize',fsize','fontweight','bold','Parent',p1)
        p2 = axes('Position',[0.45 0.3 0.52 1]);
        frame = sweidata.IQ(:,:,1);
        imagesc(linspace(0,range(sweidata.lat(:))*nacqT,size(sweidata.IQ,2)),sweidata.IQaxial,db(frame/max(frame(:))),options.display.IQrange)
        hold on
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(1)*ones(1,nacqT),'b','Linewidth',2)
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),sweidata.axial(end)*ones(1,nacqT),'b','Linewidth',2)
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2)
        plot(linspace(0,range(sweidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2)
        hold off
        axis image
        ylabel('Axial (mm)','fontsize',fsize','fontweight','bold')
        set(gca,'xTickLabel',[])
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d',par.isHarmonic),'fontsize',fsize,'fontweight','bold')
        ylim([max(edge(1)-15,sweidata.IQaxial(1)) min(edge(2)+15,sweidata.IQaxial(end))])
        
        % Compute Axially Averaged Data
        switch options.display.sw_display
            case 'disp'
                raw = squeeze(nanmean(sweidata.disp(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                raw = permute(raw, [1 3 2]);
                if options.motionFilter.enable
                    mf = squeeze(nanmean(sweidata_mf_push.disp(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                    mf = permute(mf, [1 3 2]);
                else
                    mf = nan(size(raw));
                end
                rng = options.display.disprange*0.5;
            case 'vel'
                raw = squeeze(nanmean(sweidata.vel(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                raw = permute(raw, [1 3 2]);
                if options.motionFilter.enable
                    mf = squeeze(nanmean(sweidata_mf_push.vel(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                    mf = permute(mf, [1 3 2]);
                else
                    mf = nan(size(raw));
                end
                rng = options.display.velrange;
        end
        
        figure(101)
        if isunix
            set(gcf,'Position',[1200 409 950 1050])
        elseif ispc
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
        end
        for i=1:nacqT
            subplot(nacqT/5,5,i)
            p = imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx(1,:))),:),raw(:,(par.nref+1:end),i),rng);
            set(p,'alphadata',~isnan(raw(:,(par.nref+1:end),i)))
            set(gca,'color',[0.4 0.4 0.4])
            title(sweidata.acqTime(i),'fontsize',fsize','fontweight','bold')
            xlim([0 7])
        end
        
        figure(102)
        if isunix
            set(gcf,'Position',[2150 409 950 1050])
        elseif ispc
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
        end
        for i=1:nacqT
            subplot(nacqT/5,5,i)
            p = imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx(1,:))),:),mf(:,(par.nref+1:end),i),rng);
            set(p,'alphadata',~isnan(raw(:,(par.nref+1:end),i)))
            set(gca,'color',[0.4 0.4 0.4])
            title(sweidata.acqTime(i),'fontsize',fsize','fontweight','bold')
            xlim([0 7])
        end
        pause
    end
    close 101 102
    clc
end

%% Plot Shear Wave Movies
if options.display.sw_movie
    figure(101)
    if isunix
        set(gcf,'Position',[69 435 414 906])
    elseif ispc
        set(gcf,'Position',[1002 16 423 526])
    end
    for i=1:nacqT
        for j=par.nref+1:ntrackT-1
            subplot(221)
            p1 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.disp(:,2:end,i,j),options.display.disprange);axis image; title(['Disp: Push #',num2str(i)]);
            set(p1,'alphadata',~isnan(sweidata.disp(:,2:end,i,j)))
            set(gca,'color',[0.4 0.4 0.4])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            if options.motionFilter.enable
                subplot(222)
                p2 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata_mf_push.disp(:,2:end,i,j),options.display.disprange);axis image; title(['Disp MF: Push #',num2str(i)]);
                set(p2,'alphadata',~isnan(sweidata_mf_push.disp(:,2:end,i,j)))
                set(gca,'color',[0.4 0.4 0.4])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
            end
            
            subplot(223)
            p3 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.vel(:,2:end,i,j),options.display.velrange);axis image; title(['Vel: Push #',num2str(i)]);
            set(p3,'alphadata',~isnan(sweidata.vel(:,2:end,i,j)))
            set(gca,'color',[0.4 0.4 0.4])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            if options.motionFilter.enable
                subplot(224)
                p4 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata_mf_push.vel(:,2:end,i,j),options.display.velrange);axis image; title(['Vel MF: Push #',num2str(i)]);
                set(p4,'alphadata',~isnan(sweidata_mf_push.vel(:,2:end,i,j)))
                set(gca,'color',[0.4 0.4 0.4])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
            end
            pause(0.01)
        end
        pause
    end
    close(gcf)
end

%% Plot disp vs. time sw traces
if options.display.dvt_plots
    % Calculate indices for disp. vs. time plots
    for i=1:nacqT
        idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts));
    end
    
    line_idx = 1+ [1  7  13];
    for i=1:nacqT
        raw(:,:,i,:) = squeeze(sweidata.disp(idx(i,:),line_idx,i,par.nref:end));
    end
    if options.motionFilter.enable
        for i=1:nacqT
            mf(:,:,i,:) = squeeze(sweidata_mf_push.disp(idx(i,:),line_idx,i,par.nref:end));
        end
    else
        mf = nan(size(raw));
    end
    figure(101)
    if isunix
        set(gcf,'Position',[3 284 1196 574])
    elseif ispc
        set(gcf,'Position',[1002 16 423 526])
    end
    for p = 1:nacqT
        a1 = subplot(121);hold on;grid on;title(['Raw: ',num2str(p)]);
        set(a1,'yTickLabel',[])
        a2 = subplot(122);hold on;grid on;title(['MF: ',num2str(p)]);
        set(a2,'yTickLabel',[])
        offset = options.display.disprange(2);
        for i=1:5
            plot(sweidata.trackTime(par.nref:end),-offset*(i-1) + squeeze(raw(i,:,p,:))','-x','Parent',a1)
            plot(sweidata.trackTime(par.nref:end),-offset*(i-1) + squeeze(mf(i,:,p,:))','-x','Parent',a2)
        end
        set(gcf,'currentAxes',a1)
        axis tight
        set(gcf,'currentAxes',a2)
        axis tight
        pause
        cla(a1);cla(a2);
    end
    close(gcf)
end

% Compute Axially Averaged SpatiotemporalData
switch options.display.sw_display
    case 'disp'
        for i=1:nacqT
            raw(:,i,:) = squeeze(nanmean(sweidata.disp(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
        end
        raw = permute(raw, [1 3 2]);
        if options.motionFilter.enable
            for i=1:nacqT
                mf(:,i,:) = squeeze(nanmean(sweidata_mf_push.disp(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
            end
            mf = permute(mf, [1 3 2]);
        else
            mf = nan(size(raw));
        end
        rng = options.display.disprange;
    case 'vel'
        for i=1:nacqT
            raw(:,i,:) = squeeze(nanmean(sweidata.vel(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
        end
        raw = permute(raw, [1 3 2]);
        if options.motionFilter.enable
            for i=1:nacqT
                mf(:,i,:) = squeeze(nanmean(sweidata_mf_push.vel(gate_idx(i,1):gate_idx(i,2),2:end,i,:),1));
            end
            mf = permute(mf, [1 3 2]);
        else
            mf = nan(size(raw));
        end
        rng = options.display.velrange;
end

idx = ceil(mean(gate_idx'));

int_factor = 4;
min_lat = 1.5; 
max_lat = 6.5; 
t_int = interp(double(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end)),int_factor);

% Interpolate Spatiotemporal data
for i=1:nacqT
    raw_int(:,:,i) = interp1(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),raw(:,par.nref+par.npush+par.nreverb+1:end,i)',t_int,'spline');
    if options.motionFilter.enable
        mf_int(:,:,i) = interp1(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),mf(:,par.nref+par.npush+par.nreverb+1:end,i)',t_int,'spline');
    end
end
raw_int = permute(raw_int,[2 1 3]);
sws_raw = zeros(1,nacqT);
if options.motionFilter.enable
    mf_int = permute(mf_int,[2 1 3]);
    sws_mf = zeros(1,nacqT);
end

tic
for i=1:nacqT
    temp = CalcSWSfromLatsums(raw_int(:,:,i),sweidata.lat(idx(i),:),t_int,1000,min_lat,max_lat,300);
    sws_raw(1,i) = temp.speed;
    if options.motionFilter.enable
        temp = CalcSWSfromLatsums(mf_int(:,par.nref+par.npush+par.nreverb+1:end,i),sweidata.lat(idx(i),:),t_int,1000,min_lat,max_lat,400);
        sws_mf(1,i) = temp.speed;
    end
    toc
end

figure;plot(sws_raw,'.--');ylim([0 8]);

keyboard
end

% % Calculate Shear Wave Speed
% if options.calcSWS.enable
%     switch options.calcSWS.method
%         case 'LinReg'
%             switch options.calcSWS.metric
%                 case 'TTP'
%                     if options.motionFilter.enable
%                         [pk, tpk] = subsamplepeak(sweidata_mf_push.trackTime(par.nref+par.npush+par.nreverb+1:end),sweidata_mf_push.disp(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     else
%                         [pk, tpk] = subsamplepeak(sweidata.trackTime(par.nref+par.npush+par.nreverb+1:end),sweidata.disp(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     end
%                 case 'TTPS'
%                     if options.motionFilter.enable
%                         [pk, tpk] = subsamplepeak(sweidata_mf_push.t_vel(par.nref+par.npush+par.nreverb+1:end),sweidata_mf_push.vel(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     else
%                         [pk, tpk] = subsamplepeak(sweidata.t_vel(par.nref+par.npush+par.nreverb+1:end),sweidata.vel(:,2:end,:,par.nref+par.npush+par.nreverb+1:end),4);
%                     end
%             end
%             
%             
%             figure(101);for i=1:40;imagesc(linspace(min(sweidata.lat(:)),max(sweidata.lat(:)),14),sweidata.axial,tpk(:,:,i),[min(min((tpk(:,:,i)))) max(max((tpk(:,:,i))))]);axis image;colorbar;title(i);pause;end
%             %             figure(101);for i=1:nacqT;errorbar(linspace(min(sweidata.lat(:)),max(sweidata.lat(:)),14),mean(tpk(gate_idx(1):gate_idx(2),:,i)),std(tpk(gate_idx(1):gate_idx(2),:,i)));
%             %             hold all;title(i);grid on;xlabel('Lateral (mm)');ylabel('Time to Peak (ms)');pause;end;close(gcf)
%             
%             
%             keyboard
%             klen = 5;
%             dxdi= linreg(sweidata.lat,klen,2);
%             dxdt = nan(size(tpk,1),size(tpk,3));
%             for i=1:size(tpk,3)
%                 [dtdi, r2] = linreg(tpk(:,:,i),klen,2);
%                 temp =  dxdi./dtdi;
%                 figure(102);for j=1:14;plot(temp(:,j));ylim([0 5]);grid on;title(j);pause;end
%                 dxdt(:,i) = temp(:,9).*(r2(:,9)>options.calcSWS.r2_threshold);
%                 %                 figure(102);hold all;plot(temp(:,9));ylim([0 5]);grid on;title(i);pause
%             end
%             close(gcf)
%             clear temp
%             
%             temp = dxdt;
%             temp(temp==0) = nan;
%             figure;hist(temp(:),100);title([nanmean(temp(:)), nanmedian(temp(:))])
%             pause
%             close(gcf)
%             
%             % Indices corresponding to median filter parameters
%             nax = double(ceil(options.display.medfilt(1)/(sweidata.axial(2) - sweidata.axial(1))));
%             nt = double(ceil(options.display.medfilt(2)/(sweidata.acqTime(2) - sweidata.acqTime(1))));
%             
%             axes('Position',[0.5 0.5 0.4 0.2]);
%             imagesc(sweidata.acqTime,sweidata.axial,medfilt2(dxdt,[nax nt]),options.calcSWS.SWSrange)
%             ylabel('Axial (mm)','fontsize',10','fontweight','bold')
%             title('SWS over DOF: Post Push','fontsize',10','fontweight','bold')
%             hcb = colorbar;
%             set(hcb,'Position',[0.91 0.5 0.0187/2 0.20])
%             ylabel(hcb,'SWS (m/s)','fontsize',10','fontweight','bold')
%             colormap(hot)
%             if isempty(ecgdata)
%                 xlabel('Acquisition Time (s)','fontsize',10','fontweight','bold')
%             end
%             
%         case 'LatSum'
%             keyboard
%             
%     end
% end
% 
% h2a = subplot('Position',[0.05 0.35 0.37 0.15]);
% h2b = subplot('Position',[0.05 0.1 0.37 0.15]);
% 
% for rep =1:1
%     for i=1:size(sweidata.disp,3)
%         if ~isempty(ecgdata)
%             set(pt,'Visible','off')
%         end
%         colormap(hot)
%         axes(h2a)
%         imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx)),:),raw(:,par.nref+1:end,i),rng/3)
%         xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
%         ylabel('Lateral (mm)','fontsize',10','fontweight','bold')
%         grid on
%         if strcmpi(options.display.sw_display,'disp')
%             title(sprintf('SWEI Profiles averaged over %2.2f - %2.2f mm\n Raw Displacements: Push # %d (t = %2.2fs)',gate(1),gate(2),i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%         elseif strcmpi(options.display.sw_display,'vel')
%             title(sprintf('SWEI Profiles averaged over %2.2f - %2.2f mm\n Raw Velocities: Push # %d (t = %2.2fs)',gate(1),gate(2),i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%         end
%         axes(h2b)
%         imagesc(sweidata.trackTime(par.nref+1:end),sweidata.lat(round(mean(gate_idx)),:),mf(:,par.nref+1:end,i),rng/3)
%         xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
%         ylabel('Lateral (mm)','fontsize',10','fontweight','bold')
%         grid on
%         
%         hcb = colorbar;
%         set(hcb,'Position',[0.43 0.1 0.0187/2 0.4])
%         
%         if strcmpi(options.display.sw_display,'disp')
%             title(sprintf('Motion Filtered Displacements: Push # %d (t = %2.2f s)',i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%             ylabel(hcb,'Displacement (\mum)','fontsize',10','fontweight','bold')
%         elseif strcmpi(options.display.sw_display,'vel')
%             title(sprintf('Motion Filtered Velocities: Push # %d (t = %2.2f s)',i,sweidata.acqTime(i)),'fontsize',10','fontweight','bold')
%             ylabel(hcb,'Velocity (mm/s)','fontsize',10','fontweight','bold')
%         end
%         
%         if ~isempty(ecgdata)
%             pt = plot(sweidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
%         end
%         
%         %     frame = getframe(1);
%         %     im = frame2im(frame);
%         %     [imind,cm] = rgb2ind(im,256);
%         %     if i == 1;
%         %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
%         %     else
%         %         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
%         %
%         %     end
%         pause(0.05)
%     end
% end
>>>>>>> afa77557451dedcbdaaf013db29dc3d36822037b
