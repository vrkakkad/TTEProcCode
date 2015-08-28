function [traced_gate] = dispSWEI(bdata,sweidata,options,par)
%% Check for existance of traced gate
[~,temp] = fileparts(pwd);
if isunix
    gate_path = strcat('/emfd/vrk4/Transthoracic_Clinical/traced_gates/',temp); %% Change this path!!
elseif ispc
    gate_path = strcat('E:\ClinicalDataArchive\Traced_Gates\',temp);
end
fname = strcat(gate_path,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_swei.mat');
if exist(fname,'file')
    prior_trace = 1;
    load(fname);
else
    prior_trace = 0;
    traced_gate = [];
end
dims = size(sweidata.disp);
ndepth = dims(1); nBeams = dims(2); nacqT = dims(3); ntrackT = dims(4);
edge = [sweidata.axial(1) sweidata.axial(end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Gate Parameters
if prior_trace
    gate = traced_gate;
elseif ~prior_trace
    gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
    gate = repmat(gate,[nacqT 1]);
end

if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
end
if options.display.axial_scan
    ngates = round(2*(edge(2)-edge(1))/options.display.gateWidth);
    offsets = linspace(edge(1)+options.display.gateWidth/2,edge(2)-options.display.gateWidth/2,ngates) - par.pushFocalDepth;
    options.display.gateOffset = offsets(1);
end
% Set Display Parameters
fig2 = figure;
set(fig2,'Name','TTE SWEI','UserData','main_fig','InvertHardCopy', 'off')
if isunix
    set(fig2,'Position',[1201 386 1920 1070])
    dispPar.fsize = 12;
elseif ispc
    set(fig2,'units','normalized','outerposition',[0 0 1 1])
    dispPar.fsize = 10;
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
set(fig2,'Color',dispPar.fig)
trackPRF = 1000/par.priusec(1); % kHz
rng = options.display.dispRange;
dispPar.cmap = colormap(parula);
dispPar.cmap(end,:) = dispPar.ax;
dispPar.corder = winter(options.display.n_pts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Bmode Cine
if ~isempty(bdata.ecg)
    bsamples = zeros(1,nacqT);
    for i=1:size(bdata.bimg,3)
        bsamples(i) = bdata.ecg(find(bdata.ecg(:,1)>bdata.t(i),1,'first')-1,2);
        bdata.bimg(:,:,i) = fliplr(bdata.bimg(:,:,i));
    end
    ax112 = axes('Position',[0.125 0.5 0.25 0.1]);
    plot(bdata.ecg(:,1),bdata.ecg(:,2),'LineWidth',2,'Parent',ax112);
    hold(ax112,'on')
    plot(bdata.t,bsamples,'gx','MarkerSize',6,'Parent',ax112)
    pt = plot(bdata.t(1),bsamples(1),'ro','MarkerSize',6,'MarkerFaceColor','r','Parent',ax112);
    title(sprintf('ECG (Bmode): HR = %2.0f bpm',bdata.hr),'FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    xlabel('Acquisition Time (s)','FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    xlim([0 max(bdata.t)])
    ylim([min(bdata.ecg(:,2)) max(bdata.ecg(:,2))])
    set(ax112,'Color',dispPar.ax,'XColor',dispPar.txt,'YColor',dispPar.txt,'YTickLabel',[],'FontWeight','Bold','XGrid','On','GridLineStyle','--','UserData','b_ecg_ax')
end
for i=1:size(bdata.bimg,3);
    ax11 = subplot('Position',[0.1 0.65 0.3 0.3]);
    imagesc(bdata.blat,bdata.bax,fliplr(bdata.bimg(:,:,i)));
    colormap(gray);axis image; freezeColors;
    hold(ax11,'on')
    plot(0,par.pushFocalDepth,'o','MarkerSize',6,'MarkerFaceColor','c')
    r1 = rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2,'Parent',ax11);
    if prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linestyle','--','Linewidth',2,'Parent',ax11);
    elseif ~prior_trace
        r2 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:)))  options.display.gateWidth],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    end
    %     xlabel('Lateral (mm)','FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.fig)
    %     ylabel('Axial (mm)','FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.fig)
    %     title(sprintf('TimeStamp = %s \nB-Mode Cine: Frame %d (t = %1.1f s)',num2str(options.timeStamp),i,bdata.t(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    title(sprintf('B-Mode Cine: Frame %d (t = %1.1f s)',i,bdata.t(i)),'FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    set(ax11,'XColor',dispPar.fig,'YColor',dispPar.fig,'FontWeight','Bold','UserData','bmodecine_ax')
    hold(ax11,'off')
    if ~isempty(bdata.ecg)
        delete(pt)
        pt = plot(bdata.t(i),bsamples(i),'ro','MarkerSize',6,'MarkerFaceColor','r','Parent',ax112);
    end
    pause(median(diff(bdata.t)))
    if i==size(bdata.bimg,3)
        % Display M-mode IQ
        ax12 = axes('Position',[0.45 0.7 0.50 0.2]);
        frame = abs(sweidata.IQ(:,:,1)); % Display first frame only
        frame = db(frame/max(frame(:)));
        temp = find(sweidata.IQaxial>edge(1)-1,1);
        if isempty(temp); idx(1) = 1; else idx(1) = temp; end; clear temp
        temp = find(sweidata.IQaxial>edge(2)+1,1);
        if isempty(temp); idx(2) = length(sweidata.IQaxial); else idx(2) = temp; end; clear temp
        imagesc(linspace(0,sweidata.acqTime(end),size(sweidata.IQ,2)),sweidata.IQaxial(idx(1):idx(2)),frame(idx(1):idx(2),:),options.display.IQrange)
        clear idx
        hold(ax12,'on')
        l1 = plot(linspace(0,sweidata.acqTime(end),nacqT),sweidata.axial(1)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l2 = plot(linspace(0,sweidata.acqTime(end),nacqT),sweidata.axial(end)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l3 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
        l4 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
        plot(0,par.pushFocalDepth,'>','Markersize',10,'MarkerFaceColor','c')
        hold(ax12,'off')
        xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d\n Track PRF = %1.2f kHz, Push = %d x %d \\mus',par.isHarmonic,trackPRF,par.npush,par.pushDurationusec),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        colormap(gray); freezeColors;
        set(ax12,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','UserData','mmodeIQ_ax')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute velocity data from displacement data
[vel_temp, sweidata.t_vel] = differentiateDisplacements(permute(sweidata.disp,[1 2 4 3]),sweidata.trackTime,options.motionFilter.Cutoff(2));
sweidata.vel = 1*permute(vel_temp,[1 2 4 3]);
clear vel_temp
if ~strcmpi(options.motionFilter.method,'Off')
    [vel_temp, sweidata.t_vel] = differentiateDisplacements(permute(sweidata.disp_mf_pre,[1 2 4 3]),sweidata.trackTime,options.motionFilter.Cutoff(2));
    sweidata.vel_mf_pre = 1*permute(vel_temp,[1 2 4 3]);
    %     sweidata = interpPushReverb(sweidata,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
    [vel_temp, sweidata.t_vel] = differentiateDisplacements(permute(sweidata.disp_mf_push,[1 2 4 3]),sweidata.trackTime,options.motionFilter.Cutoff(2));
    sweidata.vel_mf_push = 1*permute(vel_temp,[1 2 4 3]);
    %     sweidata= interpPushReverb(sweidata,options,par,'nan'); % NaN out push and reverb frames
    clear vel_temp
end
sweidata = interpPushReverb(sweidata,options,par,'nan'); % NaN out push and reverb frames
% Coorelation mask filter
if options.display.cc_filt
    mask = sweidata.cc>options.display.cc_thresh;
    sweidata.disp(mask==0) = nan;
    sweidata.vel(mask(:,:,:,1:end-1)==0) = nan;
    if ~strcmpi(options.motionFilter.method,'Off')
        sweidata.disp_mf_pre(mask==0) = nan;
        sweidata.vel_mf_pre(mask(:,:,:,1:end-1)==0) = nan;
        sweidata.disp_mf_push(mask==0) = nan;
        sweidata.vel_mf_push(mask(:,:,:,1:end-1)==0) = nan;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incorporate ECG Data into this figure
if ~isempty(sweidata.ecg)
    samples = zeros(1,size(sweidata.disp,3));
    for i=1:nacqT
        samples(i) = sweidata.ecg(find(sweidata.ecg(:,1)>sweidata.acqTime(i),1,'first'),2);
    end
    ax15 = axes('Position',[0.5 0.1 0.4 0.2]);
    plot(sweidata.ecg(:,1),sweidata.ecg(:,2),'Linewidth',2);
    hold(ax15,'on')
    plot(sweidata.acqTime,samples,'gx','MarkerSize',8)
    pt = plot(sweidata.acqTime(1),samples(1),'ro','Parent',ax15,'Markersize',10,'Markerfacecolor','r');
    title('ECG Trace','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    xlim([0 max(sweidata.acqTime)])
    ylim([min(sweidata.ecg(:,2)) max(sweidata.ecg(:,2))])
    set(ax15,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--','Userdata','ecg_ax')
    hold(ax15,'off')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% figure;
% for j=1:nacqT;
%     for i=2:nBeams;
%         subplot(231), imagesc(squeeze(sweidata.disp(:,i,j,:)),[-5 20]), title([i j]);
%         subplot(232), imagesc(squeeze(sweidata.disp_mf_pre(:,i,j,:)),[-5 20]), title([i j]);
%         subplot(233), imagesc(squeeze(sweidata.disp_mf_push(:,i,j,:)),[-5 20]), title([i j]);
%         subplot(234), imagesc(squeeze(sweidata.vel(:,i,j,:)),[-5 20]), title([i j]);
%         subplot(235), imagesc(squeeze(sweidata.vel_mf_pre(:,i,j,:)),[-5 10]), title([i j]);
%         subplot(236), imagesc(squeeze(sweidata.vel_mf_push(:,i,j,:)),[-5 10]), title([i j]);
%         pause
%     end
% end
%%
% figure;
% for j=1:nacqT;
%     for i=1:ntrackT-1;
%         subplot(231), imagesc(squeeze(sweidata.disp(:,[2:end],j,i)),[-5 20]), title([i j]);
%         subplot(232), imagesc(squeeze(sweidata.disp_mf_pre(:,[2:end],j,i)),[-5 20]), title([i j]);
%         subplot(233), imagesc(squeeze(sweidata.disp_mf_push(:,[2:end],j,i)),[-5 20]), title([i j]);
%         subplot(234), imagesc(squeeze(sweidata.vel(:,[2:end],j,i)),[-5 20]), title([i j]);
%         subplot(235), imagesc(squeeze(sweidata.vel_mf_pre(:,[2:end],j,i)),[-5 10]), title([i j]);
%         subplot(236), imagesc(squeeze(sweidata.vel_mf_push(:,[2:end],j,i)),[-5 10]), title([i j]);
%         pause
%     end
% end
%% Plot Axial Scan Data
if options.display.axial_scan
    for i=1:ngates
        set(0,'currentFigure',fig2)
        options.display.gateOffset = offsets(i);
        gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        gate = repmat(gate,[nacqT 1]);
        if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
            warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
        end
        for j=1:nacqT
            gate_idx(j,:) = [find(sweidata.axial>gate(j,1),1,'first') find(sweidata.axial<gate(j,2),1,'last')];
        end
        delete(r2)
        r2 = rectangle('Position',[min(sweidata.lat(:)) min(gate(:)) (max(sweidata.lat(:))-min(sweidata.lat(:))) options.display.gateWidth ],'EdgeColor','g','Linewidth',2,'Parent',ax11);
        title(sprintf('B-Mode Cine: Frame %d (t = %1.1f s)\n Gate Offset = %2.2f mm',i,bdata.t(i),offsets(i)),'FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
        ax12 = axes('Position',[0.45 0.7 0.50 0.2]);
        frame = abs(sweidata.IQ(:,:,1)); % Display first frame only
        frame = db(frame/max(frame(:)));
        temp = find(sweidata.IQaxial>edge(1)-1,1);
        if isempty(temp); idx(1) = 1; else idx(1) = temp; end; clear temp
        temp = find(sweidata.IQaxial>edge(2)+1,1);
        if isempty(temp); idx(2) = length(sweidata.IQaxial); else idx(2) = temp; end; clear temp
        imagesc(linspace(0,sweidata.acqTime(end),size(sweidata.IQ,2)),sweidata.IQaxial(idx(1):idx(2)),frame(idx(1):idx(2),:),options.display.IQrange)
        clear idx
        hold(ax12,'on')
        l1 = plot(linspace(0,sweidata.acqTime(end),nacqT),sweidata.axial(1)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l2 = plot(linspace(0,sweidata.acqTime(end),nacqT),sweidata.axial(end)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l3 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
        l4 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
        hold(ax12,'off')
        xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d\n Track PRF = %1.2f kHz, Push = %d x %d \\mus',par.isHarmonic,trackPRF,par.npush,par.pushDurationusec),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        colormap(gray); freezeColors;
        set(ax12,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','UserData','mmodeIQ_ax')
        % Compute Axially Averaged Data
        switch options.display.sw_display
            case 'disp'
                raw = squeeze(nanmean(sweidata.disp(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                raw = permute(raw, [1 3 2]);
                if ~strcmpi(options.motionFilter.method,'Off')
                    mf = squeeze(nanmean(sweidata.disp_mf_push(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                    mf = permute(mf, [1 3 2]);
                else
                    mf = nan(size(raw));
                end
                rng = options.display.dispRange;
            case 'vel'
                raw = squeeze(nanmean(sweidata.vel(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
                raw = permute(raw, [1 3 2]);
                if options.motionFilter.enable
                    mf = squeeze(nanmean(sweidata.vel_mf_push(min(gate_idx(:)):max(gate_idx(:)),2:end,:,:),1));
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
            title(sweidata.acqTime(i),'fontsize',dispPar.fsize','fontweight','bold')
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
            title(sweidata.acqTime(i),'fontsize',dispPar.fsize','fontweight','bold')
            xlim([0 7])
        end
        pause
    end
    close 101 102
    clc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trace out center of depth gate
if options.display.extras == -1  % no breaks
    trace_input = 'n';
    flatgate_input = 'y';
end
if prior_trace
    trace_input = input('\nDo you want to retrace the gate? (y/n) [n]: ','s');
    flatgate_input = input('\nDo you want to use a flat gate? (y/n) [n]: ','s');
elseif (~prior_trace && options.display.extras ~= -1)
    trace_input = input('\nDo you want to trace a depth gate? (y/n) [n]: ','s');
    flatgate_input = 'y';
end
if strcmpi(flatgate_input,'y')
    gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
    gate = repmat(gate,[nacqT 1]);
    delete(r2);delete(l3);delete(l4);
    hold(ax11,'on')
    r2 = rectangle('Position',[-2 min(gate(:)) 4 options.display.gateWidth],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    hold(ax12,'on')
    l3 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
    l4 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
    if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
    end
end
if strcmpi(trace_input,'y')
    delete(r2);delete(l3);delete(l4);
    fprintf(1,'\nReady to trace gate (GateWidth = %2.1f mm)...\nClick to define points, hit space to end tracing\n',options.display.gateWidth)
    fprintf(1,'Trace Top...Hit enter when done\n');
    i=1;
    while i<=nacqT
        [a,b,status] = ginput_alt(1,'y');
        if isempty(status)
            break
        else
            hold on
            x(i) = a; y(i) = b;
            top_mark(i) = plot(x(i),y(i),'ys','MarkerSize',6);
            i=i+1;
        end
    end
    top = interp1(x,y,sweidata.acqTime);
    clear x y status
    fprintf(1,'Trace Bottom...Hit enter when done\n');
    i=1;
    while i<=nacqT
        [a,b,status] = ginput_alt(1,'r');
        if isempty(status)
            break
        else
            hold on
            x(i) = a; y(i) = b;
            bot_mark(i) = plot(x(i),y(i),'rs','MarkerSize',6);
            i=i+1;
        end
    end
    bot = interp1(x,y,sweidata.acqTime);
    clear x y status
    delete(top_mark);delete(bot_mark)
    gate = [smooth(top',5) smooth(bot',5)];
    hold(ax11,'on')
    r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    hold(ax12,'on')
    l3 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
    l4 = plot(linspace(0,sweidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
    fprintf(1,'Gate Traced.\n')
    traced_gate = gate;
    if (min(gate(:))<sweidata.axial(1) || max(gate(:))>sweidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),sweidata.axial(1),sweidata.axial(end));
    end
    % Save traced gate
    fprintf(1,'Saving Traced Gate...\n');
    [~,temp] = fileparts(pwd);
    if isunix
        gate_path = strcat('/emfd/vrk4/Transthoracic_Clinical/traced_gates/',temp); %% Change this path!!
    elseif ispc
        gate_path = strcat('E:\ClinicalDataArchive\Traced_Gates\',temp);
    end
    if ~exist(gate_path,'dir'); mkdir(gate_path); end
    fname = strcat(gate_path,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_swei.mat');
    save(fname,'traced_gate');
end
for i=1:nacqT
    gate_idx(i,:) = [find(sweidata.axial>gate(i,1),1,'first') find(sweidata.axial<gate(i,2),1,'last')];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Shear Wave Movies
if options.display.sw_movie
    figure(101)
    if isunix
        set(gcf,'Position',[69 435 414 906])
    elseif ispc
        set(gcf,'Position',[265 28 423 526])
    end
    i=1;
    while i<=nacqT
        for j=par.nref+1:ntrackT-1
            subplot(221)
            p1 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.disp(:,2:end,i,j),options.display.dispRange);axis image;
            hold on
            d1 = plot(sweidata.lat(end,:),gate(i,1)*ones(1,nBeams-1),'g','Linewidth',2);
            d2 = plot(sweidata.lat(end,:),gate(i,2)*ones(1,nBeams-1),'g','Linewidth',2);
            hold off
            if j==1,colorbar,end
            set(p1,'alphadata',~isnan(sweidata.disp(:,2:end,i,j)))
            set(gca,'color',[0.4 0.4 0.4])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            title(sprintf('Raw Disp, t = %2.2f s',sweidata.acqTime(i)))
            if ~strcmpi(options.motionFilter.method,'Off')
                subplot(222)
                p2 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.disp_mf_push(:,2:end,i,j),options.display.dispRange);axis image;
                hold on
                d3 = plot(sweidata.lat(end,:),gate(i,1)*ones(1,nBeams-1),'g','Linewidth',2);
                d4 = plot(sweidata.lat(end,:),gate(i,2)*ones(1,nBeams-1),'g','Linewidth',2);
                hold off
                if j==1,colorbar,end
                set(p2,'alphadata',~isnan(sweidata.disp_mf_push(:,2:end,i,j)))
                set(gca,'color',[0.4 0.4 0.4])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
                title(sprintf('MF Disp, t = %2.2f s',sweidata.acqTime(i)))
            end
            
            subplot(223)
            p3 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.vel(:,2:end,i,j),options.display.velRange);axis image;
            hold on
            d5 = plot(sweidata.lat(end,:),gate(i,1)*ones(1,nBeams-1),'g','Linewidth',2);
            d6 = plot(sweidata.lat(end,:),gate(i,2)*ones(1,nBeams-1),'g','Linewidth',2);
            hold off
            if j==1,colorbar,end
            set(p3,'alphadata',~isnan(sweidata.vel(:,2:end,i,j)))
            set(gca,'color',[0.4 0.4 0.4])
            xlabel('Lateral (mm)')
            ylabel('Axial (mm)')
            title(sprintf('Raw Vel, t = %2.2f s',sweidata.acqTime(i)))
            if ~strcmpi(options.motionFilter.method,'Off')
                subplot(224)
                p4 = imagesc(sweidata.lat(end,:),sweidata.axial,sweidata.vel_mf_push(:,2:end,i,j),options.display.velRange);axis image;
                hold on
                d7 = plot(sweidata.lat(end,:),gate(i,1)*ones(1,nBeams-1),'g','Linewidth',2);
                d8 = plot(sweidata.lat(end,:),gate(i,2)*ones(1,nBeams-1),'g','Linewidth',2);
                hold off
                if j==1,colorbar,end
                set(p4,'alphadata',~isnan(sweidata.vel_mf_push(:,2:end,i,j)))
                set(gca,'color',[0.4 0.4 0.4])
                xlabel('Lateral (mm)')
                ylabel('Axial (mm)')
                title(sprintf('MF Vel, t = %2.2f s',sweidata.acqTime(i)))
            end
            pause(0.005)
        end
        waitforbuttonpress;
        val = double(get(gcf,'CurrentCharacter'));
        delete(d1),delete(d2),delete(d5),delete(d6)
        if ~strcmpi(options.motionFilter.method,'Off')
            delete(d3),delete(d4),,delete(d7),delete(d8)
        end
        if val==28, i=i-1; if i==0, i=1;end;end % Previous
        if val==31;% Same
            % Change dispRange and/or velRange
        end
        if val==29, i=i+1; end % Next
    end
    close(gcf)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        for i=1:options.display.n_pts
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
clear raw mf

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
