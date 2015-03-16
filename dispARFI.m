<<<<<<< HEAD
function [traced_gate] =  dispARFI(ecgdata,bdata,arfidata,arfidata_mf_pre,arfidata_mf_push,options,par)

%% Check for existance of traced gate
if isfield(arfidata,'traced_gate') && ~isempty(arfidata.traced_gate)
    prior_trace = 1;
    traced_gate = arfidata.traced_gate;
else
    prior_trace = 0;
    traced_gate = [];
end

dims = size(arfidata.disp);
ndepth = dims(1); nacqT = dims(2); ntrackT = dims(3);

edge = [arfidata.axial(1) arfidata.axial(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Gate Parameters
if prior_trace
    gate = repmat(arfidata.traced_gate,1,2) + repmat(options.display.gateOffset,length(arfidata.traced_gate),2) + repmat([-options.display.gateWidth/2 options.display.gateWidth/2],length(arfidata.traced_gate),1);
elseif ~prior_trace
    gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
    gate = repmat(gate,[nacqT 1]);
    traced_gate = [];
    if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

% Set Display Parameters
fig1 = figure(1);
set(fig1,'Name','TTE ARFI','UserData','main_fig')
if isunix
    set(fig1,'Position',[1201 386 1920 1070])
    dispPar.fsize = 16;
elseif ispc
    set(fig1,'units','normalized','outerposition',[0 0 1 1])
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

set(fig1,'Color',dispPar.fig)
trackPRF = 1000/par.priusec(1); % kHz

rng = options.display.disprange;

dispPar.cmap = colormap(parula);
dispPar.cmap(end,:) = dispPar.ax;

dispPar.corder = winter(options.display.n_pts);

% Auto displacement range based on quantiles
if isempty(options.display.disprange)
    if (options.motionFilter.enable && (strcmpi(options.motionFilter.method,'Polynomial') || strcmpi(options.motionFilter.method,'Both')))
        rng = quantile(arfidata_mf_push.disp(:),[0.25 0.75]);
    else
        rng = quantile(arfidata.disp(:),[0.05 0.95]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Bmode Cine
for i=1:size(bdata.bimg,3);
    ax11 = subplot('Position',[0.1 0.6 0.3 0.3]);
    imagesc(bdata.blat,bdata.bax,fliplr(bdata.bimg(:,:,i)));
    colormap(gray);axis image; freezeColors;
    hold(ax11,'on')
    r1 = rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2,'Parent',ax11);
    if prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linestyle','--','Linewidth',2,'Parent',ax11);
    elseif ~prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 options.display.gateWidth],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    end
    xlabel('Lateral (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    title(sprintf('TimeStamp = %s \nB-Mode Cine: Frame %d (t = %1.1f s)',num2str(options.timeStamp),i,bdata.t(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    set(ax11,'xcolor',dispPar.fig,'ycolor',dispPar.fig,'fontweight','bold','UserData','bmodecine_ax')
    hold(ax11,'off')
    
    %     xlim([-25 25]);ylim([max(edge(1)-15,arfidata.IQaxial(1)) min(edge(2)+15,arfidata.IQaxial(end))])
    
    pause(0.025)
    if i==size(bdata.bimg,3)
        % Display M-mode IQ
        ax12 = axes('Position',[0.45 0.7 0.52 0.2]);
        frame = abs(arfidata.IQ(:,:,1)); % Display first frame only
        frame = db(frame/max(frame(:)));
        temp = find(arfidata.IQaxial>edge(1)-7.5,1);
        if isempty(temp); idx(1) = 1; else idx(1) = temp; end; clear temp
        temp = find(arfidata.IQaxial>edge(2)+7.5,1);
        if isempty(temp); idx(2) = length(arfidata.IQaxial); else idx(2) = temp; end; clear temp
        imagesc(linspace(0,arfidata.acqTime(end),size(arfidata.IQ,2)),arfidata.IQaxial(idx(1):idx(2)),frame(idx(1):idx(2),:),options.display.IQrange)
        clear idx
        hold(ax12,'on')
        l1 = plot(linspace(0,arfidata.acqTime(end),nacqT),arfidata.axial(1)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l2 = plot(linspace(0,arfidata.acqTime(end),nacqT),arfidata.axial(end)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l3 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
        l4 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
        hold(ax12,'off')
        xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d\n Track PRF = %1.2f kHz, Push = %d x %d \\mus',par.isHarmonic,trackPRF,par.npush,par.pushDurationusec),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        colormap(gray); freezeColors;
        set(ax12,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','UserData','mmodeIQ_ax')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NaN out push reverb frames
arfidata = interpPushReverb(arfidata,options,par,'nan'); % NaN out push and reverb frames
if options.motionFilter.enable
    arfidata_mf_pre = interpPushReverb(arfidata_mf_pre,options,par,'nan'); % NaN out push and reverb frames
    arfidata_mf_push = interpPushReverb(arfidata_mf_push,options,par,'nan'); % NaN out push and reverb frames
end

% Coorelation mask filter
if options.display.cc_filt
    mask = arfidata.cc>options.display.cc_thresh;
else
    mask = [];
end

% Indices corresponding to median filter parameters
nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display M-mode ARFI
if options.motionFilter.enable
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata_mf_pre,arfidata_mf_push,par,gate,mask,dispPar,rng);
else % Motion Filter Disabled
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata,arfidata,par,gate,mask,dispPar,rng);
end

% NaN out displacements filtered out by cc_thresh
pre(pre==inf) = nan;
push(push==inf) = nan;
% arfidata.disp(mask==0) = nan;
% if options.motionFilter.enable
%     arfidata_mf_pre.disp(mask==0) = nan;
%     arfidata_mf_push.disp(mask==0) = nan;
% end

if isempty(ecgdata)
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incorporate ECG Data into this figure
if ~isempty(ecgdata)
    samples = zeros(1,nacqT);
    for i=1:nacqT
        samples(i) = ecgdata.arfi(find(ecgdata.arfi(:,1)>arfidata.acqTime(i),1,'first'),2);
    end
    ecgdata.arfi(:,2) = ecgdata.arfi(:,2)/max(ecgdata.arfi(:,2));
    
    ax15 = axes('Position',[0.5 0.1 0.4 0.2]);
    plot(ecgdata.arfi(:,1),ecgdata.arfi(:,2),'Linewidth',2);
    hold(ax15,'on')
    plot(arfidata.acqTime,samples,'gx','MarkerSize',8)
    pt = plot(arfidata.acqTime(1),samples(1),'ro','Parent',ax15,'Markersize',10,'Markerfacecolor','r');
    title('ECG Trace','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    axis tight
    set(ax15,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--','Userdata','ecg_ax')
    hold(ax15,'off')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate

for i=1:nacqT
    gate_idx(i,:) = [find(arfidata.axial>gate(i,1),1,'first') find(arfidata.axial<gate(i,2),1,'last')];
    pre_trace(i) = nanmean(pre(gate_idx(i,1):gate_idx(i,2),i));
    push_trace(i) = nanmean(push(gate_idx(i,1):gate_idx(i,2),i));
    idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts)); % Calculate indices for disp. vs. time plots
end

ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
plot(arfidata.acqTime,pre_trace,'-yo','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
hold(ax16,'on')
plot(arfidata.acqTime,push_trace,'-ro','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
grid on
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
ylim(rng)
xlim([0 max(arfidata.acqTime)])
set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
hold(ax16,'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trace Gate
if (options.dataflow.stream == 1 || options.dataflow.stream == 3)  % realTime or cluster
    trace_input = 'n';
    flatgate_input = 'y';
end

if prior_trace
    trace_input = input('\nDo you want to retrace the gate? (y/n) [n]: ','s');
    flatgate_input = input('\nDo you want to use a flat gate? (y/n) [n]: ','s');
elseif (~prior_trace && options.dataflow.stream == 2) % review
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
    l3 = plot(linspace(0,range(arfidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
    l4 = plot(linspace(0,range(arfidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
    if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

if strcmpi(trace_input,'y')
    delete(r2);delete(l3);delete(l4);
    fprintf(1,'\nReady to trace gate (GateWidth = %2.1f mm)...\nClick to define points, hit space to end tracing\n',options.display.gateWidth)
    
    % Change this to not require nacqT clicks!!
    for i=1:nacqT
        [x,y] = ginput(1);
        hold on
        mark(i) = plot(x,y,'yo','MarkerSize',10);
        gate(i,:) = (y + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        clear x y
    end
    delete(mark)
    hold(ax11,'on')
    r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    hold(ax12,'on')
    l3 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
    l4 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
    fprintf(1,'Gate Traced.\n')
    traced_gate = mean(gate,2);
    
    if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

% Display M-mode ARFI
if options.motionFilter.enable
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata_mf_pre,arfidata_mf_push,par,gate,mask,dispPar,rng);
else % Motion Filter Disabled
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata,arfidata,par,gate,mask,dispPar,rng);
end

% NaN out displacements filtered out by cc_thresh
pre(pre==inf) = nan;
push(push==inf) = nan;
% arfidata.disp(mask==0) = nan;
% if options.motionFilter.enable
%     arfidata_mf_pre.disp(mask==0) = nan;
%     arfidata_mf_push.disp(mask==0) = nan;
% end

% Compute Pre and Push Traces averaged axially over the gate
for i=1:nacqT
    gate_idx(i,:) = [find(arfidata.axial>gate(i,1),1,'first') find(arfidata.axial<gate(i,2),1,'last')];
    pre_trace(i) = nanmean(pre(gate_idx(i,1):gate_idx(i,2),i));
    push_trace(i) = nanmean(push(gate_idx(i,1):gate_idx(i,2),i));
    idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts)); % Calculate indices for disp. vs. time plots
end

delete(ax16)
ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
plot(arfidata.acqTime,pre_trace,'-yo','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
hold(ax16,'on')
plot(arfidata.acqTime,push_trace,'-ro','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
ylim(rng)
xlim([0 max(arfidata.acqTime)])
grid on
set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
hold(ax16,'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop to go through Acquisition Time

switch options.display.extras
    
    case 0
        extra1_input = input('\nDo you want to look at displacement vs. time plots? (y/n) [n]: ','s');
        extra2_input = 'n';
        if strcmpi(extra1_input,'y')
            options.display.extras = 1;
            extra2_input = input('\nDo you want to look at decorrelation/ raw motion plots? (y/n) [n]: ','s');
        end
        
        if strcmpi(extra2_input,'y')
            options.display.extras = 2;
        end
        
    case 1
        extra2_input = input('\nDo you want to look at decorrelation/ raw motion plots? (y/n) [n]: ','s');
        if strcmpi(extra2_input,'y')
            options.display.extras = 2;
        end
        
end

% Displacement vs. Time Plots
if options.display.extras > 0;
    
    i=1; skip = 0;
    fprintf(1,'\n\nPress Left/Right to move Back/Forward and Space to play through\n')
    
    while i<nacqT
        
        if ~isempty(ecgdata)
            delete(pt); %set(pt,'Visible','off')
        end
        
        if i==1
            delete(ax16)
            ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
            set(ax16,'ColorOrder',dispPar.corder)
        else
            cla(ax16)
            set(ax16,'ColorOrder',dispPar.corder)
        end
        
        % Extra Fig 1
        if options.motionFilter.enable
            plot(arfidata.trackTime(1:par.nref),squeeze(arfidata_mf_pre.disp(idx(i,:),i,1:par.nref)),'.--','Parent',ax16)
            hold(ax16,'on')
            plot(arfidata.trackTime(par.nref+1:end),squeeze(arfidata_mf_push.disp(idx(i,:),i,par.nref+1:end)),'.--','Parent',ax16)
            ylim(ax16,rng)
        else
            plot(arfidata.trackTime,squeeze(arfidata.disp(idx(i,:),i,:)),'.--','Parent',ax16)
            hold(ax16,'on')
            ylim(ax16,rng)
        end
        plot(arfidata.trackTime(idx_pre(i))*ones(1,10),linspace(-300,300,10),'y','linewidth',2,'Parent',ax16)
        plot(arfidata.trackTime(idx_push(i))*ones(1,10),linspace(-300,300,10),'g','linewidth',2,'Parent',ax16)
        title(sprintf('ARFI Displacement Profiles (within Depth Gate)\nPush # %d (t = %2.2f s)\nMotion Filter = %s',i,arfidata.acqTime(i),options.motionFilter.method*options.motionFilter.enable),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
        xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        xlim([arfidata.trackTime(1) arfidata.trackTime(end)])
        grid on
        set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
        hold(ax16,'off')
        
        if ~isempty(ecgdata)
            hold(ax15,'on')
            pt = plot(arfidata.acqTime(i),samples(i),'ro','Parent',ax15,'Markersize',10,'Markerfacecolor','r');
            hold(ax15,'off')
        end
        
        % Extra Fig 2
        if options.display.extras > 1
            
            fig2 = figure(2);
            set(fig2,'Name','IQ Traces & Correlation Map','UserData','extra_fig1')
            if isunix
                set(fig2,'Position',[-1198 580 1198 893])
            elseif ispc
                set(fig2,'units','normalized','outerposition',[0 0 1 1])
            end
            set(fig2,'Color',dispPar.fig)
            
            ax21=subplot(121); cla(ax21);
            temp = squeeze(arfidata.IQ(:,1+(i-1)*par.nBeams,:));
            I = real(temp); Q = imag(temp);
            factor = 5;
            D = size(I);
            D(1) = D(1).*factor;
            [Iup, Qup] = computeUpsampledIQdata(I,Q,factor);
            Iup = reshape(Iup,D); Qup = reshape(Qup,D);
            temp = db(abs(complex(Iup,Qup)));
            
            temp(:,par.nref+1:par.nref+par.npush+par.nreverb) = nan;
            
            axial_up = interp(arfidata.IQaxial,factor);
            
            offset = 2.5;
            for j=1:ntrackT; plot(axial_up,offset*(j-1)-temp(:,j)','b'); hold on; end; view(90,90); hold on;
            plot(gate(i,1)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
            plot(gate(i,2)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
            xlim(edge);ylim([-50 250])
            title(sprintf('Raw IQ: %d (t = %2.2f s)',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','color',dispPar.txt);
            xlabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Tracks','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);set(gca,'YTickLabel',[])
            set(ax21,'color',dispPar.ax + 0.25,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','Userdata','IQtraces_ax')
            
            ax22= subplot(122); cla(ax22);
            imagesc(arfidata.trackTime,arfidata.axial,squeeze(arfidata.cc(:,i,:)),[options.display.cc_thresh 1]);cb = colorbar; set(cb,'Color',dispPar.txt,'FontWeight','bold','UserData','cc_cb')
            title(sprintf('%s Correlation Coefficients',options.dispEst.ref_type),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);grid on;colormap(jet)
            hold on;plot(linspace(-25,25,100),gate(i,1)*ones(1,100),'g','linewidth',3);plot(linspace(-25,25,100),gate(i,2)*ones(1,100),'g','linewidth',3)
            xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold')
            set(ax22,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','ccmap_ax')
            
            raw = squeeze(arfidata.disp(:,i,:));
            mf = zeros(size(raw));
            if options.motionFilter.enable
                mf(:,1:par.nref) = squeeze(arfidata_mf_pre.disp(:,i,1:par.nref));
                mf(:,par.nref+1:end) = squeeze(arfidata_mf_push.disp(:,i,par.nref+1:end));
            end
            
            fig3 = figure(3);
            set(fig3,'Name','Displacement Data','UserData','extra_fig2')
            if isunix
                set(fig3,'Position',[-1195 -270 1198 848])
            elseif ispc
                set(fig3,'units','normalized','outerposition',[0 0 1 1])
            end
            set(fig3,'Color',dispPar.fig)
            
            ax31 = subplot(121); cla(ax31);
            imagesc(arfidata.trackTime,arfidata.axial,raw,[-150 150]);cb2 = colorbar; set(cb2,'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb_extra')
            xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            title(sprintf('Raw Displacement: Push %d\n Time = %1.2f s',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            hold on
            plot(options.display.t_disp_pre*ones(1,length(arfidata.axial)),arfidata.axial,'y','linewidth',2)
            plot(options.display.t_disp_push*ones(1,length(arfidata.axial)),arfidata.axial,'g','linewidth',2)
            l1 = plot(arfidata.trackTime,gate(i,1)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            l2 = plot(arfidata.trackTime,gate(i,2)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            set(ax31,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','rawdisp_ax')
            
            ax32 = subplot(122); cla(ax32)
            imagesc(arfidata.trackTime,arfidata.axial,mf,rng);cb3 = colorbar; set(cb3,'Color',dispPar.txt,'FontWeight','bold')
            xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            title(sprintf('MF Displacement: Push %d\n Time = %1.2f s',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            hold on
            plot(options.display.t_disp_pre*ones(1,length(arfidata.axial)),arfidata.axial,'y','linewidth',2)
            plot(options.display.t_disp_push*ones(1,length(arfidata.axial)),arfidata.axial,'g','linewidth',2)
            l3 = plot(arfidata.trackTime,gate(i,1)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            l4 = plot(arfidata.trackTime,gate(i,2)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            set(ax32,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','mfdisp_ax')
        end
        
        set(0,'CurrentFigure',fig1)
        
        if ~skip
            w = waitforbuttonpress;
            val = double(get(gcf,'CurrentCharacter'));
            if val==28
                i=i-1;
            elseif val==29
                i=i+1;
            elseif val==32
                skip = 1;
            else
                i=i+1;
            end
            if i<1;i=1;end
        end
        
        if skip
            pause(0.05)
            i=i+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate

set(0,'CurrentFigure',fig1)
for i=1:nacqT
    gate_idx(i,:) = [find(arfidata.axial>gate(i,1),1,'first') find(arfidata.axial<gate(i,2),1,'last')];
    pre_trace(i) = nanmean(pre(gate_idx(i,1):gate_idx(i,2),i));
    push_trace(i) = nanmean(push(gate_idx(i,1):gate_idx(i,2),i));
    idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts)); % Calculate indices for disp. vs. time plots
end

delete(ax16)
ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
plot(arfidata.acqTime,pre_trace,'-yo','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
hold(ax16,'on')
plot(arfidata.acqTime,push_trace,'-ro','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
grid on
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
ylim(rng)
xlim([0 max(arfidata.acqTime)])
set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
hold(ax16,'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Panel of Axial Displacement and Axial Velocity Data vs. ECG

extra3_input = 'n';

if options.display.extras ~= -1
    extra3_input = input('\nDo you want to look at temporally aligned panel? (y/n) [n]: ','s');
elseif options.display.extras >2
    extra_3_input = 'y';
end

if strcmpi(extra3_input,'y')
    track_dt = median(diff(arfidata.trackTime))*1e-3;
    acq_dt = median(diff(arfidata.acqTime));
    tt = [0:track_dt:acq_dt+arfidata.acqTime(end)];
    
    raw = arfidata.disp;
    raw_vel = nan(size(raw));
    raw_vel(:,:,1:size(raw,3)-1) = 100.*(1e-6.*diff(raw,1,3))./(track_dt);
    for i=1:nacqT;raw_vel(:,i,:) = medfilt2(squeeze(raw_vel(:,i,:)),[1 5]);end
    raw(mask==0) = nan;
    raw_vel(mask==0) = nan;
    
    raw = permute(raw,[1 3 2]);
    raw_vel = permute(raw_vel,[1 3 2]);
    
    % mf = nan(size(arfidata.disp));
    % mf(:,:,1:par.nref) = arfidata_mf_pre.disp(:,:,1:par.nref);
    % mf(:,:,par.nref+1:end) = arfidata_mf_push.disp(:,:,par.nref+1:end);
    % mf(mask==0) = nan;
    % mf = permute(mf,[1 3 2]);
    
    for i=1:nacqT
        top(i) = find(arfidata.axial>gate(i,1),1);
        bot(i) = find(arfidata.axial>gate(i,2),1);
        raw([top(i)-2:top(i)],:,i) = inf;
        raw([bot(i):bot(i)+2],:,i) = inf;
        
        raw_vel([top(i)-1:top(i)],:,i) = inf;
        raw_vel([bot(i):bot(i)+1],:,i) = inf;
        
        %     mf([top(i)-1:top(i)],:,i) = inf;
        %     mf([bot(i):bot(i)+1],:,i) = inf;
    end
    
    raw_panel = nan(size(arfidata.disp,1),length(tt));
    raw_vel_panel = nan(size(arfidata.disp,1),length(tt));
    % mf_panel = nan(size(arfidata.disp,1),length(tt));
    
    for i=1:nacqT
        start_idx = find(tt>arfidata.acqTime(i),1);
        raw_panel(:,start_idx:start_idx+size(raw,2)-1) = raw(:,:,i);
        raw_vel_panel(:,start_idx:start_idx+size(raw_vel,2)-1) = raw_vel(:,:,i);
        %     mf_panel(:,start_idx:start_idx+size(mf,2)-1) = mf(:,:,i);
    end
    
    raw_panel = raw_panel(min(top)-10:max(bot)+10,:);
    raw_vel_panel = raw_vel_panel(min(top)-10:max(bot)+10,:);
    % mf_panel = mf_panel(min(top)-10:max(bot)+10,:);
    
    width = 500; % ms
    delta = 50; % ms
    win = [0 width];
    
    close all
    fig4 = figure(4);
    set(fig4,'units','normalized','outerposition',[0 0 1 1])
    set(fig4,'Color',dispPar.fig)
    
    ax41 = axes('Position',[0.030 0.67 0.95 0.25]);
    pn1 = imagesc(1000*tt,arfidata.axial(min(top)-10:max(bot)+10),raw_panel,[-200 200]);
    xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    set(ax41,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','disppanel_ax')
    set(pn1,'AlphaData',~isnan(raw_panel))
    ylim([min(gate(:,1)) max(gate(:,2))])
    xlim(win)
    colormap(jet)
    colorbar
    
    ax42 = axes('Position',[0.030 0.07 0.95 0.25]);
    % pn2 = imagesc(1000*tt,arfidata.axial(min(top)-10:max(bot)+10),mf_panel,[-2 15]);
    pn2 = imagesc(1000*tt,arfidata.axial(min(top)-10:max(bot)+10),raw_vel_panel,[-5 5]);
    xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    set(ax42,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','velpanel_ax')
    % set(pn2,'AlphaData',~isnan(mf_panel))
    set(pn2,'AlphaData',~isnan(raw_vel_panel))
    ylim([min(gate(:,1)) max(gate(:,2))])
    xlim(win)
    colormap(jet)
    colorbar
    
    if ~isempty(ecgdata)
        ax43 = axes('Position',[0.030 0.40 0.90 0.20]);
        plot(1000*ecgdata.arfi(:,1),ecgdata.arfi(:,2),'linewidth',5);
        xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
        set(ax43,'color',[0.15 0.15 0.15],'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'yTickLabel',[],'xgrid','on','UserData','ecgpanel_ax')
        xlim(win)
    end
    val = 0;
    while val~=32
        if val==28
            win = [win(1)-delta win(2)-delta];
        end
        if val==29
            win = [win(1)+delta win(2)+delta];
        end
        if win(1)<0;win = [0 width];end
        if win(2)>arfidata.acqTime*1000;win = [arfidata.acqTime*1000-width arfidata.acqTime*1000];end
        set(ax41,'xlim',win); set(ax42,'xlim',win);
        if ~isempty(ecgdata); set(ax43,'xlim',win); end
        w = waitforbuttonpress;
        val = double(get(gcf,'CurrentCharacter'));
    end
end

%%
% % Extra Fig
% if options.display.extras > 0
%     figure(101)
%     if isunix
%         set(101,'Position',[1210 750 1909 624]);
%     elseif ispc
%         set(101,'units','normalized','outerposition',[0 0 1 1])
%     end
%     set(101,'Color',dispPar.fig)
%
%     temp_gate = zeros(ntrackT*nacqT,2);
%     for i=1:nacqT
%         temp_gate(1+ntrackT*(i-1):ntrackT*i,:) = repmat(gate(i,:),ntrackT,1);
%     end
%     disps = reshape(permute(arfidata.disp,[1 3 2]),size(arfidata.disp,1),[]);
%     vels = reshape(permute(diff(arfidata.disp,1,3),[1 3 2]),size(arfidata.disp,1),[]);
%
%     rng = quantile(arfidata.disp(:),[0.01 0.99]);
%
%     im1 = imagesc(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),arfidata.axial,disps,rng);colorbar
%     hold on
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,1),dispPar.txt,'linewidth',5)
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,2),dispPar.txt,'linewidth',5)
%     set(gca,'color',dispPar.ax);
%     set(im1,'alphaData',~isnan(disps))
%
%     figure(102)
%     if isunix
%         set(102,'Position',[1210 750 1909 624]);
%     elseif ispc
%         set(102,'units','normalized','outerposition',[0 0 1 1])
%     end
%     set(102,'Color',dispPar.fig)
%
%     im2 = imagesc(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),arfidata.axial,vels,rng/25);colorbar
%     hold on
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,1),dispPar.txt,'linewidth',5)
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,2),dispPar.txt,'linewidth',5)
%     set(gca,'color',dispPar.ax);
%     set(im2,'alphaData',~isnan(vels))
% end
%
% figure(103);set(103,'Position',[1210 390 1909 624])
% imagesc([],arfidata.axial,reshape(permute(diff(diff(arfidata.disp,1,3),1,3),[1 3 2]),size(arfidata.disp,1),[]),rng/100);colorbar

=======
function [traced_gate] =  dispARFI(ecgdata,bdata,arfidata,arfidata_mf_pre,arfidata_mf_push,options,par)

%% Check for existance of traced gate
if isfield(arfidata,'traced_gate') && ~isempty(arfidata.traced_gate)
    prior_trace = 1;
    traced_gate = arfidata.traced_gate;
else
    prior_trace = 0;
    traced_gate = [];
end

dims = size(arfidata.disp);
ndepth = dims(1); nacqT = dims(2); ntrackT = dims(3);

edge = [arfidata.axial(1) arfidata.axial(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Gate Parameters
if prior_trace
    gate = repmat(arfidata.traced_gate,1,2) + repmat(options.display.gateOffset,length(arfidata.traced_gate),2) + repmat([-options.display.gateWidth/2 options.display.gateWidth/2],length(arfidata.traced_gate),1);
elseif ~prior_trace
    gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
    gate = repmat(gate,[nacqT 1]);
    traced_gate = [];
    if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

% Set Display Parameters
fig1 = figure(1);
set(fig1,'Name','TTE ARFI','UserData','main_fig')
if isunix
    set(fig1,'Position',[1201 386 1920 1070])
    dispPar.fsize = 16;
elseif ispc
    set(fig1,'units','normalized','outerposition',[0 0 1 1])
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

set(fig1,'Color',dispPar.fig)
trackPRF = 1000/par.priusec(1); % kHz

rng = options.display.disprange;

dispPar.cmap = colormap(parula);
dispPar.cmap(end,:) = dispPar.ax;

dispPar.corder = winter(options.display.n_pts);

% Auto displacement range based on quantiles
if isempty(options.display.disprange)
    if (options.motionFilter.enable && (strcmpi(options.motionFilter.method,'Polynomial') || strcmpi(options.motionFilter.method,'Both')))
        rng = quantile(arfidata_mf_push.disp(:),[0.25 0.75]);
    else
        rng = quantile(arfidata.disp(:),[0.05 0.95]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Bmode Cine
for i=1:size(bdata.bimg,3);
    ax11 = subplot('Position',[0.1 0.6 0.3 0.3]);
    imagesc(bdata.blat,bdata.bax,fliplr(bdata.bimg(:,:,i)));
    colormap(gray);axis image; freezeColors;
    hold(ax11,'on')
    r1 = rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2,'Parent',ax11);
    if prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linestyle','--','Linewidth',2,'Parent',ax11);
    elseif ~prior_trace
        r2 = rectangle('Position',[-2 min(gate(:)) 4 options.display.gateWidth],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    end
    xlabel('Lateral (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    title(sprintf('TimeStamp = %s \nB-Mode Cine: Frame %d (t = %1.1f s)',num2str(options.timeStamp),i,bdata.t(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    set(ax11,'xcolor',dispPar.fig,'ycolor',dispPar.fig,'fontweight','bold','UserData','bmodecine_ax')
    hold(ax11,'off')
    
    %     xlim([-25 25]);ylim([max(edge(1)-15,arfidata.IQaxial(1)) min(edge(2)+15,arfidata.IQaxial(end))])
    
    pause(0.025)
    if i==size(bdata.bimg,3)
        % Display M-mode IQ
        ax12 = axes('Position',[0.45 0.7 0.52 0.2]);
        frame = abs(arfidata.IQ(:,:,1)); % Display first frame only
        frame = db(frame/max(frame(:)));
        temp = find(arfidata.IQaxial>edge(1)-7.5,1);
        if isempty(temp); idx(1) = 1; else idx(1) = temp; end; clear temp
        temp = find(arfidata.IQaxial>edge(2)+7.5,1);
        if isempty(temp); idx(2) = length(arfidata.IQaxial); else idx(2) = temp; end; clear temp
        imagesc(linspace(0,arfidata.acqTime(end),size(arfidata.IQ,2)),arfidata.IQaxial(idx(1):idx(2)),frame(idx(1):idx(2),:),options.display.IQrange)
        clear idx
        hold(ax12,'on')
        l1 = plot(linspace(0,arfidata.acqTime(end),nacqT),arfidata.axial(1)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l2 = plot(linspace(0,arfidata.acqTime(end),nacqT),arfidata.axial(end)*ones(1,nacqT),'b','Linewidth',2,'Parent',ax12);
        l3 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
        l4 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
        hold(ax12,'off')
        xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d\n Track PRF = %1.2f kHz, Push = %d x %d \\mus',par.isHarmonic,trackPRF,par.npush,par.pushDurationusec),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        colormap(gray); freezeColors;
        set(ax12,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','UserData','mmodeIQ_ax')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NaN out push reverb frames
arfidata = interpPushReverb(arfidata,options,par,'nan'); % NaN out push and reverb frames
if options.motionFilter.enable
    arfidata_mf_pre = interpPushReverb(arfidata_mf_pre,options,par,'nan'); % NaN out push and reverb frames
    arfidata_mf_push = interpPushReverb(arfidata_mf_push,options,par,'nan'); % NaN out push and reverb frames
end

% Coorelation mask filter
if options.display.cc_filt
    mask = arfidata.cc>options.display.cc_thresh;
else
    mask = [];
end

% Indices corresponding to median filter parameters
nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display M-mode ARFI
if options.motionFilter.enable
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata_mf_pre,arfidata_mf_push,par,gate,mask,dispPar,rng);
else % Motion Filter Disabled
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata,arfidata,par,gate,mask,dispPar,rng);
end

% NaN out displacements filtered out by cc_thresh
pre(pre==inf) = nan;
push(push==inf) = nan;
% arfidata.disp(mask==0) = nan;
% if options.motionFilter.enable
%     arfidata_mf_pre.disp(mask==0) = nan;
%     arfidata_mf_push.disp(mask==0) = nan;
% end

if isempty(ecgdata)
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incorporate ECG Data into this figure
if ~isempty(ecgdata)
    samples = zeros(1,nacqT);
    for i=1:nacqT
        samples(i) = ecgdata.arfi(find(ecgdata.arfi(:,1)>arfidata.acqTime(i),1,'first'),2);
    end
    ecgdata.arfi(:,2) = ecgdata.arfi(:,2)/max(ecgdata.arfi(:,2));
    
    ax15 = axes('Position',[0.5 0.1 0.4 0.2]);
    plot(ecgdata.arfi(:,1),ecgdata.arfi(:,2),'Linewidth',2);
    hold(ax15,'on')
    plot(arfidata.acqTime,samples,'gx','MarkerSize',8)
    pt = plot(arfidata.acqTime(1),samples(1),'ro','Parent',ax15,'Markersize',10,'Markerfacecolor','r');
    title('ECG Trace','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    axis tight
    set(ax15,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--','Userdata','ecg_ax')
    hold(ax15,'off')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate

for i=1:nacqT
    gate_idx(i,:) = [find(arfidata.axial>gate(i,1),1,'first') find(arfidata.axial<gate(i,2),1,'last')];
    pre_trace(i) = nanmean(pre(gate_idx(i,1):gate_idx(i,2),i));
    push_trace(i) = nanmean(push(gate_idx(i,1):gate_idx(i,2),i));
    idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts)); % Calculate indices for disp. vs. time plots
end

ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
plot(arfidata.acqTime,pre_trace,'-yo','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
hold(ax16,'on')
plot(arfidata.acqTime,push_trace,'-ro','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
grid on
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
ylim(rng)
xlim([0 max(arfidata.acqTime)])
set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
hold(ax16,'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trace Gate
if (options.dataflow.stream == 1 || options.dataflow.stream == 3)  % realTime or cluster
    trace_input = 'n';
    flatgate_input = 'y';
end

if prior_trace
    trace_input = input('\nDo you want to retrace the gate? (y/n) [n]: ','s');
    flatgate_input = input('\nDo you want to use a flat gate? (y/n) [n]: ','s');
elseif (~prior_trace && options.dataflow.stream == 2) % review
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
    l3 = plot(linspace(0,range(arfidata.lat(:))*nacqT,nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
    l4 = plot(linspace(0,range(arfidata.lat(:))*nacqT,nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
    if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

if strcmpi(trace_input,'y')
    delete(r2);delete(l3);delete(l4);
    fprintf(1,'\nReady to trace gate (GateWidth = %2.1f mm)...\nClick to define points, hit space to end tracing\n',options.display.gateWidth)
    
    % Change this to not require nacqT clicks!!
    for i=1:nacqT
        [x,y] = ginput(1);
        hold on
        mark(i) = plot(x,y,'yo','MarkerSize',10);
        gate(i,:) = (y + [-options.display.gateWidth/2 options.display.gateWidth/2]);
        clear x y
    end
    delete(mark)
    hold(ax11,'on')
    r2 = rectangle('Position',[-2 min(gate(:)) 4 range(gate(:))],'EdgeColor','g','Linewidth',2,'Parent',ax11);
    hold(ax12,'on')
    l3 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,1),'g','Linewidth',2,'Parent',ax12);
    l4 = plot(linspace(0,arfidata.acqTime(end),nacqT),gate(:,2),'g','Linewidth',2,'Parent',ax12);
    fprintf(1,'Gate Traced.\n')
    traced_gate = mean(gate,2);
    
    if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
end

% Display M-mode ARFI
if options.motionFilter.enable
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata_mf_pre,arfidata_mf_push,par,gate,mask,dispPar,rng);
else % Motion Filter Disabled
    [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,arfidata,arfidata,par,gate,mask,dispPar,rng);
end

% NaN out displacements filtered out by cc_thresh
pre(pre==inf) = nan;
push(push==inf) = nan;
% arfidata.disp(mask==0) = nan;
% if options.motionFilter.enable
%     arfidata_mf_pre.disp(mask==0) = nan;
%     arfidata_mf_push.disp(mask==0) = nan;
% end

% Compute Pre and Push Traces averaged axially over the gate
for i=1:nacqT
    gate_idx(i,:) = [find(arfidata.axial>gate(i,1),1,'first') find(arfidata.axial<gate(i,2),1,'last')];
    pre_trace(i) = nanmean(pre(gate_idx(i,1):gate_idx(i,2),i));
    push_trace(i) = nanmean(push(gate_idx(i,1):gate_idx(i,2),i));
    idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts)); % Calculate indices for disp. vs. time plots
end

delete(ax16)
ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
plot(arfidata.acqTime,pre_trace,'-yo','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
hold(ax16,'on')
plot(arfidata.acqTime,push_trace,'-ro','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
ylim(rng)
xlim([0 max(arfidata.acqTime)])
grid on
set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
hold(ax16,'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop to go through Acquisition Time

switch options.display.extras
    
    case 0
        extra1_input = input('\nDo you want to look at displacement vs. time plots? (y/n) [n]: ','s');
        extra2_input = 'n';
        if strcmpi(extra1_input,'y')
            options.display.extras = 1;
            extra2_input = input('\nDo you want to look at decorrelation/ raw motion plots? (y/n) [n]: ','s');
        end
        
        if strcmpi(extra2_input,'y')
            options.display.extras = 2;
        end
        
    case 1
        extra2_input = input('\nDo you want to look at decorrelation/ raw motion plots? (y/n) [n]: ','s');
        if strcmpi(extra2_input,'y')
            options.display.extras = 2;
        end
        
end

% Displacement vs. Time Plots
if options.display.extras > 0;
    
    i=1; skip = 0;
    fprintf(1,'\n\nPress Left/Right to move Back/Forward and Space to play through\n')
    
    while i<nacqT
        
        if ~isempty(ecgdata)
            delete(pt); %set(pt,'Visible','off')
        end
        
        if i==1
            delete(ax16)
            ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
            set(ax16,'ColorOrder',dispPar.corder)
        else
            cla(ax16)
            set(ax16,'ColorOrder',dispPar.corder)
        end
        
        % Extra Fig 1
        if options.motionFilter.enable
            plot(arfidata.trackTime(1:par.nref),squeeze(arfidata_mf_pre.disp(idx(i,:),i,1:par.nref)),'.--','Parent',ax16)
            hold(ax16,'on')
            plot(arfidata.trackTime(par.nref+1:end),squeeze(arfidata_mf_push.disp(idx(i,:),i,par.nref+1:end)),'.--','Parent',ax16)
            ylim(ax16,rng)
        else
            plot(arfidata.trackTime,squeeze(arfidata.disp(idx(i,:),i,:)),'.--','Parent',ax16)
            hold(ax16,'on')
            ylim(ax16,rng)
        end
        plot(arfidata.trackTime(idx_pre(i))*ones(1,10),linspace(-300,300,10),'y','linewidth',2,'Parent',ax16)
        plot(arfidata.trackTime(idx_push(i))*ones(1,10),linspace(-300,300,10),'g','linewidth',2,'Parent',ax16)
        title(sprintf('ARFI Displacement Profiles (within Depth Gate)\nPush # %d (t = %2.2f s)\nMotion Filter = %s',i,arfidata.acqTime(i),options.motionFilter.method*options.motionFilter.enable),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
        xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        xlim([arfidata.trackTime(1) arfidata.trackTime(end)])
        grid on
        set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
        hold(ax16,'off')
        
        if ~isempty(ecgdata)
            hold(ax15,'on')
            pt = plot(arfidata.acqTime(i),samples(i),'ro','Parent',ax15,'Markersize',10,'Markerfacecolor','r');
            hold(ax15,'off')
        end
        
        % Extra Fig 2
        if options.display.extras > 1
            
            fig2 = figure(2);
            set(fig2,'Name','IQ Traces & Correlation Map','UserData','extra_fig1')
            if isunix
                set(fig2,'Position',[-1198 580 1198 893])
            elseif ispc
                set(fig2,'units','normalized','outerposition',[0 0 1 1])
            end
            set(fig2,'Color',dispPar.fig)
            
            ax21=subplot(121); cla(ax21);
            temp = squeeze(arfidata.IQ(:,1+(i-1)*par.nBeams,:));
            I = real(temp); Q = imag(temp);
            factor = 5;
            D = size(I);
            D(1) = D(1).*factor;
            [Iup, Qup] = computeUpsampledIQdata(I,Q,factor);
            Iup = reshape(Iup,D); Qup = reshape(Qup,D);
            temp = db(abs(complex(Iup,Qup)));
            
            temp(:,par.nref+1:par.nref+par.npush+par.nreverb) = nan;
            
            axial_up = interp(arfidata.IQaxial,factor);
            
            offset = 2.5;
            for j=1:ntrackT; plot(axial_up,offset*(j-1)-temp(:,j)','b'); hold on; end; view(90,90); hold on;
            plot(gate(i,1)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
            plot(gate(i,2)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
            xlim(edge);ylim([-50 250])
            title(sprintf('Raw IQ: %d (t = %2.2f s)',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','color',dispPar.txt);
            xlabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Tracks','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);set(gca,'YTickLabel',[])
            set(ax21,'color',dispPar.ax + 0.25,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','Userdata','IQtraces_ax')
            
            ax22= subplot(122); cla(ax22);
            imagesc(arfidata.trackTime,arfidata.axial,squeeze(arfidata.cc(:,i,:)),[options.display.cc_thresh 1]);cb = colorbar; set(cb,'Color',dispPar.txt,'FontWeight','bold','UserData','cc_cb')
            title(sprintf('%s Correlation Coefficients',options.dispEst.ref_type),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);grid on;colormap(jet)
            hold on;plot(linspace(-25,25,100),gate(i,1)*ones(1,100),'g','linewidth',3);plot(linspace(-25,25,100),gate(i,2)*ones(1,100),'g','linewidth',3)
            xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold')
            set(ax22,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','ccmap_ax')
            
            raw = squeeze(arfidata.disp(:,i,:));
            mf = zeros(size(raw));
            if options.motionFilter.enable
                mf(:,1:par.nref) = squeeze(arfidata_mf_pre.disp(:,i,1:par.nref));
                mf(:,par.nref+1:end) = squeeze(arfidata_mf_push.disp(:,i,par.nref+1:end));
            end
            
            fig3 = figure(3);
            set(fig3,'Name','Displacement Data','UserData','extra_fig2')
            if isunix
                set(fig3,'Position',[-1195 -270 1198 848])
            elseif ispc
                set(fig3,'units','normalized','outerposition',[0 0 1 1])
            end
            set(fig3,'Color',dispPar.fig)
            
            ax31 = subplot(121); cla(ax31);
            imagesc(arfidata.trackTime,arfidata.axial,raw,[-150 150]);cb2 = colorbar; set(cb2,'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb_extra')
            xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            title(sprintf('Raw Displacement: Push %d\n Time = %1.2f s',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            hold on
            plot(options.display.t_disp_pre*ones(1,length(arfidata.axial)),arfidata.axial,'y','linewidth',2)
            plot(options.display.t_disp_push*ones(1,length(arfidata.axial)),arfidata.axial,'g','linewidth',2)
            l1 = plot(arfidata.trackTime,gate(i,1)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            l2 = plot(arfidata.trackTime,gate(i,2)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            set(ax31,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','rawdisp_ax')
            
            ax32 = subplot(122); cla(ax32)
            imagesc(arfidata.trackTime,arfidata.axial,mf,rng);cb3 = colorbar; set(cb3,'Color',dispPar.txt,'FontWeight','bold')
            xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            title(sprintf('MF Displacement: Push %d\n Time = %1.2f s',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
            hold on
            plot(options.display.t_disp_pre*ones(1,length(arfidata.axial)),arfidata.axial,'y','linewidth',2)
            plot(options.display.t_disp_push*ones(1,length(arfidata.axial)),arfidata.axial,'g','linewidth',2)
            l3 = plot(arfidata.trackTime,gate(i,1)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            l4 = plot(arfidata.trackTime,gate(i,2)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
            set(ax32,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','mfdisp_ax')
        end
        
        set(0,'CurrentFigure',fig1)
        
        if ~skip
            w = waitforbuttonpress;
            val = double(get(gcf,'CurrentCharacter'));
            if val==28
                i=i-1;
            elseif val==29
                i=i+1;
            elseif val==32
                skip = 1;
            else
                i=i+1;
            end
            if i<1;i=1;end
        end
        
        if skip
            pause(0.05)
            i=i+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate

set(0,'CurrentFigure',fig1)
for i=1:nacqT
    gate_idx(i,:) = [find(arfidata.axial>gate(i,1),1,'first') find(arfidata.axial<gate(i,2),1,'last')];
    pre_trace(i) = nanmean(pre(gate_idx(i,1):gate_idx(i,2),i));
    push_trace(i) = nanmean(push(gate_idx(i,1):gate_idx(i,2),i));
    idx(i,:) = ceil(linspace(gate_idx(i,1),gate_idx(i,2),options.display.n_pts)); % Calculate indices for disp. vs. time plots
end

delete(ax16)
ax16 = axes('Position',[0.05 0.15 0.4 0.3]);
plot(arfidata.acqTime,pre_trace,'-yo','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
hold(ax16,'on')
plot(arfidata.acqTime,push_trace,'-ro','Parent',ax16,'linewidth',3,'MarkerFaceColor','k');
grid on
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Parent',ax16,'Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',ax16)
ylim(rng)
xlim([0 max(arfidata.acqTime)])
set(ax16,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');
hold(ax16,'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Panel of Axial Displacement and Axial Velocity Data vs. ECG

extra3_input = 'n';

if options.display.extras ~= -1
    extra3_input = input('\nDo you want to look at temporally aligned panel? (y/n) [n]: ','s');
elseif options.display.extras >2
    extra_3_input = 'y';
end

if strcmpi(extra3_input,'y')
    track_dt = median(diff(arfidata.trackTime))*1e-3;
    acq_dt = median(diff(arfidata.acqTime));
    tt = [0:track_dt:acq_dt+arfidata.acqTime(end)];
    
    raw = arfidata.disp;
    raw_vel = nan(size(raw));
    raw_vel(:,:,1:size(raw,3)-1) = 100.*(1e-6.*diff(raw,1,3))./(track_dt);
    for i=1:nacqT;raw_vel(:,i,:) = medfilt2(squeeze(raw_vel(:,i,:)),[1 5]);end
    raw(mask==0) = nan;
    raw_vel(mask==0) = nan;
    
    raw = permute(raw,[1 3 2]);
    raw_vel = permute(raw_vel,[1 3 2]);
    
    % mf = nan(size(arfidata.disp));
    % mf(:,:,1:par.nref) = arfidata_mf_pre.disp(:,:,1:par.nref);
    % mf(:,:,par.nref+1:end) = arfidata_mf_push.disp(:,:,par.nref+1:end);
    % mf(mask==0) = nan;
    % mf = permute(mf,[1 3 2]);
    
    for i=1:nacqT
        top(i) = find(arfidata.axial>gate(i,1),1);
        bot(i) = find(arfidata.axial>gate(i,2),1);
        raw([top(i)-2:top(i)],:,i) = inf;
        raw([bot(i):bot(i)+2],:,i) = inf;
        
        raw_vel([top(i)-1:top(i)],:,i) = inf;
        raw_vel([bot(i):bot(i)+1],:,i) = inf;
        
        %     mf([top(i)-1:top(i)],:,i) = inf;
        %     mf([bot(i):bot(i)+1],:,i) = inf;
    end
    
    raw_panel = nan(size(arfidata.disp,1),length(tt));
    raw_vel_panel = nan(size(arfidata.disp,1),length(tt));
    % mf_panel = nan(size(arfidata.disp,1),length(tt));
    
    for i=1:nacqT
        start_idx = find(tt>arfidata.acqTime(i),1);
        raw_panel(:,start_idx:start_idx+size(raw,2)-1) = raw(:,:,i);
        raw_vel_panel(:,start_idx:start_idx+size(raw_vel,2)-1) = raw_vel(:,:,i);
        %     mf_panel(:,start_idx:start_idx+size(mf,2)-1) = mf(:,:,i);
    end
    
    raw_panel = raw_panel(min(top)-10:max(bot)+10,:);
    raw_vel_panel = raw_vel_panel(min(top)-10:max(bot)+10,:);
    % mf_panel = mf_panel(min(top)-10:max(bot)+10,:);
    
    width = 500; % ms
    delta = 50; % ms
    win = [0 width];
    
    close all
    fig4 = figure(4);
    set(fig4,'units','normalized','outerposition',[0 0 1 1])
    set(fig4,'Color',dispPar.fig)
    
    ax41 = axes('Position',[0.030 0.67 0.95 0.25]);
    pn1 = imagesc(1000*tt,arfidata.axial(min(top)-10:max(bot)+10),raw_panel,[-200 200]);
    xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    set(ax41,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','disppanel_ax')
    set(pn1,'AlphaData',~isnan(raw_panel))
    ylim([min(gate(:,1)) max(gate(:,2))])
    xlim(win)
    colormap(jet)
    colorbar
    
    ax42 = axes('Position',[0.030 0.07 0.95 0.25]);
    % pn2 = imagesc(1000*tt,arfidata.axial(min(top)-10:max(bot)+10),mf_panel,[-2 15]);
    pn2 = imagesc(1000*tt,arfidata.axial(min(top)-10:max(bot)+10),raw_vel_panel,[-5 5]);
    xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
    set(ax42,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','velpanel_ax')
    % set(pn2,'AlphaData',~isnan(mf_panel))
    set(pn2,'AlphaData',~isnan(raw_vel_panel))
    ylim([min(gate(:,1)) max(gate(:,2))])
    xlim(win)
    colormap(jet)
    colorbar
    
    if ~isempty(ecgdata)
        ax43 = axes('Position',[0.030 0.40 0.90 0.20]);
        plot(1000*ecgdata.arfi(:,1),ecgdata.arfi(:,2),'linewidth',5);
        xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.fig)
        set(ax43,'color',[0.15 0.15 0.15],'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'yTickLabel',[],'xgrid','on','UserData','ecgpanel_ax')
        xlim(win)
    end
    val = 0;
    while val~=32
        if val==28
            win = [win(1)-delta win(2)-delta];
        end
        if val==29
            win = [win(1)+delta win(2)+delta];
        end
        if win(1)<0;win = [0 width];end
        if win(2)>arfidata.acqTime*1000;win = [arfidata.acqTime*1000-width arfidata.acqTime*1000];end
        set(ax41,'xlim',win); set(ax42,'xlim',win);
        if ~isempty(ecgdata); set(ax43,'xlim',win); end
        w = waitforbuttonpress;
        val = double(get(gcf,'CurrentCharacter'));
    end
end

%%
% % Extra Fig
% if options.display.extras > 0
%     figure(101)
%     if isunix
%         set(101,'Position',[1210 750 1909 624]);
%     elseif ispc
%         set(101,'units','normalized','outerposition',[0 0 1 1])
%     end
%     set(101,'Color',dispPar.fig)
%
%     temp_gate = zeros(ntrackT*nacqT,2);
%     for i=1:nacqT
%         temp_gate(1+ntrackT*(i-1):ntrackT*i,:) = repmat(gate(i,:),ntrackT,1);
%     end
%     disps = reshape(permute(arfidata.disp,[1 3 2]),size(arfidata.disp,1),[]);
%     vels = reshape(permute(diff(arfidata.disp,1,3),[1 3 2]),size(arfidata.disp,1),[]);
%
%     rng = quantile(arfidata.disp(:),[0.01 0.99]);
%
%     im1 = imagesc(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),arfidata.axial,disps,rng);colorbar
%     hold on
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,1),dispPar.txt,'linewidth',5)
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,2),dispPar.txt,'linewidth',5)
%     set(gca,'color',dispPar.ax);
%     set(im1,'alphaData',~isnan(disps))
%
%     figure(102)
%     if isunix
%         set(102,'Position',[1210 750 1909 624]);
%     elseif ispc
%         set(102,'units','normalized','outerposition',[0 0 1 1])
%     end
%     set(102,'Color',dispPar.fig)
%
%     im2 = imagesc(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),arfidata.axial,vels,rng/25);colorbar
%     hold on
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,1),dispPar.txt,'linewidth',5)
%     plot(linspace(0,arfidata.acqTime(end),nacqT*ntrackT),temp_gate(:,2),dispPar.txt,'linewidth',5)
%     set(gca,'color',dispPar.ax);
%     set(im2,'alphaData',~isnan(vels))
% end
%
% figure(103);set(103,'Position',[1210 390 1909 624])
% imagesc([],arfidata.axial,reshape(permute(diff(diff(arfidata.disp,1,3),1,3),[1 3 2]),size(arfidata.disp,1),[]),rng/100);colorbar

>>>>>>> afa77557451dedcbdaaf013db29dc3d36822037b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%