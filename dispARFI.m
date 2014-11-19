function dispARFI(ecgdata,bdata,arfidata,arfidata_mf_pre,arfidata_mf_push,options,par)

[ndepth nacqT ntrackT] = size(arfidata.disp);

edge = [arfidata.axial(1) arfidata.axial(end)];

gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
if (gate(1)<arfidata.axial(1) || gate(2)>arfidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',gate(1),gate(2),arfidata.axial(1),arfidata.axial(end));
end

gate_idx = [find(arfidata.axial>gate(1),1,'first') find(arfidata.axial<gate(2),1,'last')];

% Display HQ Bmode Frames
figure
if isunix
    set(gcf,'Position',[203 286 1196 1170])
elseif ispc
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

for i=1:size(bdata.bimg,3);
    p1 = subplot('Position',[0.1 0.6 0.3 0.3]);
    imagesc(bdata.blat,bdata.bax,bdata.bimg(:,:,i));
    colormap(gray);axis image; freezeColors;
    title(i);
    hold on
    rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2)
    rectangle('Position',[-2 gate(1) 4 gate(2)-gate(1)],'EdgeColor','g','Linewidth',2)
    hold off
    xlabel('Lateral (mm)','fontsize',16,'fontweight','bold')
    ylabel('Axial (mm)','fontsize',16,'fontweight','bold')
    title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)',i,bdata.t(i)),'fontsize',16,'fontweight','bold')
    xlim([-25 25]);ylim(edge + [-15 15])
    pause(0.025)
    if i==size(bdata.bimg,3)
        % Display M-mode IQ
        p2 = axes('Position',[0.45 0.3 0.52 1]);
        frame = arfidata.IQ(:,:,1);
        imagesc(linspace(0,range(arfidata.lat(:))*nacqT,size(arfidata.IQ,2)),arfidata.IQaxial,db(frame/max(frame(:))),options.display.IQrange)
        hold on
        plot(linspace(0,range(arfidata.lat(:))*nacqT,size(arfidata.IQ,2)),arfidata.axial(1)*ones(size(arfidata.IQ,2)),'b','Linewidth',2)
        plot(linspace(0,range(arfidata.lat(:))*nacqT,size(arfidata.IQ,2)),arfidata.axial(end)*ones(size(arfidata.IQ,2)),'b','Linewidth',2)
        plot(linspace(0,range(arfidata.lat(:))*nacqT,size(arfidata.IQ,2)),gate(1)*ones(size(arfidata.IQ,2)),'g','Linewidth',2)
        plot(linspace(0,range(arfidata.lat(:))*nacqT,size(arfidata.IQ,2)),gate(2)*ones(size(arfidata.IQ,2)),'g','Linewidth',2)
        hold off
        axis image
        ylabel('Axial (mm)','fontsize',16,'fontweight','bold')
        set(gca,'xTickLabel',[])
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d',par.isHarmonic),'fontsize',16,'fontweight','bold')
        ylim(edge + [-15 15])
        colormap(gray); freezeColors;
        %         grid on
    end
end

% NaN out push reverb frames
arfidata = interpPushReverb(arfidata,options,par,'nan'); % NaN out push and reverb frames
if options.motionFilter.enable
    arfidata_mf_pre = interpPushReverb(arfidata_mf_pre,options,par,'nan'); % NaN out push and reverb frames
    arfidata_mf_push = interpPushReverb(arfidata_mf_push,options,par,'nan'); % NaN out push and reverb frames
end

% Coorelation mask filter
if options.display.cc_filt
    mask = arfidata.cc>options.display.cc_thresh;
end

% Indices corresponding to median filter parameters
nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));
cmap = colormap(hot);
cmap(65,:) = [0.5 0.5 0.5];

% Display M-mode ARFI
if options.motionFilter.enable
    
    if strcmpi(options.display.t_disp_pre,'max')
        [pk idx] = max(arfidata_mf_pre.disp(:,:,1:par.nref),[],3);
        idx_pre = round(median(idx)); % to get a single representative number through all depth
        pre = medfilt2(double(pk),[nax nt]);
        % figure out way to cc_filt these properly
        clear pk idx
    else
        idx_pre = find(arfidata.trackTime>options.display.t_disp_pre,1);
        pre = medfilt2(double(arfidata_mf_pre.disp(:,:,idx_pre)),[nax nt]);
        if options.display.cc_filt; pre(mask(:,:,idx_pre)==0) = inf; end
        idx_pre = repmat(idx_pre,[1 size(pre,2)]);
    end
    
    if strcmpi(options.display.t_disp_push,'max')
        [pk idx] = max(arfidata_mf_push.disp(:,:,par.nref:end),[],3);
        idx_push = par.nref -1 + round(median(idx)); % to get a single representative number through all depth
        push = medfilt2(double(pk),[nax nt]);
        % figure out way to cc_filt these properly
        clear pk idx
    else
        idx_push = find(arfidata.trackTime>options.display.t_disp_push,1);
        push = medfilt2(double(arfidata_mf_push.disp(:,:,idx_push)),[nax nt]);
        if options.display.cc_filt; push(mask(:,:,idx_push)==0) = inf; end
        idx_push = repmat(idx_push,[1 size(push,2)]);
    end
    
    rng = options.display.disprange;
    
    if options.display.normalize
        m1 = 1/range(pre(:)); b1 = 1 - max(pre(:))/range(pre(:));
        m2 = 1/range(push(:)); b2 = 1 - max(push(:))/range(push(:));
        pre = m1.*pre + b1; push = m2.*push + b2;
        rng = [0 1];
    end
        
    p3 = axes('Position',[0.5 0.53 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,pre,rng);
    hold on
    plot(arfidata.acqTime,gate(1)*ones(nacqT),'g','Linewidth',2)
    plot(arfidata.acqTime,gate(2)*ones(nacqT),'g','Linewidth',2)
    cp3 = copyobj(p3,gcf);
    set(cp3,'box','on','linewidth',3,'xcolor','y','ycolor','y','xticklabel',[],'yticklabel',[]);
    ylabel('Axial (mm)','Parent',p3,'fontsize',16,'fontweight','bold')
    if strcmpi(options.display.t_disp_pre,'max')
        title(sprintf('ARFI Displacements over DOF at t_m_a_x (pre push)'),'parent',p3,'fontsize',16,'fontweight','bold')
    else
        title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (pre push)',arfidata.trackTime(idx_pre(i))),'parent',p3,'fontsize',16,'fontweight','bold')
    end
    p4 = axes('Position',[0.5 0.37 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,push,rng);
    hold on
    plot(arfidata.acqTime,gate(1)*ones(nacqT),'g','Linewidth',2)
    plot(arfidata.acqTime,gate(2)*ones(nacqT),'g','Linewidth',2)
    cp4 = copyobj(p4,gcf);
    set(cp4,'box','on','linewidth',3,'xcolor','g','ycolor','g','xticklabel',[],'yticklabel',[]);
    cb = colorbar;
    set(cb,'Position',[0.91 0.37 0.0187/2 0.26])
    if options.display.normalize
        ylabel(cb,'Normalized Displacement','fontsize',16,'fontweight','bold')
    else
        ylabel(cb,'Displacement (\mum)','fontsize',16,'fontweight','bold')
    end
    ylabel('Axial (mm)','parent',p4,'fontsize',16,'fontweight','bold')
    if strcmpi(options.display.t_disp_pre,'max')
        title(sprintf('ARFI Displacements over DOF at t_m_a_x (at push)'),'parent',p4,'fontsize',16,'fontweight','bold')
    else
        title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (at push)',arfidata.trackTime(idx_push(i))),'parent',p4,'fontsize',16,'fontweight','bold')
    end
    colormap(cmap)
    set(cb,'yLim',[options.display.disprange]+[0 -1])    
else
    
    if strcmpi(options.display.t_disp_pre,'max')
        [pk idx] = max(arfidata.disp(:,:,1:par.nref),[],3);
        idx_pre = round(median(idx)); % to get a single representative number through all depth
        pre = medfilt2(double(pk),[nax nt]);
        % figure out way to cc_filt these properly
        clear pk idx
    else
        idx_pre = find(arfidata.trackTime>options.display.t_disp_pre,1);
        pre = medfilt2(double(arfidata.disp(:,:,idx_pre)),[nax nt]);
        idx_pre = repmat(idx_pre,[1 size(pre,2)]);
    end
    
    if strcmpi(options.display.t_disp_push,'max')
        [pk idx] = max(arfidata.disp(:,:,par.nref:end),[],3);
        idx_push = par.nref -1 + round(median(idx)); % to get a single representative number through all depth
        push = medfilt2(double(pk),[nax nt]);
        % figure out way to cc_filt these properly
        clear pk idx
    else
        idx_push = find(arfidata.trackTime>options.display.t_disp_push,1);
        push = medfilt2(double(arfidata.disp(:,:,idx_push)),[nax nt]);
        idx_push = repmat(idx_push,[1 size(push,2)]);
    end
    
    rng = options.display.disprange;
    
    if options.display.normalize
        m1 = 1/range(pre(:)); b1 = max(pre(:))/range(pre(:));
        m2 = 1/range(push(:)); b2 = max(push(:))/range(push(:));
        pre = m1.*pre + b1; push = m2.*push + b2;
        rng = [0 1];
    end
    
    if options.display.cc_filt
        pre(mask(:,:,idx_pre)==0) = nan;
        push(mask(:,:,idx_push)==0) = nan;
    end
    
    p3 = axes('Position',[0.5 0.53 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,pre,rng);
    hold on
    plot(arfidata.acqTime,gate(1)*ones(nacqT),'g','Linewidth',2)
    plot(arfidata.acqTime,gate(2)*ones(nacqT),'g','Linewidth',2)
    cp3 = copyobj(p3,gcf);
    set(cp3,'box','on','linewidth',3,'xcolor','y','ycolor','y','xticklabel',[],'yticklabel',[]);
    ylabel('Axial (mm)','Parent',im1,'fontsize',16,'fontweight','bold')
    if strcmpi(options.display.t_disp_pre,'max')
        title(sprintf('ARFI Displacements over DOF at t_m_a_x (pre push)',arfidata.trackTime(idx_pre(i))),'parent',p3,'fontsize',16,'fontweight','bold')
    else
        title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (pre push)',arfidata.trackTime(idx_pre(i))),'parent',p3,'fontsize',16,'fontweight','bold')
    end
    p4 = axes('Position',[0.5 0.37 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,push,rng);
    hold on
    plot(arfidata.acqTime,gate(1)*ones(nacqT),'g','Linewidth',2)
    plot(arfidata.acqTime,gate(2)*ones(nacqT),'g','Linewidth',2)
    cp4 = copyobj(p4,gcf);
    set(cp4,'box','on','linewidth',3,'xcolor','g','ycolor','g','xticklabel',[],'yticklabel',[]);
    cb = colorbar;
    set(cb,'Position',[0.91 0.37 0.0187/2 0.26])
    if options.display.normalize
        ylabel(cb,'Normalized Displacement','fontsize',16,'fontweight','bold')
    else
        ylabel(cb,'Displacement (\mum)','fontsize',16,'fontweight','bold')
    end
    ylabel('Axial (mm)','parent',im2,'fontsize',16,'fontweight','bold')
    if strcmpi(options.display.t_disp_pre,'max')
        title(sprintf('ARFI Displacements over DOF at t_m_a_x (at push)',arfidata.trackTime(idx_push(i))),'parent',p4,'fontsize',16,'fontweight','bold')
    else
        title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (at push)',arfidata.trackTime(idx_push(i))),'parent',p4,'fontsize',16,'fontweight','bold')
    end
    colormap(cmap)
    set(cb,'yLim',[options.display.disprange]+[0 -1])
end

if isempty(ecgdata)
    xlabel('Acquisition Time (s)','fontsize',16,'fontweight','bold')
end

% NaN out displacements filtered out by cc_thresh
pre(pre==inf) = nan;
push(push==inf) = nan;
arfidata.disp(mask==0) = nan;
if options.motionFilter.enable
    arfidata_mf_pre.disp(mask==0) = nan;
    arfidata_mf_push.disp(mask==0) = nan;
end

% Compute axially averaged Pre and Push Traces
pre_trace = nanmean(pre(gate_idx(1):gate_idx(2),:));
push_trace = nanmean(push(gate_idx(1):gate_idx(2),:));

% Incorporate ECG Data into this figure
if ~isempty(ecgdata)
    samples = zeros(1,nacqT);
    for i=1:nacqT
        samples(i) = ecgdata.arfi(find(ecgdata.arfi(:,1)>arfidata.acqTime(i),1,'first'),2);
    end
    ecgdata.arfi(:,2) = ecgdata.arfi(:,2)/max(ecgdata.arfi(:,2));
    
    h1 = axes('Position',[0.5 0.1 0.4 0.2]);
    plot(ecgdata.arfi(:,1),ecgdata.arfi(:,2),'Linewidth',2);
    hold on
    plot(arfidata.acqTime,samples,'kx','MarkerSize',8)
    pt = plot(arfidata.acqTime(1),samples(1),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    hold off
    grid on
    title('ECG Trace','fontsize',16,'fontweight','bold')
    xlabel('Acquisition Time (s)','fontsize',16,'fontweight','bold')
    axis tight
    hold(h1)
end

h2 = axes('Position',[0.05 0.15 0.4 0.3]);
set(h2,'Color',[0.5 0.5 0.5]);

if (options.motionFilter.enable && (strcmpi(options.motionFilter.method,'Polynomial') || strcmpi(options.motionFilter.method,'Both')))
    rng = options.display.disprange*2;
else
    rng = [-100 100];
end

% filename = 'test.gif';
for i=1:nacqT
    cla(h2)
    if ~isempty(ecgdata)
        set(pt,'Visible','off')
        set(h1,'Color',[0.5 0.5 0.5]);
    end
    if options.motionFilter.enable
        plot(arfidata.trackTime(1:par.nref),squeeze(arfidata_mf_pre.disp(ceil(linspace(gate_idx(1),gate_idx(2),5)),i,1:par.nref)),'.--','Parent',h2)
        hold on
        plot(arfidata.trackTime(par.nref+1:end),squeeze(arfidata_mf_push.disp(ceil(linspace(gate_idx(1),gate_idx(2),5)),i,par.nref+1:end)),'.--','Parent',h2)
        set(h2,'Color',[0.5 0.5 0.5]);
        ylim(h2,rng)
    else
        plot(arfidata.trackTime,squeeze(arfidata.disp(ceil(linspace(gate_idx(1),gate_idx(2),5)),i,:)),'.--','Parent',h2)
        set(h2,'Color',[0.5 0.5 0.5]);
        ylim(h2,rng)
    end
    hold on
    plot(arfidata.trackTime(idx_pre(i))*ones(1,10),linspace(-300,300,10),'y','linewidth',2,'Parent',h2)
    plot(arfidata.trackTime(idx_push(i))*ones(1,10),linspace(-300,300,10),'g','linewidth',2,'Parent',h2)
    title(sprintf('ARFI Displacement Profiles:\nDepth Gate = %2.2f - %2.2f mm\nPush # %d (t = %2.2f s)\nMotion Filter = %s',gate(1),gate(2),i,arfidata.acqTime(i),options.motionFilter.method*options.motionFilter.enable),'fontsize',16,'fontweight','bold','Parent',h2)
    if i==1
        xlabel('Track Time (ms)','fontsize',16,'fontweight','bold')
        ylabel('Displacement (\mum)','fontsize',16,'fontweight','bold')
        xlim([arfidata.trackTime(1) arfidata.trackTime(end)])
        grid on
    end
    if ~isempty(ecgdata)
        pt = plot(arfidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
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
    if options.display.IQtraces
    temp = db(squeeze(arfidata.IQ(:,1+(i-1)*par.nBeams,:)));
    temp(:,par.nref+1:par.nref+par.npush+par.nreverb) = nan;
    offset = 5;
    figure(101);
    if isunix
        set(101,'Position',[1203 390 1916 767])
    elseif ispc
        set(101,'units','normalized','outerposition',[0 0 1 1])
    end
    hh=subplot(121);cla(hh);
    for j=1:ntrackT;plot(arfidata.IQaxial,offset*(j-1)-temp(:,j)');hold on;end;view(90,90);
    hold on;plot(gate(1)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);plot(gate(2)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3)
    xlim(edge);title(sprintf('Raw IQ: %d (t = %2.2f s)',i,arfidata.acqTime(i)),'fontsize',16,'fontweight','bold');
    set(hh,'Color',[0.5 0.5 0.5]);
    xlabel('Axial (mm)','fontsize',16,'fontweight','bold');ylabel('Tracks','fontsize',16,'fontweight','bold');set(gca,'YTickLabel',[])
    
    subplot(122);imagesc(arfidata.trackTime,arfidata.axial,squeeze(arfidata.cc(:,i,:)),[options.display.cc_thresh 1]);colorbar;
    title(sprintf('%s Correlation Coefficients',options.dispEst.ref_type),'fontsize',16,'fontweight','bold');grid on;colormap(jet)
    hold on;plot(linspace(-8,8,100),gate(1)*ones(1,100),'g','linewidth',3);plot(linspace(-8,8,100),gate(2)*ones(1,100),'g','linewidth',3)
    xlabel('Track Time (ms)','fontsize',16,'fontweight','bold');ylabel('Axial (mm)','fontsize',16,'fontweight','bold')
    end
    pause
end

cla(h2);
plot(arfidata.acqTime,pre_trace,'y.--','Parent',h2);hold all;plot(arfidata.acqTime,push_trace,'gx--','Parent',h2);
set(h2,'Color',[0.5 0.5 0.5]);
xlabel('Acquisition Time (s)','fontsize',16,'fontweight','bold','Parent',h2); 
title(sprintf('Axially Averaged ARFI Displacements:\nDepth Gate = %2.2f - %2.2f mm\n',gate(1),gate(2)),'fontsize',16,'fontweight','bold','Parent',h2)
axis(h2,'tight')
