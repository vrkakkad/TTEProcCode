function dispARFI(ecgdata,bdata,arfidata,arfidata_mf_pre,arfidata_mf_push,options,par)

% Not currently compatible with ref_type = 'independent'
% figure;plot(ecgdata.bmode(:,2))

dof = 7.22*1.540/par.pushFreq*(par.pushFnum)^2;
edge = (par.pushFocalDepth + [-dof/2 dof/2]);
idx_pre = find(arfidata.trackTime>options.display.t_disp_pre,1);
idx_push = find(arfidata.trackTime>options.display.t_disp_push,1);

gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
if (gate(1)<arfidata.axial(1) || gate(2)>arfidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',gate(1),gate(2),arfidata.axial(1),arfidata.axial(end));
end

gate_idx = [find(arfidata.axial>gate(1),1,'first') find(arfidata.axial<gate(2),1,'last')];

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
    rectangle('Position',[-2 gate(1) 4 gate(2)-gate(1)],'EdgeColor','r','Linewidth',2)
    hold off
    xlabel('Lateral (mm)','fontsize',10','fontweight','bold')
    ylabel('Axial (mm)','fontsize',10','fontweight','bold')
    title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)',i,bdata.t(i)),'fontsize',10','fontweight','bold')
    xlim([-10 10]);ylim(edge + [-1 1])
    pause(0.025)
    if i==1
        % Display M-mode IQ
        subplot('Position',[0.5 0.7 0.4 0.2])
        imagesc(arfidata.acqTime,arfidata.IQaxial,abs(db(arfidata.IQ(:,:,1))),options.display.IQrange)
        hold on
        plot(arfidata.acqTime,gate(1)*ones(length(arfidata.acqTime)),'r','Linewidth',2)
        plot(arfidata.acqTime,gate(2)*ones(length(arfidata.acqTime)),'r','Linewidth',2)
        plot(arfidata.acqTime,arfidata.axial(1)*ones(length(arfidata.acqTime)),'b','Linewidth',2)
        plot(arfidata.acqTime,arfidata.axial(end)*ones(length(arfidata.acqTime)),'b','Linewidth',2)
        hold off
        ylabel('Axial (mm)','fontsize',10','fontweight','bold')
        title(sprintf('M-Mode Frames\n Harmonic Tracking = %d',par.isHarmonic),'fontsize',10','fontweight','bold')
        ylim(edge + [-1 1])
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
        mask = arfidata.ccout>options.display.cc_thresh;
%         arfidata.disp = arfidata.disp.*mask;
        if options.motionFilter.enable
%             arfidata_mf_pre.disp = arfidata_mf_pre.disp.*mask;
%             arfidata_mf_push.disp = arfidata_mf_push.disp.*mask;
        end
end

% Indices corresponding to median filter parameters
nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));

% Display M-mode ARFI
if options.motionFilter.enable
    pre = medfilt2(double(arfidata_mf_pre.disp(:,:,idx_pre)),[nax nt]);
    push = medfilt2(double(arfidata_mf_push.disp(:,:,idx_push)),[nax nt]);
    %     pre = medfilt2(double(arfidata_mf_pre.disp(gate_idx(1):gate_idx(2),:,idx_pre)),[nax nt]);
    %     push = medfilt2(double(arfidata_mf_push.disp(gate_idx(1):gate_idx(2),:,idx_push)),[nax nt]);
    [r1,c1] = find(mask(:,:,idx_pre)==0);%find(pre==0); 
    [r2,c2] = find(mask(:,:,idx_push)==0);%find(push==0);
    rng = options.display.disprange;
    if options.display.normalize
        m1 = 1/range(pre(:)); b1 = 1 - max(pre(:))/range(pre(:));
        m2 = 1/range(push(:)); b2 = 1 - max(push(:))/range(push(:));
        pre = m1.*pre + b1; push = m2.*push + b2;
        rng = [0 1];
    end
    im1 = subplot('Position',[0.5 0.53 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,pre,rng);
    hold on
    plot(arfidata.acqTime(c1), arfidata.axial(r1), 'bx') ;
    %     imagesc(arfidata.acqTime,arfidata.axial(gate_idx(1):gate_idx(2)),pre,rng);
    cim1 = copyobj(im1,gcf);
    set(cim1,'box','on','linewidth',3,'xcolor','y','ycolor','y','xticklabel',[],'yticklabel',[]);
    %     set(cim1,'xgrid','on','ygrid','on')
    ylabel('Axial (mm)','Parent',im1,'fontsize',10','fontweight','bold')
    title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (pre push)',arfidata.trackTime(idx_pre)),'parent',im1,'fontsize',10','fontweight','bold')
    im2 = subplot('Position',[0.5 0.37 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,push,rng);
    hold on
    plot(arfidata.acqTime(c2), arfidata.axial(r2), 'bx') ;
%     imagesc(arfidata.acqTime,arfidata.axial(gate_idx(1):gate_idx(2)),push,rng);
    cim2 = copyobj(im2,gcf);
    set(cim2,'box','on','linewidth',3,'xcolor','g','ycolor','g','xticklabel',[],'yticklabel',[]);
    %     set(cim1,'xgrid','on','ygrid','on')
    hcb = colorbar;
    set(hcb,'Position',[0.91 0.37 0.0187/2 0.26])
    if options.display.normalize
        ylabel(hcb,'Normalized Displacement','fontsize',10','fontweight','bold')
    else
        ylabel(hcb,'Displacement (\mum)','fontsize',10','fontweight','bold')
    end
    ylabel('Axial (mm)','parent',im2,'fontsize',10','fontweight','bold')
    title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (at push)',arfidata.trackTime(idx_push)),'parent',im2,'fontsize',10','fontweight','bold')
    colormap(hot)
else
    pre = medfilt2(double(arfidata.disp(:,:,idx_pre)),[nax nt]);
    push = medfilt2(double(arfidata.disp(:,:,idx_push)),[nax nt]);
%     pre = medfilt2(double(arfidata.disp(gate_idx(1):gate_idx(2),:,idx_pre)),[nax nt]);
%     push = medfilt2(double(arfidata.disp(gate_idx(1):gate_idx(2),:,idx_push)),[nax nt]);
    [r1,c1] = find(mask(:,:,idx_pre)==0);%find(pre==0); 
    [r2,c2] = find(mask(:,:,idx_push)==0);%find(push==0);
    rng = options.display.disprange;
    if options.display.normalize
        m1 = 1/range(pre(:)); b1 = max(pre(:))/range(pre(:));
        m2 = 1/range(push(:)); b2 = max(push(:))/range(push(:));
        pre = m1.*pre + b1; push = m2.*push + b2;
        rng = [0 1];
    end
    im1 = subplot('Position',[0.5 0.53 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,pre,rng);
    hold on
    plot(arfidata.acqTime(c1), arfidata.axial(r1), 'bx') ;
%     imagesc(arfidata.acqTime,arfidata.axial(gate_idx(1):gate_idx(2)),pre,rng);
    cim1 = copyobj(im1,gcf);
    set(cim1,'box','on','linewidth',3,'xcolor','y','ycolor','y','xticklabel',[],'yticklabel',[]);
    %     set(cim1,'xgrid','on','ygrid','on')
    ylabel('Axial (mm)','Parent',im1,'fontsize',10','fontweight','bold')
    title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (pre push)',arfidata.trackTime(idx_pre)),'Parent',im1,'fontsize',10','fontweight','bold')
    im2 = subplot('Position',[0.5 0.37 0.4 0.1]);
    imagesc(arfidata.acqTime,arfidata.axial,push,rng);
    hold on
    plot(arfidata.acqTime(c2), arfidata.axial(r2),'bx') ;
%     imagesc(arfidata.acqTime,arfidata.axial(gate_idx(1):gate_idx(2)),push,rng);
    cim2 = copyobj(im2,gcf);
    set(cim2,'box','on','linewidth',3,'xcolor','g','ycolor','g','xticklabel',[],'yticklabel',[]);
    %     set(cim1,'xgrid','on','ygrid','on')
    hcb = colorbar;
    set(hcb,'Position',[0.91 0.37 0.0187/2 0.26])
    if options.display.normalize
        ylabel(hcb,'Normalized Displacement','fontsize',10','fontweight','bold')
    else
        ylabel(hcb,'Displacement (\mum)','fontsize',10','fontweight','bold')
    end
    ylabel('Axial (mm)','parent',im2,'fontsize',10','fontweight','bold')
    title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (post push)',arfidata.trackTime(idx_push)),'parent',im2,'fontsize',10','fontweight','bold')
    colormap(hot)
end
if isempty(ecgdata)
    xlabel('Acquisition Time (s)','fontsize',10','fontweight','bold')
end

% NaN out zero displacements
arfidata.disp(mask==0) = nan;
if options.motionFilter.enable
    arfidata_mf_pre.disp(mask==0) = nan;
    arfidata_mf_push.disp(mask==0) = nan;
end

% % Display averaged M-mode ARFI
% subplot('Position',[0.05 0.37 0.4 0.17])
% if options.motionFilter.enable
%     temp1 = medfilt2(double(arfidata_mf_push.disp(:,:,idx_pre)));
%     temp2 = medfilt2(double(arfidata_mf_push.disp(:,:,idx_push)));
% else
%     temp1 = medfilt2(double(arfidata.disp(:,:,idx_pre)));
%     temp2 = medfilt2(double(arfidata.disp(:,:,idx_push)));
% end
%
% avg_pre = mean(temp1(edge_idx(1):edge_idx(2),:),1);
% avg_push = mean(temp2(edge_idx(1):edge_idx(2),:),1);
% plot(arfidata.acqTime,avg_pre,'b*--','linewidth',2);
% hold all
% plot(arfidata.acqTime,avg_push,'r*--','linewidth',2);
% ylim(options.display.disprange*2)
% grid on
% xlabel('Acquisition Time (s)')
% ylabel('\mum')
% legend('Pre','Push')

% Incorporate ECG Data into this plot
if ~isempty(ecgdata)
    samples = zeros(1,size(arfidata.disp,2));
    for i=1:size(arfidata.disp,2)
        samples(i) = ecgdata.arfi(find(ecgdata.arfi(:,1)>arfidata.acqTime(i),1,'first'),2);
    end
    ecgdata.arfi(:,2) = ecgdata.arfi(:,2)/max(ecgdata.arfi(:,2));
    
    h1 = subplot('Position',[0.5 0.1 0.4 0.2]);
    plot(ecgdata.arfi(:,1),ecgdata.arfi(:,2),'Linewidth',2);
    hold on
    plot(arfidata.acqTime,samples,'kx','MarkerSize',8)
    pt = plot(arfidata.acqTime(1),samples(1),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    hold off
    grid on
    title('ECG Trace','fontsize',10','fontweight','bold')
    xlabel('Acquisition Time (s)','fontsize',10','fontweight','bold')
    axis tight
    hold(h1)
end

h2 = subplot('Position',[0.05 0.15 0.4 0.3]);

if (options.motionFilter.enable && (strcmpi(options.motionFilter.method,'Polynomial') || strcmpi(options.motionFilter.method,'Both')))
    rng = options.display.disprange*2;
else
    rng = [-100 100];
%     rng = options.display.disprange*2;
end

% filename = 'test.gif';
for i=1:size(arfidata.disp,2)
    cla(h2)
    if ~isempty(ecgdata)
        set(pt,'Visible','off')
    end
    if options.motionFilter.enable
        plot(arfidata.trackTime(1:par.nref),squeeze(arfidata_mf_pre.disp(ceil(linspace(gate_idx(1),gate_idx(2),5)),i,1:par.nref)),'.--','Parent',h2)
        hold on
        plot(arfidata.trackTime(par.nref+1:end),squeeze(arfidata_mf_push.disp(ceil(linspace(gate_idx(1),gate_idx(2),5)),i,par.nref+1:end)),'.--','Parent',h2)
        ylim(h2,rng)
    else
        plot(arfidata.trackTime,squeeze(arfidata.disp(ceil(linspace(gate_idx(1),gate_idx(2),5)),i,:)),'.--','Parent',h2)
        ylim(h2,rng)
    end
    hold on
    plot(options.display.t_disp_pre*ones(1,10),linspace(-100,100,10),'y','linewidth',2,'Parent',h2)
    plot(options.display.t_disp_push*ones(1,10),linspace(-100,100,10),'g','linewidth',2,'Parent',h2)
    title(sprintf('ARFI Displacement Profiles:\nDepth Gate = %2.2f - %2.2f mm\nPush # %d (t = %2.2f s)\nMotion Filter = %s',gate(1),gate(2),i,arfidata.acqTime(i),options.motionFilter.method*options.motionFilter.enable),'fontsize',10','fontweight','bold','Parent',h2)
    if i==1
        xlabel('Track Time (ms)','fontsize',10','fontweight','bold')
        ylabel('Displacement (\mum)','fontsize',10','fontweight','bold')
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
    if options.display.debug
    temp = db(abs(squeeze(arfidata.IQ(:,i,:))));
    temp(:,par.nref+1:par.nref+par.npush+par.nreverb) = nan;
    offset = 5;
    figure(101);
    if isunix
        set(101,'Position',[1201 386 1920 1070])
    elseif ispc
        set(101,'units','normalized','outerposition',[0 0 1 1])
    end
    hh=subplot(121);cla(hh);for j=1:size(arfidata.disp,3);plot(arfidata.IQaxial,offset*(j-1)-temp(:,j)');hold on;end;view(90,90);
    hold on;plot(gate(1)*ones(1,100),linspace(-offset*10,offset*size(arfidata.disp,3),100),'r','linewidth',3);plot(gate(2)*ones(1,100),linspace(-offset*10,offset*size(arfidata.disp,3),100),'r','linewidth',3)
    xlim(edge);title(sprintf('Raw IQ: %d',i));xlabel('Axial (mm)');ylabel('Tracks');set(gca,'YTickLabel',[])
    subplot(122);imagesc(arfidata.trackTime,arfidata.axial,squeeze(arfidata.ccout(:,i,:)),[options.display.cc_thresh 1]);colorbar;title(sprintf('%s Correlation Coefficients',options.dispEst.ref_type));grid on;colormap(jet)
    hold on;plot(linspace(-8,8,100),gate(1)*ones(1,100),'r','linewidth',3);plot(linspace(-8,8,100),gate(2)*ones(1,100),'r','linewidth',3)
    xlabel('Track Time (ms)');ylabel('Axial (mm)')
    end
    pause
end

% test1 = arfidata.disp(ceil(linspace(edge_idx(1),edge_idx(2),5)),:,:);
% test2 = max(test1,[],3);
% figure;stem(mean(test2,1));grid on;ylim([0 30])
% title([mean(test2(:)),median(test2(:))])
% figure
% for i=1:40
% subplot(211);plot(arfidata.trackTime,squeeze(arfidata.disp(ceil(linspace(edge_idx(1),edge_idx(2),5)),i,:))','*--');grid on;,ylim([-100 100]);title(sprintf('Raw: %d',i))
% subplot(212);plot(arfidata.trackTime,squeeze(arfidata_mf.disp(ceil(linspace(edge_idx(1),edge_idx(2),5)),i,:))','*--');grid on;ylim([-15 30]);title(sprintf('MF: %d',i));
% pause
% end

