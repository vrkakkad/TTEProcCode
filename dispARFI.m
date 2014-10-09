function dispARFI(ecgdata,bdata,arfidata,arfidata_mf,options,par)

% Not currently compatible with ref_type = 'independent'

dof = 7.22*1.540/par.pushFreq*(par.pushFnum)^2;
edge = (par.pushFocalDepth + [-dof/2 dof/2]);

edge2 = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
if (edge2(1)<arfidata.axial(1) || edge2(2)>arfidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',edge2(1),edge2(2),arfidata.axial(1),arfidata.axial(end));
end
edge_idx = [find(arfidata.axial>edge2(1),1,'first') find(arfidata.axial<edge2(2),1,'last')];

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
    colormap(gray);axis image;
    title(i);
    hold on
    rectangle('Position',[-10 edge(1) 20 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2)
    hold off
    xlabel('Lateral (mm)')
    ylabel('Axial (mm)')
    title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)',i,bdata.t(i)))
    pause(0.025)
    
    if i==1
        % Display M-mode Data
        subplot('Position',[0.5 0.7 0.4 0.2])
        imagesc(arfidata.acqTime,arfidata.IQaxial,abs(db(arfidata.IQ(:,:,1))))
        hold on
        plot(arfidata.acqTime,edge2(1)*ones(length(arfidata.acqTime)),'r','Linewidth',2)
        plot(arfidata.acqTime,edge2(2)*ones(length(arfidata.acqTime)),'r','Linewidth',2)
        plot(arfidata.acqTime,arfidata.axial(1)*ones(length(arfidata.acqTime)),'b--','Linewidth',2)
        plot(arfidata.acqTime,arfidata.axial(end)*ones(length(arfidata.acqTime)),'b--','Linewidth',2)
        hold off
%         xlabel('Acquisition Time (s)')
        ylabel('Axial (mm)')
        title('M-Mode Frames')
        colormap(gray)
        grid on
    end
end

gate_avg = squeeze(mean(arfidata.disp(edge_idx(1):edge_idx(2),:,:),1))';
if options.motionFilter.enable
    gate_avg_mf = squeeze(mean(arfidata_mf.disp(edge_idx(1):edge_idx(2),:,:),1))';
end

subplot('Position',[0.5 0.4 0.4 0.2])
if options.motionFilter.enable
    imagesc(arfidata.acqTime,arfidata.trackTime,gate_avg_mf,options.display.disprange)
else
    imagesc(arfidata.acqTime,arfidata.trackTime,gate_avg,options.display.disprange)
end
if isempty(ecgdata)
xlabel('Acquisition Time (s)')
end
ylabel('Track Time (ms)')
title(sprintf('ARFI Displacement Profiles:\nDepth Gate = %2.2f - %2.2f mm\nHarmonic Tracking = %d',edge2(1),edge2(2),par.isHarmonic))
grid on
% hcb = colorbar;
% set(hcb,'location','EastOutside')
% xlabel(hcb,'Displacement (\mum)')
% Incorporate ECG Data into this plot

if ~isempty(ecgdata)
    samples = zeros(1,size(gate_avg,2));
    for i=1:size(gate_avg,2)
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
    title('ECG Trace')
    xlabel('Acquisition Time (s)')
    axis tight
    hold(h1)
end

h2 = subplot('Position',[0.05 0.15 0.35 0.3]);
for i=1:size(gate_avg,2)
    cla
    if ~isempty(ecgdata)
        set(pt,'Visible','off')
    end
    plot(arfidata.trackTime,gate_avg(:,i),'.--','Parent',h2)
    hold on
    if options.motionFilter.enable
        plot(arfidata.trackTime,gate_avg_mf(:,i),'r*-','linewidth',2,'Parent',h2)
    end
    title(sprintf('Profile %d (t = %2.2f s)',i,arfidata.acqTime(i)))
    if i==1
        xlabel('Track Time (ms)')
        ylabel('Displacement (\mum)')
        if options.motionFilter.enable
        legend('Raw','Motion Filtered','location','southeast')
        else
        legend('Raw','location','southeast')
        end
        ylim([-45 45])
        xlim([arfidata.trackTime(1) arfidata.trackTime(end)])
        grid on
    end
    if ~isempty(ecgdata)
        pt = plot(arfidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10,'Markerfacecolor','r');
    end
    pause(0.1)
end

if options.motionFilter.enable
    if strcmpi(options.display.t_disp,'max')
        [disp idx] = max(gate_avg_mf);
    else
        idx = find(arfidata.trackTime>options.display.t_disp,1)*ones(1,size(gate_avg,2));
        disp = gate_avg_mf(idx(1),:);
    end
end

% temp = colormap(jet);
% temp = interp1([1:64],temp,linspace(1,64,length(arfidata.trackTime)));

% pause
% subplot('Position',[0.5 0.4 0.4 0.2])
% plot(arfidata.acqTime,disp,'*-','Linewidth',2)
% grid on
% ylim([0 50])
% ylabel('Displacement (\mum)')
% if strcmpi(options.display.t_disp,'max')
%     title('Max Displacement vs. Acquisition Time')
% else
%     title(sprintf('Displacement at t=%1.2f ms vs. Acquisition Time',options.display.t_disp))
% end


% for i=1:size(gate_avg,2)
%     stem(arfidata.acqTime(i),disp(i),'o','markersize',10,'markerfacecolor',temp(idx(i),:))
%     hold on
% end

% figure;errorbar(mean(gate_avg'), std(gate_avg'));grid on

