function dispSWEI(ecgdata,bdata,sweidata,sweidata_mf,options,par)

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
set(gcf,'units','normalized','outerposition',[0 0 1 1])
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
    title(sprintf('HQ B-Mode: Frame %d (t = %1.1f s)\n',i,bdata.t(i)))
    pause(0.025)
    
    if i==1
        % Display M-mode Data
        subplot('Position',[0.5 0.7 0.4 0.2])
        imagesc(sweidata.acqTime,sweidata.IQaxial,abs(db(sweidata.IQ(:,ceil(par.nBeams/2):par.nBeams:end,1))))
        hold on
        plot(sweidata.acqTime,edge2(1)*ones(length(sweidata.acqTime)),'r','Linewidth',2)
        plot(sweidata.acqTime,edge2(2)*ones(length(sweidata.acqTime)),'r','Linewidth',2)
        plot(sweidata.acqTime,sweidata.axial(1)*ones(length(sweidata.acqTime)),'b--','Linewidth',2)
        plot(sweidata.acqTime,sweidata.axial(end)*ones(length(sweidata.acqTime)),'b--','Linewidth',2)
        hold off
%         xlabel('Acquisition Time (s)')
        ylabel('Axial (mm)')
        title('M-Mode Frames')
        colormap(gray)
        grid on
    end
end

vel = diff(sweidata.disp,1,4);
gate_avg = squeeze(mean(sweidata.disp(edge_idx(1):edge_idx(2),:,:,:),1));
gate_avg = permute(gate_avg,[1 3 2]);
vel_avg = 1000*diff(gate_avg,1,2)./par.priusec(1);
if options.motionFilter.enable
    vel_mf = diff(sweidata_mf.disp,1,4);
    gate_avg_mf = squeeze(mean(sweidata_mf.disp(edge_idx(1):edge_idx(2),:,:,:),1));
    gate_avg_mf = permute(gate_avg_mf,[1 3 2]);
    vel_avg_mf = 1000*diff(gate_avg_mf,1,2)./par.priusec(1);
end
lat_idx = round(mean(edge_idx));
lat = sweidata.lat(lat_idx,:);

% Interpolate Shear Wave Plot in lateral dimension by 5
lat_int = linspace(lat(1),lat(end),length(lat)*5);
gate_avg_int = interp1(lat,gate_avg,lat_int);
vel_avg_int = interp1(lat,vel_avg,lat_int);
if options.motionFilter.enable
    gate_avg_mf_int = interp1(lat,gate_avg_mf,lat_int);
    vel_avg_mf_int = interp1(lat,vel_avg_mf,lat_int);
end

if options.display.calcSWS
    % p1 = vel_avg_int(find(lat_int==0):end,find(sweidata.trackTime==0):end,:);
    % p2 = vel_avg_int(1:find(lat_int==0),find(sweidata.trackTime==0):end,:);
    p1 = vel_avg(find(lat==0):end,find(sweidata.trackTime==0):end,:);
    p2 = vel_avg(1:find(lat==0),find(sweidata.trackTime==0):end,:);
    
    t = sweidata.trackTime(find(sweidata.trackTime==0):end-1);
    t_int = linspace(t(1),t(end),length(t)*4); % Interpolate Track Time by a factor of 4 to 25 kHz
    
    % l = lat_int(find(lat_int==0):end);
    l = lat(find(lat==0):end);
    
    p1 = permute(p1,[2 1 3]); p2 = permute(p2,[2 1 3]);
    p1 = interp1(t,p1,t_int); p2 = interp1(t,p2,t_int);
    p1 = permute(p1,[2 1 3]); p2 = permute(p2,[2 1 3]);
    
    sws1 = zeros(1,size(p1,3));
    sws2 = zeros(1,size(p1,3));
    for i=1:size(p1,3)
        tic
        temp = CalcSWSfromLatsums(p1(:,:,i),l,t_int,10,1,8,0);
        sws1(i) = temp.speed;
        temp = CalcSWSfromLatsums(flipud(p2(:,:,i)),l,t_int,10,1,8,0);
        sws2(i) = temp.speed;
        toc
    end
    
    sws = mean([sws1;sws2]);
    
    h0 = subplot('Position',[0.5 0.4 0.4 0.2]);
    plot(sweidata.acqTime,sws,'*--')
    hold on
    pt1 = plot(sweidata.acqTime(1),sws(1),'ro','Parent',h0,'Markersize',10);
    hold off
    grid on
    ylabel('SWS (m/s)')
    ylim([0 10])
    title('Estimated Shear Wave Speed')
    hold(h0)
end

if ~isempty(ecgdata)
    samples = zeros(1,size(gate_avg,3));
    for i=1:size(gate_avg,3)
        samples(i) = ecgdata.swei(find(ecgdata.swei(:,1)>sweidata.acqTime(i),1,'first'),2);
    end
    ecgdata.swei(:,2) = ecgdata.swei(:,2)/max(ecgdata.swei(:,2));
    
    if options.display.calcSWS
        h1 = subplot('Position',[0.5 0.1 0.4 0.2]);
    else
        h1 = subplot('Position',[0.5 0.25 0.4 0.2]);
    end
    plot(ecgdata.swei(:,1),ecgdata.swei(:,2),'Linewidth',2);
    hold on
    plot(sweidata.acqTime,samples,'kx','MarkerSize',8)
    pt2 = plot(sweidata.acqTime(1),samples(1),'ro','Parent',h1,'Markersize',10);
    hold off
    grid on
    title('ECG Trace')
    xlabel('Acquisition Time (s)')
    axis tight
    hold(h1)
end

h2 = subplot('Position',[0.1 0.12 0.35 0.35]);
% Generate displacement or Velocity Plots
if strcmpi(options.display.sw_display,'disp')
    for i=1:size(gate_avg,3)
        if options.display.calcSWS
            set(pt1,'Visible','off')
        end
        if ~isempty(ecgdata)
            set(pt2,'Visible','off')
        end
        
        if options.motionFilter.enable
            imagesc(lat_int,sweidata.trackTime,gate_avg_mf_int(:,:,i)',options.display.disprange/2)
        else
            imagesc(lat_int,sweidata.trackTime,gate_avg_int(:,:,i)',options.display.disprange/2)
        end
        xlabel('Lateral (mm)')
        ylabel('Track Time (ms)')
        title(sprintf('SWEI Disp Profile: Push # %d  (t = %2.2f s)\nDepth Gate = %2.2f - %2.2f mm\nHarmonic Tracking = %d',i,sweidata.acqTime(i),edge2(1),edge2(2),par.isHarmonic))
        colorbar
        %         for j=1:size(sweidata.disp,4)
        %             subplot(224)
        %             if options.motionFilter.enable
        %                 imagesc(sweidata.lat(lat_idx,:),sweidata.axial,sweidata_mf.disp(:,:,i,j),options.display.disprange/2)
        %             else
        %                 imagesc(sweidata.lat(lat_idx,:),sweidata.axial,sweidata.disp(:,:,i,j),options.display.disprange/2)
        %             end
        %             xlabel('Lateral (mm)')
        %             ylabel('Axial (mm)')
        %             title(sprintf('SWEI Disp: Push # %d\nHarmonic Tracking = %d',i,par.isHarmonic));
        %             axis image
        %             pause(0.001)
        %         end
        if options.display.calcSWS
            pt1 = plot(sweidata.acqTime(i),sws(i),'ro','Parent',h0,'Markersize',10);
        end
        if ~isempty(ecgdata)
            pt2 = plot(sweidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10);
        end
        pause
    end
elseif strcmpi(options.display.sw_display,'vel')
    for i=1:size(gate_avg,3)
        if options.display.calcSWS
            set(pt1,'Visible','off')
        end
        if ~isempty(ecgdata)
            set(pt2,'Visible','off')
        end
        if options.motionFilter.enable
            imagesc(lat_int,sweidata.trackTime,vel_avg_mf_int(:,:,i)',options.display.velrange)
        else
            imagesc(lat_int,sweidata.trackTime,vel_avg_int(:,:,i)',options.display.velrange)
        end
        xlabel('Lateral (mm)')
        ylabel('Track Time (ms)')
        title(sprintf('SWEI Vel Profile: Push # %d  (t = %2.2f s)\nDepth Gate = %2.2f - %2.2f mm\nHarmonic Tracking = %d',i,sweidata.acqTime(i),edge2(1),edge2(2),par.isHarmonic))
        colorbar
        %         for j=1:size(vel,4)
        %             subplot(224)
        %             if options.motionFilter.enable
        %                 imagesc(sweidata.lat(lat_idx,:),sweidata.axial,vel_mf(:,:,i,j),options.display.velrange/2)
        %             else
        %                 imagesc(sweidata.lat(lat_idx,:),sweidata.axial,vel(:,:,i,j),options.display.velrange/2)
        %             end
        %             xlabel('Lateral (mm)')
        %             ylabel('Axial (mm)')
        %             title(sprintf('SWEI Vel: Push # %d\nHarmonic Tracking = %d',i,par.isHarmonic));
        %             axis image
        %             pause(0.001)
        %         end
        if options.display.calcSWS
            pt1 = plot(sweidata.acqTime(i),sws(i),'ro','Parent',h0,'Markersize',10);
        end
        if ~isempty(ecgdata)
            pt2 = plot(sweidata.acqTime(i),samples(i),'ro','Parent',h1,'Markersize',10);
        end
        pause
    end
end

