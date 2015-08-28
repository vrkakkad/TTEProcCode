function dispARFI(bdata,arfidata,options,par)
%% Check for existance of traced gate
[~,temp] = fileparts(pwd);
if isunix
    bordersPath = strcat('/emfd/vrk4/Transthoracic_Clinical/TracedBorders/',temp); %% Change this path!!
elseif ispc
    bordersPath = strcat('E:\ClinicalDataArchive\TracedBorders\',temp);
end

if ~arfidata.df_flag
    fname = strcat(bordersPath,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_arfi.mat');
else
    fname = strcat(bordersPath,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_swei.mat');
end

if exist(fname,'file')
    prior_trace = 1;
    load(fname); % loads variable called traced_borders
else
    prior_trace = 0;
    traced_borders = [];
end

dims = size(arfidata.disp);
ndepth = dims(1); nacqT = dims(2); ntrackT = dims(3);
edge = [arfidata.axial(1) arfidata.axial(end)];
trackPRF = 1000/par.priusec(1); % kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Gate Parameters
if prior_trace
    borders = traced_borders;
elseif ~prior_trace
    borders = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2])';
    borders = repmat(borders,[1 nacqT]);
end

if (min(borders(:))<arfidata.axial(1) || max(borders(:))>arfidata.axial(end))
    warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(borders(:)),max(borders(:)),arfidata.axial(1),arfidata.axial(end));
end

% Set Display Parameters
fig = figure;
set(fig,'Name','TTE ARFI','UserData','main_fig','InvertHardCopy', 'off')
if isunix
    set(fig,'Position',[1201 386 1920 1070])
    dispPar.fsize = 12;
elseif ispc
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    dispPar.fsize = 10;
end

if strcmpi(options.display.theme,'light')
    dispPar.fig = [1 1 1]; dispPar.txt = 'k'; dispPar.ax = [1 1 1];
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

set(fig,'Color',dispPar.fig)
dispPar.cmap = colormap(parula);
dispPar.cmap(end,:) = dispPar.ax;
dispPar.corder = jet(options.display.n_pts);
dispPar.trace_cols = prism(size(borders,1));
if size(borders,1)>1, dispPar.trace_cols(end,:) = dispPar.trace_cols(end-1,:); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Info Boxes
% Create Data Set Info Box
[~,name] = fileparts(pwd);
if par.pushDurationusec>5
    seq_type = 'Clinical';
elseif trackPRF<2.5
    seq_type = 'Continuous';
else
    seq_type = 'No Push';
end

ln1 = sprintf('Subject ID: %s',name);
ln2 = sprintf('Set Info: %d - %s',options.dataflow.setID,seq_type);
ln3 = sprintf('TimeStamp: %s',options.timeStamp);
info_box = annotation('textbox',[0.01 0.9, 0.085, 0.085],'EdgeColor',dispPar.txt,'string',{ln1,ln2,'',ln3},'Color',dispPar.txt,'Interpreter','none','FontSize',dispPar.fsize,'FontWeight','Bold');

% Create Parameters Info Box
ln1 = 'Motion Filter Settings:';
ln2 = sprintf('Method - %s',options.motionFilter.method);
if strfind(options.motionFilter.method,'PF')
    ln3 = sprintf('Cutoff = %s Hz',mat2str(options.motionFilter.Cutoff));
else
    ln3 = '';
end

if strfind(options.motionFilter.method,'Poly')
    ln4 = sprintf('TimeRange (pre) = \n%s ms', mat2str(options.motionFilter.timeRange_pre));
    ln5 = sprintf('TimeRange (push) = \n%s ms', mat2str(options.motionFilter.timeRange_push));
    ln6 = sprintf('Order = %d',options.motionFilter.order);
else
    ln4 = ''; ln5 = ''; ln6 = '';
end

par_box1 = annotation('textbox',[0.01 0.70, 0.085, 0.18],'EdgeColor',dispPar.txt,'string',{ln1,ln2,'',ln3,'',ln4,ln5,ln6},'Color',dispPar.txt,'Interpreter','none','FontSize',dispPar.fsize,'FontWeight','Bold');
ln1 = 'Display Settings:';
if options.display.cc_filt
    ln2 = sprintf('CC Threshhold = %1.3f',options.display.cc_thresh);
else
    ln2 = '';
end

ln3 = sprintf('Median Filter = [%1.1f mm %1.2f s]',options.display.medfilt(1),options.display.medfilt(2));
ln4 = sprintf('Systole Range =\n%s/100',mat2str(100*options.display.sysRange));
ln5 = sprintf('Diastole Range =\n%s/100',mat2str(100*options.display.diaRange));
ln6 = '(as fraction of normalized duration of heartbeat)';
par_box2 = annotation('textbox',[0.01 0.45, 0.085, 0.2],'EdgeColor',dispPar.txt,'string',{ln1,ln2,ln3,ln4,ln5,ln6},'Color',dispPar.txt,'Interpreter','none','FontSize',dispPar.fsize,'FontWeight','Bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Bmode Cine
[recBors,bimg_ax] = drawBmode(fig,bdata,edge,borders,dispPar,par,options.display.playCine,options);
% for i=1:100;drawBmode(fig,bdata,edge,borders,dispPar,par,options.display.playCine,options);end
% Display M-mode IQ
[IQBors,IQ_ax] = drawIQ(fig,arfidata.IQ,arfidata.IQaxial,arfidata.acqTime,edge,borders,dispPar,par,options.display.IQrange,arfidata.df_flag,trackPRF);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trace Tissue border on Bmode data
% threshhold = 0.75;
% [pks, idx] = findpeaks(double(bdata.ecg(:,2)),'MinPeakHeight',threshhold);
% t_pk = bdata.ecg(idx(1),1);
% tnorm = (bdata.t-t_pk)*(bdata.hr/60);
% mod_tnorm = mod(tnorm,1);
% range(tnorm)
% frame_idx = [1:find(tnorm>tnorm(1)+1,1,'first')+1];
% if isempty(frame_idx) frame_idx = [1:length(bdata.t)];end
%
% nbms = 51; mid = ceil(size(bdata.bimg,2)/2);
% temp = bdata.bimg(:,mid-(nbms-1)/2:mid+(nbms-1)/2,frame_idx);
% temp = reshape(temp,size(temp,1),size(temp,2)*size(temp,3));
% xvec = linspace(0,length(frame_idx)*range(bdata.blat(mid-(nbms-1)/2:mid+(nbms-1)/2)),size(temp,2));
%
% figure;
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% imagesc(xvec,bdata.bax,temp);colormap(gray);axis image;
%
% for i=1:size(frame_idx,2)
%     [x,y] = ginput(1);
%     hold on
%     plot(x,y,'kx','MarkerSize',6);
%     top(i) = y;
%     clear x y
% end
% for i=1:size(frame_idx,2)
%     [x,y] = ginput(1);
%     hold on
%     plot(x,y,'kx','MarkerSize',6);
%     bot(i) = y;
%     clear x y
% end
% close gcf
%
% [pks2, idx2] = findpeaks(double(arfidata.ecg(:,2)),'MinPeakHeight',threshhold);
% t_pk2 = arfidata.ecg(idx2(1),1);
% tnorm_arfi = (arfidata.acqTime-t_pk2)*(arfidata.hr/60);
%
%
% mod_tnorm_arfi = mod(tnorm_arfi,1);
%
% for i=1:nacqT
%     temp = find(mod_tnorm(frame_idx)>mod_tnorm_arfi(i),1,'first');
%     if isempty(temp); temp = find(mod_tnorm>mod_tnorm_arfi(i)-median(diff(mod_tnorm_arfi)),1,'first');end
%     index(i) = temp;
% end
%
% top_trace = smooth(top(index));  bot_trace = smooth(bot(index));
%
% hold(ax12,'on')
%
% l3_alt = plot(linspace(0,arfidata.acqTime(end),nacqT),top_trace,'r','Linewidth',2,'Parent',ax12);
% l4_alt = plot(linspace(0,arfidata.acqTime(end),nacqT),bot_trace,'r','Linewidth',2,'Parent',ax12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display M-mode ARFI
[pre,push,dt,options,preBors,pushBors,pre_ax,push_ax] = drawMmode(fig,arfidata,edge,borders,dispPar,par,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display ECG
if ~isempty(arfidata.ecg)
    [samples,ecg_ax] = drawECG(fig,arfidata.ecg,arfidata.acqTime,dispPar);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate
[pre_trace,push_trace,layer_idx,dispTrace_ax] = drawDispTraces(fig,pre,push,borders,arfidata.axial,arfidata.acqTime,dispPar,options);
% if ~isempty(arfidata.ecg)
%     [pre_DR,push_DR] = computeRatios(ecg_ax,dispTrace_ax,pre_trace,push_trace,arfidata.acqTime,arfidata.ecg,arfidata.hr,samples,dispPar,options);
%     for i=1:length(whos('rr_*'));eval(sprintf('delete(rr_%d)',i)), eval(sprintf('clear rr_%d',i)), end
%     ratio_box = annotation('textbox',[0.05 0.07 0.4 0.025],'String','Ratios: ','Color','w','EdgeColor','w','FontWeight','Bold');
%     for i=1:length(layer_idx)
%         str = num2str(push_DR(layer_idx(i)),'%1.2f');
%         eval(sprintf('rr_%d = annotation(''textbox'',[0.1+0.04*(i-1) 0.07 0.04 0.025],''LineStyle'',''none'',''Color'',dispPar.trace_cols(layer_idx(i),:),''FontWeight'',''bold'',''string'',str);',i));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot IQ Traces
% keyboard
%
% %%
% if options.dataflow.oneSided
%     temp = squeeze(arfidata.IQ(:,1:par.nBeams:end,:));
% else
%     temp = squeeze(arfidata.IQ(:,ceil(par.nBeams/2):par.nBeams:end,:));
% end
% I = real(temp); Q = imag(temp);
% factor = 5;
% D = size(I);
% D(1) = D(1).*factor;
% [Iup, Qup] = computeUpsampledIQdata(I,Q,factor);
% Iup = reshape(Iup,D); Qup = reshape(Qup,D);
% test = unwrap(angle(complex(Iup,Qup)));
% axial = interp(arfidata.IQaxial,factor);
%
%
% test = reshape(permute(test,[1 3 2]),size(test,1),[]);
% figure(200);for i=1:size(test,2);plot(axial,test(:,i));title(mod(i,size(test,2)/nacqT));ylim([-100 100]);grid on;pause;end

%% Trace Borders
if options.display.extras == -1  % no breaks
    trace_input = 'n';
    %     flatgate_input = 'y';
end
if prior_trace
    trace_input = input('\nDo you want to retrace tissue borders? (y/n) [n]: ','s');
    %     flatgate_input = input('\nDo you want to use a single flat gate? (y/n) [n]: ','s');
elseif (~prior_trace && options.display.extras ~= -1)
    trace_input = input('\n>>>>> Do you want to trace tissue borders? (y/n) [n]: ','s');
    %     flatgate_input = 'y';
end
% if strcmpi(flatgate_input,'y')
%     gate = (par.pushFocalDepth + options.display.gateOffset + [-options.display.gateWidth/2 options.display.gateWidth/2]);
%     gate = repmat(gate,[nacqT 1]);
%     keyboard
%     hold(ax11,'on')
%     r2 = rectangle('Position',[-2 min(gate(:)) 4 options.display.gateWidth],'EdgeColor','g','Linewidth',1,'Parent',ax11);
%     hold(ax12,'on')
%     for i=1:size(borders,1)
%         eval(sprintf('IQBor_%d = plot(linspace(0,arfidata.acqTime(end),nacqT),borders(i,:),''Color'',%s,''Linewidth'',1,''Parent'',ax12);',i,mat2str(dispPar.trace_cols(i,:))))
%     end
%     if (min(gate(:))<arfidata.axial(1) || max(gate(:))>arfidata.axial(end))
%         warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
%     end
% end
if strcmpi(trace_input,'y')
    % Delete existing borders
    for i=1:size(borders,1)-1; delete(recBors{i});end
    for i=1:size(borders,1), delete(IQBors{i}), if ~isempty(preBors),delete(preBors{i}),end, delete(pushBors{i}),end
    n = input('How many borders do you want to trace? ');
    borders = nan(n,nacqT);
    dispPar.trace_cols = prism(size(borders,1));
    if size(borders,1)>1, dispPar.trace_cols(end,:) = dispPar.trace_cols(end-1,:); end
    fprintf(1,'\nReady to trace %d borders (shallow to deep)...\nClick to define points, hit space to end tracing\n',n)
    for i=1:n
        fprintf(1,sprintf('Trace border %d...Hit enter when done\n',i));
        [borders(i,:),mark{i}] = traceBorder(nacqT,arfidata.acqTime,IQ_ax);
    end
    for i=1:n; delete(mark{i}); end
    % Draw new borders
    [recBors,bimg_ax] = drawBmode(fig,bdata,edge,borders,dispPar,par,options.display.playCine,options);
    [IQBors,IQ_ax] = drawIQ(fig,arfidata.IQ,arfidata.IQaxial,arfidata.acqTime,edge,borders,dispPar,par,options.display.IQrange,arfidata.df_flag,trackPRF);
    
    %     hold(bimg_ax,'on')
    %     r2 = rectangle('Position',[-2 min(borders(:)) 4 range(borders(:))],'EdgeColor','g','Linewidth',1,'Parent',ax11);
    %     hold(IQ_ax,'on')
    %     for i=1:size(borders,1)
    %         eval(sprintf('IQBor_%d = plot(linspace(0,arfidata.acqTime(end),nacqT),borders(i,:),''Color'',%s,''Linewidth'',1,''Parent'',ax12);',i,mat2str(dispPar.trace_cols(i,:))))
    %     end
    fprintf(1,'Borders Traced.\n')
    traced_borders = borders;
    if (min(borders(:))<arfidata.axial(1) || max(borders(:))>arfidata.axial(end))
        warning('Depth gate requested (%2.2f-%2.2f mm) falls outside the range over which displacements are computed (%2.2f-%2.2f mm)',min(gate(:)),max(gate(:)),arfidata.axial(1),arfidata.axial(end));
    end
    % Save traced gate
    fprintf(1,'>>> Saving Traced Borders...\n');
    [~,temp] = fileparts(pwd);
    if isunix
        bordersPath = strcat('/emfd/vrk4/Transthoracic_Clinical/TracedBorders/',temp); %% Change this path!!
    elseif ispc
        bordersPath = strcat('E:\ClinicalDataArchive\TracedBorders\',temp);
    end
    if ~exist(bordersPath,'dir'); mkdir(bordersPath); end
    if ~arfidata.df_flag
        fname = strcat(bordersPath,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_arfi.mat');
    else
        fname = strcat(bordersPath,filesep,sprintf('%02d',options.dataflow.setID),'_',options.timeStamp,'_swei.mat');
    end
    save(fname,'traced_borders');
end
border_idx = size(borders);
for i=1:size(borders,1)
    for j=1:size(borders,2)
        border_idx(i,j) = find(arfidata.axial>borders(i,j),1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display M-mode ARFI
[pre,push] = drawMmode(fig,arfidata,edge,borders,dispPar,par,options);% NaN out displacements filtered out by cc_thresh
% dispARFIMovie(arfidata,edge,borders,dispPar,par,options);
% arfidata.disp(mask==0) = nan;
% if ~strcmpi(options.motionFilter.method,'Off')
%     arfidata_mf_pre.disp(mask==0) = nan;
%     arfidata_mf_push.disp(mask==0) = nan;
% end
% arfidata.disp(mask==0) = nan;
% if ~strcmpi(options.motionFilter.method,'Off')
%     arfidata_mf_pre.disp(mask==0) = nan;
%     arfidata_mf_push.disp(mask==0) = nan;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate
[pre_trace,push_trace,layer_idx,dispTrace_ax] = drawDispTraces(fig,pre,push,borders,arfidata.axial,arfidata.acqTime,dispPar,options);
if ~isempty(arfidata.ecg)
    [pre_DR,push_DR] = computeRatios(ecg_ax,dispTrace_ax,pre_trace,push_trace,arfidata.acqTime,arfidata.ecg,arfidata.hr,samples,dispPar,options);
    for i=1:length(whos('rr_*'));eval(sprintf('delete(rr_%d)',i)), eval(sprintf('clear rr_%d',i)), end
    ratio_box = annotation('textbox',[0.05 0.05 0.4 0.045],'String','Ratios: ','Color','w','EdgeColor','w','FontWeight','Bold');
    for i=1:length(layer_idx)
        str = sprintf(' %1.2f\n(%1.2f)',push_DR(layer_idx(i)),pre_DR(layer_idx(i)));
        eval(sprintf('rr_%d = annotation(''textbox'',[0.1+0.04*(i-1) 0.07 0.04 0.025],''LineStyle'',''none'',''Color'',dispPar.trace_cols(layer_idx(i),:),''FontWeight'',''bold'',''string'',str);',i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% keyboard
% return
%% Contrast and Raw Motion Analysis

% Contrast Computation
nbeams = 1;
detdata = abs(arfidata.IQ);
detdata = detdata(:,:,[1:options.dispEst.dims(1)+1,options.dispEst.dims(1)+options.dispEst.dims(2)+options.dispEst.dims(3)+2:options.dispEst.dims(end)]); % Remove push/reverb frames
detdata = reshape(detdata,size(detdata,1),par.nBeams,nacqT,[]);
if ~options.dataflow.oneSided
detdata = detdata(:,ceil(par.nBeams/2):ceil(par.nBeams/2)+nbeams,:,:); % Remove far off lateral rx beams
else
detdata = detdata(:,1:1+nbeams,:,:);
end
    for j=1:nacqT
    for i=1:size(borders,1)
        IQ_idx(i,j) = find(arfidata.IQaxial>borders(i,j),1,'first');
    end
end
RV_idx = 1; %input(sprintf('Enter layer number for the chamber: Shallowest(1) to Deepest(%d): ',size(borders,1)-1));
sep_idx = 2; %input(sprintf('Enter layer number for the background: Shallowest(1) to Deepest(%d): ',size(borders,1)-1));
% LV_idx = 3; %input(sprintf('Enter layer number for the chamber: Shallowest(1) to Deepest(%d): ',size(borders,1)-1));

RV_mean = nan(nacqT,size(detdata,4)); RV_std = nan(size(RV_mean));
sep_mean = nan(nacqT,size(detdata,4)); sep_std = nan(size(sep_mean));
% LV_mean = nan(nacqT,size(detdata,4)); LV_std = nan(size(LV_mean));

for j=1:size(detdata,4)
    for i=1:nacqT
        RV = detdata(IQ_idx(RV_idx,i):IQ_idx(RV_idx+1,i),:,i,j);
        sep = detdata(IQ_idx(sep_idx,i):IQ_idx(sep_idx+1,i),:,i,j);
%         LV = detdata(IQ_idx(LV_idx,i):IQ_idx(LV_idx+1,i),:,i,j);

        RV_mean(i,j) = nanmean(RV(:)); RV_std(i,j) = nanstd(RV(:));
        sep_mean(i,j) = nanmean(sep(:)); sep_std(i,j) = nanstd(sep(:));
%         LV_mean(i,j) = nanmean(LV(:)); LV_std(i,j) = nanstd(LV(:));
    end
end

SNR(:,:,1) = RV_mean./RV_std; SNR(:,:,2) = sep_mean./sep_std; %SNR(:,:,3) = LV_mean./LV_std;
contrast(:,:,1) = db(sep_mean./RV_mean); %contrast(:,:,2) = sep_mean./LV_mean;
s2c(:,:,1) = db((sep_mean - RV_mean)./RV_mean); %s2c(:,:,2) = 20*log10((sep_mean - LV_mean)./LV_mean);

% Raw Motion RMS Computation
t_range = [-7 -2];

raw = arfidata.disp;
filt = arfidata.disp_mf_pre;

for ii=1:size(t_range,1)
    
    disp(sprintf('Calculating Residuals for t_range = [%2.2f %2.2f] ms',t_range(ii,1),t_range(ii,2)))
    raw_rms = sqrt(mean(raw(:,:,find(arfidata.trackTime>t_range(ii,1),1)-1:find(arfidata.trackTime>t_range(ii,2),1)-1).^2,3));
    filt_rms = sqrt(mean(filt(:,:,find(arfidata.trackTime>t_range(ii,1),1)-1:find(arfidata.trackTime>t_range(ii,2),1)-1).^2,3));
    for i=1:size(borders,1)-1
        for j=1:nacqT
            raw_rms_seg{i}(1:border_idx(i+1,j)-border_idx(i,j)+1,j,ii) = raw_rms(border_idx(i,j):border_idx(i+1,j),j);
            filt_rms_seg{i}(1:border_idx(i+1,j)-border_idx(i,j)+1,j,ii) = filt_rms(border_idx(i,j):border_idx(i+1,j),j);
        end
    end
end

figure
set(gcf,'color',dispPar.fig,'Position',[969    49   944   948])
ax1 = subplot(411);
shadedErrorBar(arfidata.acqTime,mean(contrast(:,:,1),2),std(contrast(:,:,1),[],2),'-*r',1)
% errorbar(arfidata.acqTime,mean(contrast(:,:,1),2),std(contrast(:,:,1),[],2),'Linewidth',2,'Color','r');
hold on
shadedErrorBar(arfidata.acqTime,mean(s2c(:,:,1),2),std(s2c(:,:,1),[],2),'-*g',1)
% errorbar(arfidata.acqTime,mean(s2c(:,:,1),2),std(s2c(:,:,1),[],2),'Linewidth',2,'Color','g');
title(sprintf('Septum vs. RV Chamber\n Mean Contrast = %2.2f dB     Mean Signal to Clutter = %2.2f dB',mean(mean(contrast(:,:,1))),mean(mean(s2c(:,:,1)))),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
xlim([0 max(arfidata.acqTime)]); ylim([-6 31]);
ylabel('dB')
set(ax1,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','ygrid','on','gridLineStyle','--');
% leg = legend([],'Contrast',[],'Signal-to-Clutter');
% set(leg,'TextColor',dispPar.txt,'EdgeColor',dispPar.txt)

if ~isempty(arfidata.ecg)
    ax2 = subplot(412);
    plot(arfidata.ecg(:,1),arfidata.ecg(:,2),'Linewidth',2);
    set(ax2,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--')
end

ax3 = subplot(413);
shadedErrorBar(arfidata.acqTime,nanmean(raw_rms_seg{2}),nanstd(raw_rms_seg{2}),'-*r',1)
hold all
shadedErrorBar(arfidata.acqTime,nanmean(filt_rms_seg{2}),nanstd(filt_rms_seg{2}),'-*g',1)
% plot(arfidata.acqTime,mean_rms,'Linewidth',2)
title(sprintf('RMS of Raw Pre ARF Motion Profiles: [%1.0f %1.0f] ms',t_range(1),t_range(2)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
xlim([0 max(arfidata.acqTime)]); ylim([-40 360]);
ylabel('RMS Disp (\mum)')
set(ax3,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','ygrid','on','gridLineStyle','--');
% leg = legend(sprintf('[%1.0f   %1.0f] ms',t_range(1,1),t_range(1,2)),sprintf('[%1.0f   %1.0f] ms',t_range(2,1),t_range(2,2)));
% set(leg,'TextColor',dispPar.txt,'EdgeColor',dispPar.txt)

ax4 = subplot(414);
shadedErrorBar(arfidata.acqTime,nanmean(filt_rms_seg{2}),nanstd(filt_rms_seg{2}),'-*g',1)
% plot(arfidata.acqTime,mean_rms,'Linewidth',2)
title(sprintf('RMS of Filtered Pre ARF Motion Profiles: [%1.0f %1.0f] ms',t_range(1),t_range(2)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
xlim([0 max(arfidata.acqTime)]); ylim([-4 14]);
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)

ylabel('RMS Disp (\mum)')
set(ax4,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--','ygrid','on','gridLineStyle','--');
% leg = legend(sprintf('[%1.0f   %1.0f] ms',t_range(1,1),t_range(1,2)),sprintf('[%1.0f   %1.0f] ms',t_range(2,1),t_range(2,2)));
% set(leg,'TextColor',dispPar.txt,'EdgeColor',dispPar.txt)
%%
% keyboard
% return
%% Gross Axial Wall Thickness/Strain

% test = diff(borders);
% figure;
% plot(arfidata.ecg(:,1),arfidata.ecg(:,2))
% hold all
% plot(arfidata.acqTime,test(2,:)./test(2,1))
% % plot(arfidata.acqTime,push_trace)
% grid on
% 
% return
%% Plot Surfaces
% figure;
% for i=1:nacqT
%
%     plot(arfidata.trackTime,squeeze(arfidata.disp(gate_idx(i,1):gate_idx(i,2),i,:))')
%     grid on;set(gca,'ylim',[-300 300])
%
% %     surf(arfidata.trackTime,arfidata.axial(gate_idx(i,1):gate_idx(i,2)),double(squeeze(arfidata.disp(gate_idx(i,1):gate_idx(i,2),i,:))))
% %     set(gca,'zlim',[-200 200],'ydir','reverse')
%
%     title(num2str(arfidata.acqTime(i)))
%     pause
% end
%
%
% % Plot multiple acqTime points
%
% clear tt tt_idx
% tt = [0.2 0.95 1.7];
% for i=1:length(tt)
%     tt_idx(i) = find(round(double(arfidata.acqTime),2)==tt(i));
% end
%
% figure
% plot(arfidata.trackTime,squeeze(arfidata.disp(gate_idx(tt_idx(1),1):gate_idx(tt_idx(1),2),tt_idx(1),:))','b')
% hold all
% plot(arfidata.trackTime,squeeze(arfidata.disp(gate_idx(tt_idx(2),1):gate_idx(tt_idx(2),2),tt_idx(2),:))','r')
% plot(arfidata.trackTime,squeeze(arfidata.disp(gate_idx(tt_idx(3),1):gate_idx(tt_idx(3),2),tt_idx(3),:))','g')
% grid on; set(gca,'ylim',[-300 300])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% keyboard
%% Analysis of Motion Filter Residuals (if BPF)/ Jitter (if LPF)

% for i=1:size(borders,1)-1
%     jitter_layer{i} = nan(max(border_idx(i+1,:))-min(border_idx(i,:))+1,nacqT);
%     res_layer{i} = nan(max(border_idx(i+1,:))-min(border_idx(i,:))+1,nacqT);
% end
% 
% disp('Calculating jitter rms...')
% jitter_rms = sqrt(mean(arfidata.jitter.^2,3));
% for i=1:size(borders,1)-1
%     for j=1:nacqT
%         jitter_layer{i}(1:border_idx(i+1,j)-border_idx(i,j)+1,j) = jitter_rms(border_idx(i,j):border_idx(i+1,j),j);
%     end
% end
% 
% t_range = [arfidata.trackTime(1) -1]
% % t_range = [-1 5];
% 
% disp('Calculating Residuals...')
% if ~isempty(arfidata.disp_mf_push)
%     res = arfidata.disp_mf_push;
% else
%     res = arfidata.disp;
% end
% res_rms = sqrt(mean(res(:,:,find(arfidata.trackTime>t_range(1),1)-1:find(arfidata.trackTime>t_range(2),1)-1).^2,3));
% for i=1:size(borders,1)-1
%     for j=1:nacqT
%         res_layer{i}(1:border_idx(i+1,j)-border_idx(i,j)+1,j) = res_rms(border_idx(i,j):border_idx(i+1,j),j);
%     end
% end
% 
%     figure
%     set(gcf,'color',dispPar.fig,'units','normalized','outerposition',[0 0 1 1])
%     jitter_ax = subplot(311);
%     for i=1:size(borders,1)-1
%     shadedErrorBar(arfidata.acqTime,nanmean(jitter_layer{i}),nanstd(jitter_layer{i}),{'Linewidth',2,'color',dispPar.trace_cols(i,:)},1)
%     hold on
%     end
%     hold off
%     xlim([0 max(arfidata.acqTime)]); ylim([-0.5 2.5]);grid on
%     ylabel('Displacement (\mum)')
%     title('RMS Jitter','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%     set(jitter_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--')
% 
%     res_ax = subplot(312);
%     for i=1:size(borders,1)-1
%     shadedErrorBar(arfidata.acqTime,nanmean(res_layer{i}),nanstd(res_layer{i}),{'Linewidth',2,'color',dispPar.trace_cols(i,:)},1)
%     hold on
%     end
%     hold off
%     xlim([0 max(arfidata.acqTime)]);grid on
%     ylabel('Displacement (\mum)')
%     title(sprintf('RMS Residual: %s - %s Hz',options.motionFilter.method,mat2str(options.motionFilter.Cutoff)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%     set(res_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--')
% 
%     if ~isempty(arfidata.ecg)
%         ecg_ax = subplot(313);
%         plot(arfidata.ecg(:,1),arfidata.ecg(:,2),'Linewidth',2);
%         xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
%         set(ecg_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--')
%     end
% 
%     figure
%     set(gcf,'Position',[950 350 800 200])
%     nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
%     nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));
%     imagesc(arfidata.acqTime,arfidata.axial,medfilt2(double(jitter_rms),[nax nt]),[0 0.75])
%     hold all
%     plot(arfidata.acqTime,borders,'Linewidth',2,'color','k');colorbar

% else
% 
%     disp('Calculating Residuals...')
%     res = arfidata.disp_mf_push;
%     res_rms = sqrt(mean(res.^2,3));
%     for i=1:size(borders,1)-1
%         for j=1:nacqT
%             res_layer{i}(1:border_idx(i+1,j)-border_idx(i,j)+1,j) = res_rms(border_idx(i,j):border_idx(i+1,j),j);
%         end
%     end
% 
%     figure
%     set(gcf,'color',dispPar.fig,'units','normalized','outerposition',[0 0 1 1])
% 
%     res_ax = subplot(211);
%     for i=1:size(borders,1)-1
%     shadedErrorBar(arfidata.acqTime,nanmean(res_layer{i}),nanstd(res_layer{i}),{'Linewidth',2,'color',dispPar.trace_cols(i,:)},1)
%     hold on
%     end
%     hold off
%     xlim([0 max(arfidata.acqTime)]);grid on
%     ylabel('Displacement (\mum)')
%     title(sprintf('RMS Residual: %s - %s Hz',options.motionFilter.method,mat2str(options.motionFilter.Cutoff)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%     set(res_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--')
% 
%     if ~isempty(arfidata.ecg)
%         ecg_ax = subplot(212);
%         plot(arfidata.ecg(:,1),arfidata.ecg(:,2),'Linewidth',2);
%         xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
%         set(ecg_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--')
%     end
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if size(borders,1) == 2
%
%     for i=1:nacqT
%         idx(:,i) = round(linspace(find(arfidata.axial>borders(1,i),1),find(arfidata.axial>borders(2,i),1),options.display.n_pts));
%     end
%     raw = nan(max(idx(end,:)-idx(1,:)+1),nacqT,ntrackT);
%     flt = nan(size(raw));
%
%     for i=1:nacqT
%         raw(1:idx(end,i)-idx(1,i)+1,i,:) = arfidata.disp(idx(1,i):idx(end,i),i,:);
%         flt(1:idx(end,i)-idx(1,i)+1,i,:) = arfidata.disp_mf_push(idx(1,i):idx(end,i),i,:);
%     end
%     if strcmpi(options.motionFilter.method,'LPF')
%         res = flt - raw;
%         jitter = nanmean(abs(res),3);
%         figure
%         set(gcf,'color',dispPar.fig,'Position',[329 347 1337 489])
%         ax1 = subplot(211);
%         shadedErrorBar(arfidata.acqTime,nanmean(jitter),nanstd(jitter)',{'Linewidth',2,'color','r'},1)
%         xlim([0 max(arfidata.acqTime)]); ylim([-2.5 12.5]); grid on
%         ylabel('Displacement (\mum)')
%         title(sprintf('Mean Jitter = %2.2f \\mum',mean(mean(jitter))),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         set(ax1,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--')
%
%         if ~isempty(arfidata.ecg)
%             ax2 = subplot(212);
%             plot(arfidata.ecg(:,1),arfidata.ecg(:,2),'Linewidth',2);
%             xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
%             set(ax2,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--')
%         end
%     else
%         res = flt;
%         E = nansumm(res.^2,3); % energy of residuals
%         % save energy trace
%     end
%
%
% else
%
%     for i=1:nacqT
%         RV_idx(:,i) = round(linspace(find(arfidata.axial>borders(1,i),1),find(arfidata.axial>borders(2,i),1),options.display.n_pts));
%         sep_idx(:,i) = round(linspace(find(arfidata.axial>borders(2,i),1),find(arfidata.axial>borders(3,i),1),options.display.n_pts));
%     end
%     RV_raw = nan(max(RV_idx(end,:)-RV_idx(1,:)+1),nacqT,ntrackT);
%     RV_flt = nan(size(RV_raw));
%     sep_raw = nan(max(sep_idx(end,:)-sep_idx(1,:)+1),nacqT,ntrackT);
%     sep_flt = nan(size(sep_raw));
%
%     for i=1:nacqT
%         RV_raw(1:RV_idx(end,i)-RV_idx(1,i)+1,i,:) = arfidata.disp(RV_idx(1,i):RV_idx(end,i),i,:);
%         RV_flt(1:RV_idx(end,i)-RV_idx(1,i)+1,i,:) = arfidata.disp_mf_push(RV_idx(1,i):RV_idx(end,i),i,:);
%
%         sep_raw(1:sep_idx(end,i)-sep_idx(1,i)+1,i,:) = arfidata.disp(sep_idx(1,i):sep_idx(end,i),i,:);
%         sep_flt(1:sep_idx(end,i)-sep_idx(1,i)+1,i,:) = arfidata.disp_mf_push(sep_idx(1,i):sep_idx(end,i),i,:);
%     end
%     if strcmpi(options.motionFilter.method,'LPF')
%         RV_res = RV_flt - RV_raw;
%         sep_res = sep_flt - sep_raw;
%         RV_jitter = nanmean(abs(RV_res),3);
%         sep_jitter = nanmean(abs(sep_res),3);
%         figure
%         set(gcf,'color',dispPar.fig,'Position',[329 347 1337 489])
%         ax1 = subplot(211);
%         shadedErrorBar(arfidata.acqTime,nanmean(RV_jitter),nanstd(RV_jitter)',{'Linewidth',2,'color','r'},1)
%         hold on
%         shadedErrorBar(arfidata.acqTime,nanmean(sep_jitter),nanstd(sep_jitter)',{'Linewidth',2,'color','g'},1)
%         ylabel('Displacement (\mum)')
%         title(sprintf('Mean RV Jitter  = %2.2f \\mum    Mean Septum Jitter = %2.2f \\mum',nanmean(nanmean(RV_jitter)),nanmean(nanmean(sep_jitter))),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         xlim([0 max(arfidata.acqTime)]); ylim([-2.5 15]); grid on
%         set(ax1,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','gridLineStyle','--')
%         if ~isempty(arfidata.ecg)
%             ax2 = subplot(212);
%             plot(arfidata.ecg(:,1),arfidata.ecg(:,2),'Linewidth',2);
%             xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
%             set(ax2,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--')
%         end
%     else
%         RV_res = RV_flt;
%         sep_res = sep_flt;
%         RV_E = nansum(RV_res.^2,3); % energy of residuals
%         sep_E = nansum(sep_res.^2,3);
%         % save energy trace
%     end
%
% end

%%
keyboard
%% Loop to go through Acquisition Time
extra1_input = input('\n>>>>> Do you want to look at correlation and disp vs. time plots? (y/n) [n]: ','s');
if strcmpi(extra1_input,'y')
    options.display.extras = 1;
    if size(borders,1)>2
        border_input = input('\n>>>>> Pick a layer number for disp vs. time plots: ');
    else
        border_input = 1;
    end
end
% Displacement vs. Time Plots
if options.display.extras > 0;
    figure(100);set(100,'Color',dispPar.fig,'units','normalized','outerposition',[0 0 1 1])
    if ~isempty(arfidata.ecg)
        [samples,ecg_ax2,pt] = drawECG(100,arfidata.ecg,arfidata.acqTime,dispPar,[0.575 0.075 0.4 0.075]);
    end
    IQ = abs(arfidata.IQ);
    IQ = db(IQ/max(max(max(IQ(:,[1:options.dispEst.dims(1),sum(options.dispEst.dims(1:3))+1:end])))));
    fprintf(1,'\n\n>>> Press Left/Right to move Back/Forward and Space to play through\n');
    if options.display.cc_filt
        mask = arfidata.cc>options.display.cc_thresh;
    else
        mask = ones(size(arfidata.disp));
    end
    idx = nan(options.display.n_pts,nacqT);
    for i=1:nacqT
        idx(:,i) = round(linspace(find(arfidata.axial>borders(border_input,i),1),find(arfidata.axial>borders(border_input+1,i),1),options.display.n_pts));
    end
    raw_rng = [-650 650];
    mf_rng = [-10 10]; %1.5*options.display.dispRange;
    if (strcmpi(options.motionFilter.method,'Off') || strcmpi(options.motionFilter.method,'LPF'))
        mf_rng = raw_rng;
    end
    i=1; skip = 0;
    %     keyboard
    %     filename = 'C:\Users\vrk4\Desktop\test.gif';
    while i<=nacqT
        if ~isempty(arfidata.ecg)
            delete(pt); %set(pt,'Visible','off')
            hold(ecg_ax2,'on')
            pt = plot(arfidata.acqTime(i),samples(i),'ro','Parent',ecg_ax2,'Markersize',10,'Markerfacecolor','r');
            hold(ecg_ax2,'off')
        end
        if i==1
            dvt_ax1 = axes('Position',[0.575 0.625 0.4 0.275],'Parent',100);
            set(dvt_ax1,'xlim',[arfidata.trackTime(1) arfidata.trackTime(end)],'ylim',raw_rng,'ColorOrder',dispPar.corder,'Color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold')
            hold(dvt_ax1,'on')
            t = title('');
            dvt_ax2 = axes('Position',[0.575 0.25 0.4 0.275],'Parent',100);
            set(dvt_ax2,'xlim',[arfidata.trackTime(1) arfidata.trackTime(end)],'ylim',mf_rng,'ColorOrder',dispPar.corder,'Color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold')
            hold(dvt_ax2,'on')
            IQ_ax2 =  axes('Position',[0.025 0.15 0.1 0.75],'color',dispPar.ax);
            corr_ax =  axes('Position',[0.1475 0.15+(max(arfidata.IQaxial)-max(arfidata.axial))*(0.75/max(arfidata.IQaxial)) 0.1 0.75*range(arfidata.axial/range(arfidata.IQaxial))],'color',dispPar.ax);
            raw_ax = axes('Position',[0.2825 0.15+(max(arfidata.IQaxial)-max(arfidata.axial))*(0.75/max(arfidata.IQaxial)) 0.1 0.75*range(arfidata.axial/range(arfidata.IQaxial))],'color',dispPar.ax);
            mf_ax = axes('Position',[0.405 0.15+(max(arfidata.IQaxial)-max(arfidata.axial))*(0.75/max(arfidata.IQaxial)) 0.1 0.75*range(arfidata.axial/range(arfidata.IQaxial))],'color',dispPar.ax);
        else
            cla(dvt_ax1); cla(dvt_ax2);delete(t);
            %             cla(motion_ax)
        end
        % Extra Fig: Correlations and disp vs. time plots
        set(100,'CurrentAxes',dvt_ax1)
        plot(arfidata.trackTime,squeeze(arfidata.disp(idx(:,i),i,:)),'--','Parent',dvt_ax1)
        plot(options.display.t_disp_pre*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax1)
        plot((options.display.t_disp_pre-dt)*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax1)
        plot(options.display.t_disp_push*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax1)
        plot(arfidata.trackTime(options.dispEst.dims(1))*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax1)
        grid(dvt_ax1,'on')
        xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        t = title(sprintf('Raw ARFI Displacement Profiles (within Depth Gate)\nPush # %d (t = %2.2f s)\nTraced Tissue Layer: %d (from top)',i,arfidata.acqTime(i),border_input),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',dvt_ax1);
        set(100,'CurrentAxes',dvt_ax2)
        if ~strcmpi(options.motionFilter.method,'Off')
            plot(arfidata.trackTime(1:options.dispEst.dims(1)),squeeze(arfidata.disp_mf_pre(idx(:,i),i,1:options.dispEst.dims(1))),'--','Parent',dvt_ax2)
            plot(arfidata.trackTime(options.dispEst.dims(1)+1:end),squeeze(arfidata.disp_mf_push(idx(:,i),i,options.dispEst.dims(1)+1:end)),'--','Parent',dvt_ax2)
            plot(options.display.t_disp_pre*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax2)
            plot((options.display.t_disp_pre-dt)*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax2)
            plot(options.display.t_disp_push*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax2)
            plot(arfidata.trackTime(options.dispEst.dims(1))*ones(1,10),linspace(-650,650,10),'r-','linewidth',2,'Parent',dvt_ax2)
        end
        str = '';
        if ~isempty(strfind(options.motionFilter.method,'Poly'))
            str = sprintf('Order = %s',num2str(options.motionFilter.order));
            %             for j=1:length(options.motionFilter.timeRange_pre)
            %                 if options.motionFilter.timeRange_pre(j) < min(arfidata.trackTime)
            %                     plot(0.9*min(arfidata.trackTime),0.9*mf_rng(2),'<','Markersize',10,'MarkerFaceColor','y','Parent',dvt_ax2)
            %                 elseif options.motionFilter.timeRange_pre(j) > max(arfidata.trackTime)
            %                     plot(0.9*max(arfidata.trackTime),0.9*mf_rng(2),'>','Markersize',10,'MarkerFaceColor','y','Parent',dvt_ax2)
            %                 else
            %                     plot(options.motionFilter.timeRange_pre(j)*ones(1,10),linspace(-650,650,10),'y-','linewidth',1,'Parent',dvt_ax2)
            %                 end
            %             end
            %             for j=1:length(options.motionFilter.timeRange_push)
            %                 if options.motionFilter.timeRange_push(j) < min(arfidata.trackTime)
            %                     plot(0.9*min(arfidata.trackTime),0.9*mf_rng(2),'<','Markersize',10,'MarkerFaceColor','g','Parent',dvt_ax2)
            %                 elseif options.motionFilter.timeRange_push(j) > max(arfidata.trackTime)
            %                     plot(0.9*max(arfidata.trackTime),0.9*mf_rng(2),'>','Markersize',10,'MarkerFaceColor','g','Parent',dvt_ax2)
            %                 else
            %                     plot(options.motionFilter.timeRange_push(j)*ones(1,10),linspace(-650,650,10),'g-','linewidth',1,'Parent',dvt_ax2)
            %                 end
            %             end
        end
        grid(dvt_ax2,'on')
        xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        title(sprintf('MF ARFI Displacement Profiles (within Depth Gate)\nMotion Filter Type: %s\n%s',options.motionFilter.method,str),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',dvt_ax2);
        %         figure(101); set(101,'Color',dispPar.fig,'Position',[177 587 818 310])
        %         motion_ax = axes;
        %         set(motion_ax,'xlim',[arfidata.trackTime(1) arfidata.trackTime(end)],'ylim',raw_rng,'ColorOrder',dispPar.corder,'Color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold')
        %         hold(motion_ax,'on')
        %         plot(arfidata.trackTime,squeeze(arfidata.motion_push(idx(:,i),i,:)),'.--','Parent',motion_ax)
        %         grid on
        %         xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        %         ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        %         title('Estimated motion profiles','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        figure(100)
        set(100,'CurrentAxes',IQ_ax2)
        imagesc(arfidata.trackTime,arfidata.IQaxial,squeeze(IQ(:,i,:)),[-50 0]);colormap(IQ_ax2,gray);
        hold on
        for j=1:size(borders,1),plot(arfidata.trackTime,borders(j,i)*ones(1,size(arfidata.disp,3)),'color',dispPar.trace_cols(j,:),'LineWidth',2);end
        plot(min(arfidata.trackTime),par.pushFocalDepth,'>','Markersize',10,'MarkerFaceColor','c')
        plot(min(arfidata.trackTime),arfidata.axial(idx(1,i)),'<','Markersize',10,'MarkerFaceColor','r')
        plot(min(arfidata.trackTime),arfidata.axial(idx(end,i)),'<','Markersize',10,'MarkerFaceColor','r')
        hold off
        xlabel('Track Time (ms)')
        ylabel('Axial (mm)')
        title('IQ Data','color',dispPar.txt,'fontweight','bold')
        set(IQ_ax2,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold');
        if strcmpi(options.dispEst.ref_type,'anchored'), corr_rng = [0.75 1]; elseif strcmpi(options.dispEst.ref_type,'progressive'), corr_rng = [0.99 1]; end
        set(100,'CurrentAxes',corr_ax)
        cc_temp = squeeze(arfidata.cc(:,i,:));
        cc_SNR = cc_temp./(ones(size(cc_temp))-cc_temp);
        %         imagesc(arfidata.trackTime,arfidata.axial,squeeze(arfidata.cc(:,i,:)),corr_rng);colormap(corr_ax,jet);
        imagesc(arfidata.trackTime,arfidata.axial,db(cc_SNR),[30 60]);colormap(corr_ax,jet);
        hold on
        for j=1:size(borders,1),plot(arfidata.trackTime,borders(j,i)*ones(1,size(arfidata.disp,3)),'k','LineWidth',2);end
        hold off
        title('CC Coeff','color',dispPar.txt,'fontweight','bold')
        set(corr_ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold');
        cb1 = colorbar('Peer',corr_ax,'Location','SouthOutside');
        set(cb1,'Position',[0.1475 0.1 0.1 0.015],'Color',dispPar.txt,'FontWeight','bold');
        % draw borders on here
        raw = squeeze(arfidata.disp(:,i,:));
        filt = squeeze(mask(:,i,:)); alphaMask = ones(size(filt)); alphaMask(filt==0) = 0;
        raw(filt==0) = nan;
        set(100,'CurrentAxes',raw_ax)
        raw_im = imagesc(arfidata.trackTime,arfidata.axial,raw,raw_rng);colormap(raw_ax,parula);
        set(raw_im,'alphaData',alphaMask)
        hold on
        for j=1:size(borders,1),plot(arfidata.trackTime,borders(j,i)*ones(1,size(arfidata.disp,3)),'k','LineWidth',2);end
        hold off
        title('Raw Disp','color',dispPar.txt,'fontweight','bold')
        set(raw_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold');
        cb2 = colorbar('Peer',raw_ax,'Location','SouthOutside');
        set(cb2,'Position',[0.2825 0.1 0.1 0.015],'Color',dispPar.txt,'FontWeight','bold');
        % draw borders on here
        mf = nan(size(raw));
        if ~strcmpi(options.motionFilter.method,'Off')
            mf(:,1:options.dispEst.dims(1)) = squeeze(arfidata.disp_mf_pre(:,i,1:options.dispEst.dims(1)));
            mf(:,options.dispEst.dims(1)+1:end) = squeeze(arfidata.disp_mf_push(:,i,options.dispEst.dims(1)+1:end));
        end
        mf(squeeze(mask(:,i,:))==0) = nan;
        set(100,'CurrentAxes',mf_ax)
        mf_im = imagesc(arfidata.trackTime,arfidata.axial,mf,mf_rng);colormap(mf_ax,parula);
        set(mf_im,'alphaData',alphaMask)
        hold on
        for j=1:size(borders,1),plot(arfidata.trackTime,borders(j,i)*ones(1,size(arfidata.disp,3)),'k','LineWidth',2);end
        hold off
        title('MF Disp','color',dispPar.txt,'fontweight','bold')
        set(mf_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold');
        cb3 = colorbar('Peer',mf_ax,'Location','SouthOutside');
        set(cb3,'Position',[0.405 0.1 0.1 0.015],'Color',dispPar.txt,'FontWeight','bold');
        % draw borders on here
        %         figure(101)
        %         frame = getframe(100);
        %         im = frame2im(frame);
        %         [imind,cm] = rgb2ind(im,256);
        %         if i == 1;
        %             imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0);
        %         else
        %             imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
        %         end
        if ~skip
            pause
            %             keyboard
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
            i=i+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%         trace_ax =  axes('Position',[0.175 0.1 0.1 0.8]);
%         temp = squeeze(arfidata.IQ(:,1+(i-1)*par.nBeams,:));
%         I = real(temp); Q = imag(temp);
%         factor = 5;
%         D = size(I);
%         D(1) = D(1).*factor;
%         [Iup, Qup] = computeUpsampledIQdata(I,Q,factor);
%         Iup = reshape(Iup,D); Qup = reshape(Qup,D);
%         temp = db(abs(complex(Iup,Qup)));
%         temp(:,par.nref+1:par.nref+par.npush+par.nreverb) = nan;
%         axial_up = interp(arfidata.IQaxial,factor);
%         offset = 2.5;
%         for j=1:ntrackT; plot(axial_up,offset*(j-1)-temp(:,j)','b'); hold on; end;
%         view(90,90); hold on;
%         plot(borders(i,1)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
%         plot(borders(i,2)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
% %         xlim(edge);ylim([-100 300])
%         title(sprintf('Raw IQ: %d (t = %2.2f s)',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','color',dispPar.txt);
%         xlabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Tracks','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);set(gca,'YTickLabel',[])
%         set(ax21,'color',dispPar.ax + 0.25,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','Userdata','IQtraces_ax')
%         ax22= subplot(122); cla(ax22);
%         im1 = imagesc(arfidata.trackTime,arfidata.axial,double(squeeze(arfidata.cc(:,i,:))),[0.99 1]);cb = colorbar; set(cb,'Color',dispPar.txt,'FontWeight','bold','UserData','cc_cb')
%         title(sprintf('%s Correlation Coefficients',options.dispEst.ref_type),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);grid on;colormap(jet)
%         hold on;plot(linspace(-25,25,100),gate(i,1)*ones(1,100),'g','linewidth',3);plot(linspace(-25,25,100),gate(i,2)*ones(1,100),'g','linewidth',3)
%         xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold')
%         opengl software


% % Calculate alphaData mask
%         filt = squeeze(mask(:,i,:));
%         alphaMask = ones(size(filt)); alphaMask(filt==0) = 0.5;
%
%         fig2 = figure(2);
%         set(fig2,'Name','IQ Traces & Correlation Map','UserData','extra_fig1')
%         if isunix
%             set(fig2,'Position',[-1198 580 1198 893])
%         elseif ispc
%             set(fig2,'units','normalized','outerposition',[0 0 1 1])
%         end
%         set(fig2,'Color',dispPar.fig)
%
%         ax21=subplot(121); cla(ax21);
%         temp = squeeze(arfidata.IQ(:,1+(i-1)*par.nBeams,:));
%         I = real(temp); Q = imag(temp);
%         factor = 5;
%         D = size(I);
%         D(1) = D(1).*factor;
%         [Iup, Qup] = computeUpsampledIQdata(I,Q,factor);
%         Iup = reshape(Iup,D); Qup = reshape(Qup,D);
%         temp = db(abs(complex(Iup,Qup)));
%         temp(:,par.nref+1:par.nref+par.npush+par.nreverb) = nan;
%         axial_up = interp(arfidata.IQaxial,factor);
%         offset = 10;
%         for j=1:ntrackT; plot(axial_up,offset*(j-1)-temp(:,j)','b'); hold on; end;
%         view(90,90); hold on;
%         plot(borders(i,1)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
%         plot(borders(i,2)*ones(1,100),linspace(-offset*10,offset*ntrackT,100),'g','linewidth',3);
% %         xlim(edge);ylim([-100 300])
%         title(sprintf('Raw IQ: %d (t = %2.2f s)',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','color',dispPar.txt);
%         xlabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Tracks','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);set(gca,'YTickLabel',[])
%         set(ax21,'color',dispPar.ax + 0.25,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','Userdata','IQtraces_ax')
%         ax22= subplot(122); cla(ax22);
%         im1 = imagesc(arfidata.trackTime,arfidata.axial,double(squeeze(arfidata.cc(:,i,:))),[0.99 1]);cb = colorbar; set(cb,'Color',dispPar.txt,'FontWeight','bold','UserData','cc_cb')
%         title(sprintf('%s Correlation Coefficients',options.dispEst.ref_type),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);grid on;colormap(jet)
%         hold on;plot(linspace(-25,25,100),gate(i,1)*ones(1,100),'g','linewidth',3);plot(linspace(-25,25,100),gate(i,2)*ones(1,100),'g','linewidth',3)
%         xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold');ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold')
%         opengl software
%         set(im1,'alphaData',alphaMask)
%         set(ax22,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','ccmap_ax')
%         raw = squeeze(arfidata.disp(:,i,:));
%         mf = nan(size(raw));
%         if ~strcmpi(options.motionFilter.method,'Off')
%             mf(:,1:par.nref) = squeeze(arfidata.disp_mf_pre(:,i,1:par.nref));
%             mf(:,par.nref+1:end) = squeeze(arfidata.disp_mf_push(:,i,par.nref+1:end));
%         end
%         fig3 = figure(3);
%         set(fig3,'Name','Displacement Data','UserData','extra_fig2')
%         if isunix
%             set(fig3,'Position',[-1195 -270 1198 848])
%         elseif ispc
%             set(fig3,'units','normalized','outerposition',[0 0 1 1])
%         end
%         set(fig3,'Color',dispPar.fig)
%         ax31 = subplot(121); cla(ax31);
%         im2 = imagesc(arfidata.trackTime,arfidata.axial,raw,[-150 150]);cb2 = colorbar; set(cb2,'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb_extra')
%         xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         title(sprintf('Raw Displacement: Push %d\n Time = %1.2f s',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         hold on
%         plot(options.display.t_disp_pre*ones(1,length(arfidata.axial)),arfidata.axial,'y','linewidth',2)
%         plot(options.display.t_disp_push*ones(1,length(arfidata.axial)),arfidata.axial,'g','linewidth',2)
%         l1 = plot(arfidata.trackTime,gate(i,1)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
%         l2 = plot(arfidata.trackTime,gate(i,2)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
%         set(im2,'alphaData',alphaMask)
%         set(ax31,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','rawdisp_ax')
%         ax32 = subplot(122); cla(ax32)
%         im3 = imagesc(arfidata.trackTime,arfidata.axial,mf,rng);cb3 = colorbar; set(cb3,'Color',dispPar.txt,'FontWeight','bold')
%         xlabel('Track Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         title(sprintf('MF Displacement: Push %d\n Time = %1.2f s',i,arfidata.acqTime(i)),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
%         hold on
%         plot(options.display.t_disp_pre*ones(1,length(arfidata.axial)),arfidata.axial,'y','linewidth',2)
%         plot(options.display.t_disp_push*ones(1,length(arfidata.axial)),arfidata.axial,'g','linewidth',2)
%         l3 = plot(arfidata.trackTime,gate(i,1)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
%         l4 = plot(arfidata.trackTime,gate(i,2)*ones(1,length(arfidata.trackTime)),'g-','linewidth',2);
%         set(im3,'alphaData',alphaMask)
%         set(ax32,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','mfdisp_ax')
%     end
%     set(0,'CurrentFigure',fig)

%     if ~skip
%         w = waitforbuttonpress;
%         val = double(get(gcf,'CurrentCharacter'));
%         if val==28
%             i=i-1;
%         elseif val==29
%             i=i+1;
%         elseif val==32
%             skip = 1;
%         else
%             i=i+1;
%         end
%         if i<1;i=1;end
%     end
%     if skip
%         pause(0.05)
%         i=i+1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Pre and Push Traces averaged axially over the gate
[pre_trace,push_trace,layer_idx,dispTrace_ax] = drawDispTraces(fig,pre,push,borders,arfidata.axial,arfidata.acqTime,dispPar,options);
if ~isempty(arfidata.ecg)
    [pre_DR,push_DR] = computeRatios(ecg_ax,dispTrace_ax,pre_trace,push_trace,arfidata.acqTime,arfidata.ecg,arfidata.hr,samples,dispPar,options);
    for i=1:length(whos('rr_*'));eval(sprintf('delete(rr_%d)',i)), eval(sprintf('clear rr_%d',i)), end
    ratio_box = annotation('textbox',[0.05 0.07 0.4 0.025],'String','Ratios: ','Color','w','EdgeColor','w','FontWeight','Bold');
    for i=1:length(layer_idx)
        str = num2str(push_DR(layer_idx(i)),'%1.2f');
        eval(sprintf('rr_%d = annotation(''textbox'',[0.1+0.04*(i-1) 0.07 0.04 0.025],''LineStyle'',''none'',''Color'',dispPar.trace_cols(layer_idx(i),:),''FontWeight'',''bold'',''string'',str);',i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
keyboard
%% Panel of Axial Displacement and Axial Velocity Data vs. ECG
extra2_input = 'n';
if options.display.extras ~= -1
    extra2_input = input('\n>>>>> Do you want to look at temporally aligned panel? (y/n) [n]: ','s');
elseif options.display.extras >1
    extra2_input = 'y';
end
if strcmpi(extra2_input,'y')
    track_dt = median(diff(arfidata.trackTime))*1e-3;
    acq_dt = median(diff(arfidata.acqTime));
    tt = [0:track_dt:acq_dt+arfidata.acqTime(end)]; %ms
    raw = arfidata.disp;
    [raw_vel,t_vel] = differentiateData(arfidata.disp,arfidata.trackTime,options.motionFilter.Cutoff(2)); % mm/s
    raw = permute(raw,[1 3 2]);
    raw_vel = permute(raw_vel,[1 3 2]);
    for i=1:nacqT
        for j=1:size(borders,1)
            idx(j,i) = find(arfidata.axial>borders(j,i),1);
            raw(idx(j,i)-1:idx(j,i)+1,:,i) = inf;
            raw_vel(idx(j,i)-1:idx(j,i)+1,:,i) = inf;
        end
    end
    raw_panel = nan(size(arfidata.disp,1),length(tt));
    raw_vel_panel = nan(size(arfidata.disp,1),length(tt));
    shift_disp = zeros(size(raw,1),size(raw,2));
    shift_vel = zeros(size(raw_vel,1),size(raw_vel,2));
    start_idx = nan(1,nacqT);
    % Initialize
    disp_rng = [-650 650]; % miscrons
    vel_rng = [-10 10]; % mm/s
    start_idx(1) = 1;
    raw_panel(:,start_idx(1):start_idx(1)+size(raw,2)-1) = raw(:,:,1) + shift_disp;
    raw_vel_panel(:,start_idx(1):start_idx(1)+size(raw_vel,2)-1) = raw_vel(:,:,1) + shift_vel;
%     shift_disp = repmat(raw_panel(:,start_idx(1)+ntrackT-1),1,size(raw,2));
%     shift_vel = repmat(raw_vel_panel(:,start_idx(1)+ntrackT-2),1,size(raw_vel,2));
    for i=2:nacqT
        start_idx(i) = find(tt>arfidata.acqTime(i),1)-1;
        raw_panel(:,start_idx(i):start_idx(i)+size(raw,2)-1) = raw(:,:,i) + shift_disp;
%         shift_disp = repmat(raw_panel(:,start_idx(i)+ntrackT-1),1,size(raw,2));
        raw_vel_panel(:,start_idx(i):start_idx(i)+size(raw_vel,2)-1) = raw_vel(:,:,i) + shift_vel;
%         shift_vel = repmat(raw_vel_panel(:,start_idx(i)+ntrackT-2),1,size(raw_vel,2));
    end
    raw_panel = raw_panel(min(idx(:))-10:max(idx(:))+10,:);
    raw_vel_panel = raw_vel_panel(min(idx(:))-10:max(idx(:))+10,:);
    width = 1500; % ms
    delta = 50; % ms
    tt = tt*1000;
    win = [0 width];
    win_idx = [find(tt>win(1),1,'first'):find(tt>win(2),1,'first')];
    %     close all
    fig200 = figure(200);
    set(fig200,'units','normalized','outerposition',[0 0 1 1],'Color',dispPar.fig)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functionize this
    ax41 = axes('Position',[0.030 0.67 0.95 0.25]);
    pn1 = imagesc(tt(win_idx),arfidata.axial(min(idx(:))-10:max(idx(:))+10),raw_panel(:,win_idx),disp_rng);
    xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    set(ax41,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','disppanel_ax')
    set(pn1,'AlphaData',~isnan(raw_panel(:,win_idx)))
    %     ylim([min(borders(:))-1 max(borders(:))+1])
    colormap(jet)
    cb = colorbar;set(cb,'Color',dispPar.txt)
    ylabel(cb,'Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    ax42 = axes('Position',[0.030 0.07 0.95 0.25]);
    pn2 = imagesc(tt(win_idx),arfidata.axial(min(idx(:))-10:max(idx(:))+10),raw_vel_panel(:,win_idx),vel_rng);
    xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    set(ax42,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','velpanel_ax')
    set(pn2,'AlphaData',~isnan(raw_vel_panel(:,win_idx)))
    %     ylim([min(gate(:,1))-1 max(gate(:,2))+1])
    colormap(jet)
    cb = colorbar; set(cb,'Color',dispPar.txt)
    ylabel(cb,'Velocity (mm/s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
    if ~isempty(arfidata.ecg)
        ax43 = axes('Position',[0.030 0.40 0.91 0.20]);
        plot(1000*arfidata.ecg(:,1),arfidata.ecg(:,2),'linewidth',5);
        xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
        set(ax43,'color',[0.15 0.15 0.15],'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'yTickLabel',[],'xgrid','on','UserData','ecgpanel_ax')
        xlim(win);ylim([-0.5 1])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functionize this
    val = 0;
    while val~=32
        if val==28
            win = [win(1)-delta win(2)-delta];
            if win(1)<0;win = [0 width];end
            if win(2)>arfidata.acqTime*1000;win = [arfidata.acqTime*1000-width arfidata.acqTime*1000];end
            win_idx = [find(tt>win(1),1,'first'):find(tt>win(2),1,'first')];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functionize this
            ax41 = axes('Position',[0.030 0.67 0.95 0.25]);
            pn1 = imagesc(tt(win_idx),arfidata.axial(min(idx(:))-10:max(idx(:))+10),raw_panel(:,win_idx),disp_rng);
            xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            set(ax41,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','disppanel_ax')
            set(pn1,'AlphaData',~isnan(raw_panel(:,win_idx)))
            %             ylim([min((:,1))-1 max(gate(:,2))+1])
            colormap(jet)
            cb = colorbar; set(cb,'Color',dispPar.txt)
            ylabel(cb,'Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            ax42 = axes('Position',[0.030 0.07 0.95 0.25]);
            pn2 = imagesc(tt(win_idx),arfidata.axial(min(idx(:))-10:max(idx(:))+10),raw_vel_panel(:,win_idx),vel_rng);
            xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            set(ax42,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','velpanel_ax')
            set(pn2,'AlphaData',~isnan(raw_vel_panel(:,win_idx)))
            %             ylim([min(gate(:,1))-2.5 max(gate(:,2))+2.5])
            colormap(jet)
            cb = colorbar; set(cb,'Color',dispPar.txt)
            ylabel(cb,'Velocity (mm/s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            if ~isempty(arfidata.ecg)
                ax43 = axes('Position',[0.030 0.40 0.90 0.20]);
                plot(1000*arfidata.ecg(:,1),arfidata.ecg(:,2),'linewidth',5);
                xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
                set(ax43,'color',[0.15 0.15 0.15],'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'yTickLabel',[],'xgrid','on','UserData','ecgpanel_ax')
                xlim(win);ylim([-0.5 1])
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functionize this
        end
        if val==29
            win = [win(1)+delta win(2)+delta];
            if win(1)<0;win = [0 width];end
            if win(2)>arfidata.acqTime(end)*1000;win = [arfidata.acqTime(end)*1000-width arfidata.acqTime(end)*1000];end
            win_idx = [find(tt>win(1),1,'first'):find(tt>win(2),1,'first')];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functionize this
            ax41 = axes('Position',[0.030 0.67 0.95 0.25]);
            pn1 = imagesc(tt(win_idx),arfidata.axial(min(idx(:))-10:max(idx(:))+10),raw_panel(:,win_idx),disp_rng);
            xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            set(ax41,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','disppanel_ax')
            set(pn1,'AlphaData',~isnan(raw_panel(:,win_idx)))
            %             ylim([min(gate(:,1))-1 max(gate(:,2))+1])
            colormap(jet)
            cb = colorbar; set(cb,'Color',dispPar.txt)
            ylabel(cb,'Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            ax42 = axes('Position',[0.030 0.07 0.95 0.25]);
            pn2 = imagesc(tt(win_idx),arfidata.axial(min(idx(:))-10:max(idx(:))+10),raw_vel_panel(:,win_idx),vel_rng);
            xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            set(ax42,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'xgrid','on','UserData','velpanel_ax')
            set(pn2,'AlphaData',~isnan(raw_vel_panel(:,win_idx)))
            %             ylim([min(gate(:,1))-1 max(gate(:,2))+1])
            colormap(jet)
            cb = colorbar; set(cb,'Color',dispPar.txt)
            ylabel(cb,'Velocity (mm/s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
            if ~isempty(arfidata.ecg)
                ax43 = axes('Position',[0.030 0.40 0.90 0.20]);
                plot(1000*arfidata.ecg(:,1),arfidata.ecg(:,2),'linewidth',5);
                xlabel('Time (ms)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
                set(ax43,'color',[0.15 0.15 0.15],'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','TickLength',[0 0.025],'yTickLabel',[],'xgrid','on','UserData','ecgpanel_ax')
                xlim(win);ylim([-0.5 1])
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functionize this
        end
        w = waitforbuttonpress;
        val = double(get(gcf,'CurrentCharacter'));
%         delete(ax41)
%         delete(ax42)
%         delete(ax43)
    end
%     close 200
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
keyboard