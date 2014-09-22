function dispRes(DataDir,fidx,varargin)

if ~exist('DataDir','var')
    DataDir = pwd;
end

if ~exist('fidx','var')
    fidx = -1;
end

cd(DataDir)

%% Extract timeStamp
if ispc
    addpath C:\users\vrk4\Documents\GitHub\SC2000\arfiProcCode\
    addpath(genpath('C:\users\vrk4\Documents\GitHub\TTEProcCode')) 
elseif isunix
    addpath /emfd/vrk4/GitHub/SC2000/arfiProcCode
    addpath(genpath('/emfd/vrk4/GitHub/TTEProcCode'))
end

list = dir('res_*'); % get timeStamp based on existance of ARFI par files

if size(list,1)<fidx
    error('Data set index requested greater than number of data sets')
end

% Reading in timestamp for data set
if length(fidx)==14
    timeStamp = fidx;
elseif fidx == -1
    timeStamp = list(end).name(end-17:end-4);
else
    timeStamp = list(fidx).name(end-17:end-4);
end

fprintf('Loading data with timeStamp = %s\n', timeStamp);

% Loading Res file
tic
load(strcat('res_',timeStamp));
par = load(strcat('arfi_par_',timeStamp,'.mat'));
par.isHarmonic
fprintf(1,'Load Time = %2.2fs\n',toc)

keyboard

%% Input Parameters
options.dispARFIParams = struct(...
    'BRange',[0 40] ...
    ,'ARange',[0 5] ...
    );


dof = 7.22*1.540/par.pushFreq*(par.pushFnum)^2;
edge = (par.pushFocalDepth + [-dof/4 dof/4]);

%% Display Bmode

figure(1)
% set(1,'Position',[0 480 1440 335]);
for i=1:size(bdata.bimg,3)
    imagesc(bdata.blat*10,bdata.bax*10,bdata.bimg(:,:,i));
    xlabel('Lateral (mm)')
    ylabel('Axial (mm)')
    title(sprintf('B-Mode: Frame %d (t = %1.1f s)\n',i,bdata.t(i)))
    colormap(gray)
    axis image
    hold on
    plot(bdata.blat*10,edge(1)*ones(length(bdata.blat)),'b')
    plot(bdata.blat*10,edge(2)*ones(length(bdata.blat)),'b')
    plot(-5*ones(length(bdata.bax)),bdata.bax*10,'b')
    plot(5*ones(length(bdata.bax)),bdata.bax*10,'b')
    hold off
    pause
end

%% Display ARFI

figure(2)
% set(1,'Position',[0 480 1440 335]);
for i=1:size(arfidata.disp_off,3)
    subplot(221)
    imagesc(arfidata.acqTime,arfidata.axial,abs(db(arfidata.IQ_off(:,:,i))))
    hold on
    plot(arfidata.acqTime,edge(1)*ones(length(arfidata.acqTime)),'b')
    plot(arfidata.acqTime,edge(2)*ones(length(arfidata.acqTime)),'b')
    hold off
    xlabel('Acquisition Time (s)')
    ylabel('Axial (mm)')
    title(sprintf('IQ Data (no push)\n Frame %d/ t = %1.2f ms\n',i,arfidata.trackTime(i)))
    colormap(gray)
    
    subplot(222)
    imagesc(arfidata.acqTime,arfidata.axial,abs(db(arfidata.IQ_on(:,:,i))))
    hold on
    plot(arfidata.acqTime,edge(1)*ones(length(arfidata.acqTime)),'b')
    plot(arfidata.acqTime,edge(2)*ones(length(arfidata.acqTime)),'b')
    hold off
    xlabel('Acquisition Time (s)')
    ylabel('Axial (mm)')
    title(sprintf('IQ Data (push)\n Frame %d/ t = %1.2f ms\n',i,arfidata.trackTime(i)))
    colormap(gray)
    
    subplot(223)
    imagesc(arfidata.acqTime,arfidata.axial,arfidata.disp_off(:,:,i),options.dispARFIParams.ARange)
    hold on
    plot(arfidata.acqTime,edge(1)*ones(length(arfidata.acqTime)),'b')
    plot(arfidata.acqTime,edge(2)*ones(length(arfidata.acqTime)),'b')
    hold off
    xlabel('Acquisition Time (s)')
    ylabel('Axial (mm)')
    title(sprintf('ARFI Data (no push)\n Frame %d/ t = %1.2f ms\n',i,arfidata.trackTime(i)))
    colormap(gray)
    
    subplot(224)
    imagesc(arfidata.acqTime,arfidata.axial,arfidata.disp_on(:,:,i),options.dispARFIParams.ARange)
    hold on
    plot(arfidata.acqTime,edge(1)*ones(length(arfidata.acqTime)),'b')
    plot(arfidata.acqTime,edge(2)*ones(length(arfidata.acqTime)),'b')
    hold off
    xlabel('Acquisition Time (s)')
    ylabel('Axial (mm)')
    title(sprintf('ARFI Data (push)\n Frame %d/ t = %1.2f ms\n',i,arfidata.trackTime(i)))
    colormap(gray)
    
    pause
end

%% 
edge(1) = input('Edge 1 = ');
edge(2) = input('Edge 2 = ');

edge_idx = [find(arfidata.axial>edge(1),1) find(arfidata.axial>edge(2),1)];

figure(3)
% set(1,'Position',[0 480 1440 335]);
subplot(211)
imagesc(arfidata.acqTime, arfidata.trackTime, squeeze(mean(arfidata.disp_off(edge_idx(1):edge_idx(2),:,:),1))',options.dispARFIParams.ARange);
xlabel('Acquisition Time (s)')
ylabel('Track Time (ms)')
title(sprintf('Mean Displacements: No push\n%1.2f - %1.2f mm',edge(1),edge(2)))

subplot(212)
imagesc(arfidata.acqTime, arfidata.trackTime, squeeze(mean(arfidata.disp_on(edge_idx(1):edge_idx(2),:,:),1))',options.dispARFIParams.ARange);
xlabel('Acquisition Time (s)')
ylabel('Track Time (ms)')
title(sprintf('Mean Displacements: Push\n%1.2f - %1.2f mm',edge(1),edge(2)))


%% Display SWEI

foc_idx = find(sweidata.axial>par.pushFocalDepth,1);

sw_off = squeeze(mean(sweidata.disp_off(edge_idx(1):edge_idx(2),:,:,:),1));
sw_off = permute(sw_off,[1 3 2]);
vel_off = 1000*diff(sw_off,2,1)./par.priusec(1);

sw_on = squeeze(mean(sweidata.disp_on(edge_idx(1):edge_idx(2),:,:,:),1));
sw_on = permute(sw_on,[1 3 2]);
vel_on = 1000*diff(sw_on,2,1)./par.priusec(1);


figure(4)
% set(1,'Position',[0 480 1440 335]);

for i=1:size(sw_off,3)
    subplot(221)
    imagesc(sweidata.lat(foc_idx,:),sweidata.trackTime,sw_off(:,:,i)');
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('Mean SW Displacements: No push\n%1.2f - %1.2f mm\nt = %1.2f s',edge(1),edge(2),sweidata.acqTime(i)))
    
    subplot(222)
    imagesc(sweidata.lat(foc_idx,:),sweidata.trackTime,sw_on(:,:,i)');
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('Mean SW Displacements: No push\n%1.2f - %1.2f mm\nt = %1.2f s',edge(1),edge(2),sweidata.acqTime(i)))
    
    subplot(223)
    imagesc(sweidata.lat(foc_idx,:),sweidata.trackTime,vel_off(:,:,i)');
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('Mean SW Velocity: No push\n%1.2f - %1.2f mm\nt = %1.2f s',edge(1),edge(2),sweidata.acqTime(i)))
    
    subplot(224)
    imagesc(sweidata.lat(foc_idx,:),sweidata.trackTime,vel_on(:,:,i)');
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('Mean SW Velocity: No push\n%1.2f - %1.2f mm\nt = %1.2f s',edge(1),edge(2),sweidata.acqTime(i)))
    pause
end




end

