clear
close all
clc

setidx = [1 1];

list = dir('res_*');

tic
fprintf(1,'Loading Data...\n')

fund = load(list(setidx(1)).name);
fund_ts = list(setidx(1)).name(end-17:end-4)
fund.apar = load(strcat('arfi_par_',fund_ts));

harm = load(list(setidx(2)).name);
harm_ts = list(setidx(2)).name(end-17:end-4)
harm.apar = load(strcat('arfi_par_',harm_ts));

fprintf(1,'Load Time = %2.2fs\n',toc)

clear fund_ts harm_ts
%% Display Bmode

figure(1)
set(1,'Position',[0 480 1440 335]);
for i=1:size(fund.bdata.bimg,3)
    subplot(121)
    imagesc(fund.bdata.blat*10,fund.bdata.bax*10,fund.bdata.bimg(:,:,i));
    xlabel('Lateral (mm)')
    ylabel('Axial (mm)')
    title(sprintf('B-Mode: Frame %d (t = %1.1f s)\n(Fundamental)',i,fund.bdata.t(i)))
    colormap(gray)
    axis image
    grid on
    
    subplot(122)
    imagesc(harm.bdata.blat*10,harm.bdata.bax*10,harm.bdata.bimg(:,:,i));
    xlabel('Lateral (mm)')
    ylabel('Axial (mm)')
    title(sprintf('B-Mode: Frame %d (t = %1.1f s)\n(Harmonic)',i,harm.bdata.t(i)))
    colormap(gray)
    axis image
    grid on
    
    pause(0.025)
end

% %% ARFI Movie
% arange = [0 10];
% 
% figure(2)
% set(2,'Position',[1 41 1440 783])
% for i=1:size(fund.arfidata.disp_off,3)
%     subplot(221)
%     imagesc(fund.arfidata.acqTime,fund.arfidata.axial,fund.arfidata.disp_off(:,:,i),arange)
%     xlabel('Acquisition Time (s)')
%     ylabel('Axial (mm)')
%     title(sprintf('"No push" M-Mode Displacements\nFrame %d (t = %1.2f ms)\n Fundamental',i,fund.arfidata.trackTime(i)))
%     colormap(gray)
%     
%     subplot(222)
%     imagesc(fund.arfidata.acqTime,fund.arfidata.axial,fund.arfidata.disp_on(:,:,i),arange)
%     xlabel('Acquisition Time (s)')
%     ylabel('Axial (mm)')
%     title(sprintf('M-Mode ARFI Displacements\nFrame %d (t = %1.2f ms)\n Fundamental',i,fund.arfidata.trackTime(i)))
%     colormap(gray)
%     
%     subplot(223)
%     imagesc(harm.arfidata.acqTime,harm.arfidata.axial,harm.arfidata.disp_off(:,:,i),arange)
%     xlabel('Acquisition Time (s)')
%     ylabel('Axial (mm)')
%     title(sprintf('"No push" M-Mode Displacements\nFrame %d (t = %1.2f ms)\n Harmonic',i,harm.arfidata.trackTime(i)))
%     colormap(gray)
%     
%     subplot(224)
%     imagesc(harm.arfidata.acqTime,harm.arfidata.axial,harm.arfidata.disp_on(:,:,i),arange)
%     xlabel('Acquisition Time (s)')
%     ylabel('Axial (mm)')
%     title(sprintf('M-Mode ARFI Displacements\nFrame %d (t = %1.2f ms)\n Harmonic',i,harm.arfidata.trackTime(i)))
%     colormap(gray)
%     
%     pause
% end

% %% SWEI Movies
% 
% disp_off = reshape(sweidata.disp_off,[],15*100,50);
% disp_on = reshape(sweidata.disp_on,[],15*100,50);
% 
% 
% figure(3)
% set(3,'units','normalized','position',[0 0 1 1]) 
% for i=1:50
%     subplot(211)
%     imagesc(sweidata.acqTime*1000,sweidata.axial,disp_off(:,:,i),[0 10])
%     xlabel('Acquisition Time (ms)')
%     ylabel('Axial (mm)')
%     title(sprintf('"No push" M-Mode Displacements\nFrame %d (t = %1.2f ms)',i,sweidata.trackTime(i)))
%     colormap(gray)
%     
%     subplot(212)
%     imagesc(sweidata.acqTime*1000,sweidata.axial,disp_on(:,:,i),[0 10])
%     xlabel('Acquisition Time (ms)')
%     ylabel('Axial (mm)')
%     title(sprintf('M-Mode SWEI Displacements\nFrame %d (t = %1.2f ms)',i,sweidata.trackTime(i)))
%     colormap(gray)
%     
%     pause
% end
%% At Focus

% ARFI
arange = [0 10];

ax_gate = fund.apar.pushFocalDepth + [-2.5 2.5];
idx_gate_fund = [find(fund.arfidata.axial>ax_gate(1),1):find(fund.arfidata.axial>ax_gate(2),1)];
idx_gate_harm = [find(harm.arfidata.axial>ax_gate(1),1):find(harm.arfidata.axial>ax_gate(2),1)];

fund.atfoc_off = squeeze(mean(fund.arfidata.disp_off(idx_gate_fund,:,:),1))';
fund.atfoc_on = squeeze(mean(fund.arfidata.disp_on(idx_gate_fund,:,:),1))';

harm.atfoc_off = squeeze(mean(harm.arfidata.disp_off(idx_gate_harm,:,:),1))';
harm.atfoc_on = squeeze(mean(harm.arfidata.disp_on(idx_gate_harm,:,:),1))';

figure(2)
set(2,'Position',[0 50 1440 335])
subplot(221)
imagesc(fund.arfidata.acqTime,fund.arfidata.trackTime,fund.atfoc_off,arange)
xlabel('Acquisition Time (s)')
ylabel('Track Time (ms)')
title(sprintf('"No push" Displacements at %d mm (Fundamental)',mean(ax_gate)))
colorbar
subplot(222)
imagesc(fund.arfidata.acqTime,fund.arfidata.trackTime,fund.atfoc_on,arange)
xlabel('Acquisition Time (s)')
ylabel('Track Time (ms)')
title(sprintf('ARFI Displacements at %d mm (Fundamental)',mean(ax_gate)))
colorbar

subplot(223)
imagesc(harm.arfidata.acqTime,harm.arfidata.trackTime,harm.atfoc_off,arange)
xlabel('Acquisition Time (s)')
ylabel('Track Time (ms)')
title(sprintf('"No push" Displacements at %d mm (Harmonic)',mean(ax_gate)))
colorbar
subplot(224)
imagesc(harm.arfidata.acqTime,harm.arfidata.trackTime,harm.atfoc_on,arange)
xlabel('Acquisition Time (s)')
ylabel('Track Time (ms)')
title(sprintf('ARFI Displacements at %d mm (Harmonic)',mean(ax_gate)))
colorbar

% SWEI
fund_lat = fund.sweidata.lat(ceil(median(idx_gate_fund)),:);
harm_lat = harm.sweidata.lat(ceil(median(idx_gate_harm)),:);

fund_temp_off = squeeze(mean(fund.sweidata.disp_off(idx_gate_fund,:,:,:)));
fund_temp_off = permute(fund_temp_off,[3 1 2]);
fund_temp_on = squeeze(mean(fund.sweidata.disp_on(idx_gate_fund,:,:,:)));
fund_temp_on = permute(fund_temp_on,[3 1 2]);
harm_temp_off = squeeze(mean(harm.sweidata.disp_off(idx_gate_harm,:,:,:)));
harm_temp_off = permute(harm_temp_off,[3 1 2]);
harm_temp_on = squeeze(mean(harm.sweidata.disp_on(idx_gate_harm,:,:,:)));
harm_temp_on = permute(harm_temp_on,[3 1 2]);

pause

figure(3)
set(3,'Position',[0 50 1440 335])
for i=1:size(fund_temp_off,3)
    subplot(221)
    imagesc(fund_lat,fund.sweidata.trackTime,fund_temp_off(:,:,i),arange)
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('"No push" Displacements at %d mm (Fundamental: %d)',mean(ax_gate),i))
    subplot(222)
    imagesc(fund_lat,fund.sweidata.trackTime,fund_temp_on(:,:,i),arange)
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('SWEI Displacements at %d mm (Fundamental: %d)',mean(ax_gate),i))
    subplot(223)
    imagesc(harm_lat,harm.sweidata.trackTime,harm_temp_off(:,:,i),arange)
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('"No push" Displacements at %d mm (Harmonic: %d)',mean(ax_gate),i))
    subplot(224)
    imagesc(harm_lat,harm.sweidata.trackTime,harm_temp_on(:,:,i),arange)
    xlabel('Lateral (mm)')
    ylabel('Track Time (ms)')
    title(sprintf('SWEI Displacements at %d mm (Harmonic: %d)',mean(ax_gate),i))
    
    pause
end

%%
keyboard

    