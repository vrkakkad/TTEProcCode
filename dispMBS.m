clear all; clc

sub_id = 'CS_05';

src = strcat('E:\TTE_Data\',sub_id);
cd(src)

%PSL
% set_idx = [2 3 1 8 9]; % CS_02
set_idx = [2:8]; % CS_05
% set_idx = [2:10]; % CS_07
% set_idx = [1:6]; % CS_10
% set_idx = [1:6 8]; % CS_11
% set_idx = [3:9]; % CS_12

%PSS
% set_idx = [12 13]; % CS_02
% set_idx = [10:17]; % CS_05
% set_idx = [11:15 17:22]; % CS_07
% set_idx = [23:27]; % CS_11
% set_idx = [14 15 18:20 22 23]; % CS_12

% Apical
% set_idx = [17 18 20:24]; % CS_12
% set_idx = [28:30 32 34 35]; % CS_12

list = dir('arfi_par_*');

cd MBS

figure
ax1 = subplot(211);
ax2 = subplot(212);

for i=1:length(set_idx)
    timeStamp = list(set_idx(i)).name(end-17:end-4);
    if exist(strcat('mbs_arfi_',timeStamp,'.mat'),'file')
        load(strcat('mbs_arfi_',timeStamp));
    else
        error(sprintf('File for Set Idx = %d (%s) does not exist',set_idx(i),timeStamp));
    end
    heartrate(i) = hr;
    set(gcf,'CurrentAxes',ax1)
    plot(acqTime,pre_trace,'-co','linewidth',1,'MarkerFaceColor','k','MarkerSize',4);
    hold(ax1,'on')
    grid on
    plot(acqTime,push_trace,'-ro','linewidth',1,'MarkerFaceColor','k','MarkerSize',4);
    
    set(gcf,'CurrentAxes',ax2)
    plot(ecg(:,1),ecg(:,2),'b','Linewidth',1);
    hold(ax2,'on')
    plot(acqTime,samples,'ro','MarkerSize',4)
    grid on
    
%     [a,b] = ginput;
%     [max(mean(b(1:2:end)),mean(b(2:2:end))),min(mean(b(1:2:end)),mean(b(2:2:end)))]
%     [max(mean(b(2:2:end))/mean(b(1:2:end)),mean(b(1:2:end))/mean(b(2:2:end)))]
%     
%     pause
%     clear a b
%     clc
%     cla(ax1); cla(ax2)
end

set(gcf,'CurrentAxes',ax1)
set(gcf,'color',[1 1 1])
set(gcf,'Position',[100 100 800 400])
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'Color',[0 0 0]);
title(sprintf('ARFI Displacements: Axially Averaged over Width of IV Septum (after motion filtering to remove cardiac motion)\n7 Repeated Acquisitions (Parasternal Long Axis View)'),'fontsize',dispPar.fsize,'fontname','trebuchet ms','fontweight','bold','Color',[0 0 0])
grid on
set(ax1,'Color',[1 1 1],'ColorOrder',dispPar.corder,'xcolor',[0 0 0],'ycolor',[0 0 0],'fontsize',6,'UserData','res_ax','xlim',[-0.75 2.5],'ylim',[-3 7],'Position',[0.13 0.51 0.775 0.34],'xticklabel',[]);
leg = legend('Passive Tracking of Myocardium (No ARF Push)',' Myocardial Response to ARF Excitation','location','SouthEast');
set(leg,'fontsize',6)

set(gcf,'CurrentAxes',ax2)
axis tight
ylabel(sprintf('ECG Trace'),'fontsize',dispPar.fsize,'Color',[0 0 0])
% set(get(gca,'YLabel'),'Rotation',0)
% temp = get(gca,'yLabel');
xlabel(sprintf('Acquisition Time\n(Normalized by Heartrate (%2.0f bpm))',hr),'fontsize',dispPar.fsize,'Color',[0 0 0])
set(ax2,'color',[1 1 1],'xcolor',[0 0 0],'ycolor',[0 0 0],'yTickLabel',[],'fontsize',6,'xgrid','on','Userdata','ecg_ax','xlim',[-0.75 2.5],'Position',[0.13 0.15 0.775 0.34])
% set(temp,'Position',get(temp,'Position') - [0.2 0 0])


clear
cd ../..