clear all; close all; clc

set_idx = [2 3 4 5 6 7 8 9 ];

list = dir('align_arfi*');

figure
ax1 = subplot(211);
ax2 = subplot(212);

for i=1:length(set_idx) 
    load(list(i).name);
    
    set(gcf,'CurrentAxes',ax1)
    plot(acqTime,pre_seg,'-yo','linewidth',2,'MarkerFaceColor','k');
    hold(ax1,'on')
    plot(acqTime,push_seg,'-ro','linewidth',2,'MarkerFaceColor','k');
    
    set(gcf,'CurrentAxes',ax2)
    plot(ecg(:,1),ecg(:,2),'Linewidth',2);
    hold(ax2,'on')
    plot(acqTime(idx_start:end)-t_start,samples(idx_start:end),'gx','MarkerSize',5)
end
set(gcf,'CurrentAxes',ax1)
set(gcf,'color',dispPar.fig)
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt);
title(sprintf('Axially Averaged ARFI Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
ylim([-2 10])
xlim([0 max(acqTime)])
grid on
set(ax1,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','res_ax');

set(gcf,'CurrentAxes',ax2)
xlim([0 max(arfidata.acqTime)])
ylim([min(ecgdata.arfi(:,2)) max(ecgdata.arfi(:,2))])
title('ECG Trace','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
set(ax2,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--','Userdata','ecg_ax')