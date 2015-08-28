function [samples,ecg_ax,pt] = drawECG(fig,ecg,acqTime,dispPar,varargin)

if ~isempty(varargin)
    pos = varargin{1};
else
    pos = [0.5 0.1 0.4 0.2];
end

ch = get(fig,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'ecg_ax'))
        delete(ch(i))
    end
end
samples = zeros(1,length(acqTime));
for i=1:length(acqTime)
    samples(i) = ecg(find(ecg(:,1)>acqTime(i),1,'first')-1,2);
end
% fig = figure(101);set(fig,'Color','k'); ecg_ax = axes;
ecg_ax = axes('Position',pos,'Parent',fig);
plot(ecg(:,1),ecg(:,2),'Linewidth',2);
hold(ecg_ax,'on')
plot(acqTime,samples,'gx','MarkerSize',8)
pt = plot(acqTime(1),samples(1),'ro','Markersize',10,'Markerfacecolor','r','Parent',ecg_ax);
title('ECG (ARFI)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
xlim([0 max(acqTime)])
ylim([min(ecg(:,2)) max(ecg(:,2))])
set(ecg_ax,'color',dispPar.ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'yTickLabel',[],'fontweight','bold','xgrid','on','gridLineStyle','--','Userdata','ecg_ax')
