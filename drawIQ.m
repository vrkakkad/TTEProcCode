function [IQBors,IQ_ax] = drawIQ(fig,IQ,axial,acqTime,edge,borders,dispPar,par,dynrange,df_flag,trackPRF)

ch = get(fig,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'IQ_ax'))
        delete(ch(i))
    end
end
%%
IQ_ax = axes('Position',[0.5 0.7 0.40 0.2],'Parent',fig);
% fig = figure(101);set(101,'Color','k');IQ_ax = axes;
frame = abs(IQ(:,:,1)); % Display first frame only
temp = IQ(:,:,par.nref-1);
frame = db(frame/max(temp(:)));
clear temp
temp = find(axial>edge(1)-1,1);
if isempty(temp); idx(1) = 1; else idx(1) = temp; end; clear temp
temp = find(axial>edge(2)+1,1);
if isempty(temp); idx(2) = length(axial); else idx(2) = temp; end; clear temp
imagesc(linspace(0,acqTime(end),size(IQ,2)),axial(idx(1):idx(2)),frame(idx(1):idx(2),:),dynrange)
clear idx
hold(IQ_ax,'on')
l1 = plot(linspace(0,acqTime(end),length(acqTime)),axial(1)*ones(1,length(acqTime)),'b','Linewidth',2,'Parent',IQ_ax);
l2 = plot(linspace(0,acqTime(end),length(acqTime)),axial(end)*ones(1,length(acqTime)),'b','Linewidth',2,'Parent',IQ_ax);
for i=1:size(borders,1)
    IQBors{i} = plot(linspace(0,acqTime(end),length(acqTime)),borders(i,:),'Color',dispPar.trace_cols(i,:),'Linewidth',2,'Parent',IQ_ax);
%     IQBors{i} = plot(linspace(0,acqTime(end),length(acqTime)),borders(i,:),'Color','g','Linewidth',2,'Parent',IQ_ax);
%     eval(sprintf('IQBor_%d = plot(linspace(0,acqTime(end),length(acqTime)),borders(i,:),''Color'',%s,''Linewidth'',1,''Parent'',IQ_ax);',i,mat2str(dispPar.trace_cols(i,:))))
%     eval(sprintf('IQBors{%d} = IQBor_%d;',i,i));
end
plot(0,par.pushFocalDepth,'>','Markersize',10,'MarkerFaceColor','c')
hold(IQ_ax,'off')
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
% title('M-Mode Frames','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
if ~df_flag
    title(sprintf('M-Mode Frames\nFoc Harmonic Tracking = %d\n Track PRF = %1.2f kHz, Push = %d x %d \\mus',par.isHarmonic,trackPRF,par.npush,par.pushDurationusec),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
else
    title(sprintf('M-Mode Frames\nDe-Foc Harmonic Tracking = %d\n Track PRF = %1.2f kHz, Push = %d x %d \\mus',par.isHarmonic,trackPRF,par.npush,par.pushDurationusec),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end
colormap(gray); freezeColors;
set(IQ_ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','off','gridLineStyle','--','UserData','IQ_ax')
