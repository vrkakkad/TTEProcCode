function [pre_trace,push_trace,layer_idx,dispTrace_ax] = drawDispTraces(fig,pre,push,borders,axial,acqTime,dispPar,options)

ch = get(fig,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'dispTrace_ax'))
        delete(ch(i))
    end
end

pre_trace = nan(size(borders,1)-1,length(acqTime));
push_trace = nan(size(pre_trace));
idx = nan(options.display.n_pts,length(acqTime),size(borders,1)-1);
for i=1:size(borders,1)-1
    gate = [borders(i,:);borders(i+1,:)];
    [pre_trace(i,:),push_trace(i,:),idx(:,:,i)] = computeDispTrace(pre,push,axial,gate,options.display.n_pts);
end

if size(borders,1)>2
    layer_input = input(sprintf('Choose layer numbers to display axially avg. displacement traces\n\tShallowest(1) to Deepest(%d): ',size(borders,1)-1));
    if ~isempty(layer_input), layer_idx = layer_input; else layer_idx = [1:size(borders,1)-1];end
else
    layer_idx = 1;
end

dispTrace_ax = axes('Position',[0.05 0.15 0.4 0.25],'Color',dispPar.fig,'Parent',fig);
hold(dispTrace_ax,'on')
for i=1:length(layer_idx)
    plot(acqTime,pre_trace(layer_idx(i),:),'o--','Color',dispPar.trace_cols(layer_idx(i),:),'linewidth',1,'Parent',dispTrace_ax);
    plot(acqTime,push_trace(layer_idx(i),:),'*-','Color',dispPar.trace_cols(layer_idx(i),:),'linewidth',1,'Parent',dispTrace_ax);
end

grid on
xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',dispTrace_ax);
ylabel('Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',dispTrace_ax);
title(sprintf('Median Tracked Displacements\n(within Depth Gate)'),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',dispTrace_ax)
ylim(options.display.dispRange)
xlim([0 max(acqTime)])
set(dispTrace_ax,'Color',dispPar.ax,'ColorOrder',dispPar.corder,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','dispTrace_ax');

% Auto Range
if options.display.autoRange
    if min(pre_trace(:))>=0;low = 0.75*min(pre_trace(:));else low = 1.25*min(pre_trace(:));end
    high = 1.25*max(push_trace(:));
    setRange([low high])
end
