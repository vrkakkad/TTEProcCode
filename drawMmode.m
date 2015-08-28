function [pre,push,dt,options,preBors,pushBors,pre_ax,push_ax] = drawMmode(fig,arfidata,edge,borders,dispPar,par,options)

ch = get(fig,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'pre_ax') || strcmpi(nm,'push_ax') || strcmpi(nm,'cp_pre_ax') || strcmpi(nm,'cp_push_ax') || strcmpi(nm,'disp_cb'))
        delete(ch(i))
    end
end
%%
if options.display.cc_filt
    mask = arfidata.cc>options.display.cc_thresh;
else
    mask = [];
end

% Indices corresponding to median filter parameters
nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));

if (isempty(arfidata.disp_mf_pre) && isempty(arfidata.disp_mf_push))
    disp_pre = arfidata.disp;
    disp_push = arfidata.disp;
else
    disp_pre = arfidata.disp_mf_pre;
    disp_push = arfidata.disp_mf_push;
end

if isempty(options.display.t_disp_push)
    idx_push = sum(options.dispEst.dims(1:3))+1;
    options.display.t_disp_push = arfidata.trackTime(idx_push);
else
    idx_push = find(arfidata.trackTime>options.display.t_disp_push,1)-1;
end

if isempty(options.display.t_disp_pre)
    options.display.t_disp_pre = options.motionFilter.timeRange_pre(1) + (options.display.t_disp_push - options.motionFilter.timeRange_push(1));
    idx_pre = find(arfidata.trackTime>options.display.t_disp_pre,1)-1;
else
    idx_pre = find(arfidata.trackTime>options.display.t_disp_pre,1)-1;
end
    
di = idx_push - options.dispEst.dims(1); dt = arfidata.trackTime(idx_push) - arfidata.trackTime(options.dispEst.dims(1));
push = medfilt2(double(abs(disp_push(:,:,idx_push)-disp_push(:,:,idx_push-di))),[nax nt]);
% push = medfilt2(double((disp_push(:,:,idx_push)-disp_push(:,:,idx_push-di))),[nax nt]);
if options.display.cc_filt; push(mask(:,:,idx_push)==0) = inf; end
idx_push = repmat(idx_push,[1 size(push,2)]);

pre = medfilt2(double(abs(disp_pre(:,:,idx_pre)-disp_pre(:,:,idx_pre-di))),[nax nt]);
% pre = medfilt2(double((disp_pre(:,:,idx_pre)-disp_pre(:,:,idx_pre-di))),[nax nt]);
if options.display.cc_filt; pre(mask(:,:,idx_pre)==0) = inf; end
idx_pre = repmat(idx_pre,[1 size(pre,2)]);

% if options.display.normalize % Check functionality
%     m1 = 1/range(pre(:)); b1 = 1 - max(pre(:))/range(pre(:));
%     m2 = 1/range(push(:)); b2 = 1 - max(push(:))/range(push(:));
%     pre = m1.*pre + b1; push = m2.*push + b2;
%     img_rng = [0 1];
% end

if options.display.showPre
    pre_ax = axes('Position',[0.5 0.51 0.4 0.1],'Parent',fig);
    imagesc(arfidata.acqTime,arfidata.axial,pre,options.display.dispRange);
    hold on
    %     cp_pre_ax = copyobj(pre_ax,fig);
    for i=1:size(borders,1)
        preBors{i} = plot(linspace(0,arfidata.acqTime(end),length(arfidata.acqTime)),borders(i,:),'Color',dispPar.trace_cols(i,:),'Linewidth',2);
        %         eval(sprintf('preBor_%d = plot(linspace(0,arfidata.acqTime(end),length(arfidata.acqTime)),borders(i,:),''Color'',%s,''Linewidth'',1);',i,mat2str(dispPar.trace_cols(i,:))))
        %         eval(sprintf('preBors{%d} = preBor_%d;',i,i));
    end
    set(pre_ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','UserData','pre_ax','xlim',[0 max(arfidata.acqTime)],'ylim',[min(edge) max(edge)])
    %     set(cp_pre_ax,'box','on','linewidth',1,'xcolor','y','ycolor','y','xticklabel',[],'yticklabel',[],'xgrid','on','UserData','cp_pre_ax','xlim',[0 max(arfidata.acqTime)],'ylim',[min(edge) max(edge)]);
    ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',pre_ax)
    title(sprintf('Tracked Displacements over DOF at t = %2.2f ms (pre push)',arfidata.trackTime(idx_pre(1))),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',pre_ax)
    set(pre_ax,'ylim',[min(borders(:))-2.5 max(borders(:))+2.5])
    push_ax = axes('Position',[0.5 0.35 0.4 0.1],'Parent',fig);
else
    pre_ax = []; preBors = [];   %cp_pre_ax = [];
    push_ax = axes('Position',[0.5 0.42 0.4 0.15],'Parent',fig);
    %     fig = figure(101);set(fig,'Color','k');push_ax = axes;
end

imagesc(arfidata.acqTime,arfidata.axial,push,options.display.dispRange);
hold on
% cp_push_ax = copyobj(push_ax,gcf);
for i=1:size(borders,1)
    pushBors{i} = plot(linspace(0,arfidata.acqTime(end),length(arfidata.acqTime)),borders(i,:),'Color',dispPar.trace_cols(i,:),'Linewidth',2);
    %     pushBors{i} = plot(linspace(0,arfidata.acqTime(end),length(arfidata.acqTime)),borders(i,:),'Color','g','Linewidth',2);
    %     eval(sprintf('pushBor_%d = plot(linspace(0,arfidata.acqTime(end),length(arfidata.acqTime)),borders(i,:),''Color'',%s,''Linewidth'',1);',i,mat2str(dispPar.trace_cols(i,:))))
    %     eval(sprintf('pushBors{%d} = pushBor_%d;',i,i));
end

set(push_ax,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','xgrid','on','UserData','push_ax','xlim',[0 max(arfidata.acqTime)],'ylim',[min(edge) max(edge)])
% set(cp_push_ax,'box','on','linewidth',1,'xcolor','r','ycolor','r','xticklabel',[],'yticklabel',[],'xgrid','on','UserData','cp_push_ax','xlim',[0 max(arfidata.acqTime)],'ylim',[min(edge) max(edge)]);
ylabel('Axial (mm)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',push_ax)
title(sprintf('Tracked Displacements over DOF at t = %2.2f ms (at push) (dt = %2.2f ms)',arfidata.trackTime(idx_push(1)),dt),'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt,'Parent',push_ax)
set(push_ax,'ylim',[min(borders(:))-2.5 max(borders(:))+2.5])
cb = colorbar;
% set(cb,'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb')

if options.display.showPre
    set(cb,'Position',[0.91 0.35 0.0187/2 0.26],'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb')
else
    set(cb,'Position',[0.91 0.42 0.0187/2 0.15],'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb')
end

if options.display.normalize
    ylabel(cb,'Normalized Displacement','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
else
    ylabel(cb,'Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end

if isempty(arfidata.ecg)
    xlabel('Acquisition Time (s)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end

colormap(dispPar.cmap)
% NaN out displacements filtered out by cc_thresh
pre(pre==inf) = nan;
push(push==inf) = nan;