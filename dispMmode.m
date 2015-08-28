function [pre,push,idx_pre,idx_push] = dispMmode(options,nax,nt,datastruct,par,borders,edge,mask,dispPar,img_rng)
if (isempty(datastruct.disp_mf_pre) && isempty(datastruct.disp_mf_push))
    disp_pre = datastruct.disp;
    disp_push = datastruct.disp;
else
    disp_pre = datastruct.disp_mf_pre;
    disp_push = datastruct.disp_mf_push;
end
if strcmpi(options.display.t_disp_pre,'max')
    [pk idx] = max(disp_pre(:,:,1:par.nref),[],3);
    idx_pre = round(median(idx)); % to get a single representative number through all depth
    pre = medfilt2(double(pk),[nax nt]);
    % figure out way to cc_filt these properly
    clear pk idx
else
    idx_pre = find(datastruct.trackTime>options.display.t_disp_pre,1);
    pre = medfilt2(double(disp_pre(:,:,idx_pre)),[nax nt]);
    if options.display.cc_filt; pre(mask(:,:,idx_pre)==0) = inf; end
    idx_pre = repmat(idx_pre,[1 size(pre,2)]);
end
if strcmpi(options.display.t_disp_push,'max')
    [pk idx] = max(disp_push(:,:,par.nref:end),[],3);
    idx_push = par.nref -1 + round(median(idx)); % to get a single representative number through all depth
    push = medfilt2(double(pk),[nax nt]);
    % figure out way to cc_filt these properly
    clear pk idx
else
    idx_push = find(datastruct.trackTime>options.display.t_disp_push,1);
    push = medfilt2(double(disp_push(:,:,idx_push)),[nax nt]);
    if options.display.cc_filt; push(mask(:,:,idx_push)==0) = inf; end
    idx_push = repmat(idx_push,[1 size(push,2)]);
end
if options.display.normalize % Check functionality!!
    m1 = 1/range(pre(:)); b1 = 1 - max(pre(:))/range(pre(:));
    m2 = 1/range(push(:)); b2 = 1 - max(push(:))/range(push(:));
    pre = m1.*pre + b1; push = m2.*push + b2;
    img_rng = [0 1];
end
ax13 = axes('Position',[0.5 0.51 0.4 0.1]);
imagesc(datastruct.acqTime,datastruct.axial,pre,img_rng);
hold on
for i=1:size(borders,1)
    eval(sprintf('preBor_%d = plot(linspace(0,datastruct.acqTime(end),length(datastruct.acqTime)),borders(i,:),''Color'',%s,''Linewidth'',1);',i,mat2str(dispPar.trace_cols(i,:))))
end
cax13 = copyobj(ax13,gcf);
set(ax13,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','pre_ax','xlim',[0 max(datastruct.acqTime)],'ylim',[min(edge)-1 max(edge)+1])
set(cax13,'box','on','linewidth',1,'xcolor','y','ycolor','y','xticklabel',[],'yticklabel',[],'xgrid','on','UserData','copy_pre_ax','xlim',[0 max(datastruct.acqTime)],'ylim',[min(edge)-1 max(edge)+1]);
ylabel('Axial (mm)','Parent',ax13,'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
if strcmpi(options.display.t_disp_pre,'max')
    title(sprintf('ARFI Displacements over DOF at t_m_a_x (pre push)'),'parent',ax13,'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
else
    title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (pre push)',datastruct.trackTime(idx_pre(1))),'parent',ax13,'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end
ax14 = axes('Position',[0.5 0.35 0.4 0.1]);
imagesc(datastruct.acqTime,datastruct.axial,push,img_rng);
hold on
for i=1:size(borders,1)
    eval(sprintf('pushBor_%d = plot(linspace(0,datastruct.acqTime(end),length(datastruct.acqTime)),borders(i,:),''Color'',%s,''Linewidth'',1);',i,mat2str(dispPar.trace_cols(i,:))))
end
cax14 = copyobj(ax14,gcf);
set(ax14,'xcolor',dispPar.txt,'ycolor',dispPar.txt,'fontweight','bold','UserData','push_ax','xlim',[0 max(datastruct.acqTime)],'ylim',[min(edge)-1 max(edge)+1])
set(cax14,'box','on','linewidth',1,'xcolor','r','ycolor','r','xticklabel',[],'yticklabel',[],'xgrid','on','UserData','copy_push_ax','xlim',[0 max(datastruct.acqTime)],'ylim',[min(edge)-1 max(edge)+1]);
ylabel('Axial (mm)','parent',ax14,'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
if strcmpi(options.display.t_disp_pre,'max')
    title(sprintf('ARFI Displacements over DOF at t_m_a_x (at push)'),'parent',ax14,'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
else
    title(sprintf('ARFI Displacements over DOF at t = %2.2f ms (at push)',datastruct.trackTime(idx_push(1))),'parent',ax14,'fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end
cb = colorbar;
set(cb,'Position',[0.91 0.35 0.0187/2 0.26],'Color',dispPar.txt,'FontWeight','bold','UserData','disp_cb')
if options.display.normalize
    ylabel(cb,'Normalized Displacement','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
else
    ylabel(cb,'Displacement (\mum)','fontsize',dispPar.fsize,'fontweight','bold','Color',dispPar.txt)
end
colormap(dispPar.cmap)