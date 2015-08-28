function [recBors,bimg_ax,becg_ax] = drawBmode(fig,bdata,edge,borders,dispPar,par,playFlag,options)

ch = get(fig,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'bimg_ax') || strcmpi(nm,'becg_ax'))
        delete(ch(i))
    end
end
% [~,temp] = fileparts(pwd);
% filename = strcat('E:\ClinicalDataArchive\TracedBorders\',temp,'\',num2str(options.dataflow.setID),'_',options.timeStamp,'.gif');
% if exist(filename,'file') saveFlag = 0; else saveFlag = 1; end
% if saveFlag = 1
%     fig = figure(101);set(fig,'Color','k');playFlag = 1;
% end
if ~isempty(bdata.ecg)
    bsamples = zeros(1,size(bdata.bimg,3));
    for i=1:size(bdata.bimg,3)
        bsamples(i) = bdata.ecg(find(bdata.ecg(:,1)>bdata.t(i),1,'first')-1,2);
    end
    becg_ax = axes('Position',[0.125 0.5 0.25 0.1],'Parent',fig);
%     becg_ax = axes('Position',[0.15 0.1 0.7 0.2],'Parent',fig);
    plot(bdata.ecg(:,1),bdata.ecg(:,2),'LineWidth',2,'Parent',becg_ax);
    hold(becg_ax,'on')
    plot(bdata.t,bsamples,'gx','MarkerSize',4,'Parent',becg_ax)
    pt = plot(bdata.t(1),bsamples(1),'ro','MarkerSize',6,'MarkerFaceColor','r','Parent',becg_ax);
    title(sprintf('ECG (Bmode): HR = %3.0f bpm',bdata.hr),'FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    xlabel('Acquisition Time (s)','FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    xlim([0 max(bdata.t)])
    ylim([min(bdata.ecg(:,2)) max(bdata.ecg(:,2))])
    set(becg_ax,'Color',dispPar.ax,'XColor',dispPar.txt,'YColor',dispPar.txt,'YTickLabel',[],'FontWeight','Bold','XGrid','On','GridLineStyle','--','UserData','becg_ax')
end
bimg_ax = axes('Position',[0.1 0.65 0.3 0.3],'Parent',fig);
% bimg_ax = axes('Position',[0.2 0.35 0.6 0.6],'Parent',fig); % with ecg
% bimg_ax = axes('Position',[0.1 0.1 0.8 0.8],'Parent',fig); % w/o ecg
if playFlag, n = size(bdata.bimg,3); else n = 1; end
for i=1:n
    delete(bimg_ax)
    bimg_ax = axes('Position',[0.1 0.65 0.3 0.3],'Parent',fig);
%     bimg_ax = axes('Position',[0.2 0.35 0.6 0.6],'Parent',fig); % with ecg
%     bimg_ax = axes('Position',[0.1 0.1 0.8 0.8],'Parent',fig); % w/o ecg
    imagesc(bdata.blat,bdata.bax,fliplr(bdata.bimg(:,:,i)));
    colormap(gray);axis image; freezeColors;
    hold(bimg_ax,'on')
    plot(0,par.pushFocalDepth,'o','MarkerSize',6,'MarkerFaceColor','c')
    rectangle('Position',[-7 edge(1) 14 edge(2)-edge(1)],'EdgeColor','b','Linewidth',2,'Parent',bimg_ax);
    for j=1:size(borders,1)-1
        recBors{j} = rectangle('Position',[-2 min(borders(j,:)) 4 max(borders(j+1,:))-min(borders(j,:))],'EdgeColor',dispPar.trace_cols(j,:),'Linestyle','-','Linewidth',2,'Parent',bimg_ax);
        %         recBors{j} = rectangle('Position',[-2 min(borders(j,:)) 4 max(borders(j+1,:))-min(borders(j,:))],'EdgeColor','g','Linestyle','-','Linewidth',2,'Parent',bimg_ax);
    end
%         xlabel('Lateral (mm)','FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
%         ylabel('Axial (mm)','FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    title(sprintf('B-Mode Cine: Frame %d (t = %1.1f s)',i,bdata.t(i)),'FontSize',dispPar.fsize,'FontWeight','Bold','Color',dispPar.txt)
    set(bimg_ax,'XColor',dispPar.txt,'YColor',dispPar.txt,'FontWeight','Bold','UserData','bimg_ax')
    grid on
    hold(bimg_ax,'off')
    if ~isempty(bdata.ecg)
        delete(pt)
        pt = plot(bdata.t(i),bsamples(i),'ro','MarkerSize',6,'MarkerFaceColor','r','Parent',becg_ax);
    end
    pause(median(diff(bdata.t)))
%         pause
    
%     frame = getframe(101);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1;
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0);
%     else
%         imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
%     end

end