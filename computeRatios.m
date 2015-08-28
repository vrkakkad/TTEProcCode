function [pre_DR,push_DR] = computeRatios(ecg_ax,dispTrace_ax,pre_trace,push_trace,acqTime,ecg,hr,samples,dispPar,options)

ch = get(ecg_ax,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'sys_markers') || strcmpi(nm,'dia_markers'))
        delete(ch(i))
    end
end
ch = get(dispTrace_ax,'Children');
for i=1:length(ch)
    nm = get(ch(i),'UserData');
    if (strcmpi(nm,'ind_plot'))
        delete(ch(i))
    end
end
[pre_DR,push_DR,idx_sys,idx_dia] = calcDR(pre_trace,push_trace,acqTime,ecg,hr,options.display.sysRange,options.display.diaRange);
indicator = zeros(1,length(acqTime));
indicator(idx_sys) = -2; indicator(idx_dia) = 2;
stem(acqTime,indicator,'Marker','.','LineWidth',3,'Color',dispPar.txt,'Parent',dispTrace_ax,'UserData','ind_plot');
plot(acqTime(idx_sys),samples(idx_sys),'ws','MarkerSize',10,'Parent',ecg_ax,'UserData','sys_markers');
plot(acqTime(idx_dia),samples(idx_dia),'wo','MarkerSize',10,'Parent',ecg_ax,'UserData','dia_markers');
