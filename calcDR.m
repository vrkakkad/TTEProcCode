function [pre_DR,push_DR,idx_sys,idx_dia] = calcDR(pre,push,acqTime,ecg,hr,sysRange,diaRange)
% Calculates Diastolic to Systolic Displacement ratios

t = ecg(:,1); sig = ecg(:,2);
threshhold = 0.75;
[pks, idx] = findpeaks(double(sig),'MinPeakHeight',threshhold);
if length(idx)>1
    new_hr = calcHR(t,sig,0);
    if isempty(new_hr)
        new_hr = calcHR(t,sig,1);
    else
        hr = new_hr;
    end
end
mod_t = mod((t - t(idx(1)))*(hr/60),1);
% figure;plot(mod_t,sig,'.');grid on
% pause
% close gcf
mod_acqTime = mod((acqTime - t(idx(1)))*(hr/60),1);
idx_sys = find(mod_acqTime>=sysRange(1) & mod_acqTime<=sysRange(2));
idx_dia = find(mod_acqTime>=diaRange(1) & mod_acqTime<=diaRange(2));
for i=1:size(pre,1)
    pre_DR(i) = mean(pre(i,idx_dia))/mean(pre(i,idx_sys));
    push_DR(i) = mean(push(i,idx_dia))/mean(push(i,idx_sys));
end
    


