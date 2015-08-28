function hr = calcHR(time,ecg,plot_flag)
% Calculates Heartrate based using time vector and ECG data

threshhold = 0.75;
[pks, idx] = findpeaks(double(ecg),'MinPeakHeight',threshhold);
if plot_flag
    figure(501)
    set(501,'Position',[9 49 704 767])
    subplot(312)
    plot(time,ecg); grid on;hold on
    plot(time,threshhold*ones(1,length(time)),'g','linewidth',2)
    plot(time(idx),pks,'v','MarkerFaceColor','r')
    xlabel('Acq Time(s)')
    axis tight
end
if length(idx)<2
    hr = [];
    return
end
for i=1:length(idx)-1
    dt(i)=(time(idx(i+1))-time(idx(i)));
end
dt = median(dt);
hr = 60/dt;

if (hr<35 || hr>135)
    hr = [];
end