function hr = calcHR(time,ecg,plot_flag)

% Calculates Heartrate based using time vector and ECG data

threshhold = 0.75;
[pks, idx] = findpeaks(double(ecg),'MinPeakHeight',threshhold);

if plot_flag
    figure(101)
    set(101,'Position',[9 49 704 767])
    subplot(312)
    plot(time,ecg); grid on;hold on
    plot(time,threshhold*ones(1,length(time)),'g','linewidth',2)
    plot(time(idx),pks,'v','MarkerFaceColor','r')
    xlabel('Acq Time(s)')
    axis tight
end


for i=1:length(idx)-1
    dt(i)=(time(idx(i+1))-time(idx(i)));
end

dt = median(dt);
hr = round(60/dt);