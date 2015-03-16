function ecgstruct = extractECG(timeStamp,plot_flag,dt)

load(strcat('ECG_data_',timeStamp));

fs = 1/(ecgdata(2,1)-ecgdata(1,1));
trig_trace_threshold = 5e-3;
high = ecgdata(:,2)>trig_trace_threshold;
ind_high = find(high==1);
t_high = find(high==1)/fs;
% ind_gap = find(ind_high);
t_gap = diff(t_high);

if plot_flag
    figure(501)
    subplot(211)
    stem(t_gap)
    grid on
    title(timeStamp)
    subplot(212)
    plot(ecgdata(:,1),ecgdata(:,2))
    grid on
    hold on
end

try
% Calculate Start and Stop Indices for HQ B-mode Grab
ind_t_gap1 = find(t_gap>0.05,1,'first');
val_t_gap1 = t_gap(ind_t_gap1);
ind_b_start = round(ind_high(ind_t_gap1) + val_t_gap1*fs);
temp = find(t_gap>0.25,2,'first');
ind_b_stop = round(ind_high(temp(end)));

% Calculate Start and Stop Indices for ARFI Grab
temp = find(t_gap>0.05,4,'last'); %4
ind_t_gap2 = temp(1);
val_t_gap2 = t_gap(ind_t_gap2);
ind_a_start = round(ind_high(ind_t_gap2) + val_t_gap2*fs);
ind_t_gap3 = temp(2);
val_t_gap3 = t_gap(ind_t_gap3);
ind_a_stop = round(ind_high(ind_t_gap3));

% Calculate Start and Stop Indices for SWEI Grab
ind_s_start = round(ind_high(ind_t_gap3) + val_t_gap3*fs);
ind_t_gap4 = temp(3); %3
val_t_gap4 = t_gap(ind_t_gap4);
ind_s_stop = round(ind_high(ind_t_gap4));

if plot_flag
    plot(ecgdata(ind_b_start,1),1.5,'o','markerfacecolor','r')
    plot(ecgdata(ind_b_stop,1),1.5,'o','markerfacecolor','r')
    plot(ecgdata(ind_a_start,1),1.5,'o','markerfacecolor','y')
    plot(ecgdata(ind_a_stop,1),1.5,'o','markerfacecolor','y')
    plot(ecgdata(ind_s_start,1),1.5,'o','markerfacecolor','g')
    plot(ecgdata(ind_s_stop,1),1.5,'o','markerfacecolor','g')
    xlabel('Acq Time (s)')
    pause
    close(501)
end

% Auto Check for ECG Data
t = ecgdata(:,1);
% dt_b = t(ind_b_stop) - t(ind_b_start);
dt_a = t(ind_a_stop) - t(ind_a_start);
dt_s = t(ind_s_stop) - t(ind_s_start);
dt_b_rng = dt(1)+[-0.1 0.1]; dt_a_rng = dt(2)+[-0.1 0.1]; dt_s_rng = dt(3)+[-0.1 0.1]; 

% if ((dt_b > dt_b_rng(1) && dt_b < dt_b_rng(2)) && (dt_a > dt_a_rng(1) && dt_a < dt_a_rng(2)) && (dt_s > dt_s_rng(1) && dt_s < dt_s_rng(2)))
if ((dt_a > dt_a_rng(1) && dt_a < dt_a_rng(2)) && (dt_s > dt_s_rng(1) && dt_s < dt_s_rng(2)))
    % Insert ECG segments into structure
    ecgstruct.bmode(:,1)= [0:1/fs:(ind_b_stop-ind_b_start)/fs];
    ecgstruct.bmode(:,2)= ecgdata(ind_b_start:ind_b_stop,3);
    
    ecgstruct.arfi(:,1)= [0:1/fs:(ind_a_stop-ind_a_start)/fs];
    ecgstruct.arfi(:,2)= ecgdata(ind_a_start:ind_a_stop,3);
    
    ecgstruct.swei(:,1)= [0:1/fs:(ind_s_stop-ind_s_start)/fs];
    ecgstruct.swei(:,2)= ecgdata(ind_s_start:ind_s_stop,3);
else
    warning('ECG Triggers are misaligned: Check manually')
    ecgstruct = [];
end

catch
    disp('Catch')
    warning('ECG Triggers are misaligned: Check manually')
    close(501)
    ecgstruct = [];
end