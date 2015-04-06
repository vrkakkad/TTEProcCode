function [bmode arfi swei hr] = extractECG(timeStamp,plot_flag,dt)

load(strcat('ECG_data_',timeStamp));

fs = 1/(ecgdata(2,1)-ecgdata(1,1));

% Filter ECG
[B A] = butter(1,[1/fs],'high');
B = double(B); A = double(A);
ecgdata(:,3) = single(filtfilt(B,A,double(ecgdata(:,3))));
ecgdata(:,3) = ecgdata(:,3)/max(ecgdata(:,3));

hr = calcHR(ecgdata(:,1),ecgdata(:,3),plot_flag);

trig_trace_threshold = 5e-3;
high = ecgdata(:,2)>trig_trace_threshold;
ind_high = find(high==1);
t_high = find(high==1)/fs;
% ind_gap = find(ind_high);
t_gap = diff(t_high);

if plot_flag
    figure(101)
    subplot(313)
    stem(t_gap)
    grid on
    subplot(311)
    plot(ecgdata(:,1),ecgdata(:,2))
    title(timeStamp)    
    grid on
    hold on
end

try
    % Calculate Start and Stop Indices for HQ B-mode Grab
    ind_t_gap1 = find(t_gap>0.075,1,'first');
    val_t_gap1 = t_gap(ind_t_gap1);
    ind_b_start = round(ind_high(ind_t_gap1) + val_t_gap1*fs);
    temp = find(t_gap>0.25,2,'first');
    ind_b_stop = round(ind_high(temp(end)));
    
    % Calculate Start and Stop Indices for ARFI Grab
    temp = find(t_gap>0.075,4,'last');
    ind_t_gap2 = temp(1);
    val_t_gap2 = t_gap(ind_t_gap2);
    ind_a_start = round(ind_high(ind_t_gap2) + val_t_gap2*fs);
    ind_t_gap3 = temp(2);
    val_t_gap3 = t_gap(ind_t_gap3);
    ind_a_stop = round(ind_high(ind_t_gap3));
    
    % Calculate Start and Stop Indices for SWEI Grab
    ind_s_start = round(ind_high(ind_t_gap3) + val_t_gap3*fs);
    ind_t_gap4 = temp(3);
    val_t_gap4 = t_gap(ind_t_gap4);
    ind_s_stop = round(ind_high(ind_t_gap4));
    
    % Auto Check for ECG Data (Stage 1)
    t = ecgdata(:,1);
    dt_b = t(ind_b_stop) - t(ind_b_start); dt_a = t(ind_a_stop) - t(ind_a_start); dt_s = t(ind_s_stop) - t(ind_s_start);
    dt_b_rng = dt(1)+[-0.1 0.1]; dt_a_rng = dt(2)+[-0.1 0.1]; dt_s_rng = dt(3)+[-0.1 0.1];
    
    % Insert ECG segments into structure
    if (dt_b > dt_b_rng(1) && dt_b < dt_b_rng(2))
        bmode(:,1) = [0:1/fs:(ind_b_stop-ind_b_start)/fs];
        bmode(:,2) = ecgdata(ind_b_start:ind_b_stop,3);
        bmode(:,2) = bmode(:,2)/max(bmode(:,2));
    else
        bmode = [];
    end
    if (dt_a > dt_a_rng(1) && dt_a < dt_a_rng(2))
        arfi(:,1) = [0:1/fs:(ind_a_stop-ind_a_start)/fs];
        arfi(:,2) = ecgdata(ind_a_start:ind_a_stop,3);
        arfi(:,2) = arfi(:,2)/max(arfi(:,2));
    else
        arfi = [];
    end
    if (dt_s > dt_s_rng(1) && dt_s < dt_s_rng(2))
        swei(:,1) = [0:1/fs:(ind_s_stop-ind_s_start)/fs];
        swei(:,2) = ecgdata(ind_s_start:ind_s_stop,3);
        swei(:,2) = swei(:,2)/max(swei(:,2));
    else
        swei = [];
    end
    
    if (isempty(bmode) || isempty(arfi) || isempty(swei))
        
        clear bmode arfi swei
        % Try again with 3 stops from end
        % Calculate Start and Stop Indices for ARFI Grab
        temp = find(t_gap>0.075,3,'last');
        ind_t_gap2 = temp(1);
        val_t_gap2 = t_gap(ind_t_gap2);
        ind_a_start = round(ind_high(ind_t_gap2) + val_t_gap2*fs);
        ind_t_gap3 = temp(2);
        val_t_gap3 = t_gap(ind_t_gap3);
        ind_a_stop = round(ind_high(ind_t_gap3));
        
        % Calculate Start and Stop Indices for SWEI Grab
        ind_s_start = round(ind_high(ind_t_gap3) + val_t_gap3*fs);
        ind_t_gap4 = temp(3);
        val_t_gap4 = t_gap(ind_t_gap4);
        ind_s_stop = round(ind_high(ind_t_gap4));
        
        % Auto Check for ECG Data (Stage 2)
        t = ecgdata(:,1);
        dt_b = t(ind_b_stop) - t(ind_b_start); dt_a = t(ind_a_stop) - t(ind_a_start); dt_s = t(ind_s_stop) - t(ind_s_start);
        dt_b_rng = dt(1)+[-0.1 0.1]; dt_a_rng = dt(2)+[-0.1 0.1]; dt_s_rng = dt(3)+[-0.1 0.1];
        
        % Insert ECG segments into structure
        if (dt_b > dt_b_rng(1) && dt_b < dt_b_rng(2))
            bmode(:,1)= [0:1/fs:(ind_b_stop-ind_b_start)/fs];
            bmode(:,2)= ecgdata(ind_b_start:ind_b_stop,3);
            bmode(:,2) = bmode(:,2)/max(bmode(:,2));
        else
            bmode = [];
        end
        if (dt_a > dt_a_rng(1) && dt_a < dt_a_rng(2))
            arfi(:,1)= [0:1/fs:(ind_a_stop-ind_a_start)/fs];
            arfi(:,2)= ecgdata(ind_a_start:ind_a_stop,3);
            arfi(:,2) = arfi(:,2)/max(arfi(:,2));
        else
            arfi = [];
        end
        if (dt_s > dt_s_rng(1) && dt_s < dt_s_rng(2))
            swei(:,1)= [0:1/fs:(ind_s_stop-ind_s_start)/fs];
            swei(:,2)= ecgdata(ind_s_start:ind_s_stop,3);
            swei(:,2) = swei(:,2)/max(swei(:,2));
        else
            swei = [];
        end
    end
    
    if (isempty(bmode) || isempty(arfi) || isempty(swei))
        warning('ECG Triggers are misaligned: Check manually')
        close gcf
    end
catch
    disp('Error!')
    warning('ECG Triggers are misaligned: Check manually')
    close gcf
    bmode = []; arfi = []; swei = [];
    return
end

if plot_flag
    plot(ecgdata(ind_b_start,1),1.5,'o','markerfacecolor','r')
    plot(ecgdata(ind_b_stop,1),1.5,'o','markerfacecolor','r')
    plot(ecgdata(ind_a_start,1),1.5,'o','markerfacecolor','y')
    plot(ecgdata(ind_a_stop,1),1.5,'o','markerfacecolor','y')
    plot(ecgdata(ind_s_start,1),1.5,'o','markerfacecolor','g')
    plot(ecgdata(ind_s_stop,1),1.5,'o','markerfacecolor','g')
    xlabel('Acq Time (s)')
    pause
    close gcf
end
