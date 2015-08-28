function [diff_data t1] = differentiateData(data,t,cutoffFreq)
diff_data = diff(data,1,3)/diff(t(1:2));
t1 = 0.5*t(1:end-1) + 0.5*t(2:end);
if exist('cutoffFreq','var') && ~isempty(cutoffFreq)
    dt = median(diff(t));
    fs = 1./dt*1e3;
    [B,A] = butter(2, cutoffFreq./(fs/2),'low');
    idx = find(isnan(squeeze(diff_data(round(size(diff_data,1)/2),1,:))));
    % Use filtfilt instead of filter
    pre = diff_data(:,:,[1:min(idx)-1]);
    D1 = size(pre);
    pre = reshape(pre,[],D1(end));
    pre = single(filtfilt(B,A,double(pre)'))';
    pre = reshape(pre,D1); diff_data(:,:,[1:min(idx)-1]) = pre;
    
    push = diff_data(:,:,[max(idx)+1:end]);
    D2 = size(push);
    push = reshape(push,[],D2(end));
    push = single(filtfilt(B,A,double(push)'))';
    push = reshape(push,D2); diff_data(:,:,[max(idx)+1:end]) = push;
end
