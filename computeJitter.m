function datastruct = computeJitter(datastruct,options)

if ~options.dispEst.interpFlag
    tempstruct = interpPushReverb(datastruct,options,''); % Interpolate through push-reverb frames
else
    tempstruct = datastruct;
end

disp = tempstruct.disp;

if size(size(disp),2)==4
    [ax,beam,push,tstep] = size(disp);
    disp = reshape(disp,ax,beam*push,tstep);
    reshape_flag = 1;
else
    reshape_flag = 0;
end

dt = median(diff(datastruct.trackTime)); fs = 1./dt*1e3;
[B,A] = butter(3, options.motionFilter.Cutoff(2)./(fs/2),'high');
jitter = temporalFilter(disp, datastruct.trackTime, datastruct.trackTime, B, A);

if reshape_flag
    jitter = reshape(jitter,ax,beam,push,tstep);
end

datastruct.jitter = jitter;

end

function data = temporalFilter(arfidata, t, tn, B, A)
% interpolate and filter
D = size(arfidata);
arfidata = reshape(arfidata, [], D(end));
data = arfidata;
% data = interp1(t(:), arfidata', tn(:), 'spline')';
D(end) = length(tn);
data = single(filtfilt(B,A,double(data)'))';
data = reshape(data,D);
end
