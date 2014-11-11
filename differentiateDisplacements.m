
function [vel t1] = differentiateDisplacements(arfidata,t,passBandHz)
vel = diff(arfidata,1,3)/diff(t(1:2));
t1 = 0.5*t(1:end-1) + 0.5*t(2:end);
if exist('passBandHz','var') && ~isempty(passBandHz)
[BB,AA] = butter(3,min(max(passBandHz*2*mean(diff(t*1e-3)),1e-12),1));
vel = flipdim(filter(BB,AA,flipdim(filter(BB,AA,vel,[],3),3),[],3),3);
end
end
