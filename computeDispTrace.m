function [pre_trace,push_trace,idx] = computeDispTrace(pre,push,axial,gate,n)

for i=1:size(pre,2)
    gate_idx(:,i) = [find(axial>gate(1,i)+1.5,1,'first') find(axial<gate(2,i)-1.5,1,'last')];
%     pre_trace(i) = nanmean(pre(gate_idx(1,i):gate_idx(2,i),i));
%     push_trace(i) = nanmean(push(gate_idx(1,i):gate_idx(2,i),i));
    pre_trace(i) = nanmedian(pre(gate_idx(1,i):gate_idx(2,i),i));
    push_trace(i) = nanmedian(push(gate_idx(1,i):gate_idx(2,i),i));
    idx(:,i) = ceil(linspace(gate_idx(1,i),gate_idx(2,i),n)); % Calculate indices for disp. vs. time plots
end
