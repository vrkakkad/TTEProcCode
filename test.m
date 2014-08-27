%%
clear
close all
clc

tic;load res_20140816090612;toc


%%  BPF Filter

temp1 = zeros(size(arfidata.disp));
temp2 = zeros(size(sweidata.disp));

[blah,temp1] = filtArfiData(arfidata.axial,arfidata.trackTime,arfidata.disp);
[blah,temp2] = filtArfiData(sweidata.axial,sweidata.trackTime,sweidata.disp);