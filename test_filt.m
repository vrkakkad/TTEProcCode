%%
clear
close all
clc

tic;load res_20140816090612;toc


%%  BPF Filter

[blah, arfidata.disp] = filtArfiData_TTE(arfidata.axial,arfidata.trackTime,arfidata.disp);
[blah, sweidata.disp] = filtArfiData_TTE(sweidata.axial,sweidata.trackTime,sweidata.disp);
clear blah