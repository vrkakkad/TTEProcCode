function dispRes(DataDir,fidx,type)

close all

if ~exist('DataDir','var')
    DataDir = pwd;
end
if ~exist('fidx','var')
    fidx = 1;
end

if ~exist('type','var')
    type = 'Both';
end
cd(DataDir)

list = dir('res_*'); % get timeStamp based on existance of ARFI par files
if size(list,1)<fidx
    error('Data set index requested greater than number of data sets')
end

timeStamp = list(fidx).name(end-17:end-4);
fprintf('Loading data with timeStamp = %s (Set # %d of %d)\n', timeStamp,fidx,size(list,1));

load(strcat('res_',timeStamp,'.mat'));

switch type
    case 'ARFI'
        options.dataflow.ARFI = 1;
        options.dataflow.SWEI = 0;
    case 'SWEI'
        options.dataflow.ARFI = 0;
        options.dataflow.SWEI = 1;
    case 'Both'
end

arfi_par = load(strcat('arfi_par_',timeStamp));
swei_par = load(strcat('swei_par_',timeStamp));
dispTTE(ecgdata,bdata,arfidata,arfi_par,sweidata,swei_par,options);