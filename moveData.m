function moveData(ecg_flag,dest,src)
if nargin<2
    ecg_flag = 0;
    error('Enter a destnation directory to move data to')
elseif nargin<3
    src = 'V:\Program Files\Siemens\syngo\Bedrock\Startup\';
end
if ~exist(dest,'dir')
    warning('Destination folder does not exist...creating it')
    mkdir(dest)
end
if ~exist(src,'dir')
    error('Scanner drives not mounted into laptop')
end
orig = pwd;
cd(src)
list = dir('par_*');
nfiles = size(list,1);
if nfiles>0
    warning('Data from previous acquisitions exists in startup, clean up source directory')
end
if strcmpi(src,'V:\Program Files\Siemens\syngo\Bedrock\Startup\')
    if ecg_flag
        fprintf(1,'Waiting for grabTTE command...\n')
        check = exist(strcat(src,'start_ecg.txt'),'file');
        while check==0
            check = exist(strcat(src,'start_ecg.txt'),'file');
        end
        fprintf(1,'ECG Acquisition Triggered...\n')
        fprintf(1,'Waiting for data files to be written to %s...\n',src);
        clear time data
        [time,data] = usbaqc();
        fprintf(1,'ECG Acquisition Complete.\n');
    end
end
if ~ecg_flag
    fprintf(1,'Waiting for data files to be written to %s...\n',src);
end
while nfiles<2
    pause(1)
    list = dir('par_*');
    nfiles = size(list,1);
end
disp('Data files detected.')
tic
timestamp_arfi = list(1).name(5:end-4);
timestamp_swei = list(2).name(5:end-4);
% Move Bmode Files
bfiles = dir('SWIF_BModeOutputImage*_20140101121212*');
for i=1:size(bfiles,1)
    name = strcat(bfiles(i).name(1:end-18),timestamp_arfi,bfiles(i).name(end-3:end));
    movefile(strcat(src,'\',bfiles(i).name),strcat(dest,'\',name));
end
% Move ARFI Files
arfifiles = dir(strcat('*',timestamp_arfi,'*'));
for i=1:size(arfifiles,1)
    name = strcat('arfi_',arfifiles(i).name);
    movefile(strcat(src,'\',arfifiles(i).name),strcat(dest,'\',name));
end
% Move SWEI Files
sweifiles = dir(strcat('*',timestamp_swei,'*'));
for i=1:size(sweifiles,1)
    name = strcat('swei_',sweifiles(i).name(1:end-18),timestamp_arfi,sweifiles(i).name(end-3:end));
    movefile(strcat(src,'\',sweifiles(i).name),strcat(dest,'\',name));
end
% Save ECG File
if strcmpi(src,'V:\Program Files\Siemens\syngo\Bedrock\Startup\')
    if ecg_flag
        ecgdata(:,1) = single(time); data = single(data);
        ecgdata(:,2) = data(:,1);
        % Filter ECG Data to remove interference with scanner triggers
        [B A] = butter(2,1/2500); 
        data(:,2) = medfilt1(double(data(:,2)),10);
        ecgdata(:,3) = single(filtfilt(B,A,double(data(:,2))));
        save(strcat(dest,'\ECG_data_',timestamp_arfi),'ecgdata');
    end
end
% Remove Dims0_old file
temp = dir('*_old*');
if ~isempty(temp)
    for i=1:size(temp,1)
        delete(temp.name);
    end
end
% Remove leftover/stray SWIF_AData* files
temp = dir('SWIF_AData*');
if ~isempty(temp)
    for i=1:size(temp,1)
        delete(temp.name);
    end
end
% Remove ECG trigger
if exist(strcat(src,'start_ecg.txt'),'file')
    delete(strcat(src,'start_ecg.txt'))
end
fprintf(1,'Data files transferred.\nTime Elapsed = %2.2f s\n',toc)
cd(orig)