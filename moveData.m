function moveData(dest,src)

if nargin<1
    error('Enter a destnation directory to move data to')
elseif nargin<2
    src = 'V:\Program Files\Siemens\syngo\Bedrock\Startup\';
end

if ~exist(dest)
    warning('Destination folder does not exist...creating it')
    mkdir(dest)
end
if ~exist(src)
    error('Scanner drives not mounted into laptop')
end

orig = pwd;
cd(src)

list = dir('par_*');
nfiles = size(list,1);

if nfiles>0
    warning('Data from previous acquisitions exists in startup, clean up source directory')
end

fprintf(1,'Waiting for data files to be written to %s...\n',src);

while nfiles<2
    pause(1)
    list = dir('par_*');
    nfiles = size(list,1);
end

disp('Data files detected...')

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

disp('Data files transferred.')
cd(orig)