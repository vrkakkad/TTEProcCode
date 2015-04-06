function bdata = extractBmode(timeStamp)

bmodeParFname = sprintf('SWIF_BModeOutputImageDims0_%s.txt', timeStamp);
bmodeFname = sprintf('SWIF_BModeOutputImage0_%s.img', timeStamp);

if ~exist(bmodeParFname, 'file')
    warning('B-mode parameters file not detected');
elseif ~exist(bmodeFname, 'file')
    warning('B-mode data file not detected');
else
    bmodePar = struct;
    bmodeParFid = fopen(bmodeParFname, 'r');
    while 1
        tline = fgetl(bmodeParFid);
        if ~ischar(tline)||isempty(tline), break, end
        k = strfind(tline, ':');
        if ~isempty(k)
            eval(sprintf('bmodePar.%s = %s;', tline(1:k-1), tline(k+1:end)));
        else
            if isempty(which('strsplit'))
                tmp = textscan(tline,'%s');
                tmp = tmp{1};
                try
                    eval(sprintf('bmodePar.%s = %s;', tmp{1}, tmp{2}));
                catch err
                    try
                        eval(sprintf('bmodePar.%s = ''%s'';', tmp{1}, tmp{2}));
                    catch
                        if length(tmp)==2
                            rethrow(err)
                        end
                    end
                end
            else
                tmp = strsplit(tline);
                try
                    eval(sprintf('bmodePar.%s = %s;', tmp{2}, tmp{3}));
                catch err
                    try
                        eval(sprintf('bmodePar.%s = ''%s'';', tmp{2}, tmp{3}));
                    catch
                        rethrow(err)
                    end
                end
            end
        end
    end
    fclose(bmodeParFid);
    bmodeFid = fopen(bmodeFname, 'rb');
    bimg = fread(bmodeFid, inf, '*uint8');
    fclose(bmodeFid);
    bimg = reshape(bimg, bmodePar.SamplesPerLine, bmodePar.LinesPerSlice, []);
    bimg = bimg(:,:,1:end-1);
    bax = (0:bmodePar.SamplesPerLine-1)./bmodePar.NumSamplesPerMm;
    blat = linspace(bmodePar.FirstLinePos, bmodePar.LastLinePos, bmodePar.LinesPerSlice);
        
    % Scan Convert B-mode (if required)
    fprintf(1,'Scan Converting B-mode Data...\n')
    if strcmpi(bmodePar.PqMethod,'EqualSineAngle')
        [bdata.bimg,bdata.bax,bdata.blat] = scan_convert('sector',single(bimg),bmodePar.FovAzimMin,bmodePar.FovAzimSpan,...
            0.1*bmodePar.VectorApexAzimZMm,2,(1540/2)/(1e-3*bax(2)-bax(1)));
    end
    
    bdata.t = 0:1/round(bmodePar.ChunkRateHz):(size(bimg,3)-1)/round(bmodePar.ChunkRateHz);
    bdata.bax = 10*bdata.bax;
    bdata.blat = 10*bdata.blat;
end