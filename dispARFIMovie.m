function dispARFIMovie(arfidata,edge,borders,dispPar,par,options)

if options.display.cc_filt
    mask = arfidata.cc>options.display.cc_thresh;
else
    mask = [];
end
% Indices corresponding to median filter parameters
nax = double(ceil(options.display.medfilt(1)/(arfidata.axial(2) - arfidata.axial(1))));
nt = double(ceil(options.display.medfilt(2)/(arfidata.acqTime(2) - arfidata.acqTime(1))));

if (isempty(arfidata.disp_mf_pre) && isempty(arfidata.disp_mf_push))
    disp_pre = arfidata.disp;
    disp_push = arfidata.disp;
else
    disp_pre = arfidata.disp_mf_pre;
    disp_push = arfidata.disp_mf_push;
end

fig = figure;

keyboard

frame = nan(size(arfidata.disp));
for i=1:par.nref+par.npush+par.nreverb
    frame(:,:,i) = medfilt2(double(disp_pre(:,:,i)),[nax nt]);
    if options.display.cc_filt; pre(mask(:,:,i)==0) = inf; end
end
for i=par.nref+par.npush+par.nreverb+1:size(arfidata.disp,3)
    frame(:,:,i) = medfilt2(double(disp_push(:,:,i)),[nax nt]);
    if options.display.cc_filt; pre(mask(:,:,i)==0) = inf; end
end
    