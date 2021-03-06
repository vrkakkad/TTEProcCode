function [u Iup Qup cc] = runPesavento_gated(I, Q, interpFactor, kernelLength,searchLength, axial, par, start_depth, end_depth)
% function [u Iup Qup] = runLoupas(I, Q, interpFactor, axial, par)
%
% Inputs: I - in-phase data
%         Q - quadrature data
%         interpFactor - upsampling factor
%         axial - axial vector (used for demodulation)
%         par - parameters structure generated by arfi_image
%
% Outputs: u - displacement matrix
%          Iup - upsampled in-phase data
%          Qup - upsampled quadrature data

start_idx = find(axial(1:interpFactor:end)>start_depth,1); 
end_idx = find(axial(1:interpFactor:end)>end_depth,1);

if isempty(start_idx);start_idx = 1;end
if isempty(end_idx);end_idx = size(I,1);end

I = I(start_idx:end_idx,:,:);
Q = Q(start_idx:end_idx,:,:);

% Setup parameters
fs = par.fs*1e6;
fc = par.fc;
c = par.c; % m/s
NumIter = 3;

D = size(I); 
D(1) = D(1).*interpFactor;
tic
if interpFactor>1
    %fprintf('Upsampling...')
    [Iup,Qup] = computeUpsampledIQdata(I,Q,interpFactor);
    %fprintf('done (%0.1fs)',toc);
else
    Iup = I;
    Qup = Q;
end
Iup = reshape(Iup, D);
Qup = reshape(Qup, D);
% idx post to upsampling data
start_idx = find(axial>start_depth,1); 
end_idx = find(axial>end_depth,1);

fs = fs*interpFactor;
KLen = ceil(kernelLength*fs/fc);
if KLen < 3
    warning('Extending Kernel Length to 3 sample (%0.1f wavelengths).',KLen*fc/fs);
    KLen = 3;
end
SrchLen = round(searchLength*fs/fc);

fdem = par.Apl3.Mod(1).DsF.data*1e6; % frequency dataset values (Hz shift)
frange = par.Apl3.Mod(1).DsF.rr;  % reference ranges (mm)
fc_vec = reshape(interp1(frange, fdem, axial(start_idx:start_idx+size(Iup,1)-1)), size(Iup,1), 1);
fdem_vec = fc_vec./fs;

tstart = tic;
IQ = complex(Iup,Qup);
IQref = repmat((IQ(:,:,par.nref)),[1 1 size(IQ,3)]);
%[Lags IQ2] = PesaventoGross(IQref,IQ,KLen,searchLength);
[Tau cc]= PesaventoParallel4(IQref,IQ,fs,fc,KLen,SrchLen,NumIter);
%[Tau0 cc0]= PesaventoParallelWholeLine2(reshape(IQref,size(IQ,1),par.nBeams,par.numBeamGroups,[]),reshape(IQ,size(IQ,1),par.nBeams,par.numBeamGroups,[]),fs,fc,SrchLen,NumIter);
%Tau0 = repmat(Tau0,[size(Tau,1),par.nBeams,1,1]);
%Tau0=  reshape(Tau0,size(Tau));
%u = (-1*Lags*(1/fs) + Tau) * c/2 * 1e6;
u = Tau*c/2*1e6;
%u0 = Tau0*c/2*1e6;
tend = toc(tstart);
%fprintf(1, 'Displacement Computation Time: %0.2fs\n', tend);
