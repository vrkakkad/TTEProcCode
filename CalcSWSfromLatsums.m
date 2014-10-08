
%==========================================================================

function out=CalcSWSfromLatsums(plane,latmm,tms,maxspeed,minlat,maxlat,fignum)

        %
        % careful, input lat in mm, tms in ms
        %

    ilatstart = min(find(latmm>=minlat));
    ilatend   = max(find(latmm<=maxlat));

    latsums=FindLatSumsUsingEndpoints(plane,ilatstart,ilatend);

%    latsums = LimitMaxSpeed(plane,latsums,latmm,tms,maxspeed);
    latsums = LimitMaxSpeed(plane,latsums,latmm,tms,maxspeed,ilatstart,ilatend);

    [itstart,itend] = FindGlobalPeak(latsums);

    if fignum
        figure(fignum)
        clf
        subplot(1,2,1)
        imagesc(tms,latmm,plane)
        hold on
        plot([tms(itstart) tms(itend)],[latmm(ilatstart) latmm(ilatend)],'w')
        hold off
        xlabel('time (ms)');
        ylabel('lateral position (mm)');
        colorbar

        subplot(1,2,2)
        imagesc(tms,tms,latsums)
        hold on
        plot([tms(itend) tms(itend)],[tms(itstart) tms(itstart)],'.w','MarkerSize',12)
        hold off
        set(gca,'DataAspectRatio',[1 1 1])
        xlabel('end time (ms)')
        ylabel('start time (ms)')
    end

    slope = (tms(itend)-tms(itstart))/(ilatend-ilatstart);
    speed = mean(diff(latmm))/slope;

    avg = nanmean(latsums(:));
    stdd = nanstd(latsums(:));
    peak = latsums(itstart,itend);

    out.avg       = avg;
    out.std       = stdd;
    out.peak      = peak;
    out.speed     = speed;
    out.itstart   = itstart;
    out.itend     = itend;
    out.ilatstart = ilatstart;
    out.ilatend   = ilatend;
    out.latsums   = latsums;
    out.plane     = plane;
    out.latmm     = latmm;
    out.tms       = tms;
end

%==========================================================================

function latsums=FindLatSumsUsingEndpoints(plane,ilatstart,ilatend)

    [nlats,ntimes]=size(plane);

    interpIndx=zeros(ntimes-1,nlats);         % get indices and fractions for
    interpFrac=zeros(ntimes-1,nlats);         %   interpolation at lat positions

    ilats=ilatstart:ilatend;

    for idiff=1:ntimes-1
        tvals=(idiff-1)/(length(ilats)-1)*(ilats-ilats(1));
        interpIndx(idiff,ilats)=floor(tvals);
        interpFrac(idiff,ilats)=1-(tvals-interpIndx(idiff,ilats));
    end

    latsums    =zeros(ntimes-1,ntimes-1);

    maxistart=ntimes-1;
    maxiend  =ntimes-1;

    startTimeIndex=1;
    endoffset = 0;


        % the routine below is must faster than the following with 3 loops
    
%    for istart=startTimeIndex:maxistart
%        for iend=istart+endoffset:maxiend             % iend >= istart ==> c > 0
%
%            idiff=iend-istart+1;            % +1 for matlab numbering
%
%            interpvals=zeros(1,nlats);              % zero array for interpolated values
%            for ilat=ilatstart:ilatend
%                idx0=interpIndx(idiff,ilat)+istart;
%                frac=interpFrac(idiff,ilat);
%                interpvals(ilat)=plane(ilat,idx0)*frac + plane(ilat,idx0+1)*(1-frac);
%            end
%            latsums1(istart,iend)=sum(interpvals);       % get Radon sum over lat positions
%        end 
%    end


    for istart=startTimeIndex:maxistart
        for iend=istart+endoffset:maxiend             % iend >= istart ==> c > 0

            idiff=iend-istart+1;            % +1 for matlab numbering

            idx0vals=interpIndx(idiff,ilats)+istart;
            fracs=interpFrac(idiff,ilats);
                
            indices1 = sub2ind(size(plane),ilats,idx0vals);
            indices2 = sub2ind(size(plane),ilats,idx0vals+1);

            latsums(istart,iend) = sum(plane(indices1).*fracs+plane(indices2).*(1-fracs));  % sum over lat pos
        end
    end
end

%==========================================================================

%function latsums = LimitMaxSpeed(plane,latsums,lat,t,maxspeed)

function latsums = LimitMaxSpeed(plane,latsums,lat,t,maxspeed,ilatstart,ilatend)

        % input lat in mm, t in ms

%    [nlats,ntimes]=size(plane);
%
%    ilatstart=1;
%    ilatend=nlats;

    dlat = mean(diff(lat));
    dt   = mean(diff(t));

    deltait = floor(dlat*(ilatend-ilatstart)/(maxspeed*dt));    % dt in ms
    
    [nsums1,nsums2]=size(latsums);
    if nsums1~=nsums2
        error('error')
    end
    
    col=(1:nsums1)';
    itstartmat=repmat(col,1,nsums2);
    itendmat=itstartmat';
    
    indices=find(itendmat<=itstartmat+deltait);
    latsums(indices)=NaN;
end

%==========================================================================

function [itstart,itend]=FindGlobalPeak(latsums)

%    dim=size(latsums,1);
%    [maxval,maxpos]=max(reshape(latsums,1,dim*dim));
%    itend=floor(maxpos/dim)+1;
%    itstart=maxpos-(itend-1)*dim;

    [maxval maxidx] = max(latsums(:));
    [itstart itend] = ind2sub(size(latsums),maxidx);
end

%==========================================================================
