% Sector scan mode:
%
%       [out,ax,lat]=scan_convert('sector',in,min_phi,span_phi,apex,flag,fs,varargin);
% 
%       where   'sector' specifies sector scan mode,
%               in is the input data
%               min is the minimum (leftmost) sector angle (degrees)
%               span is the angle spaned by the sector (degrees)
%               apex is the distance in cm from the 'ducer face to the 
%                       center of curvature of the 'ducer.  Should be negative
%               flag specifies interpolation method (1 = nearest neighbor,
%                       2 = linear interpolation laterally)
%               fs specifies the axial sampling frequency
%               varargin{1} is an optional scalar argument to set the background value
%               varargin{2} is an optional vector argument for region selection,
%                       [ax_min, ax_max, ax_inc, lat_min, lat_max, lat_inc]
%
%
% Version 0.9
% 7/23/03  Stephen McAleavey, Duke University BME
%
% Revision History:
% 02/03/04 Brian Fahey, added flag variable to allow user to choose which geometric transformation to use
%       flag = 1 uses nearest neighbor
%       flag = 2 uses linear interpolation between theta angles
% radial oversampling makes nearest neighbor approx good for most pictures, although theta artifacts 
% can be seen if look closely. Interp method about 3x slower to implement. Recommend doing bulk processing
% with nearest neighbor, then making conference/paper images w/ interp
%
% 02/06/04 sam, added optional vector argument for region selection
%          [xmin xmax ymin ymax inc] all in meters
%          selects the rectangular window to zoom in on, as well as the grid 
% 	   spacing.  You can pick these numbers based on the axes of an 
%          un-zoomed image.  If this argument isn't supplied, the default
%          is the bounding rectangle of the sector. 
%
% 2011/05/04 Dongwoon Hyun
%      Mexified code.
%
% In order to compile this code, use the following command:
%
% mex -largeArrayDims scan_convert.cpp