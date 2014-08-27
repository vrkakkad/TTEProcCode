function [arfidata axial lat] = Scan_Convert_Sector(arfidata0,sc)
[arfidata axial lat ] = scan_convert('sector',arfidata0,sc.min_phi,sc.span_phi,sc.apex,2,sc.fsiq,NaN,[sc.axialmin-sc.apex*1e-2/2 sc.axialmax-sc.apex*1e-2/2 sc.axialinc sc.latmin sc.latmax sc.latinc]);
axial = axial*10;
lat = lat*10;