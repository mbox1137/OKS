function [e,n,u] = xyz2enu(refLat, refLong, refH, X, Y, Z, ~)
  %
  % WGS-84 to ENU,  Step 2: Convert ECEF to ENU (local east, north, up)
  %
  % http://wiki.gis.com/wiki/index.php/Geodetic_system
  %
  % окт€брь-декабрь 2015, var@mail.spbstu.ru
  %
  % find reference location in ECEF coordinates
  if nargin < 7
      % WGS-84:
    [Xr,Yr,Zr] = llh2xyz(refLat,refLong, refH);
  else
      % сфера:
    [Xr,Yr,Zr] = llh2xyz(refLat,refLong, refH, 1); 
  end
  
 %% lat,long - в градусах
  e = -sind(refLong).*(X-Xr) + cosd(refLong).*(Y-Yr);
  n = -sind(refLat).*cosd(refLong).*(X-Xr) - sind(refLat).*sind(refLong).*(Y-Yr) + cosd(refLat).*(Z-Zr);
  u = cosd(refLat).*cosd(refLong).*(X-Xr) + cosd(refLat).*sind(refLong).*(Y-Yr) + sind(refLat).*(Z-Zr);
 %% lat,long - в радианах 
%   e = -sin(refLong).*(X-Xr) + cos(refLong).*(Y-Yr);
%   n = -sin(refLat).*cos(refLong).*(X-Xr) - sin(refLat).*sin(refLong).*(Y-Yr) + cos(refLat).*(Z-Zr);
%   u = cos(refLat).*cos(refLong).*(X-Xr) + cos(refLat).*sin(refLong).*(Y-Yr) + sin(refLat).*(Z-Zr);