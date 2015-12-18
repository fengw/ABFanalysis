function output = radpatv7(in)
%
% usage: output = radpatv7(in);
%
% P. Spudich 9/29/07 v7.0
%
% function to evaluate the radiation pattern according to the equations of
% Spudich & Chiou Appendix A v7, as of September 27, 2007.  Uses the
% coordinate system with u increasing along strike, d increasing downward,
% and u cross t = d.  All coordinates, stations and hypocenter, must be
% with respect to the same origin.  Note, I use d as the downward axis
% rather than z, because in Appendix A z is zero at the station elevation,
% and I want to preserve the ability to handle variable station elevations.
% 
%
% inputs:  
%   in - (structure) containing fields: 
%       .utds = (matrix ns x 3): 
%           utds(j,1:3) = [u-coord t-coord d-coord] of station j, where d is
%           the depth of the station below the origin of coords, I.E. d is 
%           NEGATIVE for stations above the origin
%       .utdh = (vector 1 x 3): 
%           utdh(1:3) = [u-coord t-coord d-coord] of hypocenter with
%           respect to same origin as stations' coords.
%       .dip = fault dip in degrees
%       .rake = fault rake in degrees
%       .watlev = (scalar) between 0 and 1 used as waterlevel on magnitude
%           of radiation pattern 
%   NOTE - utds and utdh should all contain either exclusively normal
%   cartesian coords (for the case of single segment ruptures) as
%   illustrated in Spudich & Chiou Figure A1, or exclusively generalized
%   coordinates.  
%
% outputs: 
%   output - (structure) containing fields: 
%       .radpat = (matrix ns x 2):
%           radpat(j,1:2) = [Ru Rt] for station j, WITHOUT water level
%       .radpatmag = (vector ns x 1): 
%           = max( watlev, sqrt( Ru.^2 + Rt.^2 )), magnitude of radiation 
%           pattern including water level
%
%
%--------------------------------------------------------------------------
% DISCLAIMER: 
%
% Although this program has been used by the USGS, no warranty,
% expressed or implied, is made by the USGS or the United States Government
% or any authors as to the accuracy and functioning of the program and
% related program material nor shall the fact of distribution constitute
% any such warranty, and no responsibility is assumed by the USGS in
% connection therewith.
%
% COPYRIGHT INFO
% This m-file is a U.S. Government work and is not subject to copyright. It
% may be freely distributed. 
%--------------------------------------------------------------------------
%

% ns should be number of stations
[ns,ncol] = size(in.utds); 


% sanity checks
if ncol ~= 3 ; error('radpatv7: in.utds does not have 3 columns'); end

[nr,nc] = size(in.utdh); 
if nr~=1 ; error('radpatv7: in.utdh does not have 1 row.  Transpose?'); end
if nc~=3 ; error('radpatv7: in.utdh does not have 3 columns.'); end

if in.watlev < 0 || in.watlev > 1 
    error(['radpatv7: illegal in.watlev value is ' num2str(in.watlev)])
end


% here utdprime is [uprime(:,1) tprime(:,1) -zh(:,1)] because 
% u' = us-uh, t'=ts-th, but zh = dh - ds !!
utdprime = in.utds - ones(ns,1)*in.utdh; 
% must flip sign of z component to be >0 when hypocenter is below station
utdprime(:,3) = -utdprime(:,3); 

% epicentral distance
R = sqrt( utdprime(:,1).^2 + utdprime(:,2).^2 );     

% replace zero values with NaN
Ris0 = find(R==0); 
R(Ris0) = NaN; 

% hypocentral distance
rh = sqrt( R.^2 + utdprime(:,3).^2 ); 

snf = R./rh;
cnf = utdprime(:,3)./rh;

sp = utdprime(:,2)./R;
cp = utdprime(:,1)./R;

cd = cos(in.dip * pi/180); 
sd = sin(in.dip * pi/180);
c2d = cos(2 * in.dip * pi/180); 

cl = cos(in.rake * pi/180); 
sl = sin(in.rake * pi/180);


nr = snf.*sp*sd + cnf*cd; 
nb = cnf.*sp*sd - snf*cd;
nc = cp*sd; 

sr = -cl*snf.*cp + sl*cd*snf.*sp - sl*sd*cnf;
sb = -cl*cnf.*cp + sl*cd*cnf.*sp + sl*sd*snf;
sc = cl*sp + sl*cd*cp; 

bcoef = nr.*sb + nb.*sr;
ccoef = nr.*sc + nc.*sr;

Ru = -cp.*bcoef + sp.*ccoef; 
Rt = -sp.*bcoef - cp.*ccoef; 

% take care of the rows for which R=0
Ru(Ris0) = cd*cl;
Rt(Ris0) = -sl*c2d; 

output.radpat = [Ru Rt]; 
output.radpatmag = max( in.watlev, sqrt( Ru.^2 + Rt.^2 )   ); 

return
