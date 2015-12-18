function [xc, inside] = tricl2(y, t1, t2, t3)
%
% function based on notes 4.26.07.I.quadrilateral
% tricl2 is tricl1 with inside-triangle algorithm fixed
%
% P. Spudich 5/9/07 v2.3
%
% Given a point y and a triangle with vertices at t1, t2, and t3, 
% tricl1 returns the point closest to y inside or on the 1-2 or 2-3
% boundary of the triangle.  (Note, it is anticipated that this routine
% will be called by a function seeking the point on a quadrilateral closest
% to y, and the 1-3 boundary of this triangle is an interior diagonal of
% the quadrilateral, for which the closest point is not needed.)
%
% 9/20/10:  distol added and if statement changed to ensure that points
% yperp slightly outside the t1-t3 boundary ( the bisector of the
% quadrilateral) are considered inside the triangle. 
%
% inputs:
%   t1, t2, t3 - each is a 3-component vector holding the coordinates of
%       the vertex of a triangle. The components are given in some orthogonal
%       right-handed coordinate system, preferrably N, E, Down.
%   y - a 3-component vector holding the coords of some point
%
% outputs:
%   xc - a 3-component vector which is the point inside the triangle or on
%      the 1-2 or 2-3 boundary closest to y
%   inside = 1 if yperp is inside or on the boundary of the triangle
%          = 0 otherwise
%
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

inside = 0; 
[nr,nc] = size(y); % so I remember whether it's a row or column vector

% make them column vectors 
y = y(:);
t1 = t1(:); 
t2 = t2(:); 
t3 = t3(:); 

v1 = t1 - t2;
v3 = t3 - t2;
v1xv3 = cross(v1,v3); 

% normal to triangle
nhat = v1xv3 / sqrt( v1xv3' * v1xv3  ); 

ymt2 = y - t2; 

% projection of y onto plane of triangle
yperp = ymt2 - (ymt2'*nhat)*nhat + t2; 

% determine whether yperp is inside or outside triangle

% find coefs such that (yperp-t2) = a1*v1 + a2*v3.  If a1>=0 and a2>=0 and 
% 0 <= a1+a2 <= 1, then yperp is inside the triangle.
E = [v1 v3]; 
a = E \ (yperp-t2); 

% if a(1)>=0 && a(2)>=0 && (a(1)+a(2)) <= 1 % pre-9/20/10
% next 2 lines added, 9/20/10, to ensure that points just slightly outside 
% the t1-t3 boundary are considered inside. 
distol = 1e-10; 
if a(1)>=0 && a(2)>=0 && (a(1)+a(2)) <= 1 + distol
    xc = reshape(yperp,nr,nc); % make it a row or col vector like y
    inside = 1; 
%   ======
    return
%   ======
end
    
% yperp is outside the triangle.  Find its projection onto the 1-2 and 2-3
% edges, and figure out which of them is closer to yperp

% phat1 is the normal to edge 1-2 lying in the plane of the triangle
phat1 = cross(nhat,v1) / sqrt(v1'*v1); 
% phat3 is the normal to edge 3-2 lying in the plane of the triangle
phat3 = cross(nhat,v3) / sqrt(v3'*v3); 

ypmt2 = yperp - t2;
% p1 is projection of yperp onto the line t2+scalar*v1
p1 = ypmt2 - (ypmt2'*phat1)*phat1 + t2; 
p3 = ypmt2 - (ypmt2'*phat3)*phat3 + t2; 

% p1 or p3 could be anywhere on the infinite lines along the two edges, but
% the closest point must lie within the vertices.  Find the closest point
% on both edges, and then select the closer of these two.

% r1 is vector from t2 to p1
r1 = p1 - t2;
r3 = p3 - t2; 

% closest point on edge between t2 and t1
r1v1 = r1'*v1; 
if r1v1 <= 0
    xc1 = t2; 
elseif r1v1 >= v1'*v1
    xc1 = t1;
else
    xc1 = p1;
end

% closest point on edge between t2 and t3
r3v3 = r3'*v3; 
if r3v3 <= 0
    xc3 = t2; 
elseif r3v3 >= v3'*v3
    xc3 = t3;
else
    xc3 = p3;
end

% choose the closer of the two
r1sq = sum( (y-xc1).^2 ); 
r3sq = sum( (y-xc3).^2 ); 

if r1sq<=r3sq
    xc = reshape(xc1,nr,nc); % make it a row or col vector like y
else
    xc = reshape(xc3,nr,nc); % make it a row or col vector like y
end

%=====
return
%=====
