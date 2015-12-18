function [xc, inside] = quadcl1(y, Q)
%
% function based on notes 4.26.07.I.quadrilateral
%
% P. Spudich 5/7/07 v1.5
%
% 5/9/07: modified to call tricl2
%
% Given a point y and a quadrilateral with vertices at Q = [q1 ; q2 ; 
% q3 ; q4], where the points are sequenced around the periphery, 
% quadcl1 returns the point closest to y inside or on the 
% boundary of the quadrilateral.  Note, this routine does not require that
% the four corner points lie exactly on the same plane, but it returns the
% closest point xc assuming that the reason that the four corners might not
% lie in the same plane is simple numerical inaccuracy.  The algorithm is
% the following.  The quad is divided into two triangles, T1 given by
% q1,q2,q3 and T2 given by q3,q4,q1.  quadcl1 finds the closest point on
% T1. If the projection of y (yperp) lies on T1, then yperp is xc and the
% routine returns.  If not, quadcl1 does the same for T2.  If yperp is
% inside neither T1 or T2, quadcl1 chooses the closer of the closest points
% from T1 and T2.
%     The boundary between T1 and T2 is the
% line from q1 to q3.  Inaccuracies cause T1 to be not coplanar with T2, so
% there is a fold in the quad along the q1-q3 line.  This routine will only
% locate xc on the fold if y is exactly perpendicular to T1 or T2 at xc. 
%
% inputs:
%   Q - 4x3 matrix , Q = [q1 ; q2 ; q3 ; q4], where qi 
%       is a 3-component vector holding the coordinates of
%       the vertex of the quadrilateral. The components are given in some 
%       orthogonalright-handed coordinate system, preferrably N, E, Down.
%       Quadrilateral must be simply-connected!!
%   y - a 3-component vector holding the coords of some point
%
% outputs:
%   xc - a 3-component vector which is the point inside the quad or on
%      its boundary closest
%   inside = 1 if y is normal to the surface of the quad. (yperp inside)
%          = 0 otherwise
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

% check quadrilateral has nonzero area.  If zero area, xc=any corner
v1 = Q(2,:) - Q(1,:); 
v2 = Q(3,:) - Q(2,:); 
v3 = Q(4,:) - Q(3,:); 
v4 = Q(1,:) - Q(4,:); 

if max(abs(v1'))<eps && max(abs(v2'))<eps && max(abs(v3'))<eps 
    xc = Q(1,:);
    inside = 0;
%   ======
    return
%   ======
end


% check quadrilateral is simply connected.  If cross products all point in
% the same direction, everything is okay. 
v12 = cross(v1,v2);
v23 = cross(v2,v3);
v34 = cross(v3,v4);
v41 = cross(v4,v1); 

m1 = v12*v23';
m2 = v23*v34';
m3 = v34*v41';
m4 = v41*v12';

if (sign(m1)~=sign(m2)) || (sign(m2)~=sign(m3)) || (sign(m3)~=sign(m4)) ...
        || (sign(m4)~=sign(m1)) 
%    error('quadcl1: quadrilateral is a triangle or not simply connected')
    display('WARNING quadcl1: quadrilateral is a triangle or not simply connected')
    xc = NaN*ones(1,3);
    inside = 0;
    %======
    return
    %======    
end

y = y(:)'; % make it a row vector

% triangle 1
[xc1, inside1] = tricl2(y, Q(1,:), Q(2,:), Q(3,:) ); 

% check whether yperp, projection of y onto T1, is inside or on bdry of T1 
if inside1
    xc = xc1;
    inside = 1; 
    %======
    return
    %======
    
else
    % yperp is outside triangle 1.  Is it inside triangle 2? 
    [xc2, inside2] = tricl2(y, Q(3,:), Q(4,:), Q(1,:) ); 
    if inside2
        xc = xc2;
        inside = 1; 
        %======
        return
        %======
    else
        
        % yperp is outside both triangles.  
        inside = 0;
        
        % xc is xc1 or xc2, whichever is closer to y
        d1 = sum( (y-xc1).^2 );  
        d2 = sum( (y-xc2).^2 );  
        
        if d1<d2
            xc = xc1;
        else
            xc = xc2; 
        end
        
    end
end

return
        
        
