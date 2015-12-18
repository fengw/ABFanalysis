function output = msegcl5(y, QQ, ih)
%
% function based on notes 4.26.07.I.quadrilateral
% P. Spudich 5/7/07 v5.2
%
% 9/17/10: msegcl5 is msegcl3 modified to resolve the ambiguity of segment
% number when the closest point is on the boundary between segments or when
% points on two different segments are equidistant from the site.  If
% this case arises, if ih is input, then the segment number of the 
% closest point is taken to be
% the segment number which is closer to ih (typically the hypocentral 
% segment number).  If ih is not input, then the lower segment number is 
% taken. In all cases, distance differences less than a hardwired tolerance
% are considered equal. 
% 9/15/10: msegcl3 is msegcl2 modified to resolve the ambiguity of segment
% number when the closest point is on the boundary between segments or when
% points on two different segments are equidistant from the site.  If
% this case arises, if refpt is input, then the segment number of the 
% closest point is taken to be
% the segment which is closer to an input reference point (typically the
% hypocenter).  If refpt is not input, then the lower segment number is 
% taken. In all cases, distance differences less than a hardwired tolerance
% are considered equal. 
%
% 10/3/07: msegcl2 is msegcl1 modified to output a structure containing the
% index of the closest segment as well as the outputs from msegcl1
%
% Given a matrix of points y and a complicated surface QQ consisting of
% quadrilaterals  
% msegcl1 returns the point inside or on the boundary of the surface
% closest to each y .   
%
% inputs:
%   QQ - matrix dimension 4 x (3*nsegs), specifying coords of corners of 
%       quadilaterals.  A quadrilateral is specified by a 4x3 matrix Q
%       where Q(i,1:3) contains the east, north, and up coordinates of the
%       i=1, 2, 3, 4th corner point of the quadilateral.  
%       QQ = [Q1 Q2 ... Qnsegs], where Qk is the Q matrix defining the kth
%       quadrilateral.  
%
%   y - a Nx3 matrix in which y(j,1), y(j,2), and y(j,3) are the
%   coordinates of the jth point for which the closest point is desired. 
%
%   ih - (optional) scalar - this is the segment number of a  
%       reference point (typically the hypocenter).  If
%       for some j the closest point to y(j) lies on the boundary between
%       two segments, this routine returns in iclosest the number of the
%       segment which is closer to ih.
%
% outputs:
%   output - (structure) containing fields: 
%   .xc - a Nx3 matrix in which xc(j,1), xc(j,2), and xc(j,3) are the
%       coordinates of the point on the complicated surface closest to 
%       y(:,j)  
%   .distance = a Nx1 vector of the distance from xc(k) to y(k)
%   .iclosest - segment number of the segment containing the closest
%       segment.  If the closest point is on the boundary between two segments,
%       the lower segment number is chosen. 
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

distol = 1e-10 ; % distance differences < distol treated as identical

[ny, ncy] = size(y);
[nrq, ncq] = size(QQ); 

if ncy~=3; error('msegcl3: y does not have 3 columns') ; end
if mod(ncq,3)~=0 ; error('msegcl3: number of columns of QQ not factor of 3'); end
if nrq~=4; error('msegcl3: QQ does not have 4 rows') ; end

nsegs = ncq / 3; 

distance = NaN*ones(size(y));

segdist2 = NaN*ones(nsegs,1) ; % will hold squared distances from y(:,j) to each segment
xc = NaN * ones(size(y)); 
imin = NaN * ones(ny,1); 

fid = 1;

for iy = 1:ny
        
    xcseg = NaN*ones(nsegs,3); % clear holding area for xc from each segment
    for iseg = 1:nsegs
        
        cols = (iseg-1)*3 + (1:3); 
        Q = QQ(1:4, cols); 
        
        [xx, inside] = quadcl1(y(iy,:), Q); 
        xcseg(iseg,:) = xx; 
                
        % squared distance from xcseg to y(iy)
        segdist2(iseg) = sum( (xx-y(iy,:)).^2   );
    end
   
   % pre 9/15/10, msegcl2: 
   % index of segment having min squared distance
   % imin(iy) = min( find( segdist2 == min(segdist2) )  ); 
    
   
   segdist = sqrt(segdist2); 
   distdif = abs(segdist - min(segdist));
      
   % segment numbers having distances < distol greater than min distance
   iseglttol = find(distdif < distol); 
      
   if length(iseglttol) == 1
       % unique minimum distance
       imin(iy) = min( find( segdist2 == min(segdist2) )  ); % msegcl2 result
   else
       if nargin == 2
           % choose the lower segment number
           imin(iy) = iseglttol(1); 
       else
           % find the segment number closer to ih
           jj = find( abs(iseglttol-ih) == min(abs(iseglttol-ih))); 
           if length(jj)> 1
               error(['msegcl5: for iy=' int2str(iy) ' more than one '...
                   'segment closest to hypocenter'])
           end
           imin(iy) = iseglttol(jj);
       end
   end

    xc(iy,:) = xcseg(imin(iy),:) ; % this one is the global min
    
end

distance = ( sum(  ((xc-y).^2)'  )  )' ; 
distance = sqrt(distance); 

output.xc = xc;
output.distance = distance;
output.iclosest = imin; 


return
        
        
        
        

