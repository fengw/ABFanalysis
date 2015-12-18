function out = ClosestPointToLineSegment(in)
%
% usage: out = ClosestPointToLineSegment(in)
%
% P Spudich, 9/2/10, v1.0
%
% Finds coordinates C of the point on line segment A-B closest to some 
% point X.  A, B, and X may be N x 3 matricies, in which case C(i,:) is the
% point on the line segment from A(i,:) to B(i,:) closest to X(i,:).  
%
% inputs: 
%
%   in - structure, with fields
%       .A - matrix, A(i,j), j=1:3, are the x,y,z coords of one end point
%           of the i'th line segment.  i = 1:n, n=1 or n=N allowed. 
%       .B - matrix, B(i,j), j=1:3, are the x,y,z coords of the other end point
%           point of the i'th line segment.  i = 1:n, n=1 or n=N allowed. 
%       .X - matrix, X(i,j), i=1:N, j=1:3, are the x,y,z coords of some 
%           points not (necessarily) on the line segment(s), for which the 
%           closest point on A-B are sought. NOTE N here and n in A, B. If 
%           N > 1 and n = 1, all X will be compared to the single line
%           segment A(i,:) to B(i,:).
%
%   out - structure, with fields
%       .C - matrix, C(i,j), j=1:3, are the x,y,z coords of the point on A 
%           - B closest to X(i,:)
%
%---------------------
% DISCLAIMER: 
%
% Although this program has been used by the USGS, no warranty,
% expressed or implied, is made by the USGS or the United States Government
% or any authors as to the accuracy and functioning of the program and
% related program material nor shall the fact of distribution constitute
% any such warranty, and no responsibility is assumed by the USGS in
% connection therewith.

% COPYRIGHT INFO
% This m-file is a U.S. Government work and is not subject to copyright. It
% may be freely distributed. 
%---------------------

X = in.X;
A = in.A;
B = in.B;

[N,nc] = size(X); 
if nc ~= 3 ; error('ClosestPointToLineSegment: X does not have 3 columns'); end

[na,nc] = size(A); 
if nc ~= 3 ; error('ClosestPointToLineSegment: A does not have 3 columns'); end
[nb,nc] = size(B); 
if nc ~= 3 ; error('ClosestPointToLineSegment: B does not have 3 columns'); end

if na~=nb ; error('ClosestPointToLineSegment: A and B do not have same # of rows'); end

if na~=1 && na~=N ; error('ClosestPointToLineSegment: A and B do not have same # of rows as X'); end

n = na; 

if n==1
    A = ones(N,1) * A; % make it N x 3, same dim as X
    B = ones(N,1) * B; % make it N x 3, same dim as X
end

l = sqrt( sum( (B-A).^2, 2 ) );% % l is distance from B to A, dim(l) = N x 1
r = sqrt( sum( (X-A).^2, 2 ) );% % r is distance from X to A, dim(r) = N x 1

costheta = sum( (X-A).*(B-A)    ,2 ) ./r ./l ; % dim costheta = N x 1
u = r .* costheta ; 
ul = u./l ;%
v = min(max(ul,0),1);% 

out.C = A + [v v v] .* (B-A);% % dim C = N x 3

return

