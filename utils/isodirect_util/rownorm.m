function out = rownorm(in)
%
% P. Spudich 9/2/10 v1.0
%
% rownorm calculates the norm of each row of a matrix
%
% inputs
%   in - a matrix
% 
% outputs:
%   out - a column vector, out(i) = sqrt( sum_on_j( (in(i,j).^2  )
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
int = in';
out = sqrt( sum( int.^2) );
out = out'; 

return