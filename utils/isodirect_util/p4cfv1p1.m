function out = p4cfv1p1( in )
%
%  usage:  out = p4cfv1p1( in );
%
% P. Spudich, 9/2/10 v2.1
%
% Function to calculate isochrone directivity parameters D and s for single
% and multisegment faults for a bunch of sites (with closest points on
% fault). This is based on the formulation given by B. Chiou in "P4CF
% v1.1.doc" dated 8/22/10, and the specific equations of this algorithm are
% derived in "PRCFv1.1_derivation.ai".
%
% An apology.  I tried to write the routine so that there were no loops and
% all calculations could be done as vector operations, but I failed to
% figure out how to accomplish that for the last step.  This routine would
% be much simpler if I had just used loops willy-nilly. 
%
% NOTE:  All coordinates should be given in a right handed coord system
% with the third coord being up, e.g. [east north UP]
%
% inputs: 
%
% in - (structure) with fields
%     .QQ - matrix of dimension 4 x 3*np, where np is the number of polygons
%       comprising the fault surface.  Each polygon is a quadrilateral.
%       QQ is the field pol.QQ in polygons?.mat.  
%       QQ consists of submatrices QQ1, QQ2, ... QQn, i.e. 
%       QQ = [QQ1 QQ2 ... QQnp], where size(QQi) = [4,3] and QQi(k,:) is the
%       [east north up] coord (km or m, not dg) of the k'th point of polygon
%       i.  For each quadrilateral points QQi(1,:), QQi(2,:), ... QQi(4,:)
%       are ordered sequentially around the quadrilateral, and QQi(1,:) and
%       QQi(2,:) for all i must be at the same depth, at the top of the
%       polygon.  QQi(1,:) corresponds to point Pi in Spudich and Chiou, 
%       and QQn(2,:) correponds to point Pn+1.  It is not necessary that 
%       the fault strike points in the direction from QQi(1,:) to QQi(2,:),
%       but the polygons must be ordered such that each strikes in the same
%       general direction. (NOTE 7/28/10: I think that it IS necessary that
%       the fault strike points from QQi(1,:) to QQi(2,:).)
%     .C - matrix of coordinates of points on fault closest to each site. 
%       C(1:3,j) = [east north up] coordinates of point on fault closest to
%       site j.  dim(C) = number of sites ns x 3
%     .c - vector of polygon number, dim (ns x 1);  ic(i) = polygon number
%       on which C(1:3,i) lies. 
%     .H - three-element vector of hypocenter coords in E, N, >UP< coord 
%       system ( may be row or column vector).  This routine does not check
%       to verify that your hypocenter actually lies on the fault surface.
%     .h = scalar, polygon number on which the hypocenter lies
%
% outputs: 
%   out - structure with fields
%     .D - isochrone directivity parameter D, dim(D) = ns x 1
%     .s - isochrone directivity parameter s, dim(s) = ns x 1
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
wantprint = 0;

QQ = in.QQ; 
C = in.C;%
c = in.c;%
H = in.H;%
h = in.h;% 

fid = 1; % terminal
if wantprint
    fprintf(fid, '\n')
    fprintf(fid, '\n')
    fprintf(fid,'On input to p4cfv1p1 ... \n')
    fprintf(fid,'\n     Clst pt to each site: %15.11f %15.11f %15.11f', C')
    fprintf(fid,'\n     Seg # of closest point: %3.0f', c)
    fprintf(fid,'\n     Hypocenter: %15.11f %15.11f %15.11f', H')
    fprintf(fid,'\n     Seg # of hypocenter: %3.0f', h)
end

[ns,ncomp] = size(C); % ns is number of sites (closest points)

ztop = QQ(1,3);%
zbot = QQ(4,3);% 
zh = H(3);% 

% limited sanity checks
if ncomp~=3 ; error(['p4cfv1p1: C matrix does not have 3 columns']) ; end
if ztop>0 ; error(['p4cfv1p1: ztop is positive, should be negative']) ; end
if zbot>0 ; error(['p4cfv1p1: zbot is positive, should be negative']) ; end
if zh>0 ; error(['p4cfv1p1: zh is positive, should be negative']) ; end
if zbot >= ztop
    error(['p4cfv1p1: zbot >= ztop, should be reverse'])
end
if zh < zbot || zh > ztop
    error(['p4cfv1p1: zh is not between ztop and zbot'])
end
% end sanity checks

% make Hp a row vector
Hp = (H(:))';

[nrq, ncq] = size(QQ); 
np = ncq / 3; % np is number of polygons comprising the fault

% set up matrix of points on auxiliary line at hypocentral depth 
r = (zh-ztop) / (zbot-ztop);% 
A1p = QQ(1,:) + r * ( QQ(4,:) - QQ(1,:) ) ; % dim A1p = 1 x 3*np
A2p = QQ(2,:) + r * ( QQ(3,:) - QQ(2,:) ) ; % dim A2p = 1 x 3*np
A1 = reshape( A1p, 3, np); 
A2 = reshape( A2p, 3, np);
% lengths of segments on auxiliary line
Lp = sqrt( sum( (A1-A2).^2  ) );% % dim(Lp) = 1 x np
% We need add only the last 3 columns of A2p to complete a matrix of all
% auxiliary line points
A = ( reshape( [ A1p A2p(1,(3*np-2):(3*np)) ] , 3, np+1 ) )';% 
%dim A = np+1 x 3

if wantprint
    fprintf(fid,'\n    Auxiliary point: %15.11f %15.11f %15.11f', A')
end

% sAp(i) = length along auxiliary line to point A(:,i)
sAp = cumsum( [0 Lp] );% % dim(sAp) = 1 x np+1
% remove last element; we don't need length to end of fault
sAp(np+1) = [];% % dim(sAp) = 1 x np

% along-aux-line distance to hypocenter
sh = sAp(h) + norm( Hp - A(h,:) );% % sh is a scalar

% must find Pstar, point on aux line of quadrilateral c closest to closest
% point C(i,:) for all i.

% the auxiliary line for quad c goes from A(c,:) to A(c+1,:)
in.A = A(c,:); % 
in.B = A(c+1,:); %
in.X = C; 
out = ClosestPointToLineSegment(in);

Pstar = out.C;%
sPstar = (sAp(c))' + rownorm( Pstar - A(c,:)  );% % dim(sPstar) = ns x 1 

% Determine whether the hypocenter polygon # is <= closest point polygon #
hlec = h <= c;% % dim(hlec) = ns x 1
hgtc = 1 - hlec;% 

% sG(i) and sL(i) are the s coordinates of the lesser and greater endpoints
% of the D path for site i
sG = hgtc .* (sh*ones(ns,1)) + hlec .* sPstar;%  % dim(sG) = ns x 1
sL = hlec .* (sh*ones(ns,1)) + hgtc .* sPstar;%  % dim(sL) = ns x 1

% s ratio: rs(i,j) = (sAp(j)-sL(i)) / (sG(i)-sL(i))
rs = (ones(ns,1)*sAp - sL*ones(1,np)) ./ (sG*ones(1,np) - sL*ones(1,np));%
% dim(rs) = ns x np

zC = C(:,3);% % zC(i) is z coord of ith closest point, dim(zC) = ns x 1
% zL and zG are z coords of lesser and greater end points of D path
zG = hgtc.*(zh*ones(ns,1)) + hlec.*zC;% % dim(zG) = ns x 1
zL = hlec.*(zh*ones(ns,1)) + hgtc.*zC;% % dim(zL) = ns x 1

% The 'O' point is the point where the D path crosses a polygon
% boundary. 
% zO(i,j) is the z coord of the 'O' point for the jth polygon for the ith
% site (closest point)
zO = zL*ones(1,np) + ( (zG-zL)*ones(1,np) ) .* rs;% % dim(zO) = ns x np
zr = (zO-ztop) / (zbot-ztop);% 

% To avoid the use of 3 dimensional arrays, I handle each coordinate
% separately

P1Ep = QQ(1, 1:3:3*np);% % dim(P1Ep) = 1 x np;% east coord of point P1 of polygon
P1Np = QQ(1,(1:3:3*np)+1);% % north coord of point P1 of polygon
P1Up = QQ(1,(1:3:3*np)+2);% % north coord of point P1 of polygon

P4Ep = QQ(4,1:3:3*np);% % dim(P4Ep) = 1 x np;% east coord of point P4 of polygon
P4Np = QQ(4,(1:3:3*np)+1);% % north coord of point P4 of polygon
P4Up = QQ(4,(1:3:3*np)+2);% % north coord of point P4 of polygon

% OE(i,j) is the east coord of the 'O' point for the jth polygon for the ith
% site (closest point)
OE = ones(ns,1)*P1Ep + (ones(ns,1)*(P4Ep-P1Ep)).*zr ;% % dim(OE) = ns x np
ON = ones(ns,1)*P1Np + (ones(ns,1)*(P4Np-P1Np)).*zr ;%
OU = ones(ns,1)*P1Up + (ones(ns,1)*(P4Up-P1Up)).*zr ;%

% The D path for site (closest point) i starts in polygon number l(i) and
% ends in polygon number g(i), where l(i) <= g(i). 
l = h*hlec + c.*hgtc ;% % dim(l)=ns x 1; l(i)= h if h<=c(i), = c(i) if c(i)<h
g = h*hgtc + c.*hlec ;% 

% L(i,:) contains the (e,n,u) coords of the end of the D path in segment
% l(i) for site i. 
% G(i,:) contains the (e,n,u) coords of the end of the D path in segment
% g(i) for site i. 
% The D path goes from point L to point G (although the sense of direction
% is irrelevant; it could be said to go from G to L).

L = hlec*Hp + (hgtc*ones(1,3)).*C ;% % dim(L) = ns x 3
G = hgtc*Hp + (hlec*ones(1,3)).*C ;% 

% Append point G to the end of the O points defining the D path
OE = [OE G(:,1)] ; % dim(OE) = ns x np+1
ON = [ON G(:,2)] ;
OU = [OU G(:,2)] ;

% Now I want to overwrite all O points outside the [L G] interval with L or
% G, so that a row of OE looks like [L L OE(3) OE(4) G G G], but
% regrettably I must do this with loops.  Ugh. 

for i = 1:ns
    for j = 1:l(i)
        OE(i,j) = L(i,1); 
        ON(i,j) = L(i,2); 
        OU(i,j) = L(i,3); 
    end

    for j = (g(i)+1) : (np+1)
        OE(i,j) = G(i,1); 
        ON(i,j) = G(i,2); 
        OU(i,j) = G(i,3); 
    end
end

OE;% 
ON;%
OU;%

% isochrone parameter D; dim(Dp) = 1 x ns
Dp = sum( sqrt(  diff(OE').^2 + diff(ON').^2 + diff(OU').^2  ) ) ;%

out.D = Dp';
out.s = abs( sPstar - sh);% 

out.D;%
out.s;%

return

        
        
        









