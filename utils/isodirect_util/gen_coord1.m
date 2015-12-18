function [U,T] = gen_coord1(QQ,Xm)
%
% usage: [U,T] = gen_coord1(QQ,Xm)
%
% P. Spudich 9/29/07, v1.2
%
% function for calculating generalized coordinates U and T for a set of
% stations and a fault geometry. 
%
% Based on notes 9.23.07.I.Implementing_Gen_Geom.ai
%
%
% inputs: 
%   QQ - matrix of dimension 4 x 3n, where n is the number of polygons
%       comprising the fault surface.  Each polygon is a quadrilateral.
%       QQ is the field pol.QQ in polygons?.mat.  
%       QQ consists of submatrices QQ1, QQ2, ... QQn, i.e. 
%       QQ = [QQ1 QQ2 ... QQn], where size(QQi) = [4,3] and QQi(k,:) is the
%       [east north up] coord (km or m, not dg) of the k'th point of polygon
%       i.  For each quadrilateral points QQi(1,:), QQi(2,:), ... QQi(4,:)
%       are ordered sequentially around the quadrilateral, and QQi(1,:) and
%       QQi(2,:) for all i must be at the same depth, at the top of the
%       polygon.  QQi(1,:) corresponds to point Pi in Spudich and Chiou, 
%       and QQn(2,:) correponds to point Pn+1.  It is not necessary that 
%       the fault strike points in the direction from QQi(1,:) to QQi(2,:),
%       but the polygons must be ordered such that each strikes in the same
%       general direction.
%
%   Xm - matrix of station locations, dimension Xm = 3 x nstn, where nstn
%       is the number of stations.  Xm(:,j) is the [east ; north ; up }
%       coord (in km or m, not degrees) of the jth station. 
%
% outputs: 
%   U, T - column vectors of length nstn, containing the generalized U and
%       T coords of each station.  U(i) is the gen U coord of station i.
%
% NOTE - 9/29/07: I have added tuning parameter disttol, which is a
% distance tolerance (same units as station coords, presumably km) for
% testing equality of closest distances.  The problem in this algorithm is
% that the way the Dm matrix, distance from each station to each segment,
% is set up using Mm, numerical inaccuracy can cause calculated distances 
% to a point on the intersection between two segments of the fault trace to
% yield slightly different distance, so that a numerical test of equality
% to D during the formation of the dm matrix gets an answer dependent on 
% numerical noise.  I have changed the calculation of dm to say that if 
% abs(Dm(s,p)-D) < distol, then the segment p is a 'closest' point, of
% which there may be more than one. 
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
% DISTANCES LESS THAN DISTOL ARE ASSUMED TO BE EQUAL!!
distol = 1e-10; 
% disp(['in gen_coord1, distol = ' num2str(distol) ' *******'])

[nrqq,ncqq]= size(QQ);
n = ncqq/3;   % number of polygons

% sanity checks
if nrqq~=4 ; error('gen_coord1: QQ does not have 4 rows') ; end
if rem(n,1) ~= 0 ; error('gen_coord1: QQ has wrong number of columns') ; end

[nrxm, nstn] = size(Xm); 
if nrxm ~= 3 ; error('gen_coord1: Xm does not have 3 rows') ; end

% check tops of all polygons at same elevation

acol = 3*(1:n); % column containing elevation
% if std([QQ(1,acol)' ; QQ(2,acol)']) ~= 0    % original one
if std([QQ(1,acol)' ; QQ(2,acol)']) > distol
    std([QQ(1,acol)' ; QQ(2,acol)'])
    error('gen_coord1: tops of all quadrilaterals at different depths')
end

% check that polygons are connected.  Point 2 of polygon i should equal
% point 1 of polygon i+1
for p = 1:(n-1)
    firstcol = 1 + (p-1)*3; % east coord of P vectors
    jj = firstcol: (firstcol+2); 
    if max(abs( QQ(2,jj)-QQ(1,jj+3) )) > 0
        error('gen_coord1: fault trace not connected')
    end
end % end for loop over polygons


% Now check to ensure that polygons all strike in same direction, and
% rearrange if necessary

Vm = [ QQ(2,:)-QQ(1,:) ; QQ(4,:)-QQ(1,:) ]; % eqn (1) matrix of vectors 
% pointing from corner 1 to corner2 and from corner 1 to corner 4

Cm = NaN * ones(n,3); % matrix to hold normals to all polygons

% loop over all polygons to check strike order
for p = 1:n
    
    firstcol = 1 + (p-1)*3; % east coord of v vector in this column
    jj = firstcol: (firstcol+2); 
    v2t = Vm(1,jj); 
    v4t = Vm(2,jj); 
    
    Cm(p,:) = cross(v2t,v4t) ; % eqns (2,3).  matrix of normals to polygons
    
end % end loop over polygons p

% up unit vector 
ahat = [ 0 0 1]'; 

% A(j) holds the up component of the normal to polygon j
A = Cm * ahat ; % eqn (4)

% row 1 of QQ holds coords of fault trace points P1 through Pn, but not
% Pn+1.  QQa is row 1 of QQ augmented by Pn+1, so QQa holds coords of fault
% trace in a row vector
QQa = [ QQ(1,:) QQ(2, (3*n-2):(3*n)) ] ; % eqn 4.1, dim(QQa)= 1 x 3(n+1)

if max(A) <= 0
    % normals all point down, fault trace strike points from P1 to Pn+1
    Fm = reshape( QQa, 3, n+1) ; % Fm(j,p) is j=1(e), j=2(n), j=3(up) 
    % coord of point 1 of pth polygon
elseif min(A) > 0
    % normal all point upward.  fault trace strikes from Pn to P1. 
    % don't try to fix
    error('gen_coord1: fault trace points are in anti-strike order')
else
    % strikes of individual polygons conflict
    error('gen_coord1: fault trace strikes conflict')
end

% intialize storage
L = NaN * ones(1,n); % L(j) is length of segment j trace
Lm = NaN * ones(nstn,n); 
Um = Lm; 
Tm = Um; 
Zero = zeros(nstn,n); 

% loop over polygons determining u and t unit vectors.  Fill matrices Um
% and Tm with u and t coords of each station w.r.t. each polygon

for p = 1:n
    
    Pp   = Fm(:,p) ; % column vector, coords of Pp, eqn (9)
    Ppp1 = Fm(:,p+1) ; % column vector, coords of Pp+1
    
    Ep = Ppp1 - Pp ; % top edge of polygon, from Pp to Pp+1
    
    Lp = sqrt( sum (Ep.^2) ); % length of edge p
    uhatp = Ep/Lp ; % unit vector pointing from Pp to Pp+1, eqn (10)
    thatp = cross(-ahat,uhatp); % t unit vector, eqn (11)
    
    Ppm = Pp * ones(1,nstn) ; % eqn 12, dimension 3 x nstn
    
    % (Xm-Ppm)(i,j) is the i'th coord of the vector from Pp to stn j, 
    % i=1=east, i=2=north, i=3=up
    Upmt = uhatp' * (Xm - Ppm) ; % 1 x nstn vector containing the u coord of
    % each station in the system with origin at Pp (13)
    Tpmt = thatp' * (Xm - Ppm) ; % 1 x nstn containing t coord of stations (14)
    
    L(p) = Lp; 
    Um(:,p) = Upmt'; % (15) Um(i,p) = u coord of station i w.r.t. Pp
    Tm(:,p) = Tpmt'; % (16) t coord of station i w.r.t. Pp
end % end loop over polygons

Um;%
Tm;%

% finish calculation

Lm = ones(nstn,1) * L ;% % (17), dim(Lm) = nstn x n

% Mm(s,p) is the u coord of station s w.r.t. polygon p, truncated to be in
% the range [0 Lp]
Mm = min(  max(Um,Zero) , Lm);% % (18)

Dm = sqrt( (Um-Mm).^2 + Tm.^2  );% % (19) Dm(s,p) is the horizontal distance 
% from station s to the closest point on the top edge of polygon p

D = min(Dm, [], 2);% % (20) D(s) = Shortest horizontal distance from station
% s to anyplace on the fault trace.  Equals Delta in Spudich and Chiou

% filling dm, logical matrix of closest seg number. Ultimately, 
% dm(s,p) = 1 if segment p holds the closest point to station s, 
%         = 0 otherwise

dm = abs(Dm - (D*ones(1,n))) < distol;% 
% dm(s,p) = 1 if Dm(s,p) is within distol of D(s), i.e. 

% ensure that there is only 1 closest segment for each station
nclosest = sum(dm,2); 

if min(nclosest)<1 % changed from ==0 to <1 7/29/10
    % there is a station for which no closest segment was found. Error.
    error('gen_coord1: there is a station with no closest segment')
end

if max(nclosest) > 1
    % for some stations there are more than 1 closest segment.  Arbitrarily
    % retain the first
    badrows = find(nclosest>1); % vector of rows of dm with more than one 
    % closest segement
    % loop over bad rows expunging superfluous closest segments
    for badrow = badrows'
        badcols = find(dm(badrow,:)==1); % vector of columns of dm that are 1
        badcols(1) = []; % the first one ain't bad
        dm(badrow,badcols) = 0; % fill the bad ones with zeros
        % dm(s,p) = 1 if segment p holds the closest point to station s, 
        %         = 0 otherwise
    end 
    
    % check
    mclosest = sum(dm,2); 
    if max(mclosest) > 1 
        error('gen_coord1: still a problem after purging bad cols from dm')
    end
end % end if nclosest > 1

dm;%

% This works regardless of matrix size, but we must
% ensure that all nclosest == 1 (done above)
Tc = sum( Tm.*dm, 2);% % 7/29/10

% we must use a different formula for generalized T for stations inside and
% outside the extent of the fault trace.  N(s) = 1 if station s is withing
% the fault trace extent, = 0 otherwise
N = (Um(:,1)>=0) & (Um(:,n)<=Lm(:,n)) ;% % (24), dim(N) = nstn x 1

%generalized T
T = N .* D .* sign(Tc) + ~N .* Tc ;% % (23)

% generalized U
if n==1
    U = Um(:,1); %(26)
elseif n==2
    U = ( min(Um(:,1),Lm(:,1)) + max(Um(:,2),0)  )  ; % (27)
else
    U = ( min(Um(:,1),Lm(:,1)) + max(Um(:,n),0)  )  ...
        + sum(  Mm(:,2:(n-1)), 2 );
end

U;%

return


        




