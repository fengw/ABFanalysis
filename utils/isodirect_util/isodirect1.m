function [ out ] = isodirect1( in )
%
% usage:  out = isodirect1(in); 
%
% latest modif: 5 Oct 2010 v1.0
% isodirect1 based on spudchiou2 created 13 sept 2010 - P. Spudich
%
% isodirect1 evaluates the Spudich and Chiou (Earthquake Spectra, 2008) 
% directivity model for a vector of site locations, given a fault geometry,
% earthquake magnitude, set of desired periods, and desired developer
% model. 
%
% Based on notes in 7.27.10.I-ImplementingSpud&Chiou.ai
%
% inputs: 
%
%   in - (structure) with fields     
%     .QQ - matrix of dimension 4 x 3n, where n is the number of polygons
%       comprising the fault surface.  Each polygon is a quadrilateral.
%       QQ is the field pol.QQ in polygons?.mat.  
%       QQ consists of submatrices QQ1, QQ2, ... QQn, i.e. 
%       QQ = [QQ1 QQ2 ... QQn], where size(QQi) = [4,3] and QQi(k,:) is the
%       [east north up] coord (km or m, not dg) of the k'th point of
%       polygon
%       i.  For each quadrilateral points QQi(1,:), QQi(2,:), ... QQi(4,:)
%       are ordered sequentially around the quadrilateral, and QQi(1,:) and
%       QQi(2,:) for all i must be at the same depth, at the top of the
%       polygon.  QQi(1,:) corresponds to point Pi in Spudich and Chiou, 
%       and QQn(2,:) correponds to point Pn+1.  It is necessary that 
%       the fault strike points in the direction from QQi(1,:) to QQi(2,:),
%       and the polygons must be ordered such that each strikes in the same
%       general direction. 
%
%     .Xs - matrix of site locations, dimension Xs = ns x 3, where ns
%       is the number of sites.  Xs(j,:) is the [east ; north ; up ]
%       coord (in km or m, not degrees) of the jth site. 
%
%     .dev - (string), either 'AS6', 'BA6', 'CB6', or 'CY6' - developer
%       model for which directivity is calculated
%
%     .Tlist - vector of desired periods (column or row okay)
%
%     .M - moment magnitude of event
%
%     .rake - fault rake
%
%      Hypocenter coords may be specified two ways:  
%      Method 1:  polygon-centered coord system: 
%
%     .ih - scalar, segment number hypocenter is located on 
%
%     .h - scalar, distance downdip from top of segment ih to hypocenter
%
%     .uh - scalar >=0, along-strike distance from QQ(1,ih) to hypocenter
%
%      Method 2: absolute space coords:
%
%     .H - row or column vector of hypocenter coords, length(H)=3, where 
%       H(:) is the [east ; north ; up ] coord (in km or m, not degrees)
%
%      NOTE - if you use method 1 to specify hypocenter, fields ih, h, and
%      uh should not exist in structure in, i.e. isfield(in,'h') should
%      return 0. If you use method 2 to specify hypocenter, field H
%      should not exist in structure in, i.e. isfield(in,'H') should
%      return 0.  
%
% outputs:
%
%   out - (structure) with fields
%
%     .fD - (array, dimension ns x length(Tlist) ).  fD(j,k) is the
%       predicted directivity from Spudich and Chiou equation 5 for the 
%       site at Xs(j,:) for period Tlist(k).
%
%     .Xc - matrix of closest points, dimension Xc = ns x 3, where ns
%       is the number of sites.  Xc(j,:) is the [east ; north ; up ]
%       coord of the closest point on the fault to site in.Xs(j,:)  
%
%     .Rrup - (vector, dimension ns x 1).  Rupture distance (Euclidean 
%       distance from out.Xs(j,:) to out.Xc(j,:) ).
%
%     .C - (vector, dimension ns x 1).  C from Spudich and Chiou equation
%       2.  C(j) is C at Xs(j,:).
%
%     .gencoords - (matrix, dimension ns x 2).  out.gencoords(j,1) is the U
%       generalized coordinate of site in.Xs(j,:).  out.gencoords(j,2) is 
%       the T generalized coordinate of site in.Xs(j,:).
%
%     .gencoordsh - (vector, dimension 1 x 2).  out.gencoordsh(1) is the U
%       generalized coordinate of the hypocenter.  out.gencoordsh(2) is 
%       the T generalized coordinate of the hypocenter.
%
%     .S - (vector, dimension ns x 1).  S from Spudich and Chiou equations
%       1 and 3.  out.S(j) is S corresponding to Xs(j,:).
%
%     .s - (vector, dimension ns x 1).  Lower-case s from Spudich and Chiou
%       equation 3 and Figure 1.  out.s(j) is s corresponding to Xs(j,:).
%
%     .D - (vector, dimension ns x 1).  D from Spudich and Chiou
%       equation 4 and Figure 1.  out.D(j) is s corresponding to Xs(j,:).
%
%     .ctildepr - (vector, dimension ns x 1).  c-tilde-prime from Spudich 
%       and Chiou equation 4.  out.ctildepr(j) corresponds to Xs(j,:).
%
%     .RuRt - (matrix, dimension ns x 2).  out.RuRt(j,1) is the Ru (strike-
%       parallel) radiation pattern term from Appendix A for the site 
%       in.Xs(j,:).  out.RuRt(j,2) is the Rt (strike-perpendicular) 
%       radiation pattern term from Appendix A for the site in.Xs(j,:)
%
%     .Rri - (vector, dimension ns x 1). Rri from Spudich and Chiou
%       appendix A.  Rri(j) is Rri at Xs(j,:).
%
%     .IDP - (vector, dimension ns x 1). IDP from Spudich and Chiou equation
%       1.  IDP(j) is IDP at Xs(j,:).
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
%

AS6 = [
    0.5 573 40 0.0000 0.0000 0.5414 0.2247 4.210E-02 0.5414 0
0.75 572 40 -0.0447 0.0298 0.5586 0.2460 -2.137E-02 0.5596 0.002
1 572 40 -0.0765 0.0510 0.5553 0.3138 2.305E-02 0.5598 0.008
1.5 562 38 -0.1213 0.0809 0.5174 0.3260 3.380E-02 0.5225 0.01
2 538 38 -0.1531 0.1020 0.5240 0.3711 4.692E-02 0.5341 0.019
3 465 35 -0.1979 0.1319 0.5178 0.3595 3.191E-02 0.5393 0.04
4 436 34 -0.2296 0.1530 0.5294 0.3922 6.307E-02 0.5506 0.039
5 328 30 -0.2542 0.1695 0.5285 0.4192 8.158E-02 0.5510 0.041
7.5 276 28 -0.3636 0.2411 0.5147 0.4332 -5.407E-02 0.5560 0.074
10 158 17 -0.5755 0.3489 0.5566 0.3656 -1.517E-01 0.6626 0.16];

BA6 = [
    0.5 419 27 0.0000 0.0000 0.5212 0.1945 -8.870E-03 0.5212 0
0.75 418 27 -0.0532 0.0355 0.5387 0.2618 -3.463E-03 0.5394 0.001
1 418 27 -0.0910 0.0607 0.5278 0.3068 5.466E-02 0.5327 0.009
1.5 412 26 -0.1443 0.0962 0.5052 0.3214 8.360E-02 0.5091 0.008
2 390 25 -0.1821 0.1214 0.5191 0.3815 9.181E-02 0.5301 0.021
3 371 25 -0.2353 0.1569 0.5197 0.3985 3.557E-02 0.5484 0.052
4 363 25 -0.2731 0.1821 0.5247 0.3637 3.217E-02 0.5559 0.056
5 263 20 -0.3021 0.2015 0.5513 0.3801 2.285E-02 0.5973 0.077
7.5 234 20 -0.4627 0.2727 0.5340 0.4514 4.121E-03 0.6005 0.111
10 129 12 -0.8285 0.4141 0.5171 0.3387 -1.210E-01 0.6503 0.205];

CB6 = [
  0.75 438 36 0.0000 0.0000 0.5298 0.2247 2.525E-02 0.5298 0
1 438 36 -0.0329 0.0220 0.5234 0.2666 5.890E-02 0.5243 0.002
1.5 431 34 -0.0795 0.0530 0.4889 0.2699 8.493E-02 0.4899 0.002
2 409 34 -0.1125 0.0750 0.4921 0.2507 7.915E-02 0.4972 0.01
3 387 31 -0.1590 0.1060 0.4964 0.2556 5.891E-02 0.5129 0.032
4 379 31 -0.1921 0.1280 0.5075 0.2323 8.163E-03 0.5250 0.033
5 276 27 -0.2172 0.1450 0.5206 0.2677 -4.461E-02 0.5481 0.05
7.5 248 27 -0.3227 0.2147 0.5151 0.3580 -8.447E-02 0.5613 0.082
10 129 16 -0.6419 0.3522 0.5365 0.4071 -2.105E-01   0.6497 0.174];

CY6 = [
    0.75 570 40 0.0000 0.0000 0.5428 0.3615 8.438E-02 0.5428 0
1 570 40 -0.0260 0.0200 0.5393 0.4042 8.660E-02 0.5404 0.002
1.5 560 38 -0.0627 0.0482 0.5097 0.3998 9.405E-02 0.5113 0.003
2 536 38 -0.0887 0.0682 0.5307 0.4044 1.002E-01 0.5349 0.008
3 462 35 -0.1254 0.0965 0.5311 0.4335 1.064E-01 0.5431 0.022
4 432 34 -0.1514 0.1165 0.5503 0.4274 1.389E-01 0.5626 0.022
5 324 30 -0.1715 0.1320 0.5527 0.4895 5.773E-02 0.5655 0.023
7.5 272 28 -0.2797 0.1865 0.5476 0.4693 5.731E-02 0.5713 0.041
10 154 17 -0.4847 0.2933 0.5819 0.3077 1.719E-02 0.6454 0.098];

QQ = in.QQ;
Xs = in.Xs;
Tlist = in.Tlist; 

[nrqq,ncqq]= size(QQ);
n = ncqq/3;   % number of polygons

% sanity checks of input geometry variables -----------------------
% equations in this section from notes 9.23.07.I.

if nrqq~=4 ; error('isodirect1: QQ does not have 4 rows') ; end
if rem(n,1) ~= 0 ; error('isodirect1: QQ has wrong number of columns') ; end

[ns, nrxm] = size(Xs); 
if nrxm ~= 3 ; error('isodirect1: Xs does not have 3 columns') ; end

% check tops of all polygons at same elevation

acol = 3*(1:n); % column containing elevation
%if std([QQ(1,acol)' ; QQ(2,acol)']) ~= 0
if std([QQ(1,acol)' ; QQ(2,acol)']) >= 1.e-10     % tol of zeros (float)
    error('isodirect1: tops of all quadrilaterals at different depths')
end

% check that polygons are connected.  Point 2 of polygon i should equal
% point 1 of polygon i+1
for p = 1:(n-1)
    firstcol = 1 + (p-1)*3; % east coord of P vectors
    jj = firstcol: (firstcol+2); 
    if max(abs( QQ(2,jj)-QQ(1,jj+3) )) > 0
        error('isodirect1: fault trace not connected')
    end
end % end for loop over polygons

% Now check to ensure that polygons all strike in same direction, and
% rearrange if necessary

Vm = [ QQ(2,:)-QQ(1,:) ; QQ(4,:)-QQ(1,:) ];  % eqn (1) matrix of vectors 
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
A = Cm * ahat; % eqn (4)

% row 1 of QQ holds coords of fault trace points P1 through Pn, but not
% Pn+1.  QQa is row 1 of QQ augmented by Pn+1, so QQa holds coords of fault
% trace in a row vector
QQa = [ QQ(1,:) QQ(2, (3*n-2):(3*n)) ]; % eqn 4.1, dim(QQa)= 1 x 3(n+1)

if max(A) <= 0
    % normals all point down, fault trace strike points from P1 to Pn+1
    Fm = reshape( QQa, 3, n+1) ; % Fm(j,p) is j=1(e), j=2(n), j=3(up) 
    % coord of point 1 of pth polygon
elseif min(A) > 0
    % normal all point upward.  fault trace strikes from Pn to P1. 
    % don't try to fix
    error('isodirect1: fault trace points are in anti-strike order')
else
    % strikes of individual polygons conflict
    error('isodirect1: fault trace strikes conflict')    
end

if isfield(in,'H')
    
    if isfield(in,'ih') || isfield(in,'h') || isfield(in,'uh') 
        error(['isodirect1: fields in.H and one or more of .ih, .h, or .uh exist'])
    end
    
    xhpr = in.H(:)';  % make it a row vector
    
    % check that determined hypocenter is actually on the input fault segment
    test = msegcl5(xhpr, QQ); 
    
    if test.distance > 0.001
        % original (just give the error and stop)
        %error( ['isodirect1: in.H is ' num2str(test.distance) ...
        %    ' from closest point on fault'])
        
        % discussion with Spudich 05/01/2012, you can use the closest point
        % on the fault segments as the hypocenter
        xhpr = test.xc;
       
    end
    
    % determine polygon number ih on which hypocenter lies
    ih = test.iclosest;
    
    % uhatpr is unit vector along strike of segment ih
    % jcols contains indices of columns of QQ for segment ih
    jcols = (3*ih-2):(3*ih); 
    uhatpr = (QQ(2,jcols)-QQ(1,jcols)) / norm( (QQ(2,jcols)-QQ(1,jcols)) ); 

    vpr = QQ(4,jcols) - QQ(1,jcols); % vpr points along edge of segment (not 
    %                             necessarily directly downdip)

    %dhatpr is a unit vector pointing directly downdip
    dpr = vpr - (uhatpr*vpr') * uhatpr; 
    dhatpr = dpr / norm(dpr); 
    diprad = pi/2 - acos( dhatpr * [0 0 -1]' ); 
    dipdg = diprad * 180/pi;
    
    % h is distance downdip from top of fault to hypocenter
    h = (xhpr - QQ(1,jcols)) * dhatpr'; 
    % uh is the along-strike dist of the hypocenter from the first point of
    % the fault trace of the ih'th polygon
    uh = (xhpr - QQ(1,jcols)) * uhatpr';

    
elseif isfield(in, 'ih')

    % Determining [e n u] coordinates of hypocenter
    ih = in.ih;
    h = in.h; 
    uh = in.uh; 

    % uhatpr is unit vector along strike of segment ih
    % jcols contains indices of columns of QQ for segment ih
    jcols = (3*ih-2):(3*ih); 
    uhatpr = (QQ(2,jcols)-QQ(1,jcols)) / norm( (QQ(2,jcols)-QQ(1,jcols)) ); 

    vpr = QQ(4,jcols) - QQ(1,jcols); % vpr points along edge of segment (not 
    %                             necessarily directly downdip)

    %dhatpr is a unit vector pointing directly downdip
    dpr = vpr - (uhatpr*vpr') * uhatpr; 
    dhatpr = dpr / norm(dpr); 

    diprad = pi/2 - acos( dhatpr * [0 0 -1]' ); 
    dipdg = diprad * 180/pi;

    % [e n u] coords of hypocenter
    xhpr = QQ(1,jcols) + uh*uhatpr + h*dhatpr; 

    % check that determined hypocenter is actually on the input fault segment
    test = msegcl5(xhpr, QQ); 

    if test.distance > 0.001
        error( ['isodirect1: xhpr is ' num2str(test.distance) ...
            ' from closest point on fault'])
    end

    if test.iclosest ~= ih
        error( ['isodirect1: iclosest = ' int2str(test.iclosest) ...
            ' ~= ih = ' int2str(ih) ])
    end
    
else
    % error inputting hypocenter coords
    error(['isodirect1: neither in.H or in.ih exist'])

end

% bogus exit point for testing with test_isodirect1_hypo_calc.m
%if ~isequal(ih,pi);
%    disp(['BOGUS EXIT FROM isodirect1.M'])
%    out.ih = ih;
%    out.uh = uh;
%    out.h = h;
%    out.H = xhpr; 
%    return
%end

% end sanity checks of input geometry variables ------------------


% determine quantities independent of period =====================

% determine Xc, matrix of closest points on fault to sites Xs
out1 = msegcl5(Xs, QQ, ih);
Xc = out1.xc; % dim(Xc) = ns x 3

% rupture distances from sites Xs to closest points Xc
Rrup = out1.distance ; % dim(Rrup) = ns x 1

Xh = ones(ns,1) * xhpr; % dim(Xh) = ns x 3

% distance D on fault from Xh to Xc
if n == 1
    % single segment D
    D = ( sqrt( sum( (Xc'-Xh').^2 )  ) )' ; % dim(D) = ns x 1
    % along strike distance s on fault
    s = abs( (Xc-Xh)*uhatpr' ); % dim(s) = ns x 1
else
    
    in2.QQ = in.QQ;
    in2.C = Xc;
    in2.c = out1.iclosest; 
    in2.H = xhpr; 
    in2.h = ih; 
    
    out2 = p4cfv1p1( in2 ); 
    D = out2.D;
    s = out2.s; 
end


% distances from sites Xs to hypocenter
Rh = ( sqrt( sum ( (Xs'-Xh').^2 )  )  )' ; % dim(Rh) = ns x 1

if in.vr_beta ~= 0
    vr_beta = in.vr_beta;
else
    vr_beta = 0.8;
end
ctilde = 1 ./ ( 1/vr_beta - (Rh-Rrup)./D ) ; % dim(ctilde) = ns x 1


% The way that Spudich and Chiou handles the case where D = 0
% you might want to fix this problem
ismallD = find(D<0.01);
if ~isempty(ismallD)
    ctilde(ismallD) = vr_beta;    % 0.8 = rupture velocity / shear velocity
end

% check for D = 0
%if ~isempty(find(D<0.01, 1))
%    error( ['isodirect1: some D < 0.01; code cannot yet handle'])
%end

one = ones(ns,1);

% ctilde is the only parameter that depends on vr/beta


% The following can be replaced by using different form of predictor
Ccut = 2.45;  % this is determined empirically 
C = ( min( ctilde, Ccut*one) - vr_beta ) / (Ccut-vr_beta) ; % dim(C) = ns x 1
S = log( min( 75*one, max( s, h*one )  )  ) ; % dim(S) = ns x 1

% Us and Ts are the u and t (generalized) coords of the sites
[Us, Ts] = gen_coord1( QQ, Xs');
[Uh, Th] = gen_coord1( QQ, xhpr'); 

radin.utds = [Us Ts zeros(ns,1)];
radin.utdh = [Uh Th -xhpr(3)]; 
radin.dip = dipdg;
radin.rake = in.rake; 
radin.watlev = 0.2; 

% calculate radiation pattern
radout = radpatv7(radin); 

Rri = radout.radpatmag; % magnitude of radiation pattern, dim = ns x 1
% note - radout.radpatmag has already been floored at radin.watlev

IDP = C .* S .* Rri; % dim(IDP) = ns x 1

out.C = C;
out.S = S;
out.Rri = Rri; 
out.IDP = IDP;
out.RuRt = radout.radpat; 
out.Rrup = Rrup; 
out.ctildepr = ctilde; 
out.s = s;
out.h = h;
out.D = D; 
out.Xc = Xc; 
out.gencoords = [Us Ts];
out.gencoordsh = [Uh Th];
zero = zeros(ns,1); 

% modified by Feng Wang (cutoff distance larger for CyberShake application)
% unity for 0<Rrup<cutdist; and tapers linearly to its value of zeros at Rrup >= cutdist1
% since distance doesn't affect IDP in which Rrup could be very large.
% what value you want to use depends on your problem size.
if 1
    % Use default original
    cutdist = 40;
    cutdist1 = 70;
else
    % CyberShake (use this) corrected in pynga/SC08.py
    cutdist = 200;
    cutdist1 = 250;
end
fr = max( zero,  (one - max(zero,Rrup-cutdist*one)/(cutdist1-cutdist) )); % dim(fr) = ns x 1
fm = min( 1, max(0,in.M-5.6)/0.4   ) ; 

% end determine quantities independent of period =====================

cmd = ['Ttable = ' in.dev '(:,1);'];
eval(cmd)

fD = NaN * ones(ns, length(Tlist) ); % initialize output storage

for iT = 1:length(Tlist)
    
    if Tlist(iT) < min(Ttable)
       disp(['isodirect1: requested period ' num2str(Tlist(iT)) ...
              ' less than minimum in coef table.  fD set to zero.'])
       a = 0; 
       b = 0;
       fD(:,iT) = fm * fr .* (a + b*IDP) ; % dim(fD) = ns x length(Tlist)
        
    else
     
        jt = find(Ttable == Tlist(iT));
        if isempty(jt)
           error(['isodirect1: requested period ' num2str(Tlist(iT)) ...
                  ' not in coef table'])
        end

        cmd = ['a = ' in.dev '(jt,4);'];
        eval(cmd);
        cmd = ['b = ' in.dev '(jt,5);'];
        eval(cmd);
        cmd = ['a = ' in.dev '(jt,4);'];
        eval(cmd);

        fD(:,iT) = fm * fr .* (a + b*IDP) ; % dim(fD) = ns x length(Tlist)
    
    end
    
end % end for loop over periods in Tlist

out.fD = fD; 

return


end

