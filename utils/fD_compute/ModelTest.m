function fD_new = ModelTest(M, dev, ihypo, Ttest, Tlist, filepath, ModelFlag) 

% SC08 model
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

% rupture related: 
name = ['hypo' num2str(ihypo) '_' dev '.txt'];
Ttest = 3.0; 

% default parameters:
vr_beta = 0.8; 
c_cut = 2.45;
s_cut = 75;   % km
r_cut = 0.2;  % radiation

% widowing
d_cut = [40,70];
m_cut = [5.6,6.0];

fileName = [filepath name];
A = dlmread( fileName, ',', 1, 0 );
Rrup = A(:,4);
Rfn = A(:,6);
Rfp = A(:,7);
s = A(:,8);
h = A(:,9);
ctildper = A(:,10);
IDP = A(:,11);
C = A(:,12);
S = A(:,13);

 % compute fD using the same parameter for files
cmd = ['Ttable = ' dev '(:,1);'];
eval(cmd);
jt = find(Ttable == Ttest);
if isempty(jt)
   error(['isodirect1: requested period ' num2str(Ttest) ...
          ' not in coef table'])
end
cmd = ['a = ' dev '(jt,4);'];
eval(cmd);
cmd = ['b = ' dev '(jt,5);'];
eval(cmd);

% compute windows:
fm = min(1,max(0,M-m_cut(1))/(m_cut(2)-m_cut(1)));
fr = max(0,1-max(0,Rrup-d_cut(1))/(d_cut(2)-d_cut(1)));

fD0= fm * fr .* (a + b*IDP) ; % dim(fD) = ns x length(Tlist)

% read from file fD
cmd = ['Ttable = ' dev '(:,1);'];
eval(cmd);
jt = find(Tlist == Ttest);
fD1 = A(:,end-jt);

% fD0 and fD1 should be the same (validation 1)

% ====================================================
% Original Parameters:
%vr_beta = 0.8; 
%c_cut = 2.45;
%s_cut = 75;   % km
%r_cut = 0.2;  % radiation water level
%d_cut = [40,70];
%m_cut = [5.6,6.0];

% use new set of parameters: (do the expriment to see the effect of
% vr/beta)
switch ModelFlag
    case 0
        vr_beta = 0.8;  % this should be exactly same as the above one
    case 1
        vr_beta = 0.7;  
    case 2
        vr_beta = 0.9;
    case 3
        vr_beta = 0.95;        
    case 4
        d_cut = [200,250];
    case 5
        c_cut = 4.0;
    case 6
        r_cut = 0.0;
    case 7
        vr_beta = 0.9; 
        d_cut = [200,250];    
    case 8
        c_cut = 4.0;
        d_cut = [200,250];
end

% Compute (Rrup-Rhyp)/D from original parameters (vr/beta = 0.8 in
% isodirect1.m)
geo = 1./0.8 - 1./ctildper;

% compute new c' by
ctildper_new = 1./(1./vr_beta - geo);
C_new = ( min(ctildper_new, c_cut) - vr_beta ) / (c_cut - vr_beta); 

% comptue S_new
S_new = log( min(s_cut,max(s,h)) );
Rri_new = max(sqrt(Rfn.^2+Rfp.^2), r_cut);

% compute IDP
IDP_new = C_new .* S_new .* Rri_new;


% compute new fD
cmd = ['Ttable = ' dev '(:,1);'];
eval(cmd);
jt = find(Ttable == Ttest);
if isempty(jt)
   error(['isodirect1: requested period ' num2str(Ttest) ...
          ' not in coef table'])
end

cmd = ['a = ' dev '(jt,4);'];
eval(cmd);
cmd = ['b = ' dev '(jt,5);'];
eval(cmd);

% compute windows:
fm = min(1,max(0,M-m_cut(1))/(m_cut(2)-m_cut(1)));
fr = max(0,1-max(0,Rrup-d_cut(1))/(d_cut(2)-d_cut(1)));

fD_new = fm * fr .* (a + b*IDP) ; % dim(fD) = ns x length(Tlist)


end
