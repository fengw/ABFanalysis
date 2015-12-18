function [outs, in, Pth] = fDcompute(wrkpth,id0,id1,source_info,segment_info,sites_info,NGA_info,outpth,Visual)
% Function to compute fD for given rupture,hypocenter,sites, and NGA info
% Coordinate: East-North-Zup
% 

% ======================
% Add Utility Path
% ======================
addpath(wrkpth);

% =================
% ID infos
% =================
erf_id = id0{1};
rup_var_id = id0{2};

sid = id1{1};
rid = id1{2};
ihypo = id1{3};

% ===============================
% Hypocenter and Earthquake info 
% ===============================
fid = fopen(source_info,'r'); %# open csv file for reading
il = 1; A = 'tmp';
while ~feof(fid)
    line = fgets(fid); %# read line by line
    A = [A ',' sscanf(line,'%s')]; %# sscanf can read only numeric data :(
end
fclose(fid);
B = regexp(A,',','split');
hypo = double(str2num(B{2}));
for i = 3:4
    hypo = [hypo double(str2num(B{i}))];
end

% hypocenter (east, north, up) coords
hypoi_str = ['hypo' ihypo];

in.H = hypo;
in.rake = double(str2num(B{5}));
% magnitude (assumed same for all hypocenters)
in.M = double(str2num(B{6})); 

% ===============================
% Fault segment info
% ===============================
QQ = dlmread( segment_info,' ');
in.QQ = QQ;
in.vr_beta = 0;    % if not 0, used for vr/beta effects on directivity

% check whether hypo is on the plane or not:
% iseg = 10;
% points = in.qq(:,3*iseg-2:3*iseg);
% [dist, point1] = pointplane( in.h, points);
% dist

% test fault segments and hypocenter locations
if Visual.fault == 1
    VisualFaults(in,Visual.faultFigPath,Visual.faultFigName)
end

% write the output to a formatted text file  
if ~isdir(outpth)
    mkdir(outpth);
end    
pth1 = [outpth '/ERF' erf_id '_RupVar' rup_var_id];
if ~isdir(pth1)
    mkdir(pth1);
end   
pth2 = [pth1 '/SourceID' sid '_RuptureID' rid];
if ~isdir(pth2)
    mkdir(pth2);
end   
prefix = [pth2 '/hypo' num2str(ihypo)];


% =============
% Site info
% =============
sites = dlmread( sites_info,' ', 0, 0);
Ns = sites(1,1);   % scatter points 
NiX = sites(1,2);  % regular mesh X-direction points
NiY = sites(1,3);  % regular mesh X-direction points
Nid = sites(1,4);  % regular mesh grid spacing in both direction
in.Ndims = [Ns,NiX,NiY,Nid];

%sites = dlmread( sites_info,' ',1,0 );  % header is the total # of sites
in.Xs = [sites(2:end,2) sites(2:end,3) sites(2:end,4)];
siteId = sites(2:end,1);   % site ID

% ====================
% NGA model info
% CB6, BA6, CY6, AS6
% ====================
fid = fopen(NGA_info,'r'); %# open csv file for reading
il = 1; A = 'header';
while ~feof(fid)
    line = fgets(fid); %# read line by line
    A = [A ';' sscanf(line,'%s')]; %# sscanf can read only numeric data :(
end
fclose(fid);
B = regexp(A,';','split');
B1 = regexp(B{2},',','split');
B2 = regexp(B{3},',','split');
in.Tlist = [];
for it = 1:length(B2)
    in.Tlist = [in.Tlist str2num(B2{it})];
end

for i = 1:length(B1)
    dev = B1{i};    % different NGA models
    in.dev = dev; 

    % compute (key)
    out = isodirect1(in); 

    if strcmp(in.dev, Visual.fDngaM)==1
        outs = out;
    end

    stem = [ prefix '_' in.dev];
    filename = [ stem '.txt'];
    fid = fopen(filename, 'w'); 
    outh = repmat( out.h, size(out.s) );
    x = [siteId in.Xs(:,1:2) out.Rrup out.D out.RuRt(:,2) out.RuRt(:,1) ...
        out.s outh out.ctildepr out.IDP out.C out.S out.gencoords out.fD];

    colhed = [' ID     X      Y      Rrup         D         ' ...
        '  Rfn         Rfp          s            h          ctildepr    IDP      ' ...
          '    C          S            U            T     '];
    pds = sprintf('    fD %4.2fs    ', in.Tlist); 
    fprintf( fid, [ colhed   pds ] );
    fmt = ['\n %d, %5.0f, %5.0f, %11.7f, %11.7f, %10.7f, ' ...
        '%10.7f, %11.7f, %11.7f, %8.6f, %8.6f, %10.7f, %11.7f, %11.7f, %10.7f'];
    pdfmt = '';
    for j = 1:length(in.Tlist)
        pdfmt = [pdfmt ', %12.7e'];
    end

    fprintf(fid, [ fmt pdfmt], x');
    fclose(fid);

    % visualize fD 
    Pth.fDFilePath = [pth2 '/'];
    Pth.fDFigPath = [pth2 '/fDplots/'];
    if Visual.fD == 1 && strcmp(in.dev, Visual.fDngaM)==1 && mod(ihypo,1)==0
        for ModelFlag = 8:8
            VisualfD(in, out, @ModelTest, ihypo, Visual, Pth, ModelFlag);
        end
    end   
    
end

end  % end of function
