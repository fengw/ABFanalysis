% run script
% in command line: matlab -nojvm < run.m
% 
clc;
wkd = '/Users/fengw/work/Project/CyberShake_analysis/utils/fD_compute/';
input_file = [wkd 'input_file'];

% INPUT info
fid = fopen(input_file,'r');
A = [];
while ~feof(fid)
    line = fgets(fid); %# read line by line
    A = [A sscanf(line,'%s') ';'];  % separate lines by ';'
end
fclose(fid);
B = regexp(A,';','split');
wrk = B{1};   % directory where the isodirect_util and fD_compute are located

% main utilities
wrkpth = [wrk 'isodirect_util/'];

% Rupture dependent files
erf_id = B{2};
rup_scenario_id = B{3};
id0 = {erf_id rup_scenario_id};

sids = regexp(B{4},',','split');
rids = regexp(B{5},',','split');
hypos = regexp(B{6},',','split');

% visualization (can read from file to turn on or turn off visualization)
VisualP = regexp(B{7}, ',', 'split');
Visual.fault = str2num(VisualP{1});        % Test fault geometry
Visual.fD = str2num(VisualP{2});           % Test fD calculation (shape of directivity)

Visual.fDABF = str2num(VisualP{3});        % averaging-based fD decomposition

Visual.fDPeriod = str2num(VisualP{4});
Visual.fDngaM =VisualP{5};

if Visual.fD == 1 && Visual.fault == 1
    Visual.fault = 0;
end

% PATH info
cwd = [wrk 'fD_compute/'];

% input and output directory
inputs = [cwd 'inputs/' ];
outputs = [cwd 'outputs'];

% common files
NGA_info = [inputs 'NGA_info'];

pth0 = [inputs 'ERF' erf_id '_RupVar' rup_scenario_id '/'];
for irup = 1:length(sids)
    sid=sids{irup};
    rid=rids{irup};
    Nh = str2num(hypos{irup});    
    pth1 = [pth0 'SourceID' sid '_RuptureID' rid '/'];    
    Visual.faultFigPath = [pth1 'FaultPlots/'];
    
    outs_h = cell(Nh,1); ins_h = cell(Nh,1);
    
    for ih =1:Nh    
        id1 = {sid rid ih};
        sites_info = [pth1 'sites_info'];
        source_info = [pth1 'hypo_' num2str(ih)];
        segment_info = [pth1 'rupture_info'];
        
        Visual.faultFigName = ['hypo_' num2str(ih)];
        
        [ outs,in,OutPth ] = fDcompute(wrkpth,id0,id1,source_info,segment_info,sites_info,NGA_info,outputs,Visual);
        outs_h{ih} = outs;
        ins_h{ih} = in;
    end
    
    if Visual.fDABF == 1
        % compute fD-<fD>_ih and visual
        ModelFlag = 4;  % extended distance window
        ABF_fD(outs_h, ins_h, Nh, Visual, OutPth,ModelFlag)
    end    
    
end
fprintf('Done! \n')
