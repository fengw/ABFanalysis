inpth = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups3/ERF35_SGT5_RupVar3_Vel1/Gkxmfs/RMS';
filename = 'CyberShake.NGAs.periods.Sigma0.60.RefBA.RMS.unweighted.txt';
filename = 'CyberShake.NGAs.periods.Sigma0.60.RefAS.RMS.unweighted.txt';
file0 = [inpth '/' filename];

inp = dlmread( file0, ' ' );

xt = 0:5;
yt = 0:4;

bar3(inp,0.4);