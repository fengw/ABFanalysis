clear all;


%%%%%The points in the polygon will be displayed
% CyberShake Region 
poly_lon = [-116.843,-118.732,-119.336,-117.463,-116.843];
poly_lat = [34.172,35.058,34.150,33.273,34.172]; 

%%%%%The total boundary
minlat = 33.29
maxlat = 35.02
minlong = -119.34
maxlong = -116.9

%%%%%%spaceing, which should agree with the psmask
spacing = 4/100.

latr =  minlat:spacing:maxlat; ;
longr = minlong:spacing:maxlong;
[llong,llat] = meshgrid(longr,latr);

size = length(latr)*length(longr);
plong = reshape(llong,size,1);
plat = reshape(llat,size,1);

%z = ones(size,1);
in = inpolygon(plong,plat,poly_lon,poly_lat);

z = [plong(in) plat(in)];

save more_grid.station z -ASCII
