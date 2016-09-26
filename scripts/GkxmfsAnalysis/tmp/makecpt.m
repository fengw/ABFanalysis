function makecpt( maxv, minv, Nclr, cptfile, cmapname )
% 
% Colormap Library
% input: colormap name (jet, or user extened names)
% output: cmap that can be use in changing colormap by: 
%         colormap( cmap ); colorbar()
% Feng Wang at USC
%

% Data Value
values = linspace(minv,maxv,Nclr);

defaults = {'jet','HSV','Hot','Cool','Spring','Summer','Autumn', ...
    'Winter','Gray','Bone','Copper','Pink','Lines'};
flag = 1;
for i = 1:length(defaults)
    if strcmp(cmapname,defaults{i})
        ['colormap( ' cmapname '(' num2str(Nclr) ') );']
        cmap0 = eval(['colormap( ' cmapname '(' num2str(Nclr) ') );'] );
        flag = 0;
        cmap = cmap0;    
        break;
    end
end

if flag
    switch cmapname
        % generate the colormap library with 64 as the default length
        % Note, You can use colormapeditor to generate whatever colormap
        % you want and save the cmap here as your colormap library
        % this also applies to Pylab and GMT (you need to convert to RGB
        % value first)
        
        case 'bwr'       
            cmap0 =[
                 0         0    0.5625
                 0         0    0.6250
                 0         0    0.6875
                 0         0    0.7500
                 0         0    0.8125
                 0         0    0.8750
                 0         0    0.9375
                 0         0    1.0000
            0.0235    0.0471    1.0000
            0.0471    0.0941    1.0000
            0.0706    0.1412    1.0000
            0.0941    0.1882    1.0000
            0.1176    0.2353    1.0000
            0.1412    0.2824    1.0000
            0.1647    0.3294    1.0000
            0.1882    0.3765    1.0000
            0.2118    0.4235    1.0000
            0.2353    0.4706    1.0000
            0.2588    0.5176    1.0000
            0.2824    0.5647    1.0000
            0.3059    0.6118    1.0000
            0.3294    0.6588    1.0000
            0.3529    0.7059    1.0000
            0.3765    0.7529    1.0000
            0.4000    0.8000    1.0000
            0.5000    0.8333    1.0000
            0.6000    0.8667    1.0000
            0.7000    0.9000    1.0000
            0.8000    0.9333    1.0000
            0.9000    0.9667    1.0000
            1.0000    1.0000    1.0000
            1.0000    0.9714    0.9143
            1.0000    0.9429    0.8286
            1.0000    0.9143    0.7429
            1.0000    0.8857    0.6571
            1.0000    0.8571    0.5714
            1.0000    0.8286    0.4857
            1.0000    0.8000    0.4000
            1.0000    0.7500    0.3750
            1.0000    0.7000    0.3500
            1.0000    0.6500    0.3250
            1.0000    0.6000    0.3000
            1.0000    0.5500    0.2750
            1.0000    0.5000    0.2500
            1.0000    0.4500    0.2250
            1.0000    0.4000    0.2000
            1.0000    0.3500    0.1750
            1.0000    0.3000    0.1500
            1.0000    0.2500    0.1250
            1.0000    0.2000    0.1000
            1.0000    0.1500    0.0750
            1.0000    0.1000    0.0500
            1.0000    0.0500    0.0250
            1.0000         0         0
            0.9375         0         0
            0.8750         0         0
            0.8125         0         0
            0.7500         0         0
            0.6875         0         0
            0.6250         0         0
            0.5625         0         0
            0.5000         0         0];

        %case 'white2red'
                    
    end
    value0 = linspace( minv, maxv, size(cmap0,1) );
    cmap = zeros(Nclr,3);
    
    % value (actual sampling) and interpolation
    Nc = size(cmap0,1);
    if Nc >= Nclr
        % downsampling  
        for i = 1:3
            cmap(:,i) = interp1( value0, cmap0(:,i), values );
        end
    else
        % up sampling doesn't work
        error('you should generate a colormap that has larger numer of points that your actural data points')
    end
end

% write in to file using fprintf
fid = fopen( cptfile, 'w' );
fprintf( fid, '%s\n', '#COLOR_MODEL = RGB' );
cmap = 255*cmap;
for ipoint = 1:length(values)-1
    fprintf( fid, '%-5.3f%6.0f%6.0f%6.0f   %-5.3f%6.0f%6.0f%6.0f\n', ...
    values(ipoint), cmap(ipoint,1), cmap(ipoint,2), cmap(ipoint,3), ...
    values(ipoint+1), cmap(ipoint+1,1), cmap(ipoint+1,2), cmap(ipoint+1,3) );    
end
fclose(fid);

end
