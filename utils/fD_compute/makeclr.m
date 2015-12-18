function map = makeclr( cmap )
% make colormap matrix given color option
% input: cmap = 'colormap name', e.g.: jet, cool, etc.
% output: rgb matrix

% first, set up the color library
map0.name = {'bwr',};
          %   r g b
map0.rgb = { [0,0,1;    % start color
              1,1,1;    % center color
              1,0,0];   % end color           
           };

% matlab supported colormap name:
default={'jet','hsv','hot','cool',...
    'spring','summer','autumn','winter',...
    'gray','bone','copper','pink','lines'};
N = length( default );
for il =1:N
    if strcmp( cmap, default{il} )
        flag = 1;
        break;
    else
        flag = 0;
    end
end

if flag == 1
    map = colormap( cmap );   % matlab supported colormaps
else    
    N = length( map0.name );
    for il = 1:N
        if strcmp( cmap, map0.name{il} )    
            rgbv = map0.rgb{il};
            Nc = 30;   % discretization
            
            red = [linspace(rgbv(1,1),rgbv(2,1),Nc),linspace(rgbv(2,1),rgbv(3,1),Nc)];
            green = [linspace(rgbv(1,2),rgbv(2,2),Nc),linspace(rgbv(2,2),rgbv(3,2),Nc)];
            blue = [linspace(rgbv(1,3),rgbv(2,3),Nc),linspace(rgbv(2,3),rgbv(3,3),Nc)];
            
            map = [red',green',blue']; 
            break;
        else
            err = 'pleace use valid cmap name';
            map = 0;
            break
        end
    end
end
