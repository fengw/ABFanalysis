function axs = subaxes(figH, subs, basic,FigSize)

% subs = (1,1) or (2,2) or others
% basic = (x0, y0, dx,dy) in pixel

% get screen size (un-normalized)
scrsz = get(0,'ScreenSize');

set(figH,'Units','normalized')
set(figH,'Position',FigSize);

W0 = scrsz(3); H0 = scrsz(4);

Nrow = subs(1); Ncol = subs(2);

x0 = basic(1)+100; y0 = basic(2);
dx = basic(3)+150; dy = basic(4);

width = (W0-2*x0-(Ncol-1)*dx) / Ncol;  % in pixel
height = (H0-2*y0-(Nrow-1)*dy) / Nrow; % in pixel

axs = zeros(Nrow*Ncol,4);   % should be in ratio
yi = y0 + (Nrow-1)*(height+dy);
i = 1;
for irow= 1:Nrow
    xi = x0;
    for icol = 1:Ncol
        
        axs(i,:) = [xi/W0,yi/H0,width/W0,height/H0];
        i = i + 1; 
        xi = xi + (width+dx);
    end
    yi = yi - (height + dy);
end

end


