
flag = input( 'Choose the operation: ');

switch flag
    case 1
        % test subplot:
        fig = figure(1); clf;
        subs = [2,2];
        basic = [100,100,100,100];   % in pixel;
        FigSize = [0.1,0.1,0.3,0.7]; % this has to be in ratio relative to screen size
        axs = subaxes(fig, subs, basic,FigSize);
        for isub= 1:size(axs,1)
            ax = axes('position',axs(isub,:));
            plot([1,2,3,isub],'ro')
        end
    case 2
        % test makeclr
        n = 6;
        r = (0:n)'/n;
        theta = pi*(-n:n)/n;
        X = r*cos(theta);
        Y = r*sin(theta);
        C = r*cos(2*theta);

        cmap = makeclr( 'bwr' );
        
        colormap( cmap );
        pcolor(X,Y,C); shading flat   
        colorbar
end
