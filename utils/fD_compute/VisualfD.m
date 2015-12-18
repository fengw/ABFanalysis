function VisualfD(in, out, ModelTest, ihypo, Visual, Pth, ModelFlag)

    FigPath = Pth.fDFigPath;
    FilePath = Pth.fDFilePath;    
    
    period = Visual.fDPeriod;
    fDngaM = Visual.fDngaM;
    
    % read calculated fD at irregular mesh
    X0 = in.Xs(:,1);
    Y0 = in.Xs(:,2);

    % get the fD for given period
    index = find( in.Tlist == period );   
    fD0 = out.fD(:,end-index);
 
    % compute the new fD
    fD1 = ModelTest(in.M,in.dev,ihypo,period,in.Tlist,FilePath,ModelFlag);
 
    clim = [-0.6,0.6];
    
    % Generate regular mesh and do the interpolation of the fD values
    % (keep consitent with what in python)
    Ndims = in.Ndims;
    Nx = Ndims(2); Ny = Ndims(3); dx = Ndims(4);
    X = (0:1:Nx-1)*dx; Y = (0:1:Ny-1)*dx;
    
    [YY,XX] = meshgrid( Y, X );

    fD = griddata( X0,Y0,fD0, XX, YY);
    fD11 = griddata( X0,Y0,fD1, XX, YY);
    
    fD = log10(exp(fD)); % convert to log10
    fD11 = log10(exp(fD11));
    
    if ~isdir( FigPath )
        mkdir( FigPath );
    end
    
    FigPath = [FigPath 'Model' num2str(ModelFlag) '/'];
    if ~isdir( FigPath )
        mkdir( FigPath );
    end
    
    % plot fault surface projection and location of the hypocenter (in 2D)
    FigName = ['hypo_' num2str(ihypo) '_' in.dev '_' num2str(period)];
    
    fig = figure(ihypo+1);
    clf()
    cn = 'bwr';
    cmap = makeclr(cn);
    
    if 1
        % plot fD_new and fault geometry
        axs = subaxes( fig, [1,1], [100,100,100,100], [0.1,0.1,0.3,0.7]);
        
        ax = axes( 'position', axs(1,:) );   % subplot position
        
        % plot fault
        Nseg = size(in.QQ,2)/3;
        sx = 1; sy = 2; sz = 3;
        inc = 3;

        for i=1:Nseg
            %multiple segments
            clr=[i-1,abs(i-1-Nseg/2),Nseg-i+1]/Nseg;  % random colors
            
            %clr='k';
            plot3([in.QQ(:,sx);in.QQ(1,sx)],[in.QQ(:,sy);in.QQ(1,sy)],[in.QQ(:,sz);in.QQ(1,sz)],'Color',clr,'LineWidth',1.5)
            hold on;
            sx = sx + inc; sy = sy + inc; sz = sz + inc;
        end
        plot3(in.H(1),in.H(2),in.H(3),  'r*','MarkerFaceColor','r')
        hold on;
        plot( [X(1),X(1),X(end),X(end),X(1)], [Y(1),Y(end),Y(end),Y(1),Y(1)],'k')
        grid on;
        
                        
        % plot fD
        colormap(cmap);
        pcolor(XX,YY,fD11); shading interp;
        caxis(clim)
        colorbar( 'location','southoutside' )
        hold on;
        
        
        axis('equal')
        xlabel('CS-X')
        ylabel('CS-Y')
        zlabel('Depth')
        title( ['fD of Model ' num2str(ModelFlag) ' with Fault Segments for ihypo ' num2str(ihypo)] )

    else
        axs = subaxes( fig, [2,2], [100,100,100,100], [0.1,0.1,0.3,0.7]);        
        % plot fault
        Nseg = size(in.QQ,2)/3;
        sx = 1; sy = 2; sz = 3;
        inc = 3;
        ax = axes( 'position', axs(1,:) );   % subplot position
        for i=1:Nseg
            %clr=[i-1,abs(i-1-Nseg/2),Nseg-i+1]/Nseg;  % random colors for
            %multiple segments
            clr='k';
            plot3([in.QQ(:,sx);in.QQ(1,sx)],[in.QQ(:,sy);in.QQ(1,sy)],[in.QQ(:,sz);in.QQ(1,sz)],'Color',clr,'LineWidth',1.5)
            hold on;
            sx = sx + inc; sy = sy + inc; sz = sz + inc;
        end
        plot3(in.H(1),in.H(2),in.H(3),  'r*','MarkerFaceColor','r')
        hold on;
        plot( [X(1),X(1),X(end),X(end),X(1)], [Y(1),Y(end),Y(end),Y(1),Y(1)],'k')
        grid on;
        axis('equal')
        xlabel('CS-X')
        ylabel('CS-Y')
        zlabel('Depth')
        title( ['Fault Segments for ihypo ' num2str(ihypo)] )       
        
        fDall = {fD, fD11, fD-fD11};
        tname = {'fD','fDnew','\Delta fD'};
        for isub=1:3
            ax = axes( 'position', axs(isub+1,:) );

            % plot fD
            colormap(cmap);
            pcolor(XX,YY,fDall{isub}); shading interp;
            caxis(clim)
            colorbar;
            hold on;

            % plot fault segments (map view)
            sx = 1; sy = 2;
            inc = 3;
            for i=1:Nseg
                clr='k';    % white color
                plot([in.QQ(:,sx);in.QQ(1,sx)],[in.QQ(:,sy);in.QQ(1,sy)],'Color',clr,'LineWidth',1)
                hold on;
                sx = sx + inc; sy = sy + inc;
            end
            plot(in.H(1),in.H(2), 'r*','MarkerFaceColor','r')
            hold on;
            plot( [X(1),X(1),X(end),X(end),X(1)], [Y(1),Y(end),Y(end),Y(1),Y(1)],'k')

            axis([min(X) max(X) min(Y) max(Y)])
            axis('equal')
            xlabel('CS-X')
            ylabel('CS-Y')
            title( [fDngaM ' ' tname{isub} ' Model ' num2str(ModelFlag) ' for ihypo ' num2str(ihypo) ' at ' num2str(period) ' sec']  )
        end
    end
    
    %print('-dpng','-r300',[FigPath FigName] );  
    saveas( fig, [FigPath FigName], 'png' );
    
end

