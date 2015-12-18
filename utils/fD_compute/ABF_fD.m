function ABF_fD(outs, ins, Nh, Visual, Pth, ModelFlag)

    % Path related
    FigPath = Pth.fDFigPath;
    FilePath = Pth.fDFilePath;    
    
    period = Visual.fDPeriod;
    fDngaM = Visual.fDngaM;    
        
    if ~isdir( FigPath )
        mkdir( FigPath );
    end
    
    FigPath = [FigPath 'Model' num2str(ModelFlag) ];
    if ~isdir( FigPath )
        mkdir( FigPath );
    end
    FigPath = [FigPath '/ModelABF/'];
    if ~isdir( FigPath )
        mkdir( FigPath );
    end
    in = ins{1};
    
    % Mesh related
    % read calculated fD at irregular mesh
    X0 = in.Xs(:,1);
    Y0 = in.Xs(:,2);
    Ndims = in.Ndims;
    
    % Generate regular mesh and do the interpolation of the fD values
    % (keep consitent with what in python)
    Nx = Ndims(2); Ny = Ndims(3); dx = Ndims(4);
    X = (0:1:Nx-1)*dx; Y = (0:1:Ny-1)*dx;
    
    [YY,XX] = meshgrid( Y, X );
    
    % get the fD for given period
    clim = [-0.6,0.6]; 
    index = find( in.Tlist == period );  
    fDxs = zeros( Nh, Nx, Ny );    
    for ih = 1:Nh 
        if ModelFlag == 0
            fD0 = outs{ih}.fD(:,end-index);   
        else
            fD0 = ModelTest(in.M,in.dev,ih,period,in.Tlist,FilePath,ModelFlag);
        end
        fD = griddata( X0,Y0,fD0, XX, YY);    
        %fDxs(ih,:,:) = log(exp(fD));
        fDxs = fD;
    end
    
    % compute fD-<fD>_x
    for ih = 1:Nh
        
        in = ins{ih};
        
        fDxstmp = fDxs(ih,:,:) - 1./Nh * sum(fDxs,1);
        %fDxstmp = fDxs(ih,:,:);
        
        ihypo = ih;
        fig = figure(ihypo);
        clf()
        cn = 'bwr';
        cmap = makeclr(cn);

        % plot fault surface projection and location of the hypocenter (in 2D)
        subs = [2,1];
        basic = [100,100,100,100];   % in pixel;
        FigSize = [0.1,0.1,0.2,0.8]; % this has to be in ratio relative to screen size
        axs = subaxes(fig, subs, basic,FigSize);

        % plot fault
        Nseg = size(in.QQ,2)/3;
        sx = 1; sy = 2; sz = 3;
        inc = 3;

        ax = axes( 'position', axs(1,:) );
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

        % plot fD
        ax = axes( 'position', axs(2,:) );
        colormap(cmap);
        pcolor(XX,YY, squeeze(fDxstmp)); shading interp;
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
        title( [fDngaM ' fD ABF for ihypo ' num2str(ihypo) ' at ' num2str(period) ' sec']  )

        %print('-dpng','-r300',[FigPath FigName] );  
        FigName = ['ABF_fD_hypo' num2str(ihypo) '_' in.dev '_' num2str(period)];
        saveas( fig, [FigPath FigName], 'png' );
        
    end  


end