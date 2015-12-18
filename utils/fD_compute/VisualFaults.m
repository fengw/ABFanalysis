function VisualFaults(in,FigPath, FigName)
 
    if ~isdir(FigPath)
        mkdir(FigPath);
    end 

    figure(1)
    clf()
    [nr,nc] = size(in.QQ);
    sx = 1; sy = 2; sz = 3;
    inc = 3;
    for i=1:nc/3
        %clr=[i-1,abs(i-1-4.5),nc/3-i+1]/(nc/3)   % random color
        clr='k';
        plot3([in.QQ(:,sx);in.QQ(1,sx)],[in.QQ(:,sy);in.QQ(1,sy)],[in.QQ(:,sz);in.QQ(1,sz)],'Color',clr,'LineWidth',1.5)
        hold on;
        sx = sx + inc; sy = sy + inc; sz = sz + inc;
    end
    plot3(in.H(1),in.H(2),in.H(3), 'ro')
    hold on;
    grid on;
    axis('equal')
    xlabel('CS-X')
    ylabel('CS-Y')
    zlabel('Depth')

    title( 'fault segments (3D veiw)' )
    print('-dpdf','-r200',[FigPath FigName] );    

end