filenam = 'CyberShake.NGAs.3.00.Source10.sigma';

Ny = 23; Nx = 40;

% aleatory uncertainty
data = dlmread( filenam );
gks = data(:,3);
hks = data(:,4);
iks = data(:,5);
jks = data(:,6);

figure(1)

subplot( 221 )
pcolor( reshape(gks,Nx,Ny)' ); shading interp
colorbar(); colormap(jet)

subplot( 222 )
pcolor( reshape(hks,Nx,Ny)' ); shading interp
colorbar(); colormap(jet)
subplot( 223 )
pcolor( reshape(iks,Nx,Ny)' ); shading interp
colorbar(); colormap(jet)
subplot( 224 )
pcolor( reshape(jks,Nx,Ny)' ); shading interp
colorbar(); colormap(jet)

filenam = 'CyberShake.NGAs.3.00.bs';
data = dlmread( filenam );
bs = data(:,3);
figure(2)
pcolor( reshape( bs, Nx, Ny)' );shading interp
colorbar(); colormap(jet)

