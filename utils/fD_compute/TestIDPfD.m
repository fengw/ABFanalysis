% should be the same
wkd = '/Users/fengw/work/Project/CyberShake_analysis/utils/fD_compute/';
sid=255; rid = 64; Nh = 4;
sigma='1.00_1.00';

outputs = [ wkd 'outputs/ERF35_RupVar3/' 'SourceID' num2str(sid) '_RuptureID' num2str(rid) '/'];
figure(1)

if 1
    fD = zeros(235,Nh);
    fD_ave = zeros(235,1);
    for ihypo = 1:Nh 
        filen = [outputs 'hypo' num2str(ihypo) '_BA6.txt' ];
        inputs = dlmread( filen, ',', 1, 0 );
        X = inputs(:,2); 
        Y = inputs(:,3); 
        IDP = inputs(:,11);   
        fD(:,ihypo) = inputs(:,end-2);
        fD_ave = fD_ave + fD(:,ihypo); 
        scatter(X,Y,100,fD,'filled')
        title(['fD at hypo ' num2str(ihypo)])    
        pause( 1.0 );
    end
    fD_ave = fD_ave ./ Nh; 

    figure(2) 
    scatter(X,Y,100,fD_ave,'filled')   % show the averaged pure directivity effect

    outputsPy = ['/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups54/ERF35_SGT5_RupVar3_Vel1/Gkxmfs0/Dkxs/Sigma' sigma '/' num2str(sid)];
    figure(3)

    for ihypo = 1:Nh
        filen = [outputsPy '/CyberShake.NGAs.3.00.Source' num2str(sid) '.Ih' num2str(ihypo-1) '.dkxs'];
        inputs = dlmread( filen, ' ', 0, 0);
        X = inputs(:,1); Y = inputs(:,2);
        IDP = inputs(:,10);    
        fD = inputs(:,9);  % normalized fD
        scatter(X,Y,100,fD,'filled') 
        title(['normalized fD at hypo ' num2str(ihypo)])    
        pause( 1.0 );
    end
end
