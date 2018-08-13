function []=reduceTest(fileID, nrRandZono, dim, order, orderRandZono, dist, maxL)
% check reduction performance
% Author: Anna Kopetzki (adapted from Mathias Althoff)
% Written: April-2016

mpt_init;


filterLength=[order*dim+8, order*dim+3]; %nr. of generators after filtering

% store the ratios redVol/orgVol and times in a File for each method 
fprintf(fileID, '#####################################################\n');
fprintf(fileID, '# Parameters: dim=%d, startOrder=%d, order=%d, distribution=%s, maxL=%d, nrTestedZono=%d\n', dim, orderRandZono, order, dist, maxL, nrRandZono);
fprintf(fileID, '# V-ratio=(redVol/origVol)^(1/dim)\n');


for iTry=1:nrRandZono

    %generate random zonotope using the distribution
    Z=randomZonotope(dim, orderRandZono*dim, dist, maxL);

    tic
    %%compute volume of zonotope
    %origVol=volume(Z);
    toc

    
%     %use method A---------------------------------------------------------    
%     [Zred,tA(iTry)]=reduce(Z,'methA', order);
%     %compute volume ratio
%     redVol=volume(Zred);
%     %ratioA(iTry)=ratio(redVol,origVol,dim);    
%     ratioA(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)


    %disp('girard');
    %use method girard----------------------------------------------------    
    [Zred,tGi(iTry)]=reduce(Z,'girard', order);
    %compute volume ratio
    %redVol=volume(Zred);
    %origVol=redVol;
    %ratioGi(iTry)=ratio(redVol,origVol,dim);
    %ratioGi(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ZG = Zred;
    ratioGi(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)
    
    
   
    %disp('methB');
    %use method B---------------------------------------------------------    
    [Zred,tB(iTry)]=reduce(Z,'methB', order, filterLength);
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioB(iTry)=ratio(redVol,origVol,dim);
    %ratioB(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioB(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)

    
    %disp('methC');
    %use method C---------------------------------------------------------    
    [Zred,tC(iTry)]=reduce(Z,'methC', order, filterLength);
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioC(iTry)=ratio(redVol,origVol,dim);
    %ratioC(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioC(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)
    
    
    
    %disp('pca');
    %use method reduceaPCA: PCA ----------------------------------------
    [Zred,tPCA(iTry)]=reduce(Z,'aPCA', order, filterLength);
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioPCA(iTry)=ratio(redVol,origVol,dim);
    %ratioPCA(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioPCA(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)

    %disp('lineCl');
    %use method reduceCluster: Line setting--------------------------  
    [Zred,tLine(iTry)]=reduce(Z,'cluster', order, filterLength, 4);
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioLine(iTry)=ratio(redVol,origVol,dim);
    %ratioLine(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioLine(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)

    %disp('hybrid');
    %use method reduceCluster: Hybrid setting--------------------------  
    [Zred,tHybrid(iTry)]=reduce(Z,'cluster', order, filterLength, 10);
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioHybrid(iTry)=ratio(redVol,origVol,dim);
    %ratioHybrid(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioHybrid(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)
    
    %disp('cOptdir');
    %use method reduceConstOpt: Constraint optimization ---------------
    [Zred,tCo_ip(iTry)]=reduce(Z,'constOpt', order, filterLength, 'det', 'interior-point');
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioCo_ip(iTry)=ratio(redVol,origVol,dim);
    %ratioCo_ip(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioCo_ip(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)
    
    %disp('cOptSVD');
    %use method reduceConstOpt: Constraint optimization ---------------
    [Zred,tCoSVD_ip(iTry)]=reduce(Z,'constOpt', order, filterLength, 'svd', 'interior-point');
    %compute volume ratio
    %redVol=volume(Zred);
    %ratioCoSVD_ip(iTry)=ratio(redVol,origVol,dim);
    %ratioCoSVD_ip(iTry)= volume_ratio_scaled(Z, Zred); % use scaled volume ratio (see Diss Althoff)
    ratioCoSVD_ip(iTry)= volume_ratio_scaled(ZG, Zred); % use scaled volume ratio (see Diss Althoff)
    
 


    

end




fprintf(fileID, '----------------------------------------------------------------\n');
writeResToFile(fileID, 'Girard___', ratioGi, tGi);
% fprintf(fileID, '----------------------------------------------------------------\n');
%writeResToFile(fileID, 'MethA____', ratioA, tA);
writeResToFile(fileID, 'MethB____', ratioB, tB);
writeResToFile(fileID, 'MethC____', ratioC, tC);
% fprintf(fileID, '----------------------------------------------------------------\n');
writeResToFile(fileID, 'PCA______', ratioPCA, tPCA);
writeResToFile(fileID, 'LineCl___', ratioLine, tLine);
writeResToFile(fileID, 'HybridPC_', ratioHybrid, tHybrid);
% fprintf(fileID, '----------------------------------------------------------------\n');
writeResToFile(fileID, 'CoOptdir_', ratioCo_ip, tCo_ip);
writeResToFile(fileID, 'CoOptSVD_', ratioCoSVD_ip, tCoSVD_ip);








%--------------------------------------------------------------------------
function [result]=ratio(redVol,origVol,dim)
result=(redVol/origVol)^(1/dim);


%--------------------------------------------------------------------------
function [result]=analyse(vector)

result.mean=mean(vector);
result.min=min(vector);
result.max=max(vector);
result.var=var(vector);
result.std=std(vector);
result.median=median(vector);


%--------------------------------------------------------------------------
function []=writeResToFile(idF, methodName, ratioV, tV)
ratio.ratioV=analyse(ratioV);
t.tV=analyse(tV);
fprintf(idF, '#%s:\t meanRatio=%.4f\t minRatio=%.4f\t maxRatio=%.4f\t varRatio=%.4f\t\t stdRatio=%.4f\t medianRatio=%.4f\n', methodName, ratio.ratioV.mean, ratio.ratioV.min, ratio.ratioV.max, ratio.ratioV.var, ratio.ratioV.std, ratio.ratioV.median);
fprintf(idF, '#%s:\t meanTime=%.4f\t minTime=%.4f\t\t maxTime=%.4f\t\t varTime=%.4f\t\t stdTime=%.4f\t medianTime=%.4f\n', methodName, t.tV.mean, t.tV.min, t.tV.max, t.tV.var, t.tV.std, t.tV.median);

function []=writeResToFileRatio(idF, methodName, ratioV)
ratio.ratioV=analyse(ratioV);
fprintf(idF, '#%s:\t meanRatio=%.4f\t minRatio=%.4f\t maxRatio=%.4f\t varRatio=%.4f\t\t stdRatio=%.4f\t medianRatio=%.4f\n', methodName, ratio.ratioV.mean, ratio.ratioV.min, ratio.ratioV.max, ratio.ratioV.var, ratio.ratioV.std, ratio.ratioV.median);

function []=writeResToFileTime(idF, methodName, tV)
t.tV=analyse(tV);
fprintf(idF, '#%s:\t meanTime=%.4f\t minTime=%.4f\t\t maxTime=%.4f\t\t varTime=%.4f\t\t stdTime=%.4f\t medianTime=%.4f\n', methodName, t.tV.mean, t.tV.min, t.tV.max, t.tV.var, t.tV.std, t.tV.median);


