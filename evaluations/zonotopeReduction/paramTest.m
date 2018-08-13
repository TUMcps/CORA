function []=paramTest()
% check reduction performance for several parameters
% Author: Anna Kopetzki (adapted from Mathias Althoff)
% Written: April-2016

mpt_init;

c=clock;
nameF=sprintf('Test2_%f_%d.txt',c(6), rand(1));


fileID = fopen(nameF, 'w');
fprintf(fileID, '# Test: Methods for Zonotope reduction\n');
fprintf(fileID, '# Constraint optimization: minimize volume i.e. det(C) or Frobenius norm, such that C^-1*G <= 1\n');
fprintf(fileID, '# Test Over-approximation methods\n');



% Test parameters (gamma=gaussian distr.)
maxL=100;
nrRandZono=100;
distData=['uniform    '; 'exponential'; 'gamma      ']; % uniform, exponential, gamma
distArr=cellstr(distData);
dimArr= [5 10 15]; 
orderRandZonoAll= [5 10 15]; 
orderArr=[1]; 



for ds=1:length(distArr)
    dist=char(distArr(ds));
    for di=1:length(dimArr)
        dim=dimArr(di);
        for or=1:length(orderArr)
            order=orderArr(or);
            % order of the random Zonotope must be > desired order
            orderRandZonoArr=orderRandZonoAll(orderRandZonoAll > order); 
            if length(orderRandZonoArr) > 0
                for oz=1:length(orderRandZonoArr)
                    orderRandZono=orderRandZonoArr(oz);
                    reduceTest(fileID, nrRandZono, dim, order, orderRandZono, dist, maxL);
                end
            end
        end
    end
end

