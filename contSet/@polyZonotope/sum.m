function pZsum = sum(pZ,pZlist)
% sum - computes the sum of multiple polynomial zonotopes
%
% Syntax:  
%    pZsum = sum(pZlist)
%
% Inputs:
%    pZ - polyZonotope object (first part of the sum)
%    pZlist - list of polynomial zonotope objects (cell-array)
%
% Outputs:
%    pZsum - resulting polyZonotope object representing the set of the sum
%
% Example: 
%    pZ1 = polyZonotope([1;2],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[2 0 1;0 2 1],[],[0 0 0;0 0 0;1 0 3;0 1 1]);
%
%    pZsum = sum(pZ1,{pZ2,pZ1});
%
%    figure
%    plot(pZ1,[1,2],'r','Filled',true,'EdgeColor','none');
%    figure
%    plot(pZ2,[1,2],'b','Filled',true,'EdgeColor','none');
%    figure
%    plot(pZsum,[1,2],'g','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Niklas Kochdumper
% Written:      17-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % initialize variables
    expMatList = cell(length(pZlist)+1,1);
    G = pZ.G;
    Grest = pZ.Grest;
    c = pZ.c;
    
    % bring exponent matrices to a common representation
    expMatList{1} = pZ.expMat;
    id = pZ.id;
    
    for i = 1:length(pZlist)
        [id,~,expMatList{i+1}] = mergeExpMatrix(id,pZlist{i}.id, ...
                                        expMatList{i},pZlist{i}.expMat);
    end

    % concatenate exponent and generator matrices and sum up center
    M = length(id);
    expMat = [pZ.expMat;zeros(M-size(pZ.expMat,1),size(pZ.expMat,2))];
    
    for i = 1:length(pZlist)
       G = [G,pZlist{i}.G];
       Grest = [Grest,pZlist{i}.Grest];
       eTemp = expMatList{i+1};
       expMat = [expMat,[eTemp;zeros(M-size(eTemp,1),size(eTemp,2))]];
       c = c + pZlist{1}.c;
    end

    % add up all generators that belong to identical exponents
    [ExpNew,Gnew] = removeRedundantExponents(expMat,G);
    
    % construct the resulting polynomial zonotope
    pZsum = polyZonotope(c,Gnew,Grest,ExpNew,id);

end

%------------- END OF CODE --------------