function res = test_polyZonotope_cubMap
% test_polyZonotope_cubMap - unit test function for the cubic 
%    multiplication of polynomial zonotopes with a tensor
%
% Syntax:  
%    res = test_polyZonotope_cubMap
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      17-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% TEST cubic multiplication

% define polynomial zonotope
pZ = polyZonotope([0;1],[1 -1;2 0],[],eye(2));

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
pZres = cubMap(pZ,T);

% define ground truth
temp = [2 13 -1 28 -4 21 -7 3 -1];
Z_ = [temp;temp];
c_ = Z_(:,1);
G_ = Z_(:,2:end);
expMat_ = [1 0 2 1 3 2 1 0;...
           0 1 0 1 0 1 2 3];
       
% check for correctness
if ~all(withinTol(pZres.c,c_)) || any(size(pZres.G)-size(G_)) ...
        || any(size(pZres.expMat)-size(expMat_))
    throw(CORAerror('CORA:testFailed'));
else
    for i = 1:size(expMat_,2)    

        ind = ismember(pZres.expMat',expMat_(:,i)','rows');  
        ind_ = find(ind > 0);

        if isempty(ind_)
            throw(CORAerror('CORA:testFailed'));
        elseif ~all(pZres.G(:,ind_(1)) == G_(:,i))
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

res = true;

%------------- END OF CODE --------------