function res = test_polyZonotope_split
% test_polyZonotope_split - unit test function for the splitting of a
%    polynomial zonotope object
%
% Syntax:  
%    res = test_polyZonotope_split
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
% Written:      26-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create polynomial zonotopes
c = [-1;3];
G = [2 0 1; 1 2 1];
expMat = [1 0 2; 0 1 1];
Grest = [];
pZ = polyZonotope(c,G,Grest,expMat);

% plotRandPoint(pZ,[1,2],100000,'.r');

% split the polynomial zonotope at the longest generator
pZsplit = splitLongestGen(pZ);

% define ground truth
c1 = [0; 3.5];
G1 =  [1 1/4 1/2 1/4; 1/2 9/4 1/2 1/4];
expMat1 = [1 0 1 2; 0 1 1 1];

c2 = [-2; 2.5];
G2 =  [1 1/4 -1/2 1/4; 1/2 9/4 -1/2 1/4];
expMat2 = [1 0 1 2; 0 1 1 1];

% check for correctness
if ~all(withinTol(c1,pZsplit{1}.c)) || ~all(withinTol(c2,pZsplit{2}.c))
    throw(CORAerror('CORA:testFailed'));
end

for i = 1:size(expMat1,2)    
    
    ind = ismember(pZsplit{1}.expMat',expMat1(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        throw(CORAerror('CORA:testFailed'));
    elseif ~all(pZsplit{1}.G(:,ind_(1)) == G1(:,i))
        throw(CORAerror('CORA:testFailed'));
    end
end

for i = 1:size(expMat2,2)    
    
    ind = ismember(pZsplit{2}.expMat',expMat2(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        throw(CORAerror('CORA:testFailed'));
    elseif ~all(pZsplit{2}.G(:,ind_(1)) == G2(:,i))
        throw(CORAerror('CORA:testFailed'));
    end
end


res = true;

%------------- END OF CODE --------------