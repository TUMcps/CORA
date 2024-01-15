function res = test_polyZonotope_plus
% test_polyZonotope_plus - unit test function for the addition of two
%    polynomial zonotope objects
%
% Syntax:
%    res = test_polyZonotope_plus
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

% Authors:       Niklas Kochdumper
% Written:       26-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% create polynomial zonotopes
c = [-1;3];
G = [-1 -1 3 2; 2 -1 -1 0];
E = [1 0 1 5; 0 1 3 0];
GI = [];
pZ1 = polyZonotope(c,G,GI,E);

c = [0;2];
G = [-2 -2 -1; -1 -2 -3];
E = [1 0 2; 0 0 0; 0 1 3];
GI = [];
pZ2 = polyZonotope(c,G,GI,E);

% exact addition of the two polynomial zonotopes
pZres = exactPlus(pZ1,pZ2);

% define ground truth
c = [-1; 5];
G =  [-3 -1 3 2 -2 -1; 1 -1 -1 0 -2 -3];
E = [1 0 1 5 0 2; 0 1 3 0 0 0; 0 0 0 0 1 3];

% check for correctness
res(end+1,1) = all(withinTol(c,pZres.c));

for i = 1:size(E,2)    
    
    ind = ismember(pZres.E',E(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        throw(CORAerror('CORA:testFailed'));
    elseif ~all(pZres.G(:,ind_(1)) == G(:,i))
        throw(CORAerror('CORA:testFailed'));
    end
end

% empty set
res(end+1,1) = representsa(pZ1 + polyZonotope.empty(2),'emptySet');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
