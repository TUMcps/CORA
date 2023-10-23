function res = test_polyZonotope_quadMap
% test_polyZonotope_quadMap - unit test function of quadMap
%
% Syntax:
%    res = test_polyZonotope_quadMap
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
% Written:       23-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate polynomial zonotope
c = [1;2];
G = [1 -2 1;2 3 -1];
GI = [0;0];
E = [1 0 2; 0 1 1];
pZ = polyZonotope(c,G,GI,E);

% create matrices of the quadratic map
Q{1} = [1 2;-1 2];
Q{2} = [-3 0;1 1];

% calculate quadratic map
pZres = quadMap(pZ,Q);

% define ground truth
G = [22 19 -5 11 19 -5 16 -11 2; 6 23 -9 3 23 -9 -9 11 -3];
c = [11;3];
E = [1 0 2 2 1 3 0 2 4; 0 1 1 0 1 1 2 2 2];

% check for correctness
if ~all(withinTol(pZres.c,c))
    throw(CORAerror('CORA:testFailed'));
end

for i = 1:size(E,2)
    
    ind = ismember(pZres.E',E(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        throw(CORAerror('CORA:testFailed'));
    elseif ~all(pZres.G(:,ind_(1)) == G(:,i))
        throw(CORAerror('CORA:testFailed'));
    end
    
end

res = true;

% ------------------------------ END OF CODE ------------------------------
