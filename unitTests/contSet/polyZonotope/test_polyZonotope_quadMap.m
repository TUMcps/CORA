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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% create polynomial matrices
c = [1;2];
G = [1 -2 1;2 3 -1];
Grest = [0;0];
expMat = [1 0 2; 0 1 1];
pZ = polyZonotope(c,G,Grest,expMat);

% create matrices of the quadratic map
Q{1} = [1 2;-1 2];
Q{2} = [-3 0;1 1];

% calculate quadratic map
pZres = quadMap(pZ,Q);

% define ground truth
G = [22 19 -5 11 19 -5 16 -11 2; 6 23 -9 3 23 -9 -9 11 -3];
c = [11;3];
expMat = [1 0 2 2 1 3 0 2 4; 0 1 1 0 1 1 2 2 2];

% check for correctness
if ~all(pZres.c == c)
    error('test_polyZonotope_quadMap: analytical test failed!');
end

for i = 1:size(expMat,2)
    
    ind = ismember(pZres.expMat',expMat(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        error('test_polyZonotope_quadMap: analytical test failed!');        
    elseif ~all(pZres.G(:,ind_(1)) == G(:,i))
        error('test_polyZonotope_quadMap: analytical test failed!');
    end
    
end


res = true;

%------------- END OF CODE --------------