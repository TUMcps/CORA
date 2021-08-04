function res = test_polyZonotope_reduce
% test_polyZonotope_reduce - unit test function of order reduction
%
% Syntax:  
%    res = test_polyZonotope_reduce
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
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1

% create polynomial zonotope
c = [0;0];
G = [0 2 1;3 0 1];
expMat = [1 0 1;0 1 1];
Grest = [0;0];
pZ = polyZonotope(c,G,Grest,expMat);

% reduce the polynomial zonotope
pZred = reduce(pZ,'girard',1);

% define ground truth
Grest =  [3 0; 0 4];
c = [0;0];

% check for correctness
if any(c-pZred.c) || any(any(Grest-pZred.Grest)) || ~isempty(pZred.G) || ...
   ~isempty(pZred.expMat) || ~isempty(pZred.id)
    error('test_polyZonotope_reduce: analytical test 1 failed!');
end



% TEST 2

% create polynomial zonotope
c = [0;0];
G = [0 3 1 1; 2 0 1 -3];
expMat = [1 0 1 2;0 1 1 0];
Grest = [-1 -4; -2 -1];
pZ = polyZonotope(c,G,Grest,expMat);

% reduce the polynomial zonotope
pZred = reduce(pZ,'girard',2);

% define ground truth
c = [0.5;-1.5];
G = [3;0];
Grest = [-4 2.5 0; -1 0 6.5];
expMat = 1;

% check for correctness
if any(c-pZred.c) || any(any(G-pZred.G)) || any(any(Grest-pZred.Grest)) || ...
   any(any(expMat-pZred.expMat))
    error('test_polyZonotope_reduce: analytical test 2 failed!');
end


res = true;

%------------- END OF CODE --------------