function res = test_conZonotope_interval
% test_conZonotope_interval - unit test function for the caclulation of
%                             a bounding box of a constrained zonotope object
%
% Syntax:
%    res = test_conZonotope_interval
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper
% Written:       22-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% TEST 1: Figure 1 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZono = conZonotope(Z,A,b);

% calculate interval
I = interval(cZono);

% % plot the result
% plot(cZono,[1,2],'r');
% hold on
% plot(int,[1,2],'b');

% compare with ground-truth
I_ = interval([-2.5;-1.5],[3.5;2.5]);

if ~isequal(I,I_)
    res = false;
end


% TEST 2: Figure 2 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1 0 1;0 1 2 -1];
A = [-2 1 -1];
b = 2;
cZono = conZonotope(Z,A,b);

% calculate interval
I = interval(cZono);

% % plot the result
% plot(cZono,[1,2],'r');
% hold on
% plot(int,[1,2],'b');

% compare with ground-truth
I_ = interval([-2;-2],[0;3]);

if ~isequal(I,I_)
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
