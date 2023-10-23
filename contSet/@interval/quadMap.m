function I = quadMap(varargin)
% quadMap - computes the quadratic map of an interval
%
% Syntax:
%    res = quadMap(I1,Q)
%    res = quadMap(I1,I2,Q)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    I - interval object
%
% Example: 
%    Z = zonotope([0 1 1;0 1 0]);
%    I = interval(Z);
%    Q{1} = [0.5 0.5; 0 -0.5];
%    Q{2} = [-1 0; 1 1];
% 
%    res1 = quadMap(Z,Q); % zonotope quadMap (for comparison)
%    res2 = quadMap(I,Q); % interval quadMap
%
%    figure; hold on;
%    plot(res1,[1,2],'b');
%    plot(res2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap

% Authors:       Niklas Kochdumper
% Written:       17-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute quadratic map for zonotopes
if nargin == 1
    throw(CORAerror('CORAerror:notEnoughInputArgs',2));
elseif nargin == 2
    I = quadMap(zonotope(varargin{1}),varargin{2});
elseif nargin == 3
    I = quadMap(zonotope(varargin{1}),zonotope(varargin{2}),varargin{3});
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% enclose the result with an interval
I = interval(I);

% ------------------------------ END OF CODE ------------------------------
