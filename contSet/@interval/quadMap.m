function res = quadMap(varargin)
% quadMap - computes the quadratic map of an interval
%
% Syntax:  
%    res = quadMap(Int1,Q)
%    res = quadMap(Int1,Int2,Q)
%
% Inputs:
%    Int1 - interval object
%    Int2 - interval object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    res - interval object
%
% Example: 
%    zono = zonotope([0 1 1;0 1 0]);
%    int = interval(zono);
%    Q{1} = [0.5 0.5; 0 -0.5];
%    Q{2} = [-1 0; 1 1];
%
%    res1 = quadMap(zono,Q);
%    res2 = quadMap(int,Q);
%
%    figure; hold on
%    plot(res1,[1,2],'b');
%    plot(res2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap

% Author:       Niklas Kochdumper
% Written:      17-December-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % compute quadratic map for zonotopes
    if nargin == 2
        res = quadMap(zonotope(varargin{1}),varargin{2});
    else
        res = quadMap(zonotope(varargin{1}),zonotope(varargin{2}), ...
                      varargin{3});
    end
    
    % enclose the result with an interval
    res = interval(res);

%------------- END OF CODE --------------