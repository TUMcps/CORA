function res = test_mptPolytope_plot
% test_mptPolytope_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked
%
% Syntax:  
%    res = test_mptPolytope_plot
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

% Author:       Mark Wetzlinger
% Written:      25-May-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% instantiate polytope (via conversion from zonotope)
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
P = mptPolytope(Z);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(P);
    
    % two arguments: object, dimensions
    plot(P,1);
    plot(P,[1,2]);
    plot(P,[2,3]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

%------------- END OF CODE --------------