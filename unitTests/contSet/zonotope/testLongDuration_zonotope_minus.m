function res = testLongDuration_zonotope_minus
% testLongDuration_zonotope_minus - unit test function of minus for approximating 
% the Minkowski difference of two zonotopes or a zonotope with a vector
% according to [1]. 
% A - B = C <-> B + C \subseteq A
%
% Syntax:  
%    res = testLongDuration_zonotope_minus
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Example: 
%    Z1 = zonotope([1 2 2; 0 0 2]);
%    Z2 = zonotope([0 0.5 0.5 0.3; 1 0 0.5 0.2]);
%
%    Z3 = Z1 - Z2;
%    Z4 = Z2 + Z3;
%
%    figure; hold on;
%    plot(Z1,[1 2], 'b');
%    plot(Z2,[1 2], 'r');
%    plot(Z3,[1 2], 'g');
%    plot(Z4,[1 2], 'k');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, conZonotope/minus

% Author:       Matthias Althoff
% Written:      06-May-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%% create zonotopes
Z_m = zonotope([1 1 0 1; 1 0 1 1]);
Z_m_degenerate = zonotope([1 1; 1 0]);
Z_s{1} = zonotope([0 0.5 0; 0 -0.2 0.2]); % see Fig. 2a in [1]
Z_s{2} = zonotope([0 0.5 0; 0 -0.5 0.5]); % see Fig. 2b in [1]
Z_s{3} = zonotope([0 2 0; 0 -0.5 0.5]); % see Fig. 2c in [1]
Z_s{4} = zonotope([0 0.5 0; 0 0 0.5]);
Z_s{5} = zonotope([0 1 0; 0 0 0.5]);
Z_s{6} = zonotope([0 2 0; 0 0 0.5]);
% the following Z_s seem to be incorrect for the MPT toolbox
% Z_s{7} = zonotope([-0.5 0.5 0; 1 0 0.5]); % non-zero center; seems incorrect for MPT toolbox
% Z_s{8} = zonotope([-0.5 1 0; 1 0 0.5]); % non-zero center; seems incorrect for MPT toolbox
% Z_s{9} = zonotope([-0.5 2 0; 1 0 0.5]); % non-zero center; seems incorrect for MPT toolbox

% convert Z_m to polytope
P_m = polytope(Z_m);

% initialize partial results
resPartial = [];

% define small box
smallBox = zonotope([[0;0],1e-8*eye(2)]);

%% loop through all subtrahends
for iSet = 1:length(Z_s)
    
    % compute result
    Z_res = Z_m - Z_s{iSet};
    
    % result from polytopes
    % IMPORTANT: Minkowski difference of MPT toolbox seems to be incorrect
    % for Z_s{7}-Z_s{9}
    P_res = P_m - polytope(Z_s{iSet});
    
    % check whether Minkowski difference returns the empty set
    if isempty(Z_res)
        % check if polytope solution is empty as well
        resPartial(end+1) = isempty(P_res);
    else
        % enclosure check (Definition of Minkoswki differnce)
        resPartial(end+1) = in(Z_m + smallBox, Z_res + Z_s{iSet});

        % enclosure check (comparison with polytope solution)
        resPartial(end+1) = in(P_res + smallBox, polytope(Z_res));
        
        % enclosure check (comparison with polytope solution; other direction)
        resPartial(end+1) = in(polytope(Z_res) + smallBox, P_res);
        
        % for debugging:
%         figure
%         hold on
%         plot(Z_res);
%         plot(P_res, [1 2], 'r');
    end
    
end

%% check the case when the minuend is a degenerate zonotope
% compute result
Z_res = Z_m_degenerate - Z_s{1};
    
% result from polytopes
P_res = P_m - polytope(Z_s{iSet});

% the result should be empty
if isempty(Z_res)
    % check if polytope solution is empty as well
    resPartial(end+1) = isempty(P_res);
end

%result of all tests
res = all(resPartial);

%------------- END OF CODE --------------