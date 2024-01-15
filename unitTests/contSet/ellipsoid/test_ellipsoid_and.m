function res = test_ellipsoid_and
% test_ellipsoid_and - unit test function of and
%
% Syntax:
%    res = test_ellipsoid_and
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

% Authors:       Victor Gassmann
% Written:       26-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
% empty set
load cases.mat E_c

% intersection with empty set
% res = (E_c{1}.E1 & E_e) == E_e;

for i=1:length(E_c)
    E1 = E_c{i}.E1; %non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    h1 = E_c{i}.h1;
    
    % test ellipsoid (degenerate)
    Eres = E1 & Ed1;
    Yd = randPoint(Ed1,2*i);
    for j=1:size(Yd,2)
        if contains(E1,Yd(:,j)) && ~contains(Eres,Yd(:,j))
            res = false;
            break;
        end
    end
    
    % test ellipsoid (all zero)
    Eres_0 = E1 & E0;
    if contains(E1,E0.q)
        if rank(Eres_0)~=0 || ~withinTol(Eres_0.q,E0.q,Eres_0.TOl)
            res = false;
            break;
        end
    else
        if ~representsa_(Eres_0,'emptySet',eps)
            res = false;
            break;
        end
    end
    
    % test halfspace
    Eres_h = E1 & h1;
    
    Y = randPoint(E1,2*i);
    for j=1:size(Y,2)
        if contains(h1,Y(:,j)) && ~contains(Eres_h,Y(:,j))
            res = false;
            break;
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
