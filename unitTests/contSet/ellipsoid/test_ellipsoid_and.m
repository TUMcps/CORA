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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; %non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    h1 = E_c{i}.h1;
    
    % test ellipsoid (degenerate)
    Eres = E1 & Ed1;
    Yd = randPoint(Ed1,2*i);
    for j=1:size(Yd,2)
        if in(E1,Yd(:,j)) && ~in(Eres,Yd(:,j))
            res = false;
            break;
        end
    end
    
    % test ellipsoid (all zero)
    Eres_0 = E1 & E0;
    if in(E1,E0.q)
        if rank(Eres_0)~=0 || ~withinTol(Eres_0.q,E0.q,Eres_0.TOl)
            res = false;
            break;
        end
    else
        if ~isempty(Eres_0)
            res = false;
            break;
        end
    end
    
    % test halfspace
    Eres_h = E1 & h1;
    
    Y = randPoint(E1,2*i);
    for j=1:size(Y,2)
        if in(h1,Y(:,j)) && ~in(Eres_h,Y(:,j))
            res = false;
            break;
        end
    end
    
end

if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
