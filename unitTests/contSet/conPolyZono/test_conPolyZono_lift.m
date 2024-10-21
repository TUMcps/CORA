function res = test_conPolyZono_lift
% test_conPolyZono_lift - unit test function of lift
%
% Syntax:
%    res = test_conPolyZono_lift
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

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];

cPZ = conPolyZono(c,G,E,A,b,EC);

% should always fail
try
    cPZ_lift = lift(cPZ,5,[2,1]);
    % should fail
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% except no new dimensions are created
cPZ_lift = lift(cPZ,2,[2,1]);
cPZ_proj = project(cPZ,[2,1]);
assert(isequal(cPZ_proj,cPZ_lift));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
