function res = testLong_polyZonotope_fhandle
% testLong_polyZonotope_fhandle - unit test function for computing
%    the function handle of a polynomial zonotope
%
% Syntax:
%    res = testLong_polyZonotope_fhandle
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
% Written:       24-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nTests = 50;

for i=1:nTests

    % random dimension
    n = randi(30);

    % init polynomial zonotope
    pZ = noIndep(polyZonotope.generateRandom('Dimension',n));

    % number of identifiers
    ne = length(pZ.id);

    f = fhandle(pZ);
    N = 2*n;
    B = 2*rand(ne,N)-1;
    for j=1:N
        fval_res = resolve(pZ,B(:,j));
        fval = f(B(:,j));
        f_abs = max(abs(fval_res),abs(fval));
        f_abs(f_abs==0) = 1;

        % check result
        if ~all(withinTol(fval_res./f_abs, fval./f_abs, 1e-6))
            res = false;
            return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
