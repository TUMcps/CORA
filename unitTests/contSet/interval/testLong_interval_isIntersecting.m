function res = testLong_interval_isIntersecting
% testLong_interval_isIntersecting - unit test function of
%    isIntersecting, note: only interval-to-interval tested!
%
% Syntax:
%    res = testLong_interval_isIntersecting
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

% Authors:       Mark Wetzlinger
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Random cases
res = true;
nrOfTests = 1000;

for i=1:nrOfTests

    % random dimension
    n = randi(50);

    % init random interval in [-1,1]
    lb = -rand(n,1);
    ub = rand(n,1);
    I = interval(lb,ub);

    % shift by center -> no intersection
    shift = ones(n,1);
    I_low = I - shift;
    I_high = I + shift;

    % check with correct solution
    if isIntersecting(I_low,I_high)
        res = false; break;
    end

    % interval [-1,1]
    I = interval(-ones(n,1),ones(n,1));
    % shift by max. 1
    I_low = I - rand(1);
    I_high = I + rand(1);

    % check with correct solution
    if ~isIntersecting(I_low,I_high)
        res = false; break;
    end
    
    % intervals meet at origin
    I = interval(-ones(n,1),ones(n,1));
    shift = ones(1,1);
    I_low = I - shift;
    I_high = I + shift;

    if ~isIntersecting(I_low,I_high)
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
