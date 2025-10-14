function res = test_linearSysDT_linearSys
% test_linearSysDT_linearSys - unit test for the conversion to equivalent
%   continuous-time system
%
% Syntax:
%    res = test_linearSysDT_linearSys
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
% See also: none

% Authors:       Niklas Kochdumper
% Written:       17-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % generate random linear discrete-time object
    A = [0.1 0.1; 0.3 0.4];
    B = [0.1; 1];
    c = [0;0];
    C = [1 0];
    dt = 0.4;

    sysDT = linearSysDT(A,B,c,C,dt);

    % convert to continuous-time system
    sys = linearSys(sysDT);

    % convert back to discrete-time system
    sysDT_ = linearSysDT(sys,dt);

    % check if the results are identical
    tol = 1e-5;

    assert(all(all(abs(sysDT.A - sysDT_.A) < tol)));
    assert(all(all(abs(sysDT.B - sysDT_.B) < tol)));
    assert(all(all(abs(sysDT.c - sysDT_.c) < tol)));
    assert(all(all(abs(sysDT.C - sysDT_.C) < tol)));
    assert(all(all(abs(sysDT.D - sysDT_.D) < tol)));
    assert(all(all(abs(sysDT.k - sysDT_.k) < tol)));

    % all checks ok
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
