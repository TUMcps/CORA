function res = test_polyZonotope_polyMap
% test_polyZonotope_polyMap - unit test function of polyMap
%
% Syntax:
%    res = test_polyZonotope_polyMap
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

% Authors:       Niklas Kochdumper
% Written:       23-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % polynomial zonotope
    pZ = polyZonotope(0,1,[],[2;3]);
    
    % compute polynomial map
    pZres = polyMap(pZ,1,2);

    % check correctness
    if pZres.c ~= 0 || pZres.G ~= 1 || ~all(all(pZres.E == [4;6]))
        throw(CORAerror('CORA:testFailed'));
    end

    res = true;
end

% ------------------------------ END OF CODE ------------------------------
