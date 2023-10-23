function res = testLong_polytope_center()
% testLong_polytope_center - random tests and check if the center 
%    is contained in the polytope every time
%
% Syntax:
%    res = testLong_polytope_center()
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

% Authors:       Viktor Kotsev
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 50;

for i = 1:nrTests

    % random dimension
    n = randi(5);

    %create random polytopes
    P = polytope.generateRandom('Dimension',n);

    %compute center
    c = center(P);

    %check if center is contained in polytope
    if ~contains(P,c)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
