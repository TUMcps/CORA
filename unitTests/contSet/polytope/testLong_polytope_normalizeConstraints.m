function res = testLong_polytope_normalizeConstraints
% testLong_polytope_normalizeConstraints - unit test function for
%    constraint normalization
%
% Syntax:
%    res = testLong_polytope_normalizeConstraints()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       04-December-2022
% Last update:   13-December-2022 (MW, add cases for type = 'A')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests
    
    % random dimension
    n = randi(10);

    % instantiation: generateRandom
    P = polytope.generateRandom('Dimension',n);

    % normalize offset in constraints
    P_ = normalizeConstraints(P);

    % check if normalization correct
    if ~isempty(P_.A)
        assertLoop(all( withinTol(P_.b,1) | withinTol(P_.b,0) | withinTol(P_.b,-1)),i)
    end
    % check if normalization correct
    if ~isempty(P_.Ae)
        assertLoop(all( withinTol(P_.be,1) | withinTol(P_.be,0) ),i)
    end

    % normalize vectors in constraints
    P_ = normalizeConstraints(P,'A');

    % check if normalization correct
    if ~isempty(P_.A)
        assertLoop(all(withinTol(vecnorm(P_.A',2,1),1)),i)
    end
    % check if normalization correct
    if ~isempty(P_.Ae) 
        assertLoop(all(withinTol(vecnorm(P_.Ae',2,1),1)),i)
    end

end

% ------------------------------ END OF CODE ------------------------------
