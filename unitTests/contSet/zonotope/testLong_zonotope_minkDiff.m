function res = testLong_zonotope_minkDiff
% testLong_zonotope_minkDiff - unit test function of minkDiff
%
% Syntax:
%    res = testLong_zonotope_minkDiff
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
% Written:       25-May-2023
% Last update:   24-January-2024 (MW, increase tolerance)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-10;

% quick test for all methods for syntax errors ----------------------------
% very simple tests; TODO: needs thorough testing
% can possible be improved with new contains functionality:
%     test containment, ignore if result not certified.

methods = {'inner','outer','outer:coarse', ...
    'approx','inner:conZonotope','inner:RaghuramanKoeln'};

for i=1:length(methods)
    try 
        % init
        n = 3;
        difference = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
        subtrahend = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
        minuend = difference + subtrahend;
        
        % minkDiff
        diff_comp = minkDiff(minuend,subtrahend,methods{i});

    catch ME
        assertLoop(false,ME.message,[],i);
    end
end

% method 'exact' ----------------------------------------------------------
method = 'exact';

% test 2-dimensional zonotopes
n=2; 
for i=1:10
    % init
    difference = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
    subtrahend = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
    minuend = difference + subtrahend;

    % minkDiff
    diff_comp = minkDiff(minuend,subtrahend,method);

    % diff + sub = min => min - sub = diff
    assert(isequal(difference,diff_comp,tol));

end

% test aligned
for i=1:10
    % init
    n = randi(10);
    minuend = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
    subtrahend = minuend*rand(1); % scale down

    % minkDiff
    difference = minkDiff(minuend,subtrahend,method);

    % min - sub = diff => diff + sub = min
    minuend_restored = zonotope(difference.c+subtrahend.c, difference.G+subtrahend.G);
    assert(isequal(minuend_restored,minuend,tol));
end

% test error for exact
ns = [3,5,10];
for n=ns
    minuend = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
    subtrahend = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
    assertThrowsAs(@minkDiff,'CORA:wrongValue',minuend,subtrahend,method);
end

% method 'outer:scaling' --------------------------------------------------
method='outer:scaling';

for i=1:10
    % init
    n = randi(10);
    difference = zonotope.generateRandom('Dimension',n,'NrGenerators',2*n);
    subtrahend = interval.generateRandom('Dimension',n);
    minuend = difference + subtrahend;
    
    % minkDiff
    assertThrowsAs(@minkDiff,'',minuend,subtrahend,method);
end

% test non full-dimensional zonotopes -------------------------------------
n = 2;
for i=1:10
    % init
    minuend = zonotope(interval.generateRandom('Dimension',n));
    subtrahend = minuend*rand(1); % scale down

    % project to higher dimensions
    P = [1 0; 0 1; 1 1];
    minuend_proj = P*minuend;
    subtrahend_proj = P*subtrahend;

    % minkDiff
    diff_comp_proj = minkDiff(minuend_proj,subtrahend_proj,'exact');

    % compute svd projection matrix 
    [U, S] = svd(generators(minuend_proj));
    newDim = nnz(~all(withinTol(S,0))); % nr. of new dimensions
    P_minuend = U(1:newDim, :); % projection matrix

    % project down using svd projection matrix
    minuend = P_minuend * minuend_proj;
    subtrahend = P_minuend * subtrahend_proj;
    diff_comp = P_minuend * diff_comp_proj;

    % compute exact difference in projected space
    difference = minkDiff(minuend,subtrahend,'exact');

    % test equality
    assertLoop(isequal(difference,diff_comp,tol),i);
end


% gather results ----------------------------------------------------------
res = true;

% ------------------------------ END OF CODE ------------------------------
