function res = test_polytope_randPoint
% test_polytope_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_polytope_randPoint
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
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% empty set
P = polytope();
p = randPoint(P);
res(end+1,1) = isnumeric(p) && isempty(p);


% standard polytope
A = [2 1; -1 2; -2 -3; 3 -1];
b = ones(4,1);
P = polytope(A,b);

% sample random point with different syntaxes
p = randPoint(P);
p(:,end+1) = randPoint(P,1);
p(:,end+1) = randPoint(P,1,'standard');
p(:,end+1) = randPoint(P,1,'extreme');
p(:,end+1:end+5) = randPoint(P,5,'standard');
p(:,end+1:end+5) = randPoint(P,5,'extreme');

res(end+1,1) = all(contains(P,p,'exact',1e-6));

% unbounded polytope
% P = polytope([1 0; -1 0; 0 1], [1;1;1]);
% p = randPoint(P);
% res(end+1,1) = contains(P,p);

% bounded, degenerate polytope 
P = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]);
p = randPoint(P);
res(end+1,1) = contains(P,p,'exact',1e-6);

if ~all(res)
    throw(CORAerror('CORA:testFailed'))
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
