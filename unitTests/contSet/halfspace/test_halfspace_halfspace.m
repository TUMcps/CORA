function res = test_halfspace_halfspace
% test_halfspace_halfspace - unit test function of halfspace
%
% Syntax:
%    res = test_halfspace_halfspace
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty halfspace
n = 2;
hs = halfspace.empty(n);
res(end+1,1) = representsa(hs,'emptySet') && dim(hs) == n;

% random normal vector, offset
c = [3; 2; -1];
d = 1;

% admissible initialization
hs = halfspace(c,d);
res(end+1,1) = compareMatrices(hs.c,c);
res(end+1,1) = withinTol(hs.d,d);

hs = halfspace(c',d);
res(end+1,1) = compareMatrices(hs.c,c);
res(end+1,1) = withinTol(hs.d,d);

% combine results
res = all(res);


% wrong initializations
if CHECKS_ENABLED

c_mat = [3 5 2; 0 4 1; 2 4 2];
d_vec = [5 -2];
d_mat = [2 5; -4 2];

% only one input
try
    hs = halfspace(c); % <- should throw error here
    res = false;
end

% normal vector is a matrix
try
    hs = halfspace(c_mat,d); % <- should throw error here
    res = false;
end

% offset is a vector
try
    hs = halfspace(c,d_vec); % <- should throw error here
    res = false;
end

% offset is a matrix
try
    hs = halfspace(c,d_mat); % <- should throw error here
    res = false;
end

% too many input arguments
try
    hs = halfspace(c,d,d); % <- should throw error here
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
