function res = test_conHyperplane_conHyperplane
% test_conHyperplane_conHyperplane - unit test function of conHyperplane
%
% Syntax:
%    res = test_conHyperplane_conHyperplane
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

% random normal vector, offset, constraint matrix, constraint vector
a = [-0.5 -1 0.1];
b = -1;
C = [-0.6 0.8 -1.7;...
      0.6 0.5 -0.8];
d = [1; 0.5];

% admissible initializations

% normal vector and offset
hyp = conHyperplane(a,b);
res(end+1,1) = all(withinTol(hyp.a,a)) && withinTol(hyp.b,b);

% empty constraint matrix/vector
hyp = conHyperplane(a,b,[],0);
res(end+1,1) = all(withinTol(hyp.a,a)) && withinTol(hyp.b,b) ...
    && isempty(hyp.C);

% normal vector, offset, and constraint matrix, constraint vector
hyp = conHyperplane(a,b,C,d);
res(end+1,1) = all(withinTol(hyp.a,a)) && withinTol(hyp.b,b) ...
        && compareMatrices(hyp.C,C) && compareMatrices(hyp.d,d);


% combine results
res = all(res);


% wrong initializations
if CHECKS_ENABLED

a_plus1 = [-0.5; -1; 0.1; 3];
b_vec = [-1; 1];
C_plus1 = [-0.6 0.8 -1.7  0.4;...
            0.6 0.5 -0.8 -0.6];
d_plus1 = [1; 0.5; 2];

% no input arguments
try
    hyp = conHyperplane(); % <- should throw error here
    res = false;
end

% offset as vector
try
    hyp = conHyperplane(a,b_vec); % <- should throw error here
    res = false;
end 

% a does not fit C
try
    hyp = conHyperplane(a,b,C_plus1,d); % <- should throw error here
    res = false;
end 

% a does not fit d
try
    hyp = conHyperplane(a,b,C,d_plus1); % <- should throw error here
    res = false;
end 

% too many input arguments
try
    hyp = conHyperplane(a,b,C,d,d); % <- should throw error here
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
