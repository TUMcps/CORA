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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;
tol = 1e-12;

% empty conHyperplane
hyp = conHyperplane();
if ~isempty(hyp)
    res = false;
end

% random normal vector, offset, constraint matrix, constraint vector
a = [-0.5; -1; 0.1];
b = -1;
hs = halfspace(a,b);
C = [-0.6 0.8 -1.7;...
      0.6 0.5 -0.8];
d = [1; 0.5];

% admissible initializations
% only halfspace
hyp = conHyperplane(hs);
if ~isequal(hyp.h,hs)
    res = false;
end

% normal vector and offset
hyp = conHyperplane(a,b);
if ~isequal(hyp.h,halfspace(a,b))
    res = false;
end

% halfspace and constraint matrix, constraint vector
hyp = conHyperplane(hs,C,d);
if ~isequal(hyp.h,hs) || any(any(abs(hyp.C - C) > tol)) || ...
        any(abs(hyp.d - d) > tol)
    res = false;
end

% normal vector, offset, and constraint matrix, constraint vector
hyp = conHyperplane(a,b,C,d);
if ~isequal(hyp.h,halfspace(a,b)) || ...
        any(any(abs(hyp.C - C) > tol)) || any(abs(hyp.d - d) > tol)
    res = false;
end


% wrong initializations
a_plus1 = [-0.5; -1; 0.1; 3];
b_vec = [-1; 1];
hs_plus1 = halfspace(a_plus1,b);
C_plus1 = [-0.6 0.8 -1.7  0.4;...
            0.6 0.5 -0.8 -0.6];
d_plus1 = [1; 0.5; 2];

% offset as vector
try
    hyp = conHyperplane(a,b_vec); % <- should throw error here
    res = false;
end

% C and d do not fit halfspace
try
    hyp = conHyperplane(hs_plus1,C,d); % <- should throw error here
    res = false;
end
try
    hyp = conHyperplane(hs,C_plus1,d); % <- should throw error here
    res = false;
end 

% C does not fit d
try
    hyp = conHyperplane(hs,C,d_plus1); % <- should throw error here
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


if res
    disp('test_conHyperplane successful');
else
    disp('test_conHyperplane failed');
end

%------------- END OF CODE --------------