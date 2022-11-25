function res = test_interval_interval
% test_interval_interval - unit test function of interval
%
% Syntax:  
%    res = test_interval_interval
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

% empty interval
I = interval();
res = isempty(I);


% random lower bound, random upper bound
a = [-2; -3];
b = [3; 5];
a_mat = [-2 -1; 0 -3];
b_mat = [3 2; 5 3];

% admissible initializations
I = interval(a,b);
if any(abs(I.inf - a) > tol) || any(abs(I.sup - b) > tol)
    res = false;
end

I = interval(a);
if any(abs(I.inf - a) > tol) || any(abs(I.sup - a) > tol)
    res = false;
end

I = interval(a_mat);
if any(any(abs(I.inf - a_mat) > tol)) || any(any(abs(I.sup - a_mat) > tol))
    res = false;
end

I = interval(a_mat,b_mat);
if any(any(abs(I.inf - a_mat) > tol)) || any(any(abs(I.sup - b_mat) > tol))
    res = false;
end

% wrong initializations
a_large = [10; 15];
b_small = [-20; -12];
a_plus1 = [-3; -5; -2; -8];
b_plus1 = [2; 6; 3; 9];

% lower limit larger than upper limit
try
    I = interval(a,b_small); % <- should throw error here
    res = false;
end
try
    I = interval(a_large,b); % <- should throw error here
    res = false;
end

% size of limits do not match
try
    I = interval(a_plus1,b); % <- should throw error here
    res = false;
end
try
    I = interval(a,b_plus1); % <- should throw error here
    res = false;
end
try
    I = interval(a_mat,b); % <- should throw error here
    res = false;
end
try
    I = interval(a,b_mat); % <- should throw error here
    res = false;
end

% too many input arguments
try
    I = interval(a,b,b); % <- should throw error here
    res = false;
end 



if res
    disp('test_interval successful');
else
    disp('test_interval failed');
end

%------------- END OF CODE --------------