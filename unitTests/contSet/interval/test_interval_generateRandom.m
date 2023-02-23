function res = test_interval_generateRandom
% test_interval_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_interval_generateRandom
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

% Author:       Mark Wetzlinger
% Written:      27-September-2019
% Last update:  19-May-2022 (MW,name-value pair syntax)
%               23-February-2023 (MW, more cases)
% Last revision:---

%------------- BEGIN CODE --------------

% empty call
I = interval.generateRandom();

% values for tests
n = 3;
c = [3;2;1];
r = 2;

% only dimension
I = interval.generateRandom('Dimension',n);
res(1) = dim(I) == n;

% only center
I = interval.generateRandom('Center',c);
res(2) = all(withinTol(center(I),c));

% only max radius
I = interval.generateRandom('MaxRadius',r);
res(3) = r >= max(rad(I));

% dimension and center
I = interval.generateRandom('Dimension',n,'Center',c);
res(4) = dim(I) == n && all(withinTol(center(I),c));

% dimension, center, and max radius
I = interval.generateRandom('Dimension',n,'Center',c,'MaxRadius',r);
res(5) = dim(I) == n && all(withinTol(center(I),c)) && r >= max(rad(I));

% dimension and center don't match
res(6) = true;
try
    I = interval.generateRandom('Dimension',2,'Center',ones(3,1));
    res(6) = false; % <- should not get to here
end

% dimension and center don't match
res(7) = true;
try
    I = interval.generateRandom('Dimension',3,'Center',ones(3,2));
    res(7) = false; % <- should not get to here
end

% unify results
res = all(res);

%------------- END OF CODE --------------