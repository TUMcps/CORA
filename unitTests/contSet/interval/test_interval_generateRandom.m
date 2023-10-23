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

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       27-September-2019
% Last update:   19-May-2022 (MW, name-value pair syntax)
%                23-February-2023 (MW, more cases)
%                22-May-2023 (TL, tests for interval matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty call
I = interval.generateRandom();

% values for tests
n = 3;
c = [3;2;1];
r = 2;
r_nd = [2;3;1];

% only dimension
I = interval.generateRandom('Dimension',n);
resvec(end+1) = dim(I) == n;

% only center
I = interval.generateRandom('Center',c);
resvec(end+1) = all(withinTol(center(I),c));

% only max radius
I = interval.generateRandom('MaxRadius',r);
resvec(end+1) = r >= max(rad(I));

% different max radius per dimension
I = interval.generateRandom('MaxRadius',r_nd);
resvec(end+1) = all(r_nd >= rad(I));

% dimension and center
I = interval.generateRandom('Dimension',n,'Center',c);
resvec(end+1) = dim(I) == n && all(withinTol(center(I),c));

% dimension, center, and max radius
I = interval.generateRandom('Dimension',n,'Center',c,'MaxRadius',r);
resvec(end+1) = dim(I) == n && all(withinTol(center(I),c)) && r >= max(rad(I));

% dimension and center don't match
if CHECKS_ENABLED

resvec(end+1) = true;
try
    I = interval.generateRandom('Dimension',2,'Center',ones(3,1));
    resvec(end+1) = false; % <- should not get to here
end

% dimension and center don't match
resvec(end+1) = true;
try
    I = interval.generateRandom('Dimension',3,'Center',ones(3,2));
    resvec(end+1) = false; % <- should not get to here
end

end

% test multi-dimensional --------------------------------------------------

n = [2,3];
c = [2 5 4; 1 -1 2];
r = 4;
r_nd = [1 0.5 2; 2 5 2];

% only dimension
I = interval.generateRandom('Dimension',n);
resvec(end+1) = all(dim(I) == n);

% only center
I = interval.generateRandom('Center',c);
resvec(end+1) = all(withinTol(center(I),c),'all');

% only max radius
I = interval.generateRandom('MaxRadius',r);
resvec(end+1) = r >= max(rad(I),[],'all');

% different max radius per dimension
I = interval.generateRandom('MaxRadius',r_nd);
resvec(end+1) = all(r_nd >= rad(I),'all');

% dimension and center
I = interval.generateRandom('Dimension',n,'Center',c);
resvec(end+1) = all(dim(I) == n) && all(withinTol(center(I),c),'all');

% dimension, center, and max radius
I = interval.generateRandom('Dimension',n,'Center',c,'MaxRadius',r);
resvec(end+1) = all(dim(I) == n) && all(withinTol(center(I),c),'all') && r >= max(rad(I),[],'all');


% unify results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
