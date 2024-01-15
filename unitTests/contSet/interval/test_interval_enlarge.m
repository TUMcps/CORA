function res = test_interval_enlarge
% test_interval_enlarge - unit test function of enlarge
%
% Syntax:
%    res = test_interval_enlarge
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

% Authors:       Mark Wetzlinger
% Written:       29-August-2019
% Last update:   03-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% bounded
I = interval([-2; -4; -3],[2; 3; 1]);
% ...bounded scaling factor
factor = 2;
I_enlarge = enlarge(I, factor);
I_true = interval([-4; -7.5; -5],[4; 6.5; 3]);
res(end+1,1) = isequal(I_enlarge,I_true);

% ...Inf as scaling factor
factor = Inf;
I_enlarge = enlarge(I, factor);
I_true = interval(-Inf(3,1),Inf(3,1));
res(end+1,1) = isequal(I_enlarge,I_true);

% ...0 as scaling factor
factor = 0;
I_enlarge = enlarge(I, factor);
I_true = interval(center(I));
res(end+1,1) = isequal(I_enlarge,I_true);


% unbounded
I = interval([-Inf;-2],[2;Inf]);
% ...bounded scaling factor
factor = 2;
I_enlarge = enlarge(I, factor);
I_true = interval(-Inf(2,1),Inf(2,1));
res(end+1,1) = isequal(I_enlarge,I_true);

% ...Inf as scaling factor
factor = Inf;
I_enlarge = enlarge(I, factor);
I_true = interval(-Inf(2,1),Inf(2,1));
res(end+1,1) = isequal(I_enlarge,I_true);

% ...0 as scaling factor (throws an error)
factor = 0;
try
    I_enlarge = enlarge(I, factor);
    res(end+1,1) = false;
end


% add results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
