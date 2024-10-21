function SpS_out = Inf(n)
% Inf - instantiates an n-dimensional spectrahedral shadow that is
%    equivalent to R^n
%
% Syntax:
%    SpS_out = spectraShadow.Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    SpS_out - spectrahedral shadow representing R^n
%
% Examples:
%    SpS_out = spectraShadow.Inf(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       13-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 0
    n = 0;
end
inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative'}}});

% parse input
if nargin == 0
    n = 0;
end
inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative'}}});

% the spectraShadow 1 + 0*x >= 0 is the full space
if n == 0
    SpS_out = spectraShadow([]);
    SpS_out.ESumRep = {[], []};
else
    SpS_out = spectraShadow([1 zeros([1 n])]);
    SpS_out.ESumRep = {[1 zeros([1 n])], []};
end
% assign properties
SpS_out.emptySet.val = false;
SpS_out.bounded.val = false;
SpS_out.fullDim.val = true;
SpS_out.center.val = zeros([n 1]);

% ------------------------------ END OF CODE ------------------------------
