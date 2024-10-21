function SpS_out = empty(n)
% empty - instantiates an empty n-dimensional spectraShadow
%
% Syntax:
%    SpS_out = spectraShadow.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    SpS_out - empty n-dimensional spectrahedral shadow
%
% Examples:
%    SpS_out = spectraShadow.empty(2);
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

% the spectraShadow -1 + 0*x >= 0 is empty
if n == 0
    SpS_out = spectraShadow([-1 0], zeros([0 1]), zeros([0 1]));
    SpS_out.ESumRep.val = {[-1 zeros([1 n])], []};
else
    SpS_out = spectraShadow([-1 zeros([1 n])]);
    SpS_out.ESumRep.val = {[-1 zeros([1 n])], []};
end
% assign properties
SpS_out.emptySet.val = true;
SpS_out.bounded.val = true;
SpS_out.fullDim.val = false;
SpS_out.center.val = [];


% ------------------------------ END OF CODE ------------------------------
