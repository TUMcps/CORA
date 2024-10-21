function res = isequal(trans1,trans2,varargin)
% isequal - checks if two transitions are equal by comparing the guard
%    sets, reset functions, target locations, and synchronization labels
%
% Syntax:
%    res = isequal(trans1,trans2)
%    res = isequal(trans1,trans2,tol)
%
% Inputs:
%    trans1 - transition object
%    trans2 - transition object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    % guard set
%    guard = polytope([0 1],0,[-1 0],0);
%
%    % reset function
%    reset1 = linearReset([1,0;0,-0.75],[1;0],[0;0]);
%    reset2 = linearReset([1,0;0,-0.75],[1;0],[1;0]);
%
%    % transition
%    trans1 = transition(guard,reset1,1);
%    trans2 = transition(guard,reset2,1);
%
%    % comparison
%    res = isequal(trans1,trans2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-November-2022
% Last update:   21-May-2023 (MW, extend to arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{trans1,'att','transition'};
                {trans2,'att','transition'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check array length
if any(size(trans1) ~= size(trans2))
    res = false; return
end

% loop over all individual transition objects
[r,c] = size(trans1);
res = false(r,c);
for i=1:r
    for j=1:c
        res(i,j) = aux_isequal(trans1(i,j),trans2(i,j),tol);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal(trans1,trans2,tol)

% easy checks first to facilitate quick exits
res = true;

% target location
if any(size(trans1.target) ~= size(trans2.target)) ...
        || ~all(trans1.target == trans2.target)
    res = false; return
end

% synchronization label
if ~strcmp(trans1.syncLabel,trans2.syncLabel)
    res = false; return
end

% reset function: can be set to [] by transition constructor, hence isempty
if xor(isempty(trans1.reset),isempty(trans2.reset)) ...
        || (~isempty(trans1.reset) && ~isequal(trans1.reset,trans2.reset,tol))
    res = false; return
end

% guard set
if any([isnumeric(trans1.guard),isnumeric(trans2.guard)])
    % empty transition object may have .guard = []
    if xor(isnumeric(trans1.guard),isnumeric(trans2.guard))
        res = false; return
    end
elseif ~isequal(trans1.guard,trans2.guard,tol)
    res = false; return
end

end

% ------------------------------ END OF CODE ------------------------------
