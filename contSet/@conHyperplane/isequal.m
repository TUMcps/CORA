function res = isequal(hyp1,hyp2,varargin)
% isequal - checks if two constrained hyperplanes are equal
%
% Syntax:  
%    res = isequal(hyp1,hyp2)
%    res = isequal(hyp1,hyp2,tol)
%
% Inputs:
%    hyp1 - conHyperplane object
%    hyp2 - conHyperplane object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp1 = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    hyp2 = hyp1;
%    hyp3 = conHyperplane(halfspace([1;-1],0),[1 0;-1 0],[2;2]);
%
%    isequal(hyp1,hyp2)
%    isequal(hyp1,hyp3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-September-2019
% Last update:  06-June-2022
% Last revision:---

%------------- BEGIN CODE --------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default value
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{hyp1,'att','conHyperplane'};
                {hyp2,'att','conHyperplane'}; ...
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% only implemented for two constrained hyperplanes
if ~isa(hyp1,'conHyperplane') || ~isa(hyp2,'conHyperplane')
    throw(CORAerror('CORA:noops',hyp1,hyp2));
end

% numerical check
res = isequal(hyp1.h,hyp2.h,tol) && ... % halfspaces
    all(all(abs(hyp1.C - hyp2.C) < tol)) && ... % C matrices
    all(abs(hyp1.d - hyp2.d) < tol); % distances

%------------- END OF CODE --------------
