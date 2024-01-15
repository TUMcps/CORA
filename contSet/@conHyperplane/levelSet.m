function ls = levelSet(hyp)
% levelSet - converts a constrained hyperplane to a level set; currently
%    only supported for non-constrained hyperplanes
%
% Syntax:
%    ls = levelSet(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    ls - levelSet object
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0));
%    ls = levelSet(hyp)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only supported for non-constrained hyperplane
if ~isempty(hyp.C)
    throw(CORAerror('CORA:notSupported',"The conversion from conHyperplane " ...
        + "to levelSet is only supported for non-constrained hyperplanes."));
end

% use equality representation
vars = sym('x',[dim(hyp),1]);

% read out halfspace
A = hyp.a';
b = hyp.b;

% define equality equation
eqs = A'*vars - b;
if size(A,2) == 1
    compOps = '==';
else
    compOps = repmat({'=='},size(A,2),1);
end

% instantiate level set object
ls = levelSet(eqs,vars,compOps);

% ------------------------------ END OF CODE ------------------------------
