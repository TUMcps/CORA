function ls = projectHighDim(ls,N,dims)
% projectHighDim - project a levelSet object to a higher-dimensional space
%
% Syntax:  
%    ls = projectHighDim(ls,N,dims)
%
% Inputs:
%    ls - levelSet object
%    N - dimension of the higher dimensional space
%    dims - states of the high dimensional space that correspond to the
%          states of the low dimensional level set
%
% Outputs:
%    ls - resulting levelSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, mptPolytope/projectHighDim

% Author:       Maximilian Perschl
% Written:      25-February-2022
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

oldVars = ls.vars;
newVars = sym('x',[N,1]);

newEq = ls.eq;

newEq = subs(newEq,oldVars,newVars(dims));

ls = levelSet(newEq,newVars,ls.compOp);

%------------- END CODE ----------------