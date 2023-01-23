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
% Last update:  21-January-2023 (MW, add sanity check for easier debugging)
% Last revision:---

%------------- BEGIN CODE --------------

% check input arguments
if length(ls.vars) ~= length(dims)
    throw(CORAerror('CORA:specialError',...
        'Number of projected dimensions do not match dimension of object.'));
end

% old symbolic variable names
oldVars = ls.vars;
% new symbolic variable names
newVars = sym('x',[N,1]);

% substitute variable names
newEq = subs(ls.eq,oldVars,newVars(dims));

% instantiate projected level set
ls = levelSet(newEq,newVars,ls.compOp);

%------------- END OF CODE ----------------
