function trans_lift = lift(trans,N,M,stateBind,inputBind,idStates)
% lift - project a transition object to a higher-dimensional space
%
% Syntax:
%    trans_lift = lift(trans,N,M,stateBind,inputBind,idStates)
%
% Inputs:
%    trans - transition object
%    N - dimension of the higher-dimensional state space
%    M - dimension of higher-dimensional input space
%    stateBind - states of the high-dimensional space that correspond to
%                the states of the low-dimensional reset object
%    inputBind - inputs of the high-dimensional space that correspond to
%                the inputs of the low-dimensional reset object
%    idStates - true/false whether additional states in reset function 
%               should be mapped via identity or zeros
%
% Outputs:
%    trans_lift - lifted transition
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, polytope/lift

% Authors:       Niklas Kochdumper, Johann Schoepfer, Mark Wetzlinger
% Written:       04-January-2022
% Last update:   ---
% Last revision: 10-October-2024 (MW, complete rewrite)

% ------------------------------ BEGIN CODE -------------------------------

% convert to polytope unless guard set is fullspace/levelSet/polytope
guard = copy(trans.guard);
if ~isa(guard,'fullspace') && ~isa(guard,'levelSet') && ~isa(guard,'polytope')
    guard = polytope(guard);
end
% lift guard set to the higher dimension
guard = lift(guard,N,stateBind);

% lift reset function to the higher dimension (in parallel hybrid automata,
% other states are not mapped by identity since they are synchronized later)
reset = lift(trans.reset,N,M,stateBind,inputBind,idStates);

% initialize the resulting transition (note: the lifted transition is a
% single transition and therefore has no knowledge about its integration
% into any automata, so the target and synchronization label remain)
trans_lift = transition(guard,reset,trans.target,trans.syncLabel);

% ------------------------------ END OF CODE ------------------------------
