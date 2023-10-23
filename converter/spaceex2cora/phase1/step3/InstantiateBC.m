function child = InstantiateBC(parent,template,bind)
% InstantiateBC - Apply a bind to a Base-Component template, this
%    instantiates the component using the variable names from the bind
%    in the invariant/flow/guard/reset equations; additionally, all
%    constant parameters are removed from the list of variables
%
% Syntax:
%    child = InstantiateBC(parent,template,bind)
%
% Inputs:
%    parent (struct) - instantiated parent component
%    template (struct) - template of BC to instantiate
%    bind (struct) - bind (i.e., instantion in parent network component)
%                    to be used for instantiation of the base component
%
% Outputs:
%    child (struct) - instantiated base component (child of "parent"
%                     network component)
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   13-January-2023 (MW, save constants for completeness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% start with a copy of the template
child = template;

% save global name for new instance
child.name = parent.name + "." + bind.localName;

% output on command window
fprintf("instantiating Base Component: '%s'\n",child.name);

% Apply variable mapping to all equations

% number of locations in base components
nrLocs = length(child.States);

% loop over each location
for iLoc = 1:nrLocs
    % read out location
    loc = child.States(iLoc);
    
    % Apply mapping to Flow
    child.States(iLoc).Flow = applyMappingToEquation(loc.Flow,bind);
        
    % Apply Mapping to Invariant
    child.States(iLoc).Invariant = applyMappingToCondition(loc.Invariant,bind);
    
    % number of outgoing transitions
    numTrans = length(loc.Trans);

    % loop over each outgoing transition
    for iTrans = 1:numTrans
        % read out transition
        trans = loc.Trans(iTrans);
        
        % Apply mapping to reset function
        child.States(iLoc).Trans(iTrans).reset = applyMappingToEquation(trans.reset,bind);
        
        % Apply mapping to guard set
        child.States(iLoc).Trans(iTrans).guard = applyMappingToCondition(trans.guard,bind);
    end
end

% constants are mapped in resolveBinds

% ------------------------------ END OF CODE ------------------------------
