function condition = applyMappingToCondition(condition,bind)
% applyMappingToCondition - Applies the variable mapping of a Bind struct
%    to a Invariant or Guard struct
%
% Syntax:
%    condition = applyMappingToCondition(condition,bind)
%
% Inputs:
%    condition (struct) - Invariant or Guard struct in BaseComponent format
%                         (required fields: inequalities,equalities)
%    bind (struct) - Bind struct in NetworkComponent format
%                    (required fields: keys,values,renames)
%
% Outputs:
%    condition (struct) - with variable mappings applied
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% substitute variables for mapped values in Inequations
condition.inequalities = applySymMapping(condition.inequalities,...
    bind.keys,bind.values);
% substitute variables for mapped values in Equations
condition.equalities = applySymMapping(condition.equalities,...
    bind.keys,bind.values);

% ------------------------------ END OF CODE ------------------------------
