function equation = applyMappingToEquation(equation,bind)
% applyMappingToEquation - Applies the variable mapping of a Bind struct to
%    a Flow or Reset struct
%
% Syntax:
%    equation = applyMappingToEquation(equation,bind)
%
% Inputs:
%    equation (struct) - Flow or Reset struct in BaseComponent format
%                        (required fields: varnames,exprs)
%    bind (struct) - Bind struct in NetworkComponent format
%                    (required fields: keys,values,renames)
%
% Outputs:
%    equation (struct) - with variable mappings applied
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

% rename left-side variables
equation.varNames = applyRenames(equation.varNames,bind.keys,bind.renames);
% substitute variables for mapped values in right-side expressions
equation.expressions = applySymMapping(equation.expressions,...
    bind.keys,bind.values);

% ------------------------------ END OF CODE ------------------------------
