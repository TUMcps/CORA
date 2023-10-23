function sym_out = applySymMapping(sym_in,keys,values)
% applySymMapping - substitute variables with symbolic values
%
% Syntax:
%    sym_out = applySymMapping(sym_in,keys,values)
%
% Inputs:
%    sym_in (symbolic) - symbolic expressions to be manipulated
%    keys (string) - names of variables to be substituted
%    values (symbolic) - values to be substituted for above variables
%                        keys & values must be vectors of same length
%                        keys(i) is substituted by values(i)
%
% Outputs:
%    sym_out (symbolic) - sym_in after substitutions
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

% convert parameters to cell array if necessary, (required for subs)
if iscell(keys)
    keys_c = keys;
else
    keys_c = cell(size(keys));
    for i = 1:numel(keys)
        keys_c{i} = keys(i);
    end
end

if iscell(values)
    values_c = values;
else
    values_c = cell(size(values));
    for i = 1:numel(values)
        values_c{i} = values(i);
    end
end

% use symbolic kit substitution function
sym_out = subs(sym_in,keys_c,values_c);

% ------------------------------ END OF CODE ------------------------------
