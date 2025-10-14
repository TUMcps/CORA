function out = splitByAmpersand(expr)
% splitByAmpersand - Flatten all '&' terms and return a single merged string.
%
% Syntax:
%           out = splitByAmpersand(expr)
%
% Inputs:
%    expr - string containing an algebraic SpaceEx expression (e.g. for invariants, flows, etc)
%
% Outputs:
%    out  - string containing an equivalent expression without redundant
%           parentheses
%
% Example:
%   res = splitByAmpersand('a = 5*(x + y) & (x = 3 & (y = 3))')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       09-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Accept string or char input; convert string to char to allow indexing
if isstring(expr)
    expr = char(expr);
end

expr = strtrim(expr);

% Repeatedly strip redundant outer parentheses
expr = aux_stripOuterParens(expr);

depth = 0;
splitIdx = [];

for i = 1:length(expr)
    ch = expr(i);
    if ch == '('
        depth = depth + 1;
    elseif ch == ')'
        depth = depth - 1;
    elseif ch == '&' && depth == 0
        splitIdx(end+1) = i; %#ok<AGROW>
    end
end

if isempty(splitIdx)
    % No top-level '&' -> this is a single term
    terms = {expr};
else
    % Split at all top-level '&'
    parts = {};
    prev = 1;
    for k = 1:length(splitIdx)
        idx = splitIdx(k);
        parts{end+1} = strtrim(expr(prev:idx-1)); %#ok<AGROW>
        prev = idx+1;
    end
    parts{end+1} = strtrim(expr(prev:end));

    % Recurse on each part and collect tokens
    terms = {};
    for k = 1:length(parts)
        sub = splitByAmpersand(parts{k});     % returns char
        tokens = strtrim(strsplit(sub, '&')); % cell array of char
        terms = [terms, tokens]; %#ok<AGROW>
    end
end

% Remove any empty tokens (in case of stray & or whitespace)
terms = terms(~cellfun(@isempty, terms));

% Final output: single merged string (char)
out = string(strjoin(terms, ' & '));
end


% Auxiliary functions -----------------------------------------------------


function s = aux_stripOuterParens(s)
% Remove *all* redundant outer parentheses from a char vector
if isstring(s)
    s = char(s);
end

changed = true;
while changed && numel(s) >= 2 && s(1) == '(' && s(end) == ')'
    depth = 0;
    balanced = true;
    for i = 1:length(s)
        if s(i) == '('
            depth = depth + 1;
        elseif s(i) == ')'
            depth = depth - 1;
            if depth == 0 && i < length(s)
                balanced = false;
                break;
            end
        end
    end
    if balanced
        s = strtrim(s(2:end-1)); % safe because s is a char vector
    else
        changed = false;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
