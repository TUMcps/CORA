function [NVpairs,value] = readNameValuePair(NVpairs,name,varargin)
% readNameValuePair - searches through a list of name-value pairs to return
%    the corresponding value of a given name and a redacted list (read is
%    (case-insensitive and works for chars and strings); additionally, the
%    provided value can be checked by a function handle and one may also
%    return a default value if the desired name-value pair is not given
%
% Syntax:
%    [NVpairs,value] = readNameValuePair(NVpairs,name)
%    [NVpairs,value] = readNameValuePair(NVpairs,name,check)
%    [NVpairs,value] = readNameValuePair(NVpairs,name,check,def)
%
% Inputs:
%    NVpairs - cell-array: list of name-value pairs
%    name - char: name of name-value pair
%    check - (optional) cell array of admissible attributes
%         (check function checkValueAttributes for details)
%    def - (optional) default value
%
% Outputs:
%    NVpairs - list of name-value pairs
%    value - value of name-value pair
%
% Example: 
%    NVpairs = {'Splits',8};
%    [NVpairs,value] = readNameValuePair(NVpairs,'Splits','scalar',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkValueAttributes, inputArgsCheck

% Authors:       Mark Wetzlinger, Niklas Kochdumper, Tobias Ladner
% Written:       15-July-2020
% Last update:   24-November-2021 (allow cell-array of check functions)
%                07-July-2022 (MW, case-insensitive, string compatibility)
%                03-March-2025 (TL, reworked using checkValueAttributes)
%                29-November-2023 (TL, check for CHECKS_ENABLED)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% write function into cell-array
funcs = [];
if nargin >= 3
    if ~iscell(varargin{1})
        funcs = varargin(1);
    else
        funcs = varargin{1};
    end
end

% default empty value (name not found)
value = [];
if nargin >= 4
    value = varargin{2};
end

% checks enabled
checks_enabled = CHECKS_ENABLED();

% check every second entry 
% (back to front to select last given value in case there are duplicates,
% which is also how Matlab handles it, e.g., in plot)
for i=length(NVpairs)-1:-2:1
    
    % has to be a char/string (safety check) and match NVpairs
    % note: case-insensitive!
    if (ischar(NVpairs{i}) || isstring(NVpairs{i})) && strcmpi(NVpairs{i},name)
        % name found, get value
        value = NVpairs{i+1};

        if checks_enabled && ~isempty(funcs)
            % check whether name complies with check
            res = checkValueAttributes(value,[],funcs);
            if ~res
                throw(CORAerror('CORA:specialError', ...
                    sprintf("Invalid assignment for name-value pair '%s': Must pass %s.", ...
                        name, aux_formatFuncs(funcs)) ...
                ))
            end
        end

        % empty the corresponding cells
        NVpairs(i:i+1) = [];
        break
    end
end

end


% Auxiliary functions -----------------------------------------------------

function str = aux_formatFuncs(funcs)
    % format given functions
    
    % prepare for formatting
    if iscell(funcs)
        % format each
        for i=1:numel(funcs)
            funcs{i} = aux_formatFunc(funcs{i});
        end
    
        % join and format
        str = sprintf('{%s}',strjoin(funcs,', '));
    else
        % format single
        str = aux_formatFunc(funcs);
    end
end

function str = aux_formatFunc(func)
    if isa(func,'function_handle')
        str = func2str(func);
    else
        str = sprintf('''%s''',func);
    end
end


% ------------------------------ END OF CODE ------------------------------
