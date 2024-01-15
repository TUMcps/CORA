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
%    check - (optional) function name or cell-array of function names:
%               check functions for value, e.g., 'isscalar'
%    def - (optional) default value
%
% Outputs:
%    NVpairs - list of name-value pairs
%    value - value of name-value pair
%
% Example: 
%    NVpairs = {'Splits',8};
%    [NVpairs,value] = readNameValuePair(NVpairs,'Splits','isscalar',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: readPlotOptions, polyZonotope/plot

% Authors:       Mark Wetzlinger, Niklas Kochdumper
% Written:       15-July-2020
% Last update:   24-November-2021 (allow cell-array of check functions)
%                07-July-2022 (MW, case-insensitive, string compatibility)
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
for i=1:2:length(NVpairs)-1
    
    % has to be a char/string (safety check) and match NVpairs
    % note: case-insensitive!
    if (ischar(NVpairs{i}) || isstring(NVpairs{i})) && strcmpi(NVpairs{i},name)
        % name found, get value
        value = NVpairs{i+1};

        if checks_enabled
            % check whether name complies with check
            for j=1:length(funcs)
                if ~feval(funcs{j},value)
                    throw(CORAerror('CORA:specialError', ...
                        sprintf("Invalid assignment for name-value pair '%s': Must pass '%s'.", name, funcs{j})))
                end
            end
        end

        % empty the corresponding cells
        NVpairs(i:i+1) = [];
        break
    end
end

end

% ------------------------------ END OF CODE ------------------------------
