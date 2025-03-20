function res = checkValueAttributes(value,class,attributes)
% checkValueAttributes - checks if the given value is off the correct class
%    and has the correct attributes
%
% Syntax:
%    res = checkValueAttributes(value,class,attributes)
%
% Inputs:
%    value - value to be tested
%    class - char, class of the value
%    attributes - cell array of attributes (should evaluate to logical), 
%           either
%           - function_handle - takes value as input and returns logical
%           - char, such that it can be evaluated
%                       via feval or custom auxiliary function
%
% Outputs:
%    res - logical
%
% Example: 
%    value = 1;
%    res = checkValueAttributes(value, 'numeric', {'integer','nonnan', @(value) value >= 1})
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: inputArgsCheck, readNameValuePair, feval

% Authors:       Tobias Ladner
% Written:       03-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
% internal function, expecting correct inputs
nchoosek(3,3); 

% init
resvec = false(1,numel(attributes)+1);

% check class
resvec(1) = isempty(class) || isa(value,class);

% check attributes
for i = 1:numel(attributes)
    resvec(1+i) = aux_checkAttribute(value,class,attributes{i},'all');
    if ~resvec(i)
        break
    end
end

% gather results
res = all(resvec);

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkAttribute(value,class,attribute,reduction)
    % check single attribute
    if ischar(attribute) || isstring(attribute)
        % check negation
        if startsWith(attribute,'non')
            res = ~aux_checkAttribute(value,class,extractAfter(attribute,3),'any');
            return
        end

        % evaluate
        % (custom attributes extracted from validateattributes)
        switch attribute
            % n-d arrays
            case {'2d','is2d'}
                res = aux_isnd(value,2);
            case {'3d','is3d'}
                res = aux_isnd(value,3);
            case {'square','issquare'}
                res = aux_issquare(value);

            % type
            case {'binary','isbinary'}
                res = aux_isbinary(value);
            case {'integer','isinteger'}
                res = aux_isinteger(value);

            % property
            case {'even','iseven'}
                res = aux_iseven(value);
            case {'odd','isodd'}
                res = aux_isodd(value);
            case {'negative','isnegative'}
                res = aux_isnegative(value);
            case {'positive','ispositive'}
                res = aux_ispositive(value);
            case {'zero','iszero'}
                res = aux_iszero(value);
                
            otherwise % use built-in function (starts with 'is')
                if ~startsWith(attribute,'is')
                    % prepend 'is'
                    attribute = sprintf('is%s',attribute);
                end

                % evaluate using built-in function
                res = feval(attribute,value);
        end
                
    elseif isa(attribute,'function_handle')
        % evaluate function handle
        res = attribute(value);
    else
        % unable to check attribute; unknown type
        throw(CORAerror('CORA:wrongValue','third','Unable to check attribute %s for %s', attribute, value))
    end

    % apply reduction method
    switch reduction
        case 'all'
            res = all(res,'all');
        case 'any'
            res = any(res,'all');
        case 'none'
            % res = res;
        otherwise
            throw(CORAerror('CORA:wrongValue','fourth',{'all','any','none'}))
    end
end

function res = aux_isnd(value,n)
    % check if n-d array
    res = numel(size(value)) == n;
end

function res = aux_issquare(value)
    % check if square matrix
    res = size(value,1) == size(value,2);
end

function res = aux_isbinary(value)
    % check if all values are 0 or 1
    res = value == 1 | value == 0;
end

function res = aux_isinteger(value)
    % check if all values are integers
    % see also: isinteger (only checks class)
    res = mod(value,1) == 0;
end

function res = aux_iseven(value)
    % check if all values are even
    res = mod(value,2) == 0;
end

function res = aux_isodd(value)
    % check if all values are odd
    res = mod(value,2) == 1;
end

function res = aux_ispositive(value)
    % check if all values are positive
    res = value > 0;
end

function res = aux_isnegative(value)
    % check if all values are negative
    res = value < 0;
end

function res = aux_iszero(value)
    % check if all values are zero
    res = value == 0;
end

% ------------------------------ END OF CODE ------------------------------
