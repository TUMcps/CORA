function inputArgsCheck(inputArgs)
% inputArgsCheck - checks input arguments of CORA functions for compliance
%    with semantic requirements; if any violation is detected, an error is
%    thrown, otherwise nothing happens
%
% Syntax:
%    inputArgsCheck(inputArgs)
%
% Inputs:
%    inputArgs - cell-array with structure
%                  { {input1,id,{classes},{attributes}};
%                    {input2,id,{possibilities}};
%                     ...; };
%                to describe how input arguments should be like:
%                   input1, input2, ... - value of input argument
%                   id - identifier for attribute check ('att') or string
%                        comparison ('str')
%                   classes - (only id = 'att') admissible classes as a
%                             cell-array of chars or a single char-array
%                   attributes - (only id = 'att') cell array of admissible attributes
%                                (check function checkValueAttributes for details)
%                   possibilities - (only id = 'str') admissible strings as
%                                   a cell-array of chars or single
%                                   char-array
%
% Outputs:
%    ---
%
% Example:
%    obj = capsule(zeros(2,1),ones(2,1));
%    N = 5;
%    type = 'standard';
%    inputArgs = { {obj, 'att',{'cell','capsule','other'},{'nonempty'}};
%                  {N,   'att',{'numeric'},               {'positive'}};
%                  {type,'str',{'standard','gaussian'}                } };
%    inputArgsCheck(inputArgs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkValueAttributes, readNameValuePair

% Authors:       Mingrui Wang, Mark Wetzlinger, Tobias Ladner
% Written:       30-May-2022
% Last update:   23-January-2024 (MW, exact match for strings)
%                03-March-2025 (TL, reworked using checkValueAttributes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if disabled by macro
if ~CHECKS_ENABLED
    return;
end

% turn off warnings (if no attributes are given)
warOrig = warning;
warning('off', 'all');

% number of input arguments to original function that have to be checked
nrInputArgs = length(inputArgs);

% loop over all input arguments
for i = 1:nrInputArgs
    % read information about i-th input argument
    inputArg = inputArgs{i};
    % read value of i-th input argument
    value = inputArg{1};
    % read identifier ('att' or 'str')
    identifier = inputArg{2};

    % case distinction
    switch identifier
        case 'att' % check classname (and attributes) in this case
            aux_checkAtt(i, inputArg, value);

        case 'str' % check the strings in this case
            aux_checkStr(i, inputArg, value)

        otherwise
            throw(CORAerror("CORA:wrongValue", 'second', '''att'' or ''str''.'))

    end

end

% turn warnings back on
warning(warOrig);

end


% Auxiliary functions -----------------------------------------------------

function aux_checkAtt(i, inputArg, value)
% check attribute
% inputArg: {value,'att',classes,attributes}

% read out classes
classes = inputArg{3};
if ~iscell(classes) % ensure cell array
    classes = {classes};
end
% ! should be {'class1', 'class2', ...} now

% read out attributes
attributes = {};
if numel(inputArg) >= 4
    attributes = inputArg{4};

    % ensure cell array
    if ~iscell(attributes)
        attributes = {attributes};
    end
end
% check empty
if isempty(attributes) 
    % ensure single (empty) entry
    attributes = {attributes};
end
% ensure cell of cell array
if ~iscell(attributes{1}) 
    attributes = {attributes};
end
% ! should be {{'attribute1', 'attribute2', ...}, ...} now

% ensure number of classes match number of attributes
if numel(classes) ~= numel(attributes)
    if isscalar(classes)
        classes = cellfun(@(~) classes{1}, attributes,'UniformOutput',false);
    elseif isscalar(attributes)
        attributes = cellfun(@(~) attributes{1}, classes,'UniformOutput',false);
    else
        throw(CORAerror('CORA:specialError',sprintf('Mismatch between given number of classes and antributes for %s argument.', aux_countingNumber(i))))
    end
end

% check class and attributes
resvec = false(1,numel(classes));
for j=1:numel(classes)
    resvec(j) = checkValueAttributes(value,classes{j},attributes{j});
end

% gather results
res = any(resvec);
if ~res
    % find best guess for given class
    classresvec = cellfun(@(class) isa(value,class), classes, 'UniformOutput',true);
    idx = find(classresvec);

    if isscalar(idx)
        % show specific error message for guess
        text = aux_errMsgClassAttributes(classes(idx),attributes{idx});
    else
        % show error message for all classes
        text = aux_errMsgClassAttributes(classes,attributes);
    end

    % throw error
    throw(CORAerror('CORA:wrongValue',aux_countingNumber(i), text))     
end

end

function aux_checkStr(i, inputArg, value)
% read string
if iscell(inputArg{3})
    validateStr = inputArg{3};
else
    validateStr = inputArg(3);
end

% check exact match with admissible values
if ~ismember(value, validateStr)
    % generate string of admissible values (user info)
    validrange = ['''', strjoin(validateStr, "', '"), ''''];

    % throw error
    throw(CORAerror('CORA:wrongValue', aux_countingNumber(i), ...
        validrange));
end
end

function text = aux_errMsgClassAttributes(class, varargin)

% optional argument: attributes
if ~isempty(varargin)
    attributes = varargin{1};

    if ~iscell(attributes)
        n_attribute = 1;
        varyingAttributes = false;
    else
        entryCells = cellfun(@(x) iscell(x), attributes, 'UniformOutput', true);
        if ~any(entryCells)
            varyingAttributes = false;
            % number of attributes (can be zero)
            n_attribute = length(attributes);
            % convert all numerics from attributes to string (not ideal
            % output text, but works for now...)
            for i = 1:n_attribute
                if isnumeric(attributes{i})
                    attributes{i} = num2str(attributes{i});
                end
            end
            if n_attribute == 1 && strcmp(attributes{1}, '')
                % deprecated case with '' -> counts as zero attributes
                n_attribute = 0;
            end
        else
            varyingAttributes = true;
            n_attribute = zeros(length(attributes), 1);
            for c = 1:length(attributes)
                % number of attributes (can be zero)
                n_attribute(c) = length(attributes{c});
            end
        end
    end
else
    attributes = [];
    n_attribute = 0;
    varyingAttributes = false;
end

% number of classes
n_class = length(class);

% text
text = "";

% loop over each class, add attributes
for c = 1:n_class
    % append text
    text = text + string(class{c});

    if varyingAttributes
        if n_attribute(c) > 0
            text = text + " (" + ...
                strjoin(string(attributes{c}), ', ') + ")";
        end
    else
        if n_attribute > 0
            text = text + " (" + ...
                strjoin(string(attributes), ', ') + ")";
        end
    end

    if c < n_class
        % 17 white spaces for correct alignment
        text = text + ",\n" ...
            +"                   ";
    end
end

end

function text = aux_countingNumber(i)
% returns counting number in text form
switch i
    % for early cases we write the enumeration out in full
    case 1
        text = 'first';
    case 2
        text = 'second';
    case 3
        text = 'third';
    case 4
        text = 'fourth';
    case 5
        text = 'fifth';
    case 6
        text = 'sixth';
    case 7
        text = 'seventh';
    case 8
        text = 'eighth';
    case 9
        text = 'ninth';
    otherwise
        % from here on we only write "No. "
        text = ['No. ', num2str(i)];
end

end

% ------------------------------ END OF CODE ------------------------------
