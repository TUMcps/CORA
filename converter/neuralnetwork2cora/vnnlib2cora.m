function [X0, spec] = vnnlib2cora(file)
% vnnlib2cora - import specifcations from .vnnlib files
%
% Syntax:
%    [X0,spec] = vnnlib2cora(file)
%
% Inputs:
%    file - path to a file .vnnlib file storing the specification
%
% Outputs:
%    X0 - initial set represented as an object of class interval
%    spec - specifications represented as an object of class specification
%
% Reference:
%    - https://www.vnnlib.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       23-November-2021
% Last update:   26-July-2023 (TL, speed up)
%                30-August-2023 (TL, bug fix multiple terms in and)
%                14-June-2024 (TL, major speed up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read in text from file
text = fileread(file);

% determine number of inputs and number of outputs
nrInputs = 0;
nrOutputs = 0;
lineBreaks = strfind(text, newline);
for ln = 1:(numel(lineBreaks)-1)
    % iterate through file
    i = lineBreaks(ln)+1;
    i1 = lineBreaks(ln+1);
    if startsWith(text(i:i1), '(declare-const ')
        temp = text(i+15:i1);
        ind = find(temp == ' ');
        temp = temp(1:ind(1)-1);
        if strcmp(temp(1), 'X')
            % found new input; read out index
            nrInputs = max(nrInputs, str2double(temp(3:end)));
        elseif strcmp(temp(1), 'Y')
            % found new output; read out index
            nrOutputs = max(nrOutputs, str2double(temp(3:end)));
        end
    end
end

% +1 due to 0-indexing in vnnlib files
data.nrInputs = nrInputs + 1;
data.nrOutputs = nrOutputs + 1;
data.currIn = 0;

% parse file
data.polyInput = [];
data.polyOutput = [];
while ~isempty(text)
    if startsWith(text, '(assert')
        text = strtrim(text(8:end));
        [len, data] = aux_parseAssert(text, data);
        text = strtrim(text(len+1:end));
    else
        ln = regexp(text,newline,'once');
        text = text((ln+1):end);
    end
end

% convert data to polytopes ---

% a) convert input
% potentially convert input polytopes to intervals
X0 = cell(1, length(data.polyInput));
for i = 1:length(X0)
    polyStruct = data.polyInput(i);
    P = polytope(polyStruct.C, polyStruct.d);
    [res, I] = representsa_(P, 'interval', eps);

    if res
        X0{i} = I;
    else
        throw(CORAerror('CORA:specialError', 'Input set is not an interval.'))
    end
end

% b) convert output
Y = cell(1, length(data.polyOutput));
for i = 1:length(data.polyOutput)
    Y{i} = polytope(data.polyOutput(i).C, data.polyOutput(i).d);
end

% construct specification from list of output polytopes
if isempty(Y)
    throw(CORAerror("CORA:converterIssue",sprintf('Unable to convert file: %s', file)));
elseif isscalar(Y)
    % We have a single unsafe set.
    spec = specification(Y{1}, 'unsafeSet');
else
    % We have a union of unsafe sets, i.e., any point contained in any set 
    % is a counter example.
    % We want to reduce the number of specification sets. Therefore, we 
    % try to convert to the union of unsafe sets to a union of safe sets.
    Y_safe = safeSet2unsafeSet(Y); % We can use the function; safe vs. unsafe is juts
    % Check if we could reduce the number of specifications.
    if length(Y_safe) < length(Y)
        % We could reduce the number of specifications; we add the sets 
        % as unsafe sets.
        Y = Y_safe;
        type = 'safeSet'; % Type of union, e.g. Y is a union of safe sets
    else
        % We could not reduce the number of specifications; we add the 
        % sets as unsafe sets.
        type = 'unsafeSet'; % Type of union, e.g. Y is a union of unsafe sets
    end
    % Initialize the result
    spec = [];
    % Add all the specifications together. We can only add inverse 
    % sets, because then all specification have to hold simultaneously.
    for i = 1:length(Y)
        switch type
            case 'safeSet'
                % We have to convert the safe set to an unsafe set.
                Yi_safe = specification(Y{i},'safeSet');
                % All inverted sets unsafe sets have to be avoided.
                spec = add(spec,Yi_safe);
            case 'unsafeSet'
                % We have to convert the unsafe set to an safe set.
                Yi_unsafe = specification(Y{i},'unsafeSet');
                % We add the unsafe set.
                spec = add(spec,Yi_unsafe);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [len, data] = aux_parseAssert(text, data)
% parse one assert statement

if startsWith(text, '(<=') || startsWith(text, '(>=')

    [len, data] = aux_parseLinearConstraint(text, data);

elseif startsWith(text, '(or')

    text = strtrim(text(4:end));
    data_.spec = [];
    data_.nrOutputs = data.nrOutputs;
    len = 5;

    % parse all or conditions
    while ~startsWith(text, ')')

        % parse one or condition
        data_.nrInputs = data.nrInputs;
        data_.nrOutputs = data.nrOutputs;
        data_.polyInput = [];
        data_.polyOutput = [];
        data_.currIn = 0;

        [len_, data_] = aux_parseAssert(text, data_);

        % update remaining text
        text = strtrim(text(len_:end));
        len = len + len_;

        % update input conditions
        if ~isempty(data_.polyInput)
            if ~isempty(data.polyInput)
                data.polyInput(end+1) = data_.polyInput(1);
            else
                data.polyInput = data_.polyInput;
            end
        end

        % update output conditions
        if ~isempty(data_.polyOutput)
            if ~isempty(data.polyOutput)
                data.polyOutput(end+1) = data_.polyOutput(1);
            else
                data.polyOutput = data_.polyOutput;
            end
        end
    end

elseif startsWith(text, '(and')

    text = strtrim(text(5:end));
    len = 6;

    % parse all or conditions
    while ~startsWith(text, ')')
        [len_, data] = aux_parseAssert(text, data);
        text = text(len_:end);

        % trim white spaces
        text_ = strtrim(text);
        len_ = len_ + (length(text)-length(text_));
        text = text_;

        % move overall length counter to current position
        len = len + len_;
        if startsWith(text, '(')
            % multiple terms in and; correct len
            len = len - 1;
        end
    end

else
    throw(CORAerror('CORA:converterIssue', sprintf('Failed to parse vnnlib file. Parsed up to line %i.', len)))
end
end

function S = aux_createPolytopeStruct(n)
S = struct;
S.C = zeros(2*n,n);
S.d = zeros(2*n,1);
end

function [len, data] = aux_parseLinearConstraint(text, data)
% parse a linear constraint

% extract operator
op = text(2:3);
text = text(5:end);
len = 5;

% get type of constraint (on inputs X or on output Y)
type = aux_getTypeOfConstraint(text);

% initialization
if strcmp(type, 'input')
    C = zeros(1, data.nrInputs);
    d = 0;
else
    C = zeros(1, data.nrOutputs);
    d = 0;
end

% parse first argument
[C1, d1, len_] = aux_parseArgument(text, C, d);
len = len + len_;
text = strtrim(text(len_:end));

% parse second argument
[C2, d2, len_] = aux_parseArgument(text, C, d);
len = len + len_;

% combine the two arguments
if strcmp(op, '<=')
    C = C1 - C2;
    d = d2 - d1;
else
    C = C2 - C1;
    d = d1 - d2;
end

% combine the current constrain with previous constraints
if strcmp(type, 'input')
    if isempty(data.polyInput)
        data.polyInput = aux_createPolytopeStruct(data.nrInputs);
    end

    data.currIn = data.currIn+1;

    for i = 1:length(data.polyInput)
        data.polyInput(i).C(data.currIn,:) = C;
        data.polyInput(i).d(data.currIn) = d;
    end

else % output
    if isempty(data.polyOutput)
        data.polyOutput = aux_createPolytopeStruct(0);
    end

    for i = 1:length(data.polyOutput)
        data.polyOutput(i).C = [data.polyOutput(i).C; C];
        data.polyOutput(i).d = [data.polyOutput(i).d; d];
    end
end
end

function [C, d, len] = aux_parseArgument(text, C, d)
% parse next argument

if startsWith(text, 'X') || startsWith(text, 'Y')

    len = [];
    for i = 1:length(text)
        if strcmp(text(i), ' ') || strcmp(text(i), ')')
            len = i;
            break;
        end
    end
    index = str2double(text(3:len-1)) + 1;
    C(index) = C(index) + 1;

elseif startsWith(text, '(+')

    % parse first argument
    [C1, d1, len] = aux_parseArgument(text, C, d);
    text = strtrim(text(len:end));

    % parse second argument
    [C2, d2, len_] = aux_parseArgument(text, C, d);
    len = len + len_;

    % combine both arguments
    C = C1 + C2;
    d = d1 + d2;

elseif startsWith(text, '(-')

    % parse first argument
    [C1, d1, len] = aux_parseArgument(text, C, d);
    text = strtrim(text(len:end));

    % parse second argument
    [C2, d2, len_] = aux_parseArgument(text, C, d);
    len = len + len_;

    % combine both arguments
    C = C1 - C2;
    d = d1 - d2;

else

    len = [];
    for i = 1:length(text)
        if strcmp(text(i), ' ') || strcmp(text(i), ')')
            len = i;
            break;
        end
    end
    d = d + str2double(text(1:len-1));
end
end

function type = aux_getTypeOfConstraint(text)
% check if the current constraint is on the inputs or on the outputs

indX = regexp(text, 'X', 'once');
indY = regexp(text, 'Y', 'once');

% either X or Y must be given
if isempty(indX)
    if isempty(indY)
        % none given
        throw(CORAerror('CORA:notSupported', 'File format not supported'));
    else
        % Y is not empty
        type = 'output';
    end
elseif isempty(indY)
    % X is not empty
    type = 'input';
else
    % return smaller
    if indX(1) < indY(1)
        type = 'input';
    else
        type = 'output';
    end
end
end

function spec = aux_combineSafeSets(spec)
% combine all specifications involving safe sets to a single safe set

% find all specifications that define a safe set
ind = [];
for i = 1:length(spec)
    if strcmp(spec(i).type, 'safeSet')
        ind = [ind, i];
    end
end

% check if safe sets exist
if length(ind) > 1

    % combine safe sets to a single polytope
    poly = spec(ind(1)).set;
    for i = 2:length(ind)
        poly = poly & spec(ind(i)).set;
    end
    specNew = specification(poly, 'safeSet');

    % remove old specifications
    ind_ = setdiff(1:length(spec), ind);

    if isempty(ind_)
        spec = specNew;
    else
        spec = spec(ind_);
        spec = add(spec, specNew);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
