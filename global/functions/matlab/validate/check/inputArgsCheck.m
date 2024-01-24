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
%                   attributes - (only id = 'att') admissible types (check
%                                MATLAB function validateattributes) as a
%                                cell-array of chars or single char-array
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
% See also: validateattributes

% Authors:       Mingrui Wang, Mark Wetzlinger
% Written:       30-May-2022
% Last update:   23-January-2024 (MW, exact match for strings)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if disabled by macro
if ~CHECKS_ENABLED
    return;
end

% turn off warnings (if no attributes are given)
warOrig = warning;
warning('off','all');

% number of input arguments to original function that have to be checked
nrInputArgs = length(inputArgs); 

% loop over all input arguments
for i = 1:nrInputArgs
    % read information about i-th input argument
    inputArg = inputArgs{i};
    % read value of i-th input argument
    arg_name = inputArg{1};
    % read identifier ('att' or 'str')
    identifier = inputArg{2};

    % case destinction
    switch identifier
        case 'att'    % check classname (and attributes) in this case

            % class of i-th input argument
            if iscell(inputArg{3})
                class = inputArg{3};
            else
                class = inputArg(3);
            end

            % 4-th entry represents required attributes (optional)
            if length(inputArg) == 3
                % only classes have to be checked
                found = false;
                for c=1:length(class)
                    if isa(arg_name,class{c})
                        found = true; break
                    end
                end
                if ~found
                    % throw error
                    throw(CORAerror('CORA:wrongValue',...
                        aux_countingNumber(i),aux_errMsgClassAttributes(class)));
                end

            else
                % attributes are given

                % different cases
                if ~iscell(inputArg{4})
                    % case 1: 'scalar' attribute: same for all classes
                    attributes = inputArg(4);
                    % check attribute using built-in MATLAB function
                    try
                        validateattributes(arg_name,class,attributes);
                    catch
                        % throw error
                        throw(CORAerror('CORA:wrongValue',aux_countingNumber(i),...
                            aux_errMsgClassAttributes(class,attributes)));
                    end
                else
                    % check
                    entryCells = cellfun(@(x) iscell(x),inputArg{4},'UniformOutput',true);
                    attributes = inputArg{4};

                    if ~any(entryCells)
                        % case 2: cell-array of attributes: same set of
                        % checks for all classes
                        
                        % check attributes using built-in MATLAB function
                        try
                            validateattributes(arg_name,class,attributes);
                        catch ME
                            % throw error
                            throw(CORAerror('CORA:wrongValue',aux_countingNumber(i),...
                                aux_errMsgClassAttributes(class,attributes)));
                        end
                    elseif all(entryCells)
                        % case 3: cell-array of cell-arrays: individual set
                        % of checks for all classes
                        
                        % note: only one check has to be ok
                        checkok = false(length(class),1);
                        for c=1:length(class)
                            % check attribute using built-in MATLAB function
                            try
                                validateattributes(arg_name,class(c),attributes(c));
                                checkok(c) = true;
                            catch
                                checkok(c) = false;
                            end
                        end
                        % throw error if no combination is ok
                        if ~any(checkok)
                            throw(CORAerror('CORA:wrongValue',aux_countingNumber(i),...
                                        aux_errMsgClassAttributes(class,attributes)));
                        end
                    else
                        % some are cells, but some are not -> should not happen
                        throw(CORAerror('CORA:specialError','Unknown case'));
                    end
                end
            end

        case 'str'    % check the strings in this case

            % read string
            if iscell(inputArg{3})
                validateStr = inputArg{3};
            else
                validateStr = inputArg(3);
            end

            % check exact match with admissible values
            if ~ismember(arg_name,validateStr)
                % generate string of admissible values (user info)
                validrange = ['''', strjoin(validateStr,"', '"), ''''];

                % throw error
                throw(CORAerror('CORA:wrongValue',aux_countingNumber(i),...
                    validrange));
            end

        otherwise
            throw(CORAerror("CORA:wrongValue",'second','''att'' or ''str''.'))

    end

end

% turn warnings back on
warning(warOrig);

end


% Auxiliary functions -----------------------------------------------------

function text = aux_errMsgClassAttributes(class,varargin)

% optional argument: attributes
if ~isempty(varargin)
    attributes = varargin{1};

    if ~iscell(attributes)
        n_attribute = 1;
        varyingAttributes = false;
    else
        entryCells = cellfun(@(x) iscell(x),attributes,'UniformOutput',true);
        if ~any(entryCells)
            varyingAttributes = false;
            % number of attributes (can be zero)
            n_attribute = length(attributes);
            % convert all numerics from attributes to string (not ideal
            % output text, but works for now...)
            for i=1:n_attribute
                if isnumeric(attributes{i})
                    attributes{i} = num2str(attributes{i});
                end
            end
            if n_attribute == 1 && strcmp(attributes{1},'')
                % deprecated case with '' -> counts as zero attributes
                n_attribute = 0;
            end
        else
            varyingAttributes = true;
            n_attribute = zeros(length(attributes),1);
            for c=1:length(attributes)
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
for c=1:n_class
    % append text
    text = text + string(class{c});

    if varyingAttributes
        if n_attribute(c) > 0
            text = text + " (" + ...
                strjoin(string(attributes{c}),', ') + ")";
        end
    else
        if n_attribute > 0
            text = text + " (" + ...
                strjoin(string(attributes),', ') + ")";
        end
    end

    if c < n_class
        % 17 white spaces for correct alignment
        text = text + ",\n" ...
            + "                   ";
    end
end

end

function text = aux_countingNumber(i)
% returns counting number in text form
switch i
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
        text = ['No. ', num2str(i)];
end

end

% ------------------------------ END OF CODE ------------------------------
