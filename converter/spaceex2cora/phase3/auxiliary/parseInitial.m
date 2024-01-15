function [initialSet,inputSet] = parseInitial(configLine,fileDesc,state_names,input_names)
% parseInitial - Constructs CORA sets from initial values given in a
%    SpaceEx-ConfigFile
%
% Syntax:
%    [initialSet,inputSet] = parseInitial(configLine,fileDesc)
%
% Inputs:
%    configLine  - Line (String) from the config-file containig the initial
%                  values
%    fileDesc    - File descriptor for the config-file, in case the initial
%                  values are longer than a single line
%    state_names - Struct array containing structs with field "name"
%                  specifying the names of the state variables
%                  NOTE: Names must be in order!
%    input_names - Struct array containing structs with field "name"
%                  specifying the names of the input variables
%                  NOTE: Names must be in order!
%
% Outputs:
%    initialSet - CORA interval object denoting the initial sets of all
%                 state variables as mentioned and in the same order as 
%                 given in state_names
%    inputSet   - CORA interval object denoting the initial sets of all
%                 input variables as mentioned and in the same order as 
%                 given in input_names
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Maximilian Perschl
% Written:       05-September-2021
% Last update:   16-September-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

operators = ["==" "<=" ">=" "<" ">"];
% cut irrelevant parts off the string
configLine_split = split(configLine,'"');
% if the initial values are entirely specified in a single line, there are
% two occurences of " -> the size of the returned cell array is 3
% if there are multiple lines, there is only one " in the line -> the size
% of the returned cell array is 2
if length(configLine_split) == 2
    % get new lines until the "initial" segment is complete
    while count(configLine,'"') < 2
        configLine = strcat(configLine,fgetl(fileDesc));
    end
    configLine_split = split(configLine,'"');
end

% Assume entirety of initial specifications is
% contained within "configLine_split"
if length(configLine_split) == 3
    % split string along "&" to get the initial conditions for every
    % variable
    configLine_split = split(configLine_split(2),"&");
    % Possible numerical operators: ==,>=,<=,
    
    % Since CORA offers no open set representations, options "<" and ">"
    % are replaced with their closed form counterparts "<=" and ">=" and a
    % warning is issued
    
    % Initialize map storing variable names as keys along with their
    % intervals as values
    variableMap = containers.Map();
    
    for i = 1:length(configLine_split)
        skip_term = 0;
        currVar = configLine_split{i};
        currVarName = "";
        % iterate over string letter by letter, seperating
        % limits,operators, and names
        % loop changes between states 0 and 1
        % "reading operator"-1 and "reading non-operator(= variable name/bound)"-0
        % start reading non-operator
        state = 0;
        % while staying in one state, build string
        % processes when exiting the state
        currString = "";
        currBounds = [];
        currOperators = [];
        % this char array memorizes the order of occurences
        % e.g. bnb means "bound - name - bound"
        % this is in order to maintain readability of the code
        order = [];
        for j = 1:length(currVar)
            currChar = currVar(j);
            if state
                if contains("<=>",currChar)
                    currString = strcat(currString,currChar);
                    currOperators = [currOperators;currString];
                    currString = "";
                else
                    currOperators = [currOperators;currString];
                    currString = currChar;
                end
                % all operators are at max length 2, so we force a state
                % change here
                state = 0;
            else
                if contains("<=>",currChar)
                    % operator found, force state change
                    % determine if currString is a variable or a bound
                    if isempty(str2double(currString)) || isnan(str2double(currString))
                        % variable
                        
                        % check if variable name is already assigned, in this case, we
                        % skip the term because multi-variable conditions like x < y
                        % are ignored
                        if ~strcmp(currVarName,"")
                            skip_term = 1;
                        end
                        currVarName = currString;
                        order = [order "n"];
                    else
                        % bound
                        currBounds = [currBounds str2double(currString)];
                        order = [order "b"];
                    end
                    currString = currChar;
                    state = 1;
                else
                    % continue current lexing
                    currString = strcat(currString,currChar);
                end
            end
        end
        % parse final string, must be variable or bound
        if isempty(str2double(currString)) || isnan(str2double(currString))
            % variable
            
            % check if variable name is already assigned, in this case, we
            % skip the term because multi-variable conditions like x < y
            % are ignored
            if ~strcmp(currVarName,"")
                skip_term = 1;                
            end
            currVarName = currString;
            order = [order "n"];
        else
            % bound
            currBounds = [currBounds str2double(currString)];
            order = [order "b"];
        end
        
        if skip_term
            warning("Multi-variable term was skipped! Initial conditions like 'x < y' " ...
                +"are ignored within the spaceex2cora conversion!");
            continue;
        end
        
        % parse results into CORA interval
        currOperators = string(currOperators);
        
        % sanity check
        
        % if there are multiple different operators, the config file is
        % suboptimal e.g. 5<=x==7
        if length(unique(currOperators)) > 1
            throw(CORAerror('CORA:converterIssue',...
                ['SpaceEx configuration file is bloated, '...
                'check for initial conditions like "5 <= x == 7".']));
        end
        
        % default interval values for upper and lower bounds
        lb = -inf;
        ub = inf;
        
        % if the variable already exists, we overwrite the current bounds
        % this is included in order to properly catch cases such as 
        % 'x <= 4 & x >= 3' where one variable is defined over multiple
        % terms
        if variableMap.isKey(currVarName)
            lb = variableMap(currVarName).infimum;
            ub = variableMap(currVarName).supremum;
        end
        
        % iterate over list of collected names/bounds/operators
        op_count = 1;
        for j = 1:length(order)-1
            if order(j) == 'n'
                % set array bounds depending on operator
                switch(find(operators ==currOperators(op_count)))
                    case 1
                        % ==
                        lb = currBounds(op_count);
                        ub = currBounds(op_count);
                    case 2
                        % <=
                        ub = currBounds(op_count);
                    case 3
                        % >=
                        lb = currBounds(op_count);
                    case 4
                        % <
                        ub = currBounds(op_count);
                        warning("CORA only takes closed set representations into account! " ...
                            +"Changed initial set to closed set!");
                    case 5
                        % >
                        lb = currBounds(op_count);
                        warning("CORA only takes closed set representations into account! " ...
                            +"Changed initial set to closed set!");
                end
                op_count = op_count + 1;
            else
                % set array bounds depending on operator
                switch(find(operators ==currOperators(op_count)))
                    case 1
                        % ==
                        lb = currBounds(op_count);
                        ub = currBounds(op_count);
                    case 2
                        % <=
                        lb = currBounds(op_count);
                    case 3
                        % >=
                        ub = currBounds(op_count);
                    case 4
                        % <
                        lb = currBounds(op_count);
                        warning("CORA only takes closed set representations into account! " ...
                            +"Changed initial set to closed set!");
                    case 5
                        % >
                        ub = currBounds(op_count);
                        warning("CORA only takes closed set representations into account! " ...
                            +"Changed initial set to closed set!");
                end
                op_count = op_count + 1;
            end
        end
        % all information is obtained, now we can build the intervals
        variableMap(currVarName) = interval(lb,ub);        
    end
    
    % If the size is neither 2 nor 3, the cfg file contains syntax errors
else
    throw(CORAerror('CORA:converterIssue',...
        'SpaceEx syntax of the .cfg-file is not correct for the initial set.'));
end

% check if all states and inputs have been assigned initial sets
if length(state_names)+length(input_names) ~= variableMap.Count
    throw(CORAerror('CORA:converterIssue',...
        ['The configuration file does not specify an initial set '...
        'for all state and input variables']));
end

% state_names and input_names are ordered lists of names, which we use to
% define the initial and input sets in the correct order
initialSet = interval.empty(1);
for i = 1:length(state_names)
    initialSet = vertcat(initialSet,variableMap(state_names(i).name));
end
inputSet = interval.empty(1);
for i = 1:length(input_names)
    inputSet = vertcat(inputSet,variableMap(input_names(i).name));
end

% ------------------------------ END OF CODE ------------------------------
