function ME = CORAerror(identifier,varargin)
% CORAerror - central hub for all error messages thrown by CORA functions
%
% Syntax:
%    ME = CORAerror(identifier,varargin)
%
% Inputs:
%    identifier - name of CORA error
%                 'CORA:wrongInputInConstructor'
%                 'CORA:noInputInSetConstructor'
%                 'CORA:dimensionMismatch'
%                 'CORA:emptySet'
%                 'CORA:converterIssue'
%                 'CORA:configFile'
%                 'CORA:fileNotFound'
%                 'CORA:wrongValue'
%                 'CORA:noExactAlg'
%                 'CORA:noSpecificAlg'
%                 'CORA:solverIssue'
%                 'CORA:noSuitableSolver'
%                 'CORA:PLInonConvergent'
%                 'CORA:emptyProperty'
%                 'CORA:wrongFieldValue'
%                 'CORA:plotProperties'
%                 'CORA:notEnoughInputArgs'
%                 'CORA:tooManyInputArgs'
%                 'CORA:evenNumberInputArgs'
%                 'CORA:oddNumberInputArgs'
%                 'CORA:degenerateSet'
%                 'CORA:YALMIP'
%                 'CORA:outOfDomain'
%                 'CORA:specialError'
%                 'CORA:notSupported'
%                 'CORA:notDefined'
%                 'CORA:nnLayerNotSupported'
%                 'CORA:reachSetExplosion'
%                 'CORA:outOfMemory'
%                 'CORA:testFailed'
%                 'CORA:noops'
%                 'CORA:install'
%    varargin - further information depending on specific error
%
% Outputs:
%    ME - MException object for thrown error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mingrui Wang, Mark Wetzlinger
% Written:       05-April-2022
% Last update:   ---
% Last revision: 31-July-2023 (TL, workspace call)

% ------------------------------ BEGIN CODE -------------------------------

% full error stack
st = dbstack("-completenames");

% read info from error stack
try
    [filename,classname,functionname,sprintfFilesep] = aux_readInfo(st);
catch ME
    if strcmp(ME.identifier, 'MATLAB:needMoreRhsOutputs')
        % show 'call from workspace' warning
        warning(['You are very likely attempting to call CORAerror ',...
                'from the MATLAB console; please do not do that, as ',...
                'CORAerror will then be unable to trace back from where ',...
                'the error originates.']);

        % set default values
        filename = 'filename';
        classname = 'classname';
        functionname = 'functionname';
        sprintfFilesep = '/';
    else
        rethrow(ME);
    end
end

% standard help message
if ~strcmp(classname,functionname)
    helpmsg = sprintf("  Type 'help %s%s%s' for more information.",...
        classname,sprintfFilesep,functionname);
else
    % different call for constructors
    helpmsg = sprintf("  Type 'help %s' for more information.",classname);
end

% copy to name-value pairs for processing in some errors
NVpairs = varargin;


% generate error message depending on identifier and further input args
switch identifier

    % constructor is not properly called; input args:
    % - information message about what went wrong
    case 'CORA:wrongInputInConstructor'
        errmsg_form = 'Wrong input arguments for constructor of class: %s\n  %s \n%s';
        infomsg = varargin{1};
        errmsg = sprintf(errmsg_form,classname,infomsg,helpmsg);


    % constuctor of a contSet class called without input arguments:
    % - information message about contSet.empty and contSet.Inf functions
    case 'CORA:noInputInSetConstructor'
        errmsg_form = 'No input arguments for constructor of class: %s\n  %s \n%s';
        infomsg = ['Please consider calling ' classname '.empty or ' classname '.Inf instead.'];
        errmsg = sprintf(errmsg_form,classname,infomsg,helpmsg);


    % dimension mismatch between sets (contSet objects, numeric); input
    % args: two objects
    case 'CORA:dimensionMismatch'
        obj1 = varargin{1};
        obj2 = varargin{2};

        name1 = class(obj1);
        name2 = class(obj2);
        if isa(obj1,'contSet') && ~isa(obj1,'interval') && ~isa(obj1,'taylm') ...
                && ~isa(obj1,'zoo') && ~isa(obj1,'affine')
            dimOrSize1 = dim(obj1);
        else % range-bounding, matrixSet, numeric
            dimOrSize1 = size(obj1);
        end
        if isa(obj2,'contSet') && ~isa(obj2,'interval') && ~isa(obj2,'taylm') ...
                && ~isa(obj2,'zoo') && ~isa(obj2,'affine')
            dimOrSize2 = dim(obj2);
        else % range-bounding, matrixSet, numeric
            dimOrSize2 = size(obj2);
        end

        % splice error message: part1
        errmsg = ['The first object (' name1 ') '];
        if isscalar(dimOrSize1)
            errmsg = [errmsg 'is defined in dimension ' num2str(dimOrSize1)];
        else
            errmsg = [errmsg 'has size ' num2str(dimOrSize1(1)) '-by-' num2str(dimOrSize1(2))];
        end

        % ...part2
        errmsg = [errmsg ',\nbut the second object (' name2 ') '];
        if isscalar(dimOrSize2)
            errmsg = [errmsg 'is defined in dimension ' num2str(dimOrSize2)];
        else
            errmsg = [errmsg 'has size ' num2str(dimOrSize2(1)) '-by-' num2str(dimOrSize2(2))];
        end
        errmsg = [errmsg '.'];


    % the given set is empty, thus the operation is not defined
    case 'CORA:emptySet'
        errmsg = 'Set is empty!';


    % error in a converter (the name of that converter is automatically
    % read from the call stack; input args:
    % - (char) information about error
    case 'CORA:converterIssue'
        infomsg = varargin{1};
        try
            convertername = aux_readConverter(st);
            errmsg = ['Error in converter ' convertername ':\n  ' infomsg];
        catch
            errmsg = ['Error in converter:\n  ' infomsg];
        end


    % error in generation of configuration file (usually developer mistake)
    % - (char) description of issue
    case 'CORA:configFile'
        issue = varargin{1};
        errmsg_form = 'Error in generation of configuration file:\n   %s';
        errmsg = sprintf(errmsg_form,issue);


    % provided file does not exist
    % - (char) file name
    case 'CORA:fileNotFound'
        providedfile = varargin{1};
        errmsg = ['File with name ' providedfile ' could not be found.'];


    % a given input argument receives a wrong value; input args:
    % - (char) number of input argument, e.g.,
    %              'first', 'second', 'first and/or second', etc.
    %          OR description of name-value pairs starting(!) with
    %              'name-value pair'
    % - (char/string) description of range of admissible values
    case 'CORA:wrongValue'
        % number of input argument
        arg = varargin{1};
        % description of expected value
        explains = varargin{2};

        % error message
        if contains(arg,'name-value pair')
            errmsg_form = "Wrong value for %s.\n" + ...
                "  The right value: %s\n%s";
        else
            errmsg_form = "Wrong value for the %s input argument.\n" + ...
                "  The right value: %s\n%s";
        end
        errmsg = sprintf(errmsg_form,arg,explains,helpmsg);


    % no exact algorithm implemented (mostly because there is none)
    % - (contSet) object of first/second/... contSet class in operation
    case 'CORA:noExactAlg'
        classlist = "";
        for i=1:length(varargin)-1
            classlist = classlist + string(class(varargin{i})) + ", ";
        end
        if ~isstring(varargin{end})
            classlist = classlist + string(class(varargin{end})) + ".";
            errmsg = sprintf(...
                "There is no exact algorithm for function %s with input arguments:\n  %s",...
                filename,classlist);
        else
            classlist = classlist;
            infomsg = varargin{end};
            errmsg = sprintf(...
                "There is no exact algorithm for function %s with input arguments:\n  %s %s",...
                filename,classlist,infomsg);
        end
        
    % a specific algorithm is not implemented
    case 'CORA:noSpecificAlg'
        algorithm = varargin{1};
        classlist = "";
        for i=2:length(varargin)-1
            classlist = classlist + string(class(varargin{i})) + ", ";
        end
        if ~isstring(varargin{end})
            classlist = classlist + string(class(varargin{end})) + ".";
            errmsg = sprintf(...
                "There is no algorithm '%s' for function %s with input arguments:\n  %s",...
                algorithm,filename,classlist);
        else
            classlist = classlist;
            infomsg = varargin{end};
            errmsg = sprintf(...
                "There is no algorithm '%s' for function %s with input arguments:\n  %s %s",...
                algorithm,filename,classlist,infomsg);
        end

    
    % solver used within operation is not working properly
    % - (char/string) name of solver
    case 'CORA:solverIssue'
        solver = '';
        if ~isempty(varargin)
            solver = [' (',varargin{1},')'];
        end
        errmsg = sprintf(...
            "Solver%s in %s failed due to numerical/other issues! \n %s",...
            solver,filename,helpmsg);
    

    % no suitable solver for problem was found; input args:
    % - (char/string) name of solver
    case 'CORA:noSuitableSolver'
        type = '';
        if ~isempty(varargin)
            type = varargin{1};
        end
        errmsg = sprintf("No suitable %s solver found!",type);


    % Picard-Lindeloef iteration did not converge (e.g., reachInner)
    % - (char/string) what did not converge?
    case 'CORA:notConverged'
        issue = varargin{1};
        errmsg_form = '%s did not converge.';
        errmsg = sprintf(errmsg_form,issue);


    case 'CORA:emptyProperty'
        errmsg = 'Property is empty.';


    % unexpected value of a struct field; input args:
    % - (char) name of field
    % - (cell array of char) list of admissible values for given field
    case 'CORA:wrongFieldValue'
        % name of field
        field = varargin{1};
        % allowed values for field
        allowedValues = "'" + strjoin(varargin{2},"', '") + "'";

        % error message
        errmsg_form = "Wrong value for the field '%s'.\n" + ...
            "  Admissible values: %s\n";
        errmsg = sprintf(errmsg_form,field,allowedValues);


    % error in plotting of sets; input args:
    % - (numeric) number of dimensions for which there is a problem
    case 'CORA:plotProperties'
        num = varargin{1};
        if num == 1
            errmsg = 'At least one dimension has to be specified.';
        elseif num == 2 % deprecated
            errmsg = 'At least two dimensions have to be specified.';
        elseif num == 3
            errmsg = 'Only up to three dimensions can be specified.';
        else
            errmsg = 'Incorrect number of dimensions specified.';
        end


    % function does not receive enough input arguments; input args:
    % - (numeric) number of minimally required input arguments
    case 'CORA:notEnoughInputArgs'
        input_num = varargin{1};
        errmsg_form = 'The function %s requires at least %d input argument(s).\n %s';
        errmsg = sprintf(errmsg_form,filename,input_num,helpmsg);


    % function receives too many input arguments (should be used for all
    % functions with varargin apart from name-value pairs); input args:
    % - (numeric) maximum number of input arguments
    case 'CORA:tooManyInputArgs'
        input_num = varargin{1};
        errmsg_form = 'The function %s requires no more than %d input argument(s).\n %s';
        errmsg = sprintf(errmsg_form,filename,input_num,helpmsg);


    % function requires an even number of input arguments (likely because
    % of name-value pairs, e.g., generateRandom)
    case 'CORA:evenNumberInputArgs'
        errmsg_form = 'The function %s requires an even number of input argument(s).\n %s';
        errmsg = sprintf(errmsg_form,filename,helpmsg);

    % function requires an odd number of input arguments (likely because
    % of name-value pairs, e.g., generateRandom)
    case 'CORA:oddNumberInputArgs'
        errmsg_form = 'The function %s requires an odd number of input argument(s).\n %s';
        errmsg = sprintf(errmsg_form,filename,helpmsg);

    % function takes name-value pairs, but provided pair has a name which
    % is not within the list of admissible names for all pairs; input args:
    % - (char) name of redundant field
    % - (cell array of char) list of admissible names
    case 'CORA:redundantNameValuePair'
        name = varargin{1};
        listofnames = varargin{2};
        errmsg_form = ['The name-value pair with the name ''%s'''...
            ' is redundant for the function %s.\n' ...
            'The list of admissible names is: ''%s''.\n%s'];
        errmsg = sprintf(errmsg_form,name,filename,...
            strjoin(listofnames,''', '''),helpmsg);


    % error occurring due to a degenerate set; input args:
    % - (char) further information about error
    case 'CORA:degenerateSet'
        errmsg = varargin{1};


    % YALMIP not working correctly; input args:
    % - (char) further information about error
    case 'CORA:YALMIP'
        errmsg = varargin{1};


    % function cannot be executed as the input argument is out of domain,
    % mainly used for range bounding operations; input args: name-value
    % pairs with names
    % - 'validDomain': description of valid domain
    case 'CORA:outOfDomain'
        [~,validDomain] = readNameValuePair(NVpairs,'validDomain');
        errmsg_form = ['The interval is not inside the valid domain of the function %s.'...
            '\n  Valid domain is: %s \n %s'];
        errmsg = sprintf(errmsg_form,filename,validDomain,helpmsg);


    % specific errors which strongly depend on the context of the
    % operation; input args:
    % - (char) additional information about error
    case 'CORA:specialError'
        errmsg = varargin{1};


    % functionality is not supported, mostly occurs due to intricate
    % combination of input arguments which cannot be sensibly checked
    % beforehand; input args:
    % - (char) information about unsupported functionality
    case 'CORA:notSupported'
        errmsg = varargin{1};


    % functionality is not defined, mostly due to 
    % mathematical inconsistencies; input args:
    % - (char) information about undefined functionality
    % - (char) See also:
    case 'CORA:notDefined'
        errmsg = sprintf('Undefined functionality: %s', varargin{1});
        if length(varargin) == 2
            errmsg = sprintf('%s\nSee also: %s', errmsg, varargin{2});
        end


    % propagation computation is not supported by a layer
    % input args:
    % - (nnLayer) layer object
    % - (char) unsupported computation
    case 'CORA:nnLayerNotSupported'
        layer = varargin{1};
        unsupComp = varargin{2};
        errmsg = sprintf("Computation of %s is not supported by %s (%s).", ...
            unsupComp, class(layer), layer.name);

    
    % reachable set explodes due to excessively large abstraction error
    case 'CORA:reachSetExplosion'
        errmsg = 'Abort analysis due to reachable set explosion!';


    % matlab runs out of memory while computing something within CORA,
    % e.g. if 'MATLAB:array:SizeLimitExceeded' would be thrown otherwise
    % input args:
    % - (char) additional information, e.g. how to fix it.
    % - (optional) MException object
    case 'CORA:outOfMemory'
        errmsg = sprintf('Out of memory. ');
        if length(varargin) == 2
            ME = varargin{2};
            errmsg = sprintf('%s%s\n',errmsg,ME.message);
        end
        errmsg = sprintf('%s%s',errmsg,varargin{1});


    % standard message for failed unit tests (sometimes quicker than
    % exiting a number of loops and returning false)
    case 'CORA:testFailed'
        errmsg = 'Unit test failed.';

        
    % certain set operation not implemented for given input arguments
    case 'CORA:noops'
        classlist = "";
        addinfo = "";
        punct = ", ";
        for i=1:length(varargin)
            if strcmp(string(class(varargin{i})),'double')
                % additionally print which kind of double
                if isscalar(varargin{i})
                    addinfo = "-scalar";
                elseif isvector(varargin{i})
                    if size(varargin{i},1) == 1
                        addinfo = "-vector (row)";
                    else
                        addinfo = "-vector (column)";
                    end
                else
                    addinfo = "-matrix";
                end
            elseif strcmp(string(class(varargin{i})),'levelSet')
                % additionally print comparison operator(s)
                compOps = strjoin(cellstr(unique(varargin{i}.compOp)),',');
                addinfo = " (comparison operator: '" + compOps + "')";
            end
            if i == length(varargin)
                punct = ".";
            end
            classlist = classlist + string(class(varargin{i})) + addinfo + punct;
            addinfo = "";
        end
        errmsg = sprintf(...
            "The function '%s' is not implemented for the following arguments:\n  %s \n%s",...
            filename,classlist,helpmsg);


    % standard message for failed installation
    case 'CORA:install'
        errmsg = varargin{1};


    % handle non-defined identifiers
    otherwise
        % note: this function can never be a CORA error!
        error("Bug: Identifier not specified!");

end

% instantiate MException object
ME = MException(identifier,errmsg);

end


% Auxiliary functions -----------------------------------------------------

function [filename,classname,functionname,sprintfFilesep] = aux_readInfo(st)

% stack length
stlength = length(st);

% find index in stack where error comes from
errIdx = [];
for i=1:stlength
    if ~any(strcmp(st(i).name,{'CORAerror','inputArgsCheck'}))
        errIdx = i;
        % checkNameValuePairs shifts index by one
        if strcmp(st(i).name,'checkNameValuePairs')
            errIdx = i+1;
        end
        break
    end
end

% name of file where error occurred
filename = st(errIdx).name;
if contains(filename,'.')
    % likely a constructor -> remove part until dot
    filename = filename(strfind(filename,'.')+1:end);
end
filename = [filename '.m'];

% position of all file separators
filesepPos = strfind(st(errIdx).file,filesep);

% position of @ (empty if not a class function)
atPos = strfind(st(errIdx).file,'@');
% position of . (from .m-file)
dotPos = strfind(st(errIdx).file,'.');

if ~isempty(atPos)
    % file separator after @-sign
    filesepPosAfterAt = filesepPos(filesepPos > atPos);
    
    % read out classname
    classname = extractBetween(st(errIdx).file,atPos+1,filesepPosAfterAt(1)-1);
    classname = classname{1};

    % read out functionname (includes 'private'!)
    functionname = extractBetween(st(errIdx).file,filesepPosAfterAt(1)+1,dotPos-1);
    functionname = functionname{1};
else
    % classname now name of folder (starting from within coraroot)
    fileSepAfterCORA = filesepPos(filesepPos > length(CORAROOT));

    if isscalar(fileSepAfterCORA)
        % in root directory
        classname = 'CORAROOT';
        functionname = st(errIdx).file(fileSepAfterCORA(1:end-1));
    else
        classname = extractBetween(st(errIdx).file,fileSepAfterCORA(1)+1,fileSepAfterCORA(end)-1);
        classname = classname{1};
    
        % read out functionname (includes 'private'!)
        functionname = extractBetween(st(errIdx).file,fileSepAfterCORA(end)+1,dotPos-1);
        functionname = functionname{1};
    end
end

% adapt filesep so that sprintf does not have an issue
if filesep == '\'
    sprintfFilesep = '\\';
    if contains(functionname,'\')
        functionname = strrep(functionname,filesep,'\\');
    end
    if isempty(atPos) && contains(classname,'\')
        % classname can be also include filesep in this case
        classname = strrep(classname,filesep,'\\');
    end
else
    sprintfFilesep = filesep;
end

end

function convertername = aux_readConverter(st)
% read the name of the currently implemented converters from the call stack

% calling function (where error occurs) is at index 2
temp = st(2).file;

% read out name of converter
if contains(temp,[CORAROOT filesep 'converter'])
    convPos = strfind(temp,'converter');
    temp = temp(convPos+10:end);
    filesepPos = strfind(temp,filesep);
    convertername = temp(1:filesepPos(1)-1);
else
    % string 'converter' cannot be found in the calling function
    throw(CORAerror('CORA:specialError',...
        'Converter error thrown outside of converter directory.'));
end

end

% ------------------------------ END OF CODE ------------------------------
