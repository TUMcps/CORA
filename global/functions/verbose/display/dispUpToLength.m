function dispUpToLength(strList,maxLength,varargin)
% dispUpToLength - Displays a list of strings separated by the delimiter
%    ", " up to a maximum length per line; optionally, one can add a string
%    with which the first line should start, all further lines are
%    vertically aligned with the end of the first string in the first line;
%    the maximum length is exceeded if the next chosen string is too long
%
% Syntax:
%    dispUpToLength(strList,maxLength,initString)
%
% Inputs:
%    strList - string array
%    maxLength - maximum length of characters per line
%    initString - (optional) string at the beginning of the first line
%
% Outputs:
%    ---
%
% Example: 
%    strList = ["first entry"; "second entry"; ...
%               "third entry that is veeeeeeeeeeeeeeeeeeeeeeeeeeery long"];
%    maxLength = 60;
%    initString = "List of strings: ";
%    dispUpToLength(strList,maxLength,initString);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton/display, location/display

% Authors:       Mark Wetzlinger
% Written:       03-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default value
initString = setDefaultValues({""},varargin);

% check input arguments
inputArgsCheck({{strList,'att',{'string'},{'vector'}}; ...
                {maxLength,'att',{'numeric'},{'scalar','positive'}}; ...
                {initString,'att',{'string'},{'nonempty'}}});

% print initial string
fprintf(initString);

% length of initial string
initLength = strlength(initString);
% admissible length for each line (for vertical alignment)
admLength = maxLength - initLength;

% lengths of each entry
lengths = arrayfun(@(x)strlength(x),strList,'UniformOutput',true);
% increase all lenghts but last one by 2 to account for ", " delimiter
lengths(1:end-1) = lengths(1:end-1) + 2;

line = 1;
% loop over list of strings and check how many can be printed each iteration
while ~isempty(strList)

    % init logical index for printing
    printed = false(length(strList),1);

    % index up to which we can print (if first entry is already too long,
    % allow it anyway)
    idx = max([find(cumsum(lengths) < admLength,1,'last'), 1]);
    printed(1:idx) = true;

    % print white space
    if line > 1
        fprintf(strjoin(repmat(" ",initLength,1),""));
    end

    % print entries with delimiter
    fprintf(strjoin(strList(1:idx),", "));

    % remove entries from list of strings that have just been printed
    strList = strList(~printed);
    lengths = lengths(~printed);

    % print line break
    if ~isempty(strList)
        fprintf(',\n');
    else
        fprintf('\n');
    end

    % increment line counter
    line = line + 1;
end

% ------------------------------ END OF CODE ------------------------------
