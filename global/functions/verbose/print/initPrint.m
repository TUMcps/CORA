function [fid,closefid,obj,accuracy,doCompact,clearLine] = initPrint(varargin)
% initPrint - parses the given parameters to determine if output should be
%    printed to a file
%
% Syntax:
%    [fid,closefid,obj,accuracy,doCompact,clearLine] = initPrint(obj)
%    [fid,closefid,obj,accuracy,doCompact,clearLine] = initPrint(obj,accuracy,doCompact,clearLine)
%    [fid,closefid,obj,accuracy,doCompact,clearLine] = initPrint(filename,obj,varargin)
%    [fid,closefid,obj,accuracy,doCompact,clearLine] = initPrint(fid,obj,varargin)
%
% Inputs:
%    obj - object
%    accuracy - (optional) floating-point precision
%    doCompact - (optional) whether the matrix is printed compactly
%    clearLine - (optional) whether to finish with '\n'
%    filename - char, filename to print given obj to
%    fid - char, fid to print given obj to
%
% Outputs:
%    fid - file identifier (1 for stdout)
%    closefid - whether file should be closed at the end of the function
%    obj, accuracy, doCompact, clearLine - see above
%
%
% See also: printMatrix, printStruct, printSet, printSystem

% Authors:       Tobias Ladner
% Written:       20-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,5);

% set defaults
fid = 1;
closefid = false;

% determine if filename or fid is given -----------------------------------
param1 = varargin{1};
if ischar(param1) || isstring(param1)
    % filename is given; open it
    fid = fopen(param1,'w');
    closefid = true;
    varargin = varargin(2:end);
elseif nargin >= 2 && isnumeric(param1)
    % fid might be given; check if second entry is not accuracy
    param2 = varargin{2};
    if ischar(param2) || isstring(param2)
        % param2 is accuracy, does params1 is not fid
    else
        % fid is given
        fid = param1;
        varargin = varargin(2:end);
    end
end

% expected remaining parameters: obj, accuracy, doCompact, clearLine ------

% read out object
obj = varargin{1};
% set default values for remaining parameters
[accuracy,doCompact,clearLine] = setDefaultValues({'%4.3f',false,true},varargin(2:end));

% check inputs
inputArgsCheck({ ...
    {fid,'att','numeric','scalar'}, ...
    {obj,'att',{'numeric','logical','cell','struct','contSet','contDynamics','simResult','matrixSet'}}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

% handle 'high' accuracy
if ischar(accuracy) && strcmp(accuracy,'high')
    accuracy = '%16.16f';
end

end

% ------------------------------ END OF CODE ------------------------------
