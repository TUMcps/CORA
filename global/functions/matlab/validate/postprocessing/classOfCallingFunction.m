function classname = classOfCallingFunction()
% classOfCallingFunction - returns the name of the class from which the
%    function (which class this function) has been called
%
% Syntax:
%    classname = classOfCallingFunction()
%
% Inputs:
%    -
%
% Outputs:
%    classname - char-array of a class
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call stack
st = dbstack("-completenames");

% name of calling function is at index 3
fullname = st(3).file;

% index of @
idxAt = strfind(fullname,'@');

% indices of delimiters
idxDelim = strfind(fullname,filesep);
% first index after @
idxDelim = idxDelim(idxDelim > idxAt);

% read out classname
classname = fullname(idxAt+1:idxDelim(1)-1);

% ------------------------------ END OF CODE ------------------------------
