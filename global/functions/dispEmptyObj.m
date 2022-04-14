function dispEmptyObj(obj,argname)
% dispEmptyObj - displays text for empty objects
%
% Syntax:  
%    dispEmptyObj(obj)
%
% Inputs:
%    obj - object of some class
%    argname - name of obj in workspace
%
% Outputs:
%    -
%
% Example:
%    I = interval();
%    disp(I); % calls dispEmptyObj
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: display
%
% Author:       Mark Wetzlinger
% Written:      01-May-2020
% Last update:  02-May-2020
% Last revision:---

%------------- BEGIN CODE --------------

fprintf(newline);
disp(argname + " =");
fprintf(newline);
disp("  1x1 empty " + class(obj));
fprintf(newline);


end

%------------- END OF CODE --------------