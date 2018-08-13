function display(obj)
% display - Displays the left and right limit of the interval
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    ---
%
% Example: 
%    a = interval(2,3);
%    display(a);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-June-2015
% Last update:  22-February-2016 now it displays the name (Dmitry Grebenyuk)
% Last revision:---

%------------- BEGIN CODE --------------

%determine size of interval
[rows, columns] = size(obj.inf);

name = [inputname(1), ' = '];
disp(name)
for i = 1:rows
    str = [];
    % display one row
    for j = 1:columns
        newStr = sprintf('[%0.5f,%0.5f]',obj.inf(i,j),obj.sup(i,j));
        str = [str,' ',newStr];
    end
    disp(str);
end

%------------- END OF CODE --------------