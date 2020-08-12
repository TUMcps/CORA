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
% Last update:  22-February-2016 (DG, now it displays the name)
%               01-May-2020 (MW, handling of empty case)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(obj)
    dispEmptyObj(obj,inputname(1));
    
else
    fprintf(newline);
    name = [inputname(1), ' ='];
    disp(name)
    fprintf(newline);

    %determine size of interval
    [rows, cols] = size(obj.inf);
    
    for i = 1:rows
        str = ' ';
        % display one row
        for j = 1:cols
            newStr = sprintf('[%0.4f, %0.4f]',obj.inf(i,j),obj.sup(i,j));
            str = [str,' ',newStr];
        end
        disp(str);
    end
    
    fprintf(newline);
end

%------------- END OF CODE --------------