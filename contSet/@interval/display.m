function display(I)
% display - Displays the properties of an interal object (lower and upper
%    bounds) on the command window
%
% Syntax:  
%    display(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    ---
%
% Example: 
%    I = interval(2,3);
%    display(I);
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

if isemptyobject(I)
    
    dispEmptyObj(I,inputname(1));
    
else
    fprintf(newline);
    disp([inputname(1), ' =']);
    fprintf(newline);

    %determine size of interval
    [rows, cols] = size(I.inf);
    
    for i = 1:rows
        str = ' ';
        % display one row
        for j = 1:cols
            newStr = sprintf('[%0.4f, %0.4f]',I.inf(i,j),I.sup(i,j));
            str = [str,' ',newStr];
        end
        disp(str);
    end
    
    fprintf(newline);
end

%------------- END OF CODE --------------