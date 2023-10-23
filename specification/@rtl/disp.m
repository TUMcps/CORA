function disp(obj)
% disp - override disp function to show object on command window
%
% Syntax:
%    disp(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    -
%
% Example: 
%    x = stl('x',2)
%    rtl(x(1) < 5 | x(2) > 3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    str = aux_displayRecursive(obj);
    disp(str);
end


% Auxiliary functions -----------------------------------------------------

function str = aux_displayRecursive(obj)
% recursive function to get display strings for each subformula

    switch obj.type
        
        case 'all'
            str = ['A(',char(formattedDisplayText(obj.lhs)),')'];

        case 'next'
            str1 = aux_displayRecursive(obj.lhs);
            str = ['X[',num2str(obj.time),'](',str1,')'];

        case '&'
            str1 = aux_displayRecursive(obj.lhs);
            str2 = aux_displayRecursive(obj.rhs);
            str = ['(',str1,' & ', str2,')'];

        case '|'
            str1 = aux_displayRecursive(obj.lhs);
            str2 = aux_displayRecursive(obj.rhs);
            str = ['(',str1,' | ', str2,')'];
    end
end

% ------------------------------ END OF CODE ------------------------------
