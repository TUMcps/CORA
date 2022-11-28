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
%    ---
%
% Example: 
%    x = stl('x',2)
%    until(x(1) < 5,x(2) < 3,interval(0.1,0.2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    str = displayRecursive(obj);
    disp(str);
end

% Auxiliary Functions -----------------------------------------------------

function str = displayRecursive(obj)
% recursive function to get display strings for each subformula

    switch obj.type
        
        case 'variable'
            str = obj.var;

        case 'true'
            str = 'true';

        case 'false'
            str = 'false';
            
        case 'until'
            str1 = displayRecursive(obj.lhs);
            str2 = displayRecursive(obj.rhs);
            str = ['(',str1,' U[',num2str(obj.from),',', ...
                                            num2str(obj.to),'] ',str2,')'];

        case 'release'
            str1 = displayRecursive(obj.lhs);
            str2 = displayRecursive(obj.rhs);
            str = ['(',str1,' R[',num2str(obj.from),',', ...
                                            num2str(obj.to),'] ',str2,')'];

        case 'globally'
            str1 = displayRecursive(obj.lhs);
            str = ['G[',num2str(obj.from),',',num2str(obj.to),'](',str1,')'];

        case 'finally'
            str1 = displayRecursive(obj.lhs);
            str = ['F[',num2str(obj.from),',',num2str(obj.to),'](',str1,')'];

        case 'next'
            str1 = displayRecursive(obj.lhs);
            str = ['X[',num2str(obj.from),'](',str1,')'];

        case '&'
            str1 = displayRecursive(obj.lhs);
            if ~ismember(obj.lhs.type,{'until','globally','finally', ...
                                                'next','true','false'})
                str1 = ['(',str1,')'];
            end
            str2 = displayRecursive(obj.rhs);
            if ~ismember(obj.rhs.type,{'until','globally','finally', ...
                                                'next','true','false'})
                str2 = ['(',str2,')'];
            end
            str = ['(',str1,' & ', str2,')'];

        case '|'
            str1 = displayRecursive(obj.lhs);
            if ~ismember(obj.lhs.type,{'until','globally','finally', ...
                                                'next','true','false'})
                str1 = ['(',str1,')'];
            end
            str2 = displayRecursive(obj.rhs);
            if ~ismember(obj.rhs.type,{'until','globally','finally', ...
                                                'next','true','false'})
                str2 = ['(',str2,')'];
            end
            str = ['(',str1,' | ', str2,')'];

        case '~'
            str1 = displayRecursive(obj.lhs);
            if ismember(obj.lhs.type,{'until','globally','finally', ...
                                                'next','true','false'})
                str = ['~',str1];
            else
                str = ['~(',str1,')'];
            end

        case '+'
            if isnumeric(obj.lhs)
                str1 = num2str(obj.lhs);
            else    
                str1 = displayRecursive(obj.lhs);
            end
            if isnumeric(obj.rhs)
                str2 = num2str(obj.rhs);
            else    
                str2 = displayRecursive(obj.rhs);
            end
            str = [str1, ' + ', str2];

        case '-'
            if isnumeric(obj.lhs)
                str1 = num2str(obj.lhs);
            else    
                str1 = displayRecursive(obj.lhs);
            end
            if isnumeric(obj.rhs)
                str2 = num2str(obj.rhs);
            else    
                str2 = displayRecursive(obj.rhs);
            end
            str = [str1, ' - ', str2];

        case '*'
            str1 = num2str(obj.lhs);
            str2 = displayRecursive(obj.rhs);
            str = [str1, '*', str2];
                                            
        case '<'
            str1 = displayRecursive(obj.lhs);
            if isnumeric(obj.rhs)
                str2 = num2str(obj.rhs);
            else    
                str2 = displayRecursive(obj.rhs);
            end
            str = [str1, ' < ',str2]; 

        case '<='
            str1 = displayRecursive(obj.lhs);
            if isnumeric(obj.rhs)
                str2 = num2str(obj.rhs);
            else    
                str2 = displayRecursive(obj.rhs);
            end
            str = [str1, ' <= ',str2]; 

        case '>'
            str1 = displayRecursive(obj.lhs);
            if isnumeric(obj.rhs)
                str2 = num2str(obj.rhs);
            else    
                str2 = displayRecursive(obj.rhs);
            end
            str = [str1, ' > ',str2]; 

        case '>='
            str1 = displayRecursive(obj.lhs);
            if isnumeric(obj.rhs)
                str2 = num2str(obj.rhs);
            else    
                str2 = displayRecursive(obj.rhs);
            end
            str = [str1, ' >= ',str2]; 
    end
end

%------------- END OF CODE --------------