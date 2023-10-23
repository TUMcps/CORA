function subscriptMatrix = i2s(siz,indexVector)
% i2s - extends the MATLAB function ind2sub so that it can be used
% conveniently for arbitrarily many dimensions. The function ind2sub 
% converts linear indices to subscripts
%
% Syntax:
%    subscriptMatrix = i2s(siz,indexVector)
%
% Inputs:
%    siz - vector of number of segments for each dimension as a row vector
%    ("size")
%    indexVector - vector of linear indices that should be converted into
%    subscripts
%
% Outputs:
%    subscriptMatrix - matrix of subscripts; each row corrsponds to one
%    index
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       14-September-2006
% Last update:   28-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    
%init subscriptMatrix
subscriptMatrix = [];

%subscript variable string
string = [];
for iChar = 1:length(siz)
    string = [string,'s',num2str(iChar),','];
end
string(end) = [];

%Generate command string
command=['[',string,']=ind2sub([',num2str(siz),'],[',num2str(indexVector),']);'];
eval(command);

%arrange variables in a vector
string = strrep(string,',',';');
command = ['subscriptMatrix=[',string,']'';'];
eval(command);
    
end
    
% ------------------------------ END OF CODE ------------------------------
