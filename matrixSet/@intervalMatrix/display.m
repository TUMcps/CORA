function display(intMat)
% display - displays an intervalMatrix object on the command window
%
% Syntax:
%    display(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    ---
%
% Example:
%    intMat = intervalMatrix([1 2 3; 2 3 1],[1 0 2; 0 1 1])
%    intMat = intervalMatrix()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger, Tobias Ladner
% Written:       18-June-2010
% Last update:   03-April-2023 (MW, add empty case)
%                25-April-2024 (TL, harmonized display with contSet classes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(intMat)
    dispEmptyObj(intMat,inputname(1));
    return
end

% display input variable
fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display dimension
fprintf('%s:\n', class(intMat))
disp(['- dimension: ', num2str(dim(intMat))]);
fprintf(newline);

% call helper function (also used for interval)
displayInterval(intMat.int,false);

end

% ------------------------------ END OF CODE ------------------------------
