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

% Authors:       Matthias Althoff
% Written:       19-June-2015
% Last update:   22-February-2016 (DG, now it displays the name)
%                01-May-2020 (MW, handling of empty case)
%                11-September-2023 (TL, respect output display format)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases (only vector)
if size(I,2) <= 1
    if representsa(I,'emptySet')
        dispEmptySet(I,inputname(1));
        return
    elseif representsa(I,'fullspace')
        dispRn(I,inputname(1));
        return
    end
end

% display input variable
fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension
display@contSet(I);
fprintf(newline);

% call helper function (also used for intervalMatrix)
displayInterval(I,false);

% ------------------------------ END OF CODE ------------------------------
