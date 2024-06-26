function displayInterval(I,showInputname)
% displayInterval - Displays the properties of an interal object
%
% Syntax:
%    displayInterval(I,showInputname)
%
% Inputs:
%    I - interval object
%    showInputname - whether inputname should be plot
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
% See also: interval/display

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       19-June-2015
% Last update:   22-February-2016 (DG, now it displays the name)
%                01-May-2020 (MW, handling of empty case)
%                11-September-2023 (TL, respect output display format)
%                25-April-2024 (TL, helper functions for intervalMatrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 1
    showInputname = false;
end


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

% show inputname if desired
if showInputname
    % display input variable
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    %display dimension
    display@contSet(I);
    fprintf(newline);
end

if issparse(I)

    % non-zero indices of infimum
    [inf_row,inf_col] = find(I.inf);
    % non-zero indices of supremum
    [sup_row,sup_col] = find(I.sup);

    % for sparse matrices, we plot all indices where either the infimum
    % or the supremum is non-zero
    indices = [[inf_row inf_col];[sup_row sup_col]];
    indices = unique(indices,'rows');
    
    % loop over all indices
    for i=1:size(indices,1)
        idxRow = indices(i,1);
        idxCol = indices(i,2);
        % set index
        idxStr = sprintf("(%i,%i)",idxRow,idxCol);
        % values
        lb = full(I.inf(idxRow,idxCol));
        ub = full(I.sup(idxRow,idxCol));
        % print value
        lb = strtrim(formattedDisplayText(lb));
        ub = strtrim(formattedDisplayText(ub));
        fprintf('  %-8s [%s, %s]\n',idxStr,lb,ub);
    end

else

    %determine size of interval
    [rows, cols] = size(I.inf);

    for i = 1:rows
        str = ' ';
        % display one row
        for j = 1:cols
            lb = strtrim(formattedDisplayText(I.inf(i,j)));
            ub = strtrim(formattedDisplayText(I.sup(i,j)));
            newStr = sprintf('[%s, %s]',lb,ub);
            str = [str,' ',newStr];
        end
        disp(str);
    end

end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
