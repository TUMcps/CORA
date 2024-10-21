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
%                27-September-2024 (MW, all-zero sparse case)
%                18-October-2024 (TL, n-d intervals)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases (only vector)
if numel(size(I)) <= 2 && any(size(I,2) <= 1)
    if representsa(I,'emptySet')
        dispEmptySet(I,inputname(1));
        return
    elseif representsa(I,'fullspace')
        dispRn(I,inputname(1));
        return
    end
end

% all-zero and sparse
if issparse(I) && representsa_(I,'origin',0)
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);

    % display
    disp("   All zero sparse interval: " + strjoin(string(size(I)),"x"));
    fprintf(newline);
    return
end

% display input variable
varname = inputname(1);
fprintf(newline);
disp(varname + " =");
fprintf(newline);

%display dimension
display@contSet(I);
fprintf(newline);

% check dimension
dims = size(I);
if numel(dims) <= 2
    % display 2-dimensional interval
    aux_display2D(I);
else
    % display n-dimensional interval page-wise
    aux_displayND(I,varname)
end

end


% Auxiliary functions -----------------------------------------------------

function aux_display2D(I)
    % display 2-dimensional interval
    
    % call helper function (also used for intervalMatrix)
    displayInterval(I,false);

end

function aux_displayND(I,varname)
    % display n-dimensional interval page-wise

    % determine number of pages
    dims = size(I);
    dims_pages = dims(3:end);
    numPages = prod(dims_pages);
    I_pages = reshape(I,dims(1),dims(2),numPages);

    % function to convert linear indices to matrix indices 
    % (similar to ind2sub)
    convertIndex = @(i) arrayfun(@(d) mod(floor((i-1) ./ prod(dims_pages(1:d-1))), dims_pages(d)) + 1, 1:length(dims_pages));

    for i=1:numPages
        % convert indices
        page_idx = convertIndex(i);
        
        % display input variable
        fprintf('%s(:,:,%s) = ', varname, strrep(num2str(page_idx),'  ',','));
        fprintf(newline);

        % call helper function for current page
        subs = struct('type', '()','subs', {{':',':',i}});
        displayInterval(subsref(I_pages,subs),false);
    end

end

% ------------------------------ END OF CODE ------------------------------
