function list = splitIntoConvexSets(pgon)
% splitIntoConvexSets - split a polygon into a set of convex region
%
% Syntax:
%    list = splitIntoConvexSets(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    list - cell containing convex subsets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute triangulation
list = triangulation(pgon);

% unite neighbouring polygons to obtain larger convex shapes
while true

    finished = true;

    % loop over all polygons in the current list
    for i = 1:length(list)

        stop = false;

        % try to combine it with other polygons in the list
        for j = i + 1:length(list)

            tmp = list{i} | list{j};

            if isConvex(tmp)
                list{i} = tmp;
                list{j} = [];
                stop = true;
                finished = false;
                break;
            end
        end

        if stop
            break;
        end
    end

    % remove empty list entries
    list = list(~cellfun('isempty', list));

    if finished
        break;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
