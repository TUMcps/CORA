function S_out = single(S)
% single - converts all numeric properties of a set to single
%
% Syntax:
%    S_out = single(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    S_out - contSet object with all numeric properties converted to single
%
% Example:
%    Z = zonotope([1;0],[1 0 -1; 0 1 1]);
%    Z_sgl = single(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: double

% Authors:       Tobias Ladner
% Written:       02-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get property info (throws CORA:noops if not supported)
[~,propertyOrder] = getPrintSetInfo(S);

% convert each property
args = cell(1,numel(propertyOrder));
for p = 1:numel(propertyOrder)
    args{p} = aux_convertProperty(S.(propertyOrder{p}));
end

% reconstruct set
S_out = feval(class(S),args{:});

end


% Auxiliary functions -----------------------------------------------------

function property = aux_convertProperty(property)
    if isnumeric(property)
        property = single(property);
    elseif iscell(property)
        for i = 1:numel(property)
            property{i} = aux_convertProperty(property{i});
        end
    elseif isstruct(property)
        fields = fieldnames(property);
        for i = 1:numel(fields)
            property.(fields{i}) = aux_convertProperty(property.(fields{i}));
        end
    elseif isa(property,'contSet')
        property = single(property);
    end
end


% ------------------------------ END OF CODE ------------------------------
