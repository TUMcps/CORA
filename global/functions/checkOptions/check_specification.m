function options = check_specification(obj, options)
% check_specification - checks if check_specification
%  1) takes an allowed value
% rewrites all halfspaces s.t. A x <= b comparison enabled
%
% Syntax:
%    check_specification(options, obj)
%
% Inputs:
%    obj     - linear system
%    options - options for object
%
% Outputs:
%    options - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      29-July-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';
option = 'specification';
if isfield(options,option)
    
    % output dimension -> compare with defined halfspaces
    if isscalar(obj.C)
        outputDim = obj.dim;
    else
        outputDim = size(obj.C,1);
    end
    
    % 1. options.specification is of type cell
    if iscell(options.specification)
        if all(size(options.specification) - 1)
            % options.specification has to be a cell vector
            error(printOptionOutOfRange(obj, option, strct));
        elseif ~all(cellfun(@(spec) isa(spec,'halfspace'), options.specification))
            % all have to be halfspaces
            error(printOptionOutOfRange(obj, option, strct));
        elseif any(cellfun(@length, cellfun(@(hs) hs.c, options.specification, 'UniformOutput', false)) - outputDim)
            % normal vector of halfspaces must be of output dimension
            error(printOptionOutOfRange(obj,option, strct));
        end
        % concatenate halfspaces
        normalvectors = cellfun(@(hs) hs.c', options.specification, 'UniformOutput', false); 
        hs.distances = cellfun(@(hs) hs.d, options.specification); 
        for vec=1:length(normalvectors)
            hs.normals(vec,:) = normalvectors{vec};
        end
        options.specification = hs;
        % rewrite for check in loop, since it is faster than
        % ... ~all(cellfun(@(spec) in(spec,Rout{i-1}), options.specification))
        % which uses cell structure from user input
        
    % 1. converted into matrix for normals and vector for distances
    elseif isstruct(options.specification)
        % e.g. when called from copyOptions
        
        % check for missing fields
        if ~isfield(options.specification,'normals')
            error('options.specification: field .normals missing');
        elseif ~isfield(options.specification,'distances')
            error('options.specification: field .distances missing');
        end
        
        % check fields for correct definition
        if ~isnumeric(options.specification.normals) || ...
                ~isnumeric(options.specification.distances)
            % normals and distances have to be numeric
            error(printOptionOutOfRange(obj,option,strct));
        elseif size(options.specification.normals,2) ~= outputDim
            % normals have to be of output dimension
            error(printOptionOutOfRange(obj,option,strct));
        elseif ~isvector(options.specification.distances)
            % distances has to be a vector
            error(printOptionOutOfRange(obj,option,strct));
        elseif size(options.specification.normals,1) ~= ...
                length(options.specification.distances)
            % same number of normals as distances
            error(printOptionOutOfRange(obj,option,strct));
        end
        
    end
end


end

%------------- END OF CODE --------------

