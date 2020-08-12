function options = check_pHA_inputs(options, obj)
% check_pHA_inputs - checks if inputs U, inputCompMap
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    options = check_pHA_inputs(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - parallel hybrid automaton object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Niklas Kochdumper
% Written:      15-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

numComps = length(obj.components);

% check if Uloc is specified -> warning that it is redundant
if isfield(options,'Uloc')
    warning('params.Uloc is redundant!');
end

option = 'U';
strct = 'params';

% assign internal value options.Uloc which stores the input set for each
% location
if ~isfield(options,option)
% initialize input set with zeros

    if isfield(options,'inputCompMap')
        warning('params.inputCompMap is redundant!');
    end

    uloc = cell(numComps,1);
    numLoc = length(obj.components{1}.location);
    uloc{1} = cell(numLoc,1);
    
    for i = 1:numLoc
        uloc{1}{i} = zonotope(zeros(obj.numInputs,1));
    end   
    
    options.Uloc = uloc;
    options.inputCompMap = ones(obj.numInputs,1);
    
elseif ~iscell(options.U) 
% same input set for each location    
    
    if isfield(options,'inputCompMap')
        warning('params.inputCompMap is redundant!');
    end

    if dim(options.U) ~= obj.numInputs
        error(printOptionOutOfRange(obj,option,strct));
    end

    uloc = cell(numComps,1);
    numLoc = length(obj.components{1}.location);
    uloc{1} = cell(numLoc,1);
    
    for i = 1:numLoc
        uloc{1}{i} = options.U;
    end   
    
    options.Uloc = uloc;
    options.inputCompMap = ones(obj.numInputs,1);
    
    options = rmfield(options,option);
    
else
% copy input set for each location

    % check if input component map exists and has the correct dimension
    if ~isfield(options,'inputCompMap')
        error(printOptionMissing(obj,'inputCompMap',strct));
    end
    
    if ~all(size(options.inputCompMap) == [obj.numInputs,1]) || ...
       ~all(size(options.inputCompMap) == [1,obj.numInputs]) || ...
       max(options.inputCompMap) > numComps || min(options.inputCompMap) < 1
       
        error(printOptionOutOfRange(obj,'inputCompMap',strct));
    end

    % check if length of cell is identical to the number of components
    if all(size(options.u) ~= [numComps,1]) && ...
       all(size(options.u) ~= [1,numComps])
        error(printOptionOutOfRange(obj,option,strct));
    else
        for i = 1:numComps
           % check if length of cell for each component is identical to the
           % number of locations
           numLoc = length(obj.components{i}.location); 
           if all(size(options.u{i}) ~= [numLoc,1]) && ...
              all(size(options.u{i}) ~= [1,numLoc])
                error(printOptionOutOfRange(obj,option,strct));
           end
           % check if size of input matches inputCompMap
           dims = length(find(options.inputCompMap == i));
           if dims > 0
               for j = 1:numLoc
                    if dim(options.U{i}{j}) ~= dims
                        error(printOptionOutOfRange(obj,option,strct));
                    end
               end
           end
        end
        options.Uloc = options.U;
    end
    
    options = rmfield(options,option);
end

end
%------------- END OF CODE --------------

