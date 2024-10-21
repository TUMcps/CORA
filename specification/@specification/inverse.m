function spec = inverse(spec)
% inverse - inverts a specification object
%
% Syntax:
%    spec = inverse(spec)
%
% Inputs:
%    spec - specification object
%
% Outputs:
%    spec - inverted specification object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper
% Written:       26-November-2021             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% divide into safe and unsafe sets
safe = []; unsafe = []; 
for i = 1:length(spec)
    
    if ~representsa_(spec(i).time,'emptySet',eps)
        throw(CORAerror('CORA:notSupported',...
            'Computing the inverse is not yet supported for timed specifications!'));
    end
    
    if strcmp(spec(i).type,'safeSet')
        safe = [safe,i]; 
    elseif strcmp(spec(i).type,'unsafeSet')
        unsafe = [unsafe,i];
    else
        throw(CORAerror('CORA:notSupported',...
            "Computing the inverse is not yet supported for " + ...
            "other types than 'safeSet' and 'unsafeSet'."));
    end
end

if ~isempty(safe) && ~isempty(unsafe)
    throw(CORAerror('CORA:notSupported',...
            "Computing the inverse is not yet supported for " + ...
            "mixed types 'safeSet' and 'unsafeSet'."));
elseif ~isempty(safe)
    
    % combine all safe set specifications to a single unsafe set
    set = polytope(spec(safe(1)).set);
    
    for i = 2:length(safe)
        set = set & polytope(spec(safe(i)).set);
    end
    
    spec = specification(set,'unsafeSet');
   
else
    
    % use different conversion depending on the number of sets 
    if length(unsafe) == 1
        spec = specification(spec(unsafe(1)).set,'safeSet');
    else
        % get all unsafe sets
        sets = cell(length(unsafe),1);
        for i = 1:length(unsafe)
            sets{i} = spec(unsafe(i)).set; 
        end
        
        % convert union of safe sets to an equivalent union of unsafe
        % set representation
        sets = safeSet2unsafeSet(sets);
        spec = [];
        for i = 1:length(sets)
            spec = add(spec,specification(sets{i},'unsafeSet'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
