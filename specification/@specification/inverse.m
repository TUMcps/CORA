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
    
    % convert polytopes with a single inequality constraints to halfspaces
    spec = aux_convertToHalfspaces(spec);
end


% Auxiliary functions -----------------------------------------------------

function spec = aux_convertToHalfspaces(spec)
% convert all polytopes in a specification with a single inequality
% constraint to a halfspace
    for i = 1:length(spec)
        if isa(spec(i,1).set,'polytope') && size(spec(i,1).set.A,1) == 1
            spec(i,1).set = halfspace(spec(i,1).set.A,spec(i,1).set.b); 
        end

    end
end

% ------------------------------ END OF CODE ------------------------------
