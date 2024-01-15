function [forbiddenSpecs,spec_mapping] = ...
    parseForbidden(configLine,state_names,component_names,location_names)
% parseForbidden - Constructs CORA specifications from forbidden states
%    given in a SpaceEx-ConfigFile
%
% Syntax:
%    [forbiddenSpecs,spec_mapping] = ...
%       parseForbidden(configLine,state_names,component_names,location_names)
%
% Inputs:
%    configLine  - Line (String) from the config-file containig the
%                  forbidden states
%    state_names - Struct array containing structs with field "name"
%                  specifying the names of the state variables
%    component_names - String array containing the names of all automaton
%                      components
%    location_names - Cell array where entry i contains a list of all
%                     location names of automaton component i
%
%    NOTE: All input parameters containing names must be in order!
%
% Outputs:
%    forbiddenSpecs - CORA specification object defining forbidden areas
%                     within the state-space
%    spec_mapping - cell array specifying which specifications are relevant
%                   to which locations of which subcomponent
%                   e.g.: spec_mapping =  {{[1,2]; [1]}, {[2]; [3]}} means
%                     subcomponent 1: specs 1,2 for location 1, spec 1 for
%                     location 2
%                     subcomponent 2: spec 2 for location 1, spec 3 for
%                     location 2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Maximilian Perschl
% Written:       08-September-2021
% Last update:   24-September-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% cut irrelevant parts off the string
configLine_split = split(configLine,'"');

% if the forbidden values are entirely specified in a single line, there are
% two occurences of " -> the size of the returned cell array is 3
% if there are multiple lines, there is only one " in the line -> the size
% of the returned cell array is 2
if length(configLine_split) == 2
    % get new lines until the "forbidden" segment is complete
    while count(configLine,'"') < 2
        configLine = strcat(configLine,fgetl(fileDesc));
    end
    configLine_split = split(configLine,'"');
end

% If the split array has more than 3 elements, there is a syntax error in
% the cfg file
if length(configLine_split) ~= 3
    throw(CORAerror('CORA:converterIssue',...
        'Syntax error in forbidden specification.'));
end
% Assume entirety of forbidden specifications is
% contained within "configLine_split(2)"
configLine = configLine_split(2);

% We start with an empty specification (empty unsafe set), and add onto it iteratively
forbiddenSpecs = specification(halfspace.empty(length(state_names)),'unsafeSet');

% We assume the specification of "forbidden" to be in conjuntive normal
% form, therefore we can get a set of specifications we can add on to each
% other by splitting the configuration line by "|"(= logical OR)
forbiddenTerms = split(configLine,"|");

% For each forbidden term, transform the string to a symbolic equation of
% the right format for the "eq2set" function, then obtain the set and
% generate the specifcation. Then add the specifcation to the others
% If a specification is tied to a specific location of a specific
% component, an entry in spec_mapping is generated
spec_mapping = {};
% By default, all specifications are assumed to be relevant for all
% locations and components. If at least one specification is only relevant
% for certain locations or components, the spec_mapping object is non-empty
% and contains entries for all automata and locations
no_mapping = 1;
% Iterating variable for creating spec_mapping entries
% Starts at 1 because the resulting specification object will have a dummy
% spec at index 1
spec_index = 1;


% Iteration over terms
for i = 1:length(forbiddenTerms)
    % increment specification counter
    spec_index = spec_index + 1;
    % split among "&" operators, then sort subterms into equalities and
    % non-equalities and add entries concerning location specifications in
    % the mapping
    subterms = split(forbiddenTerms(i),'&');
    equalities = [];
    inequalities = [];
    % flag for filtering out specifications which solely specify a
    % forbidden location
    forbidden_location_term = 0;
    % flag for determining if a location restriction was given
    loc_restriction = 0;
    for j = 1:length(subterms)
        % filter out location constraints
        if contains(subterms(j),"loc(")
            % first time this happens, initialize the spec mapping
            if no_mapping
                for k = 1:length(component_names)
                    spec_mapping{k} = {};
                    for l = 1:length(location_names{k})
                        spec_mapping{k}{l} = [];
                    end
                end
            end
            % set flaggs for location-specific specifications
            no_mapping = 0;
            loc_restriction = 1;
            
            % prune for arguments
            arg_strings = split(subterms(j),")=");
            % get component name
            component_name = extractAfter(arg_strings(1),"(");
            % get location name
            location_name = arg_strings(2);
            
            % check for high-level component-specification (currently not
            % supported)
            if contains(component_name,".")
                throw(CORAerror('CORA:converterIssue',...
                    'High-level specification of components is currently not supported.'));
            end
            
            % get index of component
            component_index = find(component_names==component_name);
            % get index of location
            location_index = find(location_names{component_index} == location_name);
            
            % check if the constraint is solely reaching a location
            if length(subterms) == 1
                % Set flag for specification generation further down
                forbidden_location_term = 1;
            end
            
            spec_mapping{component_index}{location_index} = [spec_mapping{component_index}{location_index} spec_index];
            continue;
        end
        % generate symbolic equation, order to 
        % "left-hand-side *operator* 0"
        % example: x1 <= x2 --> x1 - x2 <= 0
        currEqu = str2sym(subterms(j));
        currEqu = currEqu - rhs(currEqu);
        % determine type of equation 
        % only add left-hand side (format required by eq2set)
        if hasSymType(currEqu,'eq')
            equalities = [equalities lhs(currEqu)];
        else
            inequalities = [inequalities lhs(currEqu)];
        end            
    end
    
    % Create different new specifications, depending on if the current term
    % constrains either locations or states
    if ~forbidden_location_term
        % get resulting CORA set representation
        % (spaceEx models as a whole have no outputs, so they are irrelevant
        % for the forbidden section and are ignored here)
        set = eq2set(inequalities,equalities,state_names,[]);
        % expand the specifications
        forbiddenSpecs = add(forbiddenSpecs,specification(set,'unsafeSet'));
    else
        % Create a specification that is never satisfied, active
        % for the relevant location only
        forbiddenSpecs = add(forbiddenSpecs,...
            specification(interval.empty(length(state_names)),'invariant'));
    end
    
    % if no restriction on the location has been made, the specification
    % holds for all locations and components
    if ~loc_restriction && ~no_mapping
        spec_mapping = aux_createMapEntryForAll(spec_mapping,spec_index,component_names,location_names);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function resulting_map = aux_createMapEntryForAll(prev_map,constr_index,component_names,location_names)
% create an entry for a specification that is relevant to all locations and
% components
    n_components = length(component_names);
    for i = 1:n_components
        n_locations(i) = length(location_names{i});
    end
    for i = 1:n_components
        for j = 1:n_locations(i)
            prev_map{i}{j} = [prev_map{i,j} constr_index];
        end
    end
    resulting_map = prev_map;
end

% ------------------------------ END OF CODE ------------------------------
