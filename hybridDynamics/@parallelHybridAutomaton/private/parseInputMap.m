function res = parseInputMap(obj,options)
% parseInputMap - Calculate an input map that assignes each system input a
%                 subcomponent and a corresponding index of the
%                 subcomponent input
%
% Syntax:  
%    res = parseInputMap(obj,options)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    options - struct defining algorithm options
%
% Outputs:
%    res - resulting subcomponent-to-input map
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      09-July-2018  
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % get the number of component in the parallel hybrid automaton
    numComp = length(obj.components);

    % Explanation of component-wise input:
    % Input domain of PHA: U = [u1 u2 u3 u4 u5]
    % inputCompMap is a mapping: [1,numInputs] -> [0,numComponents]
    % (example: inputCompMap(1:5) = [1  1  3  2  3 ], for a 3-component PHA)
    % => Input domain of Component c: U_c = [u_i | m(i) == c]
    % => Additional global input domain: U_glob = [u_i | m(i) == 0]
    % (example: U_1 = [u1 u2]; U_2 = [u4]; U_3 = [u3 u5]; U_glob = [])
    % => for all component-locations l: uCompLoc{c}{l} must specify U_c
    % => uGlob must specify U_glob
    % (example: uCompLoc{1}{*} ∈ R^2; uCompLoc{3}{*} ∈ R^2; uCompLoc{2}{*} ∈ R
    %           uGlob ∈ R^0)
    
    if isfield(options,'inputCompMap')
        
        % map global input to component inputs
        res = zeros(obj.numInputs,2);
        compNumInputs = zeros(numComp,1);
        globNumInputs = 0;

        for i = 1:length(options.inputCompMap)
            comp = options.inputCompMap(i);

            if comp ~= 0
                
                % get local index of input
                localIdx = compNumInputs(comp) + 1;
                
                % update local index counter
                compNumInputs(comp) = localIdx;
            else
                localIdx = globNumInputs + 1;
                globNumInputs = localIdx;
            end
            
            % store component + index of input
            res(i,:) = [comp,localIdx];
        end
    else
        error("'options.inputCompMap' missing! You can pass zeros(numInputs,1) to indicate global input.");
    end
end

%------------- END OF CODE --------------