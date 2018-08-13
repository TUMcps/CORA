function res = mergeInputVector(obj,inputMap,loc,uGlob,uCompLoc)
% mergeInputVector - construct the input vector for the current location
%                    from the global input and the inputs for the
%                    subcomponents
%
% Syntax:  
%    res = mergeInputVector(obj,inputMap,loc,uGlob,uCompLoc)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    inputMap - precomputed matrix. The first row specifies the indes of
%               the subcomponent, the second row the index of the input for
%               this subcomponents
%    loc - id of the current location
%    uGlob - input vector for the globally defined inputs
%    uCompLoc - input vector for all subcomponents
%
% Outputs:
%    res - resulting vector of inputs
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      09-July-2018  
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    res = zeros(size(inputMap,1),1);

    % loop over all inputs
    for i = 1:obj.numInputs
        
        % read precomputed inputMap
        comp = inputMap(i,1);
        dim = inputMap(i,2);
        
        if comp ~= 0
            % lookup component input for current location
            c_loc = loc{comp};
            res(i) = uCompLoc{comp}{c_loc}(dim);
        else
            % Component 0 indicates global input
            res(i) = uGlob(dim);
        end
    end
end

%------------- END OF CODE --------------