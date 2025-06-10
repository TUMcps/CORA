function nnStruct = exportAsStruct(nn)
% exportAsStruct - exports this network as a struct
%
% Syntax:
%    nnStruct = nn.exportAsStruct()
%
% Inputs:
%    nn - neuralNetwork
%
% Outputs:
%    nnStruct - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/importFromStruct

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(nn.layers)
    nnStruct = [];
else
    % iterate through all layers
    for i = 1:length(nn.layers)
        % export as struct
        nnStruct(i,1) = exportAsStruct(nn.layers{i});
    end
end

end

% ------------------------------ END OF CODE ------------------------------
