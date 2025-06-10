function jsonstr = exportAsJSON(nn, file_path)
% exportAsJSON - exports this network in JSON format
%
% Syntax:
%    nn = nn.exportAsJSON()
%    nn = nn.exportAsJSON(file_path)
%
% Inputs:
%    nn - neuralNetwork
%    file_path - (optional) str, file path to store the network
%
% Outputs:
%    json - json string
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/importFromJSON

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to struct
nnStruct = exportAsStruct(nn);

% convert to json
jsonstr = jsonencode(nnStruct);

% save to file if path given
if nargin > 1 && ~isempty(file_path)
    try
        fid = fopen(file_path,'w');
        fprintf(fid,jsonstr);
        fclose(fid);
    catch ME
        fprintf('Unable to write to file: %s\n', file_path);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
