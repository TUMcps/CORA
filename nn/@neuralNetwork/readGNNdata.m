function data = readGNNdata(file_path, varargin)
% readGNNdata - reads and converts a gnn data json file
%
% Syntax:
%    obj = neuralNetwork.readGNNdata(file_path,verbose)
%
% Inputs:
%    file_path - path to data file
%    verbose - bool if information should be displayed
%
% Outputs:
%    data - data table
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/readGNNnetwork

% Authors:       Tobias Ladner
% Written:       29-January-2024
% Last update:   --- 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,2)
verbose = setDefaultValues({false}, varargin);
inputArgsCheck({ ...
    {file_path, 'att', {'char', 'string'}}; ...
    {verbose, 'att', 'logical'}; ...
})
if ~isfile(file_path)
    throw(CORAerror('CORA:fileNotFound', file_path));
end

% prepare and setup the model file containing weights and biases
fid = fopen(file_path);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
jsondata = jsondecode(str);

% first row are column headers
headers = jsondata{1};

% read data (prealocate?)
for i=2:length(jsondata)
    data_i = jsondata{i};
    for j=1:length(headers)
        data(i-1,1).(headers{j}) = data_i{j};
    end
end

% convert to table
data = struct2table(data,'AsArray',true);

end


% Auxiliary functions -----------------------------------------------------

function str = aux_getJSONLayerName(n, type)
    str = ['l_', num2str(n), '_', type];
end

% ------------------------------ END OF CODE ------------------------------
