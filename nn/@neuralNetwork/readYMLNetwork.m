function obj = readYMLNetwork(file_path)
% readYMLNetwork - reads and converts a network saved in yml format
%
% Syntax:
%    res = neuralNetwork.readYMLNetwork(file_path)
%
% Inputs:
%    file_path - path to file
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/getFromCellArray

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       12-November-2021
% Last update:   30-March-2022
%                30-November-2022 (removed neuralNetworkOld)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read text from file
text = fileread(file_path);
lines = strsplit(text,'\n');

% get activation functions
text = strtrim(lines{1}); actFun = []; cnt = 1; finished = false;

while ~finished
    % read next activation
    ind = strfind(text,[num2str(cnt),':']);
    cnt = cnt + 1;
    text = text(ind(1)+2:end);
    ind = strfind(text,',');
    if ~isempty(ind)
        actFunName = strtrim(text(1:ind(1)-1));
        text = text(ind(1)+1:end);
    else
        actFunName = strtrim(text(1:end-1));
        finished = true;
    end
    % add to previous activations
    if strcmp(actFunName,'Sigmoid')
        actFun = [actFun; {'sigmoid'}];
    elseif strcmp(actFunName,'Tanh')
        actFun = [actFun; {'tanh'}];
    elseif strcmp(actFunName,'ReLU')
        actFun = [actFun; {'ReLU'}];
    elseif strcmp(actFunName,'Linear')
        % no activation -> add identity
        actFun = [actFun; {'identity'}];
    else
        throw(CORAerror('CORA:converterIssue'));
    end
end

% split lines into offsets and weights
bias = []; weights = [];

for i = 3:length(lines)
   if startsWith(lines{i},'weights')
      bias = lines(3:i-1);
      weights = lines(i+1:end);
   end
end

if isempty(bias) || isempty(weights)
    throw(CORAerror('CORA:converterIssue'));
end

% parse the bias 
b = cell(size(actFun));
bias = [bias,{[num2str(length(actFun)+1),':']}];
cnt = 1;

for i = 1:length(b)
   biasStr = 'biasData = ';
   ind = strfind(bias{cnt},'[');
   bias{cnt} = bias{cnt}(ind(1)-1:end);
   while ~startsWith(strtrim(bias{cnt}),[num2str(i+1),':'])
       biasStr = [biasStr, strtrim(bias{cnt})];
       cnt = cnt + 1;
   end
   evalc(biasStr);
   b{i} = biasData';
end

% parse the weights
W = cell(size(actFun));
weights = [weights,{[num2str(length(actFun)+1),':']}];
cnt = 1;

for i = 1:length(W)
   cnt = cnt + 1;
   weightStr = [];
   ind = strfind(weights{cnt},'[');
   weights{cnt} = weights{cnt}(ind(1)+1:end);
   % parse string
   while ~startsWith(strtrim(weights{cnt}),[num2str(i+1),':'])
       weightStr = [weightStr, strtrim(weights{cnt})];
       cnt = cnt + 1;
       if startsWith(strtrim(weights{cnt}),'- [')
          ind = strfind(weights{cnt},'[');
          weights{cnt} = weights{cnt}(ind(1)+1:end);
          weightStr(end) = ';';
       end
   end
   evalc(['weightData = [',weightStr]);
   W{i} = weightData;
end

% construct neural network
obj = neuralNetwork.getFromCellArray(W,b,actFun);

% ------------------------------ END OF CODE ------------------------------
