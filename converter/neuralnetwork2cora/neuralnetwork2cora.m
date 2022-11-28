function NN = neuralnetwork2cora(input)
% neuralnetwork2cora - import neural networks into CORA
%
% Syntax:  
%    NN = neuralnetwork2cora(obj)
%    NN = neuralnetwork2cora(file)
%
% Inputs:
%    obj - DAGnetwork or dlnetwork objects from MATLABs build-in "Deep
%          Learning Toolbox"
%    file - path to a file storing the neural network (currently only
%           .onnx files are supported
%
% Outputs:
%    NN - object of class neuralNetworkOld
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetworkOld

% Author:       Niklas Kochdumper
% Written:      12-November-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse different types of inputs
    if ischar(input)
        
        % check if the file exists
        if exist(input,'file') ~= 2
            throw(CORAerror('CORA:fileNotFound',input));
        end
        
        % handle different file types
        [~,~,type] = fileparts(input);
       
        if strcmp(type,'.onnx')
            try 
                net = importONNXNetwork(input,'OutputLayerType', ...
                                                            'regression');
                NN = deepLearningToolbox2cora(net.Layers);
            catch ex
                try
                    net = importONNXNetwork(input,'OutputLayerType', ...
                                 'regression','TargetNetwork','dlnetwork');
                    NN = constructFromParams(net.Learnables.Value);
                catch
                    try
                        L = importONNXLayers(input,'OutputLayerType', ...
                                        'regression','ImportWeights',true);
                        NN = deepLearningToolbox2cora(L);
                    catch
                        rethrow(ex);
                    end
                end
            end
        elseif strcmp(type,'.nnet')
            NN = importNNET(input);
        elseif strcmp(type,'.yml')
            NN = importYML(input);
        elseif strcmp(type,'.sherlock')
            NN = importSherlock(input);
        else
            throw(CORAerror('CORA:notSupported',...
                'Files of this type are not supported!')); 
        end
        
    elseif isa(input,'DAGNetwork') || isa(input,'dlnetwork')
        NN = deepLearningToolbox2cora(input.Layers);
    end
end


% Auxiliary Functions -----------------------------------------------------

function NN = deepLearningToolbox2cora(Layers)
% Convert a neural network represented as an object of MATLABs build-in
% "Deep Learning Toolbox" to CORA format

    % initialization
    W = {}; b = {}; actFun = {}; W_ = 1; b_ = 0; openLayer = false;

    % loop over all layers
    for i = 1:length(Layers)

        L = Layers(i); 

        % handle different types of layers
        if isa(L,'nnet.cnn.layer.FullyConnectedLayer')
            W_ = L.Weights*W_; 
            if isscalar(b_)
                b_ = L.Weights*b_*ones(size(L.Weights,2),1) + L.Bias;
            else
                b_ = L.Weights*b_ + L.Bias; 
            end
            openLayer = true;
        elseif isa(L,'nnet.onnx.layer.ElementwiseAffineLayer')
            s = L.Scale; o = squeeze(L.Offset);
            if size(o,1) >= 1 && size(o,2) == 1
               W_ = s*W_; b_ = s*b_ + o; 
            elseif size(o,1) == 1 && size(o,2) > 1
               W_ = s*W_; b_ = s*b_ + o';
            else
                throw(CORAerror('CORA:converterIssue',...
                    ['Conversion failed due to incompatibility with ', ...
                                            '"Elementwise Affine Layer"!']));
            end
            openLayer = true;
        elseif isa(L,'nnet.cnn.layer.ReLULayer')
            W{end+1} = W_; b{end+1} = b_; actFun{end+1} = 'ReLU';
            W_ = 1; b_ = 0; openLayer = false;
        elseif isa(L,'nnet.cnn.layer.LeakyReLULayer')
            alpha = L.Scale;
            W{end+1} = W_; b{end+1} = b_; actFun{end+1} = 'LeakyReLU';
            W_ = 1; b_ = 0; openLayer = false;
        elseif isa(L,'nnet.onnx.layer.SigmoidLayer')
            W{end+1} = W_; b{end+1} = b_; actFun{end+1} = 'sigmoid';
            W_ = 1; b_ = 0; openLayer = false;
        end
    end
    
    % consider a potential last output layer
    if openLayer
       W{end+1} = W_; b{end+1} = b_; actFun{end+1} = 'identity'; 
    end
    
    % convert all entries to double
    for i = 1:length(W)
       W{i} = double(W{i}); b{i} = double(b{i}); 
    end

    % construct neural network object
    NN = neuralNetworkOld(W,b,actFun);
end

function NN = constructFromParams(params)
% construct a neural network from the learnable parameters

    % extract weights and biases from the learnable parameters
    W = {}; b = {};

    for i = 1:length(params)
        if size(params{i},2) == 1
            b{end+1,1} = double(extractdata(params{i}));
        else
            W{end+1,1} = double(extractdata(params{i}))';
        end
    end

    % try to construct a neuralNetworkOld object
    actFun = [repmat({'ReLU'},[length(W)-1,1]); {'identity'}];

    NN = neuralNetworkOld(W,b,actFun);

    warning(['Could not determine activation functions,', ...
                                    ' so ReLUs are used as a default.']);
end

function NN = importNNET(input)
% import a neural network from a .nnet file

    % read text from file
    text = fileread(input);
    lines = strsplit(text,'\n');
    
    % remove comments
    lines = lines(~cellfun(@(x) startsWith(x,'//'),lines));
    
    % get number of layers
    ind = find(lines{1} == ',');
    nrLayers = str2double(lines{1}(1:ind));
    
    % get number of inputs
    ind = find(lines{2} == ',');
    nrInputs = str2double(lines{2}(1:ind));
    
    % get number of neurons in each layer
    nrNeurons = zeros(nrLayers,1);
    text = strtrim(lines{2}(ind+1:end));
    
    for i = 1:nrLayers
        ind = find(text == ',');
        nrNeurons(i) = str2double(text(1:ind-1));
        text = strtrim(text(ind+1:end));
    end
    
    % loop over all layers and read the weights and biases
    W = cell(nrLayers,1); b = cell(nrLayers,1); 
    lines = lines(8:end); inpSize = nrInputs;
    
    for i = 1:nrLayers
        
        % read weights
        W{i} = zeros(nrNeurons(i),inpSize);
        
        for j = 1:nrNeurons(i)
            evalc(['temp = [',lines{j},'];']);
            W{i}(j,:) = temp;
        end
        lines = lines(nrNeurons(i)+1:end);
        
        % read bias
        b{i} = zeros(nrNeurons(i),1);
        
        for j = 1:nrNeurons(i)
            b{i}(j) = str2double(lines{j});
        end
        
        lines = lines(nrNeurons(i)+1:end);
        inpSize = nrNeurons(i);
    end
    
    % construct neural network
    NN = neuralNetworkOld(W,b,'ReLU');
end

function NN = importYML(input)
% import a neural network from a .yml file

    % read text from file
    text = fileread(input);
    lines = strsplit(text,'\n');
    
    % get activation functions
    text = strtrim(lines{1}); actFun = []; cnt = 1; finished = false;
    
    while ~finished
        ind = strfind(text,[num2str(cnt),':']);
        cnt = cnt + 1;
        text = text(ind(1)+2:end);
        ind = strfind(text,',');
        if ~isempty(ind)
            temp = strtrim(text(1:ind-1));
            text = text(ind+1:end);
        else
            temp = strtrim(text(1:end-1));
            finished = true;
        end
        if strcmp(temp,'Sigmoid')
            actFun = [actFun; {'sigmoid'}]; 
        elseif strcmp(temp,'Tanh')
            actFun = [actFun; {'tanh'}]; 
        elseif strcmp(temp,'ReLU')
            actFun = [actFun; {'ReLU'}];
        elseif strcmp(temp,'Linear')
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
       temp = 'temp = ';
       ind = strfind(bias{cnt},'[');
       bias{cnt} = bias{cnt}(ind(1)-1:end);
       while ~startsWith(strtrim(bias{cnt}),[num2str(i+1),':'])
           temp = [temp, strtrim(bias{cnt})];
           cnt = cnt + 1;
       end
       evalc(temp);
       b{i} = temp';
    end
    
    % parse the weights
    W = cell(size(actFun));
    weights = [weights,{[num2str(length(actFun)+1),':']}];
    cnt = 1;
    
    for i = 1:length(W)
       cnt = cnt + 1;
       temp = [];
       ind = strfind(weights{cnt},'[');
       weights{cnt} = weights{cnt}(ind(1)+1:end);
       while ~startsWith(strtrim(weights{cnt}),[num2str(i+1),':'])
           temp = [temp, strtrim(weights{cnt})];
           cnt = cnt + 1;
           if startsWith(strtrim(weights{cnt}),'- [')
              ind = strfind(weights{cnt},'[');
              weights{cnt} = weights{cnt}(ind(1)+1:end);
              temp(end) = ';';
           end
       end
       evalc(['temp = [',temp]);
       W{i} = temp;
    end
    
    % construct neural network
    NN = neuralNetworkOld(W,b,actFun);
end

function NN = importSherlock(input)
% import a neural network from a .sherlock file (neural network format used
% by the sherlock tool for reachability analysis)

    % read text from file
    text = fileread(input);
    lines = strsplit(text,'\n');
    
    % get network properties
    nrInputs = str2double(strtrim(lines{1}));
    nrOutputs = str2double(strtrim(lines{2}));
    hiddenLayers = str2double(strtrim(lines{3}));
    
    % get number of neurons in each layer
    nrNeurons = cell(hiddenLayers + 2,1);
    nrNeurons{1} = nrInputs;
    nrNeurons{end} = nrOutputs;
    
    for i = 1:hiddenLayers
        nrNeurons{i+1} = str2double(strtrim(lines{3+i}));
    end
    
    % initialization
    cnt = 3 + hiddenLayers;
    W = cell(hiddenLayers+1,1);
    b = cell(hiddenLayers+1,1);
    
    % loop over all layers
    for i = 1:length(nrNeurons)-1
        
        % initialization
        temp = zeros(nrNeurons{i+1},nrNeurons{i}+1);
       
        % read data
        for k = 1:nrNeurons{i+1}
            offset = (k-1)*(nrNeurons{i}+1);
            for j = 1:nrNeurons{i}+1
                temp(k,j) = str2double(strtrim(lines{cnt+offset+j}));
            end
        end
        cnt = cnt + (nrNeurons{i}+1)*nrNeurons{i+1};
        
        % get weight matrix and bias vector
        W{i} = temp(:,1:end-1);
        b{i} = temp(:,end);
    end
    
    % construct neural network
    NN = neuralNetworkOld(W,b,'ReLU');
end

%------------- END OF CODE --------------