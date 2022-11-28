classdef neuralNetworkOld
% neuralNetworkOld - class that stores neural networks
%
% Syntax:  
%    obj = neuralNetworkOld(W,b,actFun)
%
% Inputs:
%    W - cell-array storing the weights for all layers
%    b - cell-array storing the bias vectors for all layers
%    actFun - cell-array storing the activation functions for all layers as
%             strings ('ReLU', 'sigmoid', 'tanh', or 'identity') 
%
% Outputs:
%    obj - generated object
%
% Example:
%    W{1} = rand(100,4); b{1} = rand(100,1);
%    W{2} = rand(2,100); b{2} = rand(2,1);
%    nn = neuralNetworkOld(W,b,'ReLU');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neurNetContrSys

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      17-September-2021             
% Last update:  ---
% Last revision:10-August-2022 (TL: renamed to neuralNetworkOld)

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    W = [];                                         % weights
    b = [];                                         % bias vectors
    actFun = [];                                    % activation functions
    nrOfInputs = [];                                % number of inputs
    nrOfOutputs = [];                               % number of outputs
end
    
methods
    
    % class constructor
    function obj = neuralNetworkOld(W,b,actFun)
        
        validActFuns = {'ReLU','sigmoid','tanh','identity'};
        
        % check correctness of input arguments
        if ~iscell(W) || (size(W,1) ~= 1 && size(W,2) ~= 1) 
            throw(CORAerror('CORA:wrongValue','first',...
                "cell-array storing the weights"));
        else
            for i = 1:length(W)-1
                if size(W{i+1},2) ~= size(W{i},1)
                    throw(CORAerror('CORA:wrongValue','first',...
                        "cell-array storing the weights"));
                end
            end
        end
        
        if ~iscell(b) || ~all(size(b) == size(W))
            throw(CORAerror('CORA:wrongValue','second',...
                "cell-array storing the bias and match the size of 'W'"));
        else
            for i = 1:length(b)
               if size(b{i},1) ~= size(W{i},1) || size(b{i},2) ~= 1
                   throw(CORAerror('CORA:wrongValue','second',...
                       "cell-array storing the bias and match the size of 'W'"));
               end
            end
        end
        
        if iscell(actFun)
            if ~all(size(actFun) == size(W))
                throw(CORAerror('CORA:wrongValue','third',...
                    "'ReLU', 'sigmoid', 'tanh', or 'identity'"));
            else
               for i = 1:length(actFun)
                  if ~ismember(actFun{i},validActFuns)
                      throw(CORAerror('CORA:wrongValue','third',...
                          "'ReLU', 'sigmoid', 'tanh', or 'identity'"));
                  end
               end
            end
        elseif ischar(actFun)
            if ~ismember(actFun,validActFuns)
                throw(CORAerror('CORA:wrongValue','third',...
                    "'ReLU', 'sigmoid', 'tanh', or 'identity'"));
            else
               actFun = repmat({actFun},size(W)); 
            end
        else
            throw(CORAerror('CORA:wrongValue','third',...
                "'ReLU', 'sigmoid', 'tanh', or 'identity'"));
        end
        
        % assign object properties
        obj.W = W; obj.b = b; obj.actFun = actFun;
        obj.nrOfInputs = size(W{1},2); obj.nrOfOutputs = size(W{end},1);
    end   

    display(obj)
end

methods (Static = true)
    obj = generateRandom(varargin)
end
end

%------------- END OF CODE --------------
