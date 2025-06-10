function [idxFreedFeats,featOrder,timesPerFeat] = explain(nn, x, label, epsilon, varargin)
% explain - computes a minimal abductive explanation
%
% Syntax:
%    idxFreedFeats = explain(nn, x, label, epsilon, varargin)
%
% Inputs:
%    nn - neuralNetwork
%    x - numeric
%    label - numeric, label class
%    epsilon - numeric, noise radius
%    varargin - name-value pair
%       <'Verbose',verbose> - logical, verbose output
%       <'Method',method> - str, methods: 'standard', 'abstract+refine'
%       <'FeatOrder',featOrderMethod> - str, method to determine
%           order to process features: 'sensitivity', 'in-order'
%       <'InputSize',inputSize> - numeric, size of the input (H,W,C)
%       <'RefineSteps',refineSteps> - numeric vector
%       <'BucketType',bucketType> - char, for reduction: 'static', 'dynamic'
%       <'OutputThreshold',delta> - numeric, threshold for regression tasks
%       <'Timeout',timeout> - numeric, timeout for explanation computation
%
% Outputs:
%    idxFreedFeats - numeric, indices of freed features
%       for inputs with multiple channels, 
%       all channels are freed simultaneously
%    featOrder - order in which feature were processed
%    timesPerFeat - numeric, time to free each feature
%
% Reference:
%    [1] Bassan et al., "Explaining Fast and Slow: Abstraction and 
%        Refinement of Provable Explanations". ICML. 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       26-May-2024
% Last update:   24-June-2024 (TL, generalized function)
%                16-July-2024 (TL, improved plotting of explanation)
%                15-November-2024 (TL, added option for regression tasks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse inputs
narginchk(4,inf)
[nn,x,label,epsilon,verbose,method,featOrderMethod,refineMethod,inputSize,refineSteps,bucketType,delta,timeout] = aux_parseInputs(nn,x,label,epsilon,varargin{:});
aux_checkInputs(nn,x,label,epsilon,verbose,method,featOrderMethod,refineMethod,inputSize,refineSteps,bucketType,delta,timeout)

% init plotting
aux_initPlotting(x,verbose,inputSize)

% have logic difference as output for easier specification check 
% (for classification tasks: n_out > 1, otherwise regression)
if nn.neurons_out > 1
    nn = aux_appendLogitDifferenceLayer(nn,label);
end

% get processing order and init freed features
if isnumeric(featOrderMethod)
    featOrder = featOrderMethod;
else
    featOrder = getInputNeuronOrder(nn,featOrderMethod,x,inputSize);
end
idxFreedFeats = [];
timesPerFeat = nan(1,numel(featOrder));

% check if all features can be freed
timerVal = tic;
if verbose
    disp('Checking if all features can be freed..')
end
isVerified = aux_checkVerifiability(nn,x,epsilon,featOrder,inputSize,delta);
t = toc(timerVal);
if verbose
    fprintf('Time to compute initial check: %.4f.\n',t);
end
if isVerified
    disp('All features can be freed.')
    idxFreedFeats = featOrder;
    aux_updatePlotting(x,verbose,idxFreedFeats,inputSize)
    timesPerFeat = zeros(1,numel(featOrder));
    return
end

% init abstract network (used for some methods)
nn_abs = [];

% iteratively free features
for i=1:numel(featOrder)
    if verbose
        fprintf('---\n')
        fprintf('Freeing input feature: %i/%i (%i):\n',i,numel(featOrder),featOrder(i));
    end

    % temporarily add feature i to freed features
    idxFreedFeats_i = [idxFreedFeats,featOrder(i)];

    % run verification
    timerVal = tic;
    [isVerified,nn_abs] = aux_runVerification(nn,nn_abs,x,label,epsilon,verbose,method,refineMethod,idxFreedFeats_i,inputSize,refineSteps,bucketType,delta);
    timesPerFeat(i) = toc(timerVal);
    if verbose
        fprintf('Elapsed time: %.4f.\n',timesPerFeat(i));
    end


    if isVerified % verified
        % permanently add feature i to freed features
        idxFreedFeats = idxFreedFeats_i;

        % verbose output
        if verbose
            fprintf('Input feature %i was freed.\n',featOrder(i))
            
            % update plot
            aux_updatePlotting(x,verbose,idxFreedFeats,inputSize)
        end

    else % could not verify freed features
        if verbose
            fprintf('Input feature %i cannot be freed.\n',featOrder(i))
        end
    end

    % check if freeing next feature would exceed timeout
    if (i+1)*mean(timesPerFeat(1:i)) > timeout
        disp('Timeout!')
        break
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [nn,x,label,epsilon,verbose,method,featOrderMethod,refineMethod,inputSize,refineSteps,bucketType,delta,timeout] = aux_parseInputs(nn,x,label,epsilon,varargin)
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Verbose','Method','FeatOrder','RefineMethod','InputSize','RefinementSteps','BucketType','OutputThreshold','Timeout'});
    % read arguments
    [NVpairs,verbose] = readNameValuePair(NVpairs,'Verbose',{},false);
    [NVpairs,method] = readNameValuePair(NVpairs,'Method',{},'standard');
    [NVpairs,featOrderMethod] = readNameValuePair(NVpairs,'FeatOrder',{},'sensitivity');
    [NVpairs,refineMethod] = readNameValuePair(NVpairs,'RefineMethod',{},'all');
    [NVpairs,inputSize] = readNameValuePair(NVpairs,'InputSize',{},'all');
    [NVpairs,refineSteps] = readNameValuePair(NVpairs,'RefinementSteps',{},0.1:0.1:1);
    [NVpairs,bucketType] = readNameValuePair(NVpairs,'BucketType',{},'static');
    [NVpairs,delta] = readNameValuePair(NVpairs,'OutputThreshold',{},0);
    [NVpairs,timeout] = readNameValuePair(NVpairs,'Timeout',{},Inf);
end

function aux_checkInputs(nn,x,label,epsilon,verbose,method,featOrderMethod,refineMethod,inputSize,refineSteps,bucketType,delta,timeout)

    % validate name-value pairs
    if CHECKS_ENABLED
        % validate ordered input
        inputArgsCheck({ ...
            {nn, 'att', 'neuralNetwork'}; ...
            {x, 'att', {'numeric'}}; ...
            {label, 'att', {'numeric','scalar'}}; ...
            {epsilon, 'att', {'numeric','scalar','nonnegative'}}; ...
        });

        % verbose
        if ~isscalar(verbose) || ~islogical(verbose)
            throw(CORAerror('CORA:wrongValue','name-value pair ''Verbose''','logical, scalar'))
        end

        % method
        methods = {'standard','abstract+refine'};
        if ~ismember(method,methods)
            throw(CORAerror('CORA:wrongValue','name-value pair ''Method''',methods))
        end

        % refineMethod
        refineMethods = {'all','sensitivity','rand'};
        if ~ismember(refineMethod,refineMethods)
            throw(CORAerror('CORA:wrongValue','name-value pair ''RefineMethod''',refineMethods))
        end

        % inputSize
        if ~isnumeric(inputSize) || ~isvector(inputSize) || numel(inputSize) ~= 3
            throw(CORAerror('CORA:wrongValue','name-value pair ''InputSize''','has to be a three-dimensional numeric vector'))
        end

        % refineSteps
        if ~isnumeric(refineSteps) || isempty(refineSteps) || ~isvector(refineSteps)
            throw(CORAerror('CORA:wrongValue','name-value pair ''RefineSteps''','has to be a numeric vector'))
        end

        % bucketType
        bucketTypes = {'static','dynamic'};
        if ~ismember(bucketType,bucketTypes)
            throw(CORAerror('CORA:wrongValue','name-value pair ''BucketType''',bucketType))
        end

        % output threshold delta
        if ~isnumeric(delta) || ~isscalar(delta) || ~all(delta >= 0,"all")
            throw(CORAerror('CORA:wrongValue','name-value pair ''OutputThreshold''','positive numeric scalar'))
        end

        % output threshold delta
        if ~isnumeric(timeout) || ~isscalar(timeout) || ~all(timeout >= 0,"all")
            throw(CORAerror('CORA:wrongValue','name-value pair ''Timeout''','positive numeric scalar'))
        end
    end
end

function isVerified = aux_checkVerifiability(nn,x,epsilon,featOrder,inputSize,delta)
    X = aux_initInputSet(x,epsilon,featOrder,inputSize);
    Y = nn.evaluate(X);
    isVerified = aux_checkSpecs(Y,delta);
end

function X = aux_initInputSet(x,epsilon,idxFree,inputSize)
    X = interval(x);
    if isempty(inputSize)
        % allow any value for free features
        X(idxFree) = X(idxFree) + interval(-epsilon,epsilon); % pixel +/- epsilon
    else
        % allow any value for free pixels
        X = reshape(X,[inputSize(1)*inputSize(2),inputSize(3)]);
        X(idxFree,:) = X(idxFree,:) + interval(-epsilon,epsilon); % pixel +/- epsilon
        X = reshape(X,[],1);
    end
    X = zonotope(X); % convert to zonotope
end

function x = aux_normalizeInputForPlotting(x)
    % normalize image for plotting
    if min(x) < 0
        % [-1,1] -> [0,1]
        x = x*0.5 + 0.5;
    end

    if max(x) > 1
        % [0,255] -> [0,1]
        x = x / max(x);
    end
end

function aux_initPlotting(x,verbose,inputSize)
    % init plotting
    if verbose
        % normalize image for plotting
        x = aux_normalizeInputForPlotting(x);

        % plot
        figure; subplot(1,2,1); hold on
        imshow(reshape(x,inputSize),'InitialMagnification','fit')
        % aux_repositionOverlay(true) % to have the same scaling..
        title('Image')
        subplot(1,2,2); hold on;
        title('Explanation')
        imshow(reshape(x,inputSize),'InitialMagnification','fit')
        aux_repositionOverlay(true)
        drawnow
    end
end

function aux_repositionOverlay(initNewOverlay,x,alpha)
    % init new overlay (underneath all other overlays)
    ax = gca;
    if initNewOverlay
        % init new overlay
        h = imshow(zeros(1),'InitialMagnification','fit');
        set(h,'AlphaData',zeros(1));
    else
        % delete old overlay
        chH = ax.Children;
        delete(chH(end-1));

        % draw new overlay
        h = imshow(x,'InitialMagnification','fit');
        set(h,'AlphaData',alpha)
    end

    % reposition overlay to be underneath all other overlays
    chH = ax.Children;
    chH = [chH(2:end-1);chH(1);chH(end)];
    set(ax,'Children',chH);
end

function aux_updatePlotting(x,verbose,idxFree,inputSize)
    % update plot showing all freed features
    if verbose     
        % plot

        alsoPlotOrgImage = false;
        if alsoPlotOrgImage
            % normalize image for plotting
            x = aux_normalizeInputForPlotting(x);
    
            % add color to grayscale images
            if inputSize(3) == 1
                x = reshape(x,inputSize);
                x = repmat(x,[1 1 3]);
            end 
            % paint all pixels
            alpha = ones(inputSize(1:2));
        else
            % use white canvas
            x = ones([inputSize(1:2) 3]);
            % but only paint respective freed pixels
            alpha = zeros(inputSize(1:2));
            alpha(idxFree) = 1;
        end

        % find color color
        colorOrder = parula(12);
        ax = gca();
        color = colorOrder(numel(ax.Children)-1,:);

        % color freed pixels
        x = reshape(x,[],3);
        x(idxFree,1) = color(1); % r
        x(idxFree,2) = color(2); % g
        x(idxFree,3) = color(3); % b
        x = reshape(x,[inputSize(1:2) 3]);
        
        % update image
        aux_repositionOverlay(false,x,alpha)
        drawnow
    end
end

function nn = aux_appendLogitDifferenceLayer(nn,label)
    % appends a layer of the neural network to output the logit difference
    W_argmax = eye(nn.neurons_out);
    W_argmax(:,label) = W_argmax(:,label) - 1;
    nn.layers{end+1} = nnLinearLayer(W_argmax,0,'logit difference');

    % get neural network in normal form
    nn = nn.getNormalForm();
end

function [nn_abs,Y] = aux_initAbstractNetwork(nn,nn_abs,X,method,verbose,refineSteps,bucketType)
    % check last reduction
    if isempty(nn_abs)
        % default
        reductionRate = refineSteps(1);
    else
        reductionRate = max(nn_abs.reductionRate);
    end
        

    % create an initial abstract network
    switch method
        case 'abstract+refine'
            % check if network can be refined
            if reductionRate == 1
                % use original network
                if verbose
                    fprintf('Using original network..\n')
                end
                nn_abs = nn;
                Y = nn.evaluate(X);
            else
                % compute abstract network
                if verbose
                    fprintf('Using abstract network (%.0f%% of neurons)..\n', reductionRate*100)
                end
                [nn_abs,Y] = nn.computeReducedNetwork(X,"ReductionRate",reductionRate,'Verbose',false,'BucketType',bucketType);
            end

        otherwise
            % no special abstract network needed
            nn_abs = [];
            Y = [];
    end
end

function [res,nn_abs,Y] = aux_refineAbstractNetwork(nn,nn_abs,X,verbose,method,refineMethod,refineSteps,bucketType)

    % assume refinement was not successfull
    res = false;
    Y = [];

    switch method
        case 'abstract+refine'

            % compute refined abstract network

            % choose method
            switch refineMethod
                case 'all'
                    % check if network can be refine
                    reductionRate_old = max(nn_abs.reductionRate);
                    if reductionRate_old < max(refineSteps)
                        reductionRate_new = refineSteps(find(refineSteps > reductionRate_old,1));
                        if reductionRate_new == 1
                            % use original network
                            if verbose
                                disp('Returning original network..\n')
                            end
                            nn_abs = nn;
                            Y = nn_abs.evaluate(X);
                        else
                            % compute abstract network
                            if verbose
                                fprintf('Refining abstract network (%.0f%% of neurons)..\n', reductionRate_new*100)
                            end
                            [nn_abs,Y] = nn.computeReducedNetwork(X,"ReductionRate",reductionRate_new,'Verbose',false,'BucketType',bucketType);
                        end
                        res = true;
                    end

                case 'sensitivity'
                    throw(CORAerror('CORA:notSupported','Refinement using sensitivity.'))

                case 'rand'
                    nn_abs = nn.computeReducedNetwork(X,"ReductionRate",rand(1),'Verbose',false,'BucketType',bucketType);
                    res = true;

            end


        otherwise
            % no refinement implemented/needed
    end
end

function [isVerified,nn_abs] = aux_runVerification(nn,nn_abs,x,label,epsilon,verbose,method,refineMethod,idxFreedFeats_i,inputSize,refineSteps,bucketType,delta)

    % construct input set
    X = aux_initInputSet(x,epsilon,idxFreedFeats_i,inputSize);
    
    % run verification
    switch method

        case 'standard'
            if verbose
                fprintf('Using original network..\n')
            end

            % compute output set
            Y = nn.evaluate(X);    

            % check specification (classification)
            isVerified = aux_checkSpecs(Y,delta);

        case 'abstract+refine'
            
            isVerified = [];

            % init abstract network (also obtains Y)
            [nn_abs,Y] = aux_initAbstractNetwork(nn,nn_abs,X,method,verbose,refineSteps,bucketType);

            while isempty(isVerified)

                % 1.) check specs
                isVerified = aux_checkSpecs(Y,delta);
                
                if isVerified
                    % 2.a) current set of features can be freed
                    return

                else
                    % 2.b) can we refine the abstract network?

                    % try to find counter example
                    xs_ = aux_findCounterExample(nn,nn_abs,x,label,epsilon,delta);

                    if isempty(xs_)
                        % 3.a) no counter example found, refine abstraction (also obtains Y)
                        [resRefine,nn_abs,Y] = aux_refineAbstractNetwork(nn,nn_abs,X,verbose,method,refineMethod,refineSteps,bucketType);

                        if resRefine
                            % 4.a) refinement successful
                            % reset variable to continue in loop
                            isVerified = [];
                            aux_repositionOverlay(true);
                        else
                            % 4.b) network could not be refined 
                            % stop loop, feature i cannot be freed
                            isVerified = false;
                        end

                    else
                        % 3.b) falsified on original network
                        % no necessity to refine abstract network
                        % feature i cannot be freed
                        disp('Counterexample')
                        isVerified = false;
                    end
                end
            end
    
        otherwise
            % should not happen as it was checked before
            throw(CORAerror('CORA:specialError',sprintf('Unknown method: %s', method)))

    end

end

function [isVerified,idxSpurious] = aux_checkSpecs(Y,delta)
    % check specification, works for set and numeric input

    if isnumeric(Y)
        % numeric input

        if isscalar(Y)
            % regression 
            % (unable to determine deviation here)
            isVerified = true;
            idxSpurious = [];

        else
            % classification
            ys = Y;
            us = max(ys,[],1);
            idxSpurious = find(max(us,[],1) > 0);
            isVerified = isempty(idxSpurious);
        end

    else
        % set input
        if dim(Y) == 1
            % regression
            isVerified = rad(interval(Y)) <= delta;
        else
            % classification
            isVerified = all(interval(Y).sup <= 0, 'all');
        end
    end
end

function xs_ = aux_findCounterExample(nn,nn_abs,x,label,epsilon,delta)
    % compute attack
    xs_ = nn_abs.computeFGSMAttack(x,label,struct,epsilon);
    
    % evaluate attack on original network
    ys_ = nn.evaluate(xs_);
    [~,idxSpurious] = aux_checkSpecs(ys_,delta);

    % extract counter example
    xs_ = x(:,idxSpurious);
end

% ------------------------------ END OF CODE ------------------------------
