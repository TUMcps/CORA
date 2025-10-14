function [Zred, minimalVolume] = priv_minSubZonotope(Z, varargin)
% priv_minSubZonotope - Computes the optimal outer-approximation of Z with an order 
% reduction of 1/d by removing one generator, and rescaling d generators.
% Here d = dimension of Z. For higher order reduction, minSubZonotope is applied 
% iteratively onto Z s.t. Zred is approximated after second iteration.
%
% Syntax:
%    Zred = priv_minSubZonotope(Z);
%    [Zred, minVolume] = priv_minSubZonotope(Z, 'method', method,
%    'order', order);
%
% Inputs:
%    Z - zonotope object
%    varargin - Value pairs (optional) 
%    method - 'minVolume' (default)
%           - 'VSF' uses relation between scaling factor and volume
%           - 'TSF' uses product of scalers 
%           - 'shortestGenVol' remove shortest generator, determine Zred by
%              volume
%           - 'shortestGenTSF' remove shortest generator, determine Zred by
%              scaling factors
%           - 'scalingDet'  takes the product of scaling factors and
%              rescaled sub generator matrix (Performs significantly worse
%              than the other methods and is not recommended)
%    order - reduces zonotope to this order, order = (n-1)/d is default
%    composeVolume - calculates volume for methods that do not require the 
%                    volume during evaluation step. Default is false.
%
% Outputs:
%    Zred - zonotope object
%    varargout - volume of zonotope object (optional output, heavily
%              influences runtime of function)
%
% Example:
%    Z = zonotope([0;0],[1 -1 0; 0 1 -1]);
%    [Zred, runtime, minimalVolume] = minSubZonotope(Z, 'method', 'minVolume', 'composeVolume', true);
%
% References:
%    [1] E. Grover et al. "Determinants and the volumes of parallelotopes
%        and zonotopes", 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/minSubZonotope

% Authors:       Ngoc Toan Nguyen, Sebastian Sigl
% Written:       01-November-2024
% Last update:   11-June-2025 (SS, improved readability, removed unnecessary parts)
% Last revision: 25-April-2025 (SS, cleanup)

% ------------------------------ BEGIN CODE -------------------------------

% Set parameters
c = Z.c;
G = Z.G;
[d, NrOfGens] = size(G); % [Dimension, Number of Generators]

% Ensure number of input is Zonotope and pairs
if mod(nargin, 2) ~= 1
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'method','order','composeVolume'});
    % method given?
    [NVpairs, method] = readNameValuePair(NVpairs,'method', 'char', 'minVolume');
    % order given?
    [NVpairs, order] = readNameValuePair(NVpairs,'order', 'scalar', (NrOfGens - 1)/d);
    % composeVolume given?
    [NVpairs, composeVolume] = readNameValuePair(NVpairs,'composeVolume', 'logical', false);
end

% Catch errors in inputs
aux_checkInputArgs(Z, method, order);

batch_size = 1e7;
batches = batchCombinator(NrOfGens, NrOfGens - 1, batch_size, struct);

switch method
    case 'minVolume'
        [Zred, minimalVolume] = aux_minVolume(NrOfGens, d, batch_size, c, G, batches);
    case 'VSF'
        [Zred, minimalVolume] = aux_VSF(NrOfGens, d, batch_size, c, G, batches, composeVolume, order);
    case 'TSF'
        [Zred, minimalVolume] = aux_TSF(NrOfGens, d, batch_size, c, G, batches, composeVolume, order);
    case {'shortestGenVol','shortestGenTSF'}
        [Zred, minimalVolume] = aux_shortestGen(method, NrOfGens, d, batch_size, c, G, composeVolume, order);
    case 'scalingDet'
        [Zred, minimalVolume] = aux_scalingDet(NrOfGens, d, batch_size, c, G, batches, composeVolume, order);
end

% Iteratively calls minSubZonotope to reduce order of zonotope greater than 1/d.
% Zred will be approximated after the second iteration onwards.
if order < (NrOfGens - 1)/d

    % Calculate number of iteration needed after order of Z by 1/d
    NrOfIter = int16(NrOfGens - order * d) - 1;

    % Approximates zonotope to given order
    for iter = 1 : NrOfIter
        Zred = priv_minSubZonotope(Zred, 'method', method, 'composeVolume', false); % reduce previous Zred
    end
    
    % Outputs volume if varargout was requested
    if composeVolume == true && nargout > 2
        minimalVolume = volume(Zred);
    end

end

end


% Auxiliary functions -----------------------------------------------------

function aux_checkInputArgs(Z, method, order)
    % check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED

        inputArgsCheck({{Z, 'att', 'zonotope'}})

        % Lists of all possible methods
        methodList = {'minVolume', 'VSF', 'TSF', 'shortestGenVol', 'shortestGenTSF', 'scalingDet'};

        % Catch errors in inputs
        if ~any(strcmp(methodList, method))
            throw(CORAerror('CORA:wrongValue','name-value pair: method',"'minVolume', 'VSF', 'TSF', 'shortestGenVol', 'shortestGenTSF', 'scalingDet'"));
        end

        if order < 1
            throw(CORAerror("CORA:wrongValue", 'Order is smaller than 1. No outer approximation possible.'));
        end
        
    end
end

function [Zred, minimalVolume] = aux_minVolume(NrOfGens, d, batch_size, c, G, batches)
    % Minimum Volume method
    minVolume = 0;
    batchesVolume = batchCombinator(NrOfGens, int16(d), ...
        batch_size, struct); % Combination of all columns for sub determinants
    
    % Dictionary of all sub determinant values of the original generator matrix
    det_database = dictionary();
    for g=1:size(batchesVolume, 1)
        det_database({batchesVolume(g,:)}) = abs(det(G(:, batchesVolume(g, :))));
    end
    
    % Construct all possible zonotopes
    for i=1:size(batches, 1)
        subGens = G;
        logicalidx = true(1, NrOfGens);
        logicalidx(batches(i,:)) = false; % Set indices of non-discarded columns to false
        subGens(: , logicalidx) = zeros(d, 1); % Set discarded generator in generator matrix to 0
        discard_gen = G(:, logicalidx); % Discarded generator
    
        % Iterate through all d x d sub generator matrix
        for j=1:size(batchesVolume, 1)
    
            pre_scale_Gens = subGens(:, batchesVolume(j, :)); % Generators before scaling
    
            if sum(all(pre_scale_Gens == 0)) == 0 % Avoid zero-column
                factors = linsolve(pre_scale_Gens, discard_gen); % Calculate scaling factors
    
                dim_adjust_factors = ones(NrOfGens, 1);
                dim_adjust_factors(batchesVolume(j,:)) = ones(length(factors), 1 ...
                    ) + abs(factors); % Scaling vector to rescale G, dim: num_generators x 1
    
                accVol = 0; % Initiate volume calculation
    
                % Calculate volume by the using determinant dictionary
                for m=1:size(batchesVolume, 1)
                    preVolGens = subGens(:, batchesVolume(m, :));
                    if sum(all(preVolGens == 0)) == 0
                        currVol = det_database({batchesVolume(m,:)} ...
                            ) * prod(dim_adjust_factors(batchesVolume(m,:)));
                        accVol = accVol + currVol;
                    end
                end
                vol = 2^d * accVol;
    
                % Initiate first Zcandid as Zred
                if minVolume == 0
                    minG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                    minG = minG(:, batches(i,:)); % Scaled generator matrix without discarded generator
                    minVolume = vol; % Minimal volume
                end
    
                % Take smaller zonotope as Zred
                if vol < minVolume
                    minG = subGens * diag(dim_adjust_factors);
                    minG = minG(:, batches(i,:));
                    minVolume = vol;
                end
            end
        end
    end
    Zred = zonotope(c, minG);
    minimalVolume = minVolume; % Outputs volume if varargout was requested
end

function [Zred, minimalVolume] = aux_VSF(NrOfGens, d, batch_size, c, G, batches, composeVolume, order)
    % VSF method
    minScale = 0; % Initiate evaluation parameter phi_VSF

    batchesVolume = batchCombinator(NrOfGens, int16(d), ...
        batch_size, struct);

    % Construct all possible zonotopes
    for i=1:size(batches, 1)
        subGens = G;
        logicalidx = true(1, NrOfGens); % Create logical array
        logicalidx(batches(i, :)) = false; % Set indices of non-discarded columns to false2
        subGens(:, logicalidx) = zeros(d, 1); % Sub generator matrix with discarded column being 0
        discard_gen = G(:, logicalidx); % discarded generator

        % Iterate through all d x d sub generator matrix
        for j=1:size(batchesVolume, 1)

            pre_scale_Gens = subGens(:, batchesVolume(j, :)); % To be scaled generators

            if sum(all(pre_scale_Gens == 0)) == 0 % Avoid zero-column
                factors = linsolve(pre_scale_Gens, discard_gen); % New scaling factors

                dim_adjust_factors = ones(NrOfGens, 1);
                dim_adjust_factors(batchesVolume(j,:)) = ones(length(factors), 1 ...
                    ) + abs(factors); % Factor vector, dim: num_generators x 1

                factors = abs(factors) + ones(length(factors), 1);

                accVolScale = 0; %Initiate phi_VSF

                % Calculate phi_VSF of Zcandid
                for k=1:length(factors)
                    if (d - k) < (NrOfGens - d) % only calculates binomial coefficients C(n, k) if n > k
                        batch_factors = batchCombinator(length(factors), k, batch_size, struct);
                        for m=1:size(batch_factors,1)
                            currVolScale = prod(factors(batch_factors(m, :))) * nchoosek(NrOfGens - d, d - k);
                            accVolScale = accVolScale + currVolScale;
                        end
                    end
                end
                
                % Initiate first Zcandid as Zred and phi_VSF as optimal phi_VSF
                if minScale == 0
                    minG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                    minG = minG(:, batches(i,:));
                    minScale = accVolScale;
                end
                
                % Take smaller phi_VSF as Zred
                if accVolScale < minScale
                    minG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                    minG = minG(:, batches(i,:));
                    minScale = accVolScale;
                end
            end
        end
    end
    Zred = zonotope(c, minG);
    
    % Outputs volume if varargout was requested
    if composeVolume == true && nargout > 1 && order == (NrOfGens - 1) / d
        minimalVolume = volume(Zred);
    else
        minimalVolume = [];
    end
end

function [Zred, minimalVolume] = aux_TSF(NrOfGens, d, batch_size, c, G, batches, composeVolume, order)
    % TSF method
    minSimpleScale = 0; % Initiate evaluation parameter phi_TSF

    batchesVolume = batchCombinator(NrOfGens, int16(d), ...
        batch_size, struct); 
    
    % Construct all possible zonotopes
    for i=1:size(batches, 1)
        subGens = G;
        logicalidx = true(1, NrOfGens); % Create logical array
        logicalidx(batches(i, :)) = false; % Set indices of non-discarded columns to false2
        subGens(:, logicalidx) = zeros(d, 1); % Sub generator matrix with discarded column being 0
        discard_gen = G(:, logicalidx); % discarded generator

        % Iterate through all d x d sub generator matrix
        for j=1:size(batchesVolume, 1)

            pre_scale_Gens = subGens(:, batchesVolume(j, :)); % To be scaled generators

            if sum(all(pre_scale_Gens == 0)) == 0 % Avoid zero-column
                factors = linsolve(pre_scale_Gens, discard_gen); % New scaling factors

                dim_adjust_factors = ones(NrOfGens, 1);
                dim_adjust_factors(batchesVolume(j,:)) = ones(length(factors), 1 ...
                    ) + abs(factors); % Factor vector, dim: num_generators x 1

                simpleVolScale = prod(dim_adjust_factors); % product of scaling factors

                % Initiate first Zcandid as Zred and phi_TSF as optimal phi_TSF
                if minSimpleScale == 0 
                    minG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                    minG = minG(:, batches(i,:));
                    minSimpleScale = simpleVolScale;
                end

                % Take smaller phi_TSF as Zred
                if simpleVolScale < minSimpleScale 
                    minG = subGens * diag(dim_adjust_factors); 
                    minG = minG(:, batches(i,:));
                    minSimpleScale = simpleVolScale;
                end
            end
        end
    end
    Zred = zonotope(c, minG);

    % Outputs volume if varargout was requested
    if composeVolume == true && nargout > 1 && order == (NrOfGens - 1) / d
        minimalVolume = volume(Zred);
    else
        minimalVolume = [];
    end
end

function [Zred, minimalVolume] = aux_shortestGen(method, NrOfGens, d, batch_size, c, G, composeVolume, order)
    % shortest Generator method
    minVolume = 0; % Initiate minimal Volume
    minSimpleScale = 0; % Initiate evaluation parameter phi_TSF

    % Determine shortest generator
    genLengths = vecnorm(G, 2, 1);
    [~, minIndex] = min(genLengths);
    logicalidx = true(1, NrOfGens); % Create logical array
    logicalidx(minIndex) = false; % Set indices of non-discarded columns to false
    
    % Remove shortest generator from generator matrix
    subGens = G(:, logicalidx);
    discard_gen = G(:, minIndex);

    [~, NrOfSubGens] = size(subGens);

    batchesVolume = batchCombinator(NrOfSubGens, int16(d), ...
        batch_size, struct); 

    % Iterate through all d x d sub generator matrix
    for j=1:size(batchesVolume, 1)

        pre_scale_Gens = subGens(:, batchesVolume(j, :)); % To be scaled generators

        factors = linsolve(pre_scale_Gens, discard_gen); % New scaling factors

        dim_adjust_factors = ones(NrOfSubGens, 1);
        dim_adjust_factors(batchesVolume(j,:)) = ones(length(factors), 1 ...
            ) + abs(factors); % Factor vector, dim: num_generators x 1

        % Shortest Generator method + Minimal Volume
        if strcmp(method, "shortestGenVol")
            currG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
            currZ = zonotope(c, currG);
            currVol = volume(currZ);

            % Initiate first Zcandid as Zred
            if minVolume == 0 
                Zred = currZ;
                minVolume = currVol;
            end

            % Take smaller zonotope as Zred
            if currVol < minVolume
                Zred = currZ;
                minVolume = currVol;
            end

            % Outputs volume if varargout was requested
            if composeVolume == true && nargout > 1 && order == (NrOfGens - 1) / d
                minimalVolume = volume(Zred);
            else
                minimalVolume = [];
            end

        % Shortest Generator method + TSF
        elseif strcmp(method, "shortestGenTSF")
            simpleVolScale = prod(dim_adjust_factors); % product of scaling factors

            % Initiate first Zcandid as Zred and phi_TSF as optimal phi_TSF
            if minSimpleScale == 0
                minG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                minSimpleScale = simpleVolScale;
            end

            % Take smaller phi_TSF as Zred
            if simpleVolScale < minSimpleScale 
                minG = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                minSimpleScale = simpleVolScale;
            end
        end
    end

    % Outputs Zred from Shortest Generator method + TSF
    if strcmp(method, "shortestGenTSF")
        Zred = zonotope(c, minG);
        
        % Outputs volume if varargout was requested
        if composeVolume == true && nargout > 1 && order == (NrOfGens - 1) / d
            minimalVolume = volume(Zred);
        else
            minimalVolume = [];
        end
    end
end

function [Zred, minimalVolume] = aux_scalingDet(NrOfGens, d, batch_size, c, G, batches, composeVolume, order)
    % scaling Determinant method
    min_scale_det = 0; % Initiate minimal scaling-determinant parameter

    batchesVolume = batchCombinator(NrOfGens, int16(d), ...
        batch_size, struct);

    % Dictionary of all sub determinant values of the original generator matrix
    det_database = dictionary(); 
    for g=1:size(batchesVolume, 1)
        detG = abs(det(G(:, batchesVolume(g, :)))); % Determinant of Sub-generator matrix
        det_database({batchesVolume(g,:)}) = detG; % Add determinant to dictionary
    end
    
    % Construct all possible zonotopes
    for i=1:size(batches, 1)
        logicalidx = true(1, NrOfGens); % Create logical array
        logicalidx(batches(i, :)) = false; % Set indices of non-discarded columns to false

        subGens = G;
        subGens(:, logicalidx) = zeros(d, 1); % Sub generator matrix with discarded column being 0
        discard_gen = G(:, logicalidx); % Discarded generator

        % Iterate through all d x d sub generator matrix
        for j=1:size(batchesVolume, 1)

            pre_scale_Gens = subGens(:, batchesVolume(j, :)); % To be scaled generators

            if sum(all(pre_scale_Gens == 0)) == 0 % Avoid zero-column
                factors = linsolve(pre_scale_Gens, discard_gen); % New scaling factors
                scale_det = prod(factors) * det_database({batchesVolume(j, :)}); % Calculate scaling-determinant parameter

                % Initiate first Zcandid as Zred and minimal scaling-determinant parameter
                if min_scale_det == 0
                    dim_adjust_factors = ones(NrOfGens, 1);
                    dim_adjust_factors(batchesVolume(j,:)) = ones(length(factors), 1 ...
                        ) + abs(factors); % Factor vector, dim: num_generators x 1
                    new_G = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                    new_G = new_G(:, batches(i,:)); % New generator matrix without zero column

                    min_scale_det = scale_det;
                end

                % Take smaller scaling-determinant parameter as Zred
                if scale_det > min_scale_det 
                    dim_adjust_factors = ones(NrOfGens, 1);
                    dim_adjust_factors(batchesVolume(j,:)) = ones(length(factors), 1 ...
                        ) + abs(factors); % Factor vector, dim: num_generators x 1
                    new_G = subGens * diag(dim_adjust_factors); % new scaled generator matrix
                    new_G = new_G(:, batches(i,:)); % New generator matrix without zero column

                    min_scale_det = scale_det;
                end
            end
        end
    end
    Zred = zonotope(c, new_G);
    
    % Outputs volume if varargout was requested
    if composeVolume == true && nargout > 1 && order == (NrOfGens - 1) / d
        minimalVolume = volume(Zred);
    else
        minimalVolume = [];
    end
end

% ------------------------------ END OF CODE ------------------------------
