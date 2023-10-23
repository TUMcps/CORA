function df_i = lookupDf(layer,i)
% lookupDf - look-up table for i-th derivative of the given layer
%     global look-up table to reuse the computation of the derivative 
%
% Syntax:
%    res = nnHelper.lookupDf(layer,i)
%
% Inputs:
%    layer - nnActivationLayer
%    i - i-th derivative
%
% Outputs:
%    df_i - function handle of the i-th derivative
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       30-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init lookup dict
persistent lookupDf;
if isempty(lookupDf)
    lookupDf = containers.Map('KeyType','char','ValueType','any');
end

% find layer name
layerName = class(layer);

% check if present
if ~isKey(lookupDf,layerName)
    % init layer struct
    layerStruct = struct;

    % init first derivative
    syms x real
    sym_df = simplify(diff(layer.f(x)));
    layerStruct(1).sym_df = sym_df;

    % and function handle
    df = matlabFunction(sym_df);
    layerStruct(1).df = df;

    % store
    lookupDf(layerName) = layerStruct;
end

% read layer struct
layerStruct = lookupDf(layerName);

% returns a function handle of the i-th derivative of this layer
if i == 0
    % function
    df_i = layer.f;

elseif i <= length(layerStruct)
    % derivative is known
    df_i = layerStruct(i).df;

else
    % find derivative
    maxKnownDer = length(layerStruct);

    % find last known derivative
    sym_df_j1 = layerStruct(maxKnownDer).sym_df;

    % iterate up to desired derivative
    for j=maxKnownDer+1:i
        % compute symbolic derivative
        sym_df_j = simplify(diff(sym_df_j1));
        layerStruct(j).sym_df = sym_df_j;
        
        % compute function handle
        df_j = matlabFunction(sym_df_j);
        layerStruct(j).df = df_j;

        % prepare for next iteration
        sym_df_j1 = sym_df_j;
    end

    % retrieve result
    df_i = df_j;

    % store result
    lookupDf(layerName) = layerStruct;
end

end

% ------------------------------ END OF CODE ------------------------------
