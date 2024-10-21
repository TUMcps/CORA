function linReset_ = resolve(linReset,flowList,stateBinds,inputBinds)
% resolve - resolves inputs which are outputs from other equations
%
% Syntax:
%    linReset = resolve(linReset,flowList,stateBinds,inputBinds)
%
% Inputs:
%    linReset - linearReset object
%    flowList - list of dynamical equations containing output equations
%    stateBinds - states of the high-dimensional space that correspond to
%                 the states of the low-dimensional reset object
%    inputBinds - connections of input to global input/outputs of other
%                 components
%
% Outputs:
%    linReset_ - linearReset object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/resolve

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       04-April-2022
% Last update:   ---
% Last revision: 11-October-2024 (MW, refactored from other functions)

% ------------------------------ BEGIN CODE -------------------------------

% only supports linear output functions for input resolution
if ~all(cellfun(@(sys) isa(sys,'linearSys'),flowList,'UniformOutput',true))
    throw(CORAerror('CORA:wrongValue','second',...
        'All output functions for input resolution must be linear.'));
end

% loop over all inputs of synchronized reset, three possibilities:
%   - the input is a dummy input (same as global...?)
%   - the input is a global input
%   - the input is an output of another component
% in the last case, we have to resolve the equation, i.e., the first input
% to the first component is the second output of the second component
numComp = numel(inputBinds);

% rewrite inputs (cell = component, values = global indices of component's inputs)
idxCompInput = arrayfun(@(i) size(vertcat(inputBinds{1:i-1}),1)+1:...
    size(vertcat(inputBinds{1:i}),1),1:numComp,'UniformOutput',false);

% read out reset matrices and vector
A = linReset.A;
B = linReset.B;
c = linReset.c;

% loop over all components
for thisCompIdx=1:numComp
    % resolve input binds of i-th component

    for j=1:size(inputBinds{thisCompIdx},1)

        % index of other component (if global input, index is 0)
        otherCompIdx = inputBinds{thisCompIdx}(j,1);

        % check if ith input of synchronized reset is linked to a global
        % input (or is a dummy input)
        if otherCompIdx == 0
            % nothing to resolve...
            continue
        end

        % obtain output equation of other component
        C = flowList{otherCompIdx}.C;
        D = flowList{otherCompIdx}.D;
        q = flowList{otherCompIdx}.k;

        % current input is *l*-th output of other component
        l = inputBinds{thisCompIdx}(j,2);
        
        % convert the output matrix to a matrix for consistency below
        if isscalar(C) && C == 1
            C = eye(flowList{otherCompIdx}.nrOfStates);
        end
    
        % insert output l-th output of other component k, i.e., y_k(l), for
        % the j-th input of component i, i.e., u_i(j):
        % B_i(:,j) * u_i(j)
        % = B_i(:,j) * y_k(l)
        % = B_i(:,j) * ( C_k(l,:) * x_k + D_k(l,:) * u_k + q_k(l) )
        % = B_i(:,j) * C_k(l,:) * x_k
        %   + B_i(:,j) * D_k(l,:) * u_k
        %   + B_i(:,j) * q_k(l)

        % check for circular input sequences: the term D_k(l,:)*u_k may
        % point back to u_i(j), so we force D_k(l,jj) to be zero or u_k(jj)
        % to be a global input to prevent infinite loops
        for jj = 1:size(D,2)
            if ~(D(l,jj) == 0 || inputBinds{otherCompIdx}(jj,1) == 0)
                throw(CORAerror('CORA:notSupported',...
                    ['It is not allowed for the feedthrough matrix D '...
                    'to point to inputs that are defined by the '...
                    'output of other subsystems, since it would '...
                    'otherwise be able to construct infinite loops!']));
            end
        end

        % read out part of (full) B matrix relevant for the current mapping
        B_part = B(stateBinds{thisCompIdx},idxCompInput{thisCompIdx}(j));
    
        % compute effect: BC affects A, BD affects B, Bq affects c
        BC_j = B_part * C(l,:);
	    BD_j = B_part * D(l,:);
        Bq_j = B_part * q(l);
    
        % update state matrix
        A(stateBinds{thisCompIdx},stateBinds{otherCompIdx}) = ...
            A(stateBinds{thisCompIdx},stateBinds{otherCompIdx}) + BC_j;
        
        % update input matrix by effect of re-mapped global inputs
        B(stateBinds{thisCompIdx},idxCompInput{otherCompIdx}) = ...
            B(stateBinds{thisCompIdx},idxCompInput{otherCompIdx}) + BD_j;

        % update offset vector
        c(stateBinds{thisCompIdx}) = c(stateBinds{thisCompIdx}) + Bq_j;

    end
end

% we now remove all columns of non-global inputs, as the corresponding
% inputs have been resolved
allBinds = vertcat(inputBinds{:});
idxGlobalInput = allBinds(:,1) == 0;
B = B(:,idxGlobalInput);

% up until now, input binds to the same global input were treated as
% separate columns in B, one for each component with that same global input
% we now merge these separate columns into one per global input and
numGlobalInputs = max(allBinds(allBinds(:,1) == 0,2));
for i=1:numGlobalInputs
    idx = allBinds(allBinds(:,1) == 0,2) == i;
    B(:,i) = sum(B(:,idx),2);
end
B = B(:,1:numGlobalInputs);

% initialize reset function (note that the constructor adds a dummy input
% if there is none, so that we don't have to do this here)
linReset_ = linearReset(A,B,c);

% ------------------------------ END OF CODE ------------------------------
