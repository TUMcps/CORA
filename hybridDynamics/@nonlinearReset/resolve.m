function nonlinReset_ = resolve(nonlinReset,flowList,stateBinds,inputBinds)
% resolve - resolves inputs which are outputs from other equations
%
% Syntax:
%    nonlinReset_ = resolve(nonlinReset,flowList,stateBinds,inputBinds)
%
% Inputs:
%    nonlinReset - nonlinearReset object
%    flowList - list of dynamical equations containing output equations
%    stateBinds - states of the high-dimensional space that correspond to
%                 the states of the low-dimensional reset object
%    inputBinds - connections of inputs to global inputs or outputs of
%                 other components
%
% Outputs:
%    nonlinReset_ - nonlinearReset object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearReset/resolve

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       04-April-2022
% Last update:   ---
% Last revision: 14-October-2024 (MW, refactored from other functions)

% ------------------------------ BEGIN CODE -------------------------------

% loop over all inputs of synchronized reset, three possibilities:
%   - the input is a dummy input (same as global...?)
%   - the input is a global input
%   - the input is an output of another component
% in the last case, we have to resolve the equation, i.e., the first input
% to the first component is the second output of the second component
numComp = numel(inputBinds);

% save indices of global inputs, as these are the only ones that will
% remain after the input resolution
allBinds = vertcat(inputBinds{:});
idxGlobalInput = allBinds(:,1) == 0;
numGlobalInputs = nnz(idxGlobalInput);
idxCompInput = arrayfun(@(i) size(vertcat(inputBinds{1:i-1}),1)+1:...
    size(vertcat(inputBinds{1:i}),1),1:numComp,'UniformOutput',false);

% init symbolic variables
x_sym = sym('x',[nonlinReset.preStateDim,1],'real');
u_sym = sym('u',[nonlinReset.inputDim,1],'real');
% insert into function
f_sym = nonlinReset.f(x_sym,u_sym);

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

        % current input is *l*-th output of other component
        l = inputBinds{thisCompIdx}(j,2);
    
        % insert output l-th output of other component k, i.e., y_k(l), for
        % the j-th input of component i, i.e., u_i(j)
        
        % evaluate output equation of other compontent
        g_otherComp_sym = flowList{otherCompIdx}.out_mFile(...
            x_sym(stateBinds{otherCompIdx}),...
            u_sym(idxCompInput{otherCompIdx}));
        % use only expression for l-th output
        g_otherComp_sym_l = g_otherComp_sym(l);

        % check if resulting equation uses any inputs at all -> must be
        % global input, otherwise infinite loops possible
        if any(logical(symvar(g_otherComp_sym_l) == u_sym ...
                & ~idxGlobalInput),'all')
            throw(CORAerror('CORA:notSupported',...
                    ['Output functions from other components must not '
                     'have local, but only global inputs.']));
        end

        % substitute inputs by output function
        f_sym = subs(f_sym,...
            u_sym(idxCompInput{thisCompIdx}(:,j)),g_otherComp_sym_l);

    end
end

% remove all non-global inputs which we have resolved to states and/or
% global inputs, so that only global inputs remain
f_sym = subs(f_sym,u_sym(idxGlobalInput),u_sym(1:numGlobalInputs));

% convert to function handle
f = matlabFunction(f_sym,'Vars',{x_sym,u_sym});

% initialize reset function (note that the constructor adds a dummy input
% if there is none, so that we don't have to do this here)
nonlinReset_ = nonlinearReset(f);

% ------------------------------ END OF CODE ------------------------------
