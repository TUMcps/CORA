function [count,out] = inputArgsLength(f,varargin)
% inputArgsLength - computes the number of inputs of a function handle
%
% Syntax:
%    [count,out] = inputArgsLength(f)
%    [count,out] = inputArgsLength(f,inpArgs)
%
% Inputs:
%    f - function handle 
%    inpArgs - number of input arguments for the function (max. 26)
%
% Outputs:
%    count - vector storing the length of each input argument
%    out - output dimension of the function handle
%
% Example:
%    f = @(x,u) [x(1)*x(5)^2; sin(x(3)) + u(2)];
%    inputArgsLength(f)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       11-September-2020
% Last update:   17-June-2022 (MW, speed up when inpArgs exceeds nargin to f)
% Last revision: 20-November-2022 (MW, restucture alternative method)

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments: number of input arguments to function handle
narginf = setDefaultValues({nargin(f)},varargin);

% try fast way to determine number of inputs (does not work if
% statements like "length(x)" occur in the dynamic function
try
    % create symbolic variables (length 100 for each input)
    narginvars = sym('x',[100,narginf]);
    narginvars_cell = num2cell(narginvars,1);
    
    % evaluate function
    fsym = f(narginvars_cell{1:narginf});
    
    % output dimension of f
    out = size(fsym);
    if any(out == 1)
        out = max(out);
    end
    % used variables from array of symbolic variables (size: [100,inpArgs])
    vars = symvar(fsym);
    % logical indices for used variables
    mask = ismember(narginvars,vars);
    % required dimension of each inpArgs for evaluation of f
    count = zeros(narginf,1);
    for i=1:narginf
        % last 'true' value in logical indices is required dimension
        tmp = find(mask(:,i),1,'last');
        if isempty(tmp)
            tmp = 0;
        end
        count(i) = tmp;
    end

    % sanity check: call function with computed number of input arguments
    try
        inputs = cell(narginf,1);
        for i=1:narginf
            inputs{i} = narginvars_cell{i}(1:count(i));
        end
        f(inputs{:});
        % return only now...
        return;
    end
end

% upper bound
bound = 1000;

% special case: only one input argument
if narginf == 1
    % increment number of input arguments until correct
    count = 0;
    narginvars = sym('x',[bound,narginf]);
    while count < bound
        count = count + 1;
        try
            output = f(narginvars(1:count));
            % get length of symbolic output
            out = size(output);
            if any(out == 1)
                out = max(out);
            end
            return
        end
    end
end

maxVal = 3;
narginvars = sym('x',[bound,narginf]);
narginvars_cell = num2cell(narginvars,1);

% checked combinations
checked_comb = zeros(1,narginf);

while true
    
    comb = combinator(maxVal,narginf) - 1;
    if ~isempty(checked_comb)
        % logical index which combinations remain
        comb_logIdx = true(size(comb,1),1);
        % find which already checked combinations are in current list
        for i=1:size(comb,1)
            for j=1:size(checked_comb,1)
                if any(all(checked_comb(j,:) == comb(i,:),2))
                    comb_logIdx(i) = false;
                end
            end
        end
        % remove already checked combinations
        comb = comb(comb_logIdx,:);
    end
    
    for i = 1:size(comb,1)
        % current combination of input arguments
        curr_comb = comb(i,:);

        % pass symbolic inputs to function handle
        try
            input = cell(narginf,1);
            for j = 1:narginf
                input{j} = narginvars_cell{j}(1:curr_comb(j));
            end
            output = f(input{:});

            % function evaluation successful!

            % get required length for each input argument
            count = curr_comb;
            % get length of symbolic output
            out = size(output);
            if any(out == 1)
                out = max(out);
            end
            return
        catch
            % proceed to next number of inputs combination
            continue;
        end
    end
    
    % append to list of checked combinations
    checked_comb = [checked_comb; comb];

    % limit number of checks
    if size(checked_comb,1) > bound
        throw(CORAerror('CORA:specialError',...
           'Could not determine length of input arguments!')); 
    end

    % increment number of inputs
    maxVal = maxVal + 3;
end

% ------------------------------ END OF CODE ------------------------------
