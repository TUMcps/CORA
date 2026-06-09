function res = benchmark_conZonotope_GPU_bounding()
% benchmark_conZonotope_GPU_bounding - evaluate and compare the different
%    bounding methods for constrained zonotopes.
%
% Syntax:
%    res = benchmark_conZonotope_GPU_bounding()
%
% Inputs:
%    -
%
% Outputs:
%    res - 
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       28-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')

% Check if a GPU is available.
try
    if gpuDeviceCount('available') == 0
        % No GPU available.
        return;
    end
catch e
    return;
end

% Specify number of dimensions, number of generators, and number of batch
% sizes.
n = 32; % Number of dimensions.
q = 256; % Number of generators.
bSz = 1; % 64; % Batch size.
p = 64; % Number of constraints.

% Sample a batch of random constraint zonotope.
cZs.c = rand([n bSz]);
cZs.G = rand([n q bSz]);
cZs.dr = 0;
cZs.A = rand([p q bSz]);
cZs.b = rand([p bSz]);

% We convert the constrained zonotope constrained to equality constraints.
% For a fair comparision, we view the constraints as (i) inequality and 
% (ii) equality constraints, because the conversion changes the number of 
% generators or constraints.
cZsEq1 = aux_castConZonotope( ...
    aux_2ConZonoWithEqConst(cZs),double(1));
cZsEq2 = aux_castConZonotope(cZs,double(1));
% We construct the equivalent constrained zonotope with inequalty 
% constraints.
cZsIneq1 = aux_castConZonotope(cZs,gpuArray(single(1)));
cZsIneq2 = aux_castConZonotope( ...
    aux_2ConZonoWithIneqConst(cZs),gpuArray(single(1)));

% Initialize results cell.
results = {};

% 1. Default CORA Bounding ------------------------------------------------
fprintf('Computing bounds ("CORA") ...');

% Initialize result values for CORA.
results{1} = struct( ...
    'name','CORA', ...
    'time',0, ...
    'lowerBounds',[], ...
    'upperBounds',[] ...
);
% There is no batching possible. Do a loop.
for i=1:bSz
    % 1. Use converted equality constraints.

    % Extract the parameters.
    [c1,G1,~,A1,b1] = aux_extractConZonotope(cZsEq1,i);
    % Construct the constrained zonotope.
    cZ1 = conZonotope(c1,G1,A1,b1);

    timerval = tic; % Start measuring time.
    I1 = interval(cZ1); % Compute the bounds.
    t1 = toc(timerval); % Stop the timer.
    if representsa(I1,'emptySet')
        l1 = NaN([n 1]);
        u1 = NaN([n 1]);
    else
        l1 = I1.inf;
        u1 = I1.sup;
    end

    % 2. View the generated constraints as equality constraints.

    % Extract the parameters.
    [c2,G2,~,A2,b2] = aux_extractConZonotope(cZsEq2,i);
    % Construct the constrained zonotope.
    cZ2 = conZonotope(c2,G2,A2,b2);

    timerval = tic; % Start measuring time.
    I2 = interval(cZ2); % Compute the bounds.
    t2 = toc(timerval); % Stop the timer.
    if representsa(I2,'emptySet')
        l2 = NaN([n 1]);
        u2 = NaN([n 1]);
    else
        l2 = I2.inf;
        u2 = I2.sup;
    end

    % Add the times.
    results{1}.time = results{1}.time + t1 + t2;
    % Add the bounds.
    results{1}.lowerBounds = [results{1}.lowerBounds l1 l2];
    results{1}.upperBounds = [results{1}.upperBounds u1 u2];
end
fprintf(' done\n');

% 2. Zonotope Bounding (ignore constraints) -------------------------------

fprintf('Computing bounds ("zonotope") ...');
% Compute everything in a batch-wise fashion.

% Extract the parameters.
[c1,G1,~,~,~] = aux_extractConZonotope(cZsIneq1,[]);

timerval = tic; % Start measuring time.
r1 = reshape(sum(abs(G1),2),[n bSz]); % Compute the radius.
l1 = c1 - r1; % Compute the lower bound.
u1 = c1 + r1; % Compute the upper bound.
t1 = toc(timerval); % Stop the timer.

% Extract the parameters.
[c2,G2,~,~,~] = aux_extractConZonotope(cZsIneq2,[]);

timerval = tic; % Start measuring time.
r2 = reshape(sum(abs(G2),2),[n bSz]); % Compute the radius.
l2 = c2 - r2; % Compute the lower bound.
u2 = c2 + r2; % Compute the upper bound.
t2 = toc(timerval); % Stop the timer.

% Construct the results.
results{end+1} = struct( ...
    'name','zonotope', ...
    'time',t1 + t2, ...
    'lowerBounds',reshape(permute(cat(3,l1,l2),[1 3 2]),[n 2*bSz]), ...
    'upperBounds',reshape(permute(cat(3,u1,u2),[1 3 2]),[n 2*bSz]) ...
);
fprintf(' done\n');

% 3. Fourier-Motzkin Bounding ---------------------------------------------

fprintf('Computing bounds ("fourier-motzkin") ...');
% Specify the options.
options.nn.conzonotope_bounding_method = 'fourier-motzkin';
options.nn.polytope_bound_approx_max_iter = 8;
% Set default values.
options = nnHelper.validateNNoptions(options);

timerval = tic; % Start measuring time.
% Compute the bounds.
[l1,u1,~,~] = conZonotope.approximateBoundsWithGPU(cZsIneq1,1,options);
t1 = toc(timerval); % Stop the timer.

timerval = tic; % Start measuring time.
% Compute the bounds.
[l2,u2,~,~] = conZonotope.approximateBoundsWithGPU(cZsIneq2,1,options);
t2 = toc(timerval); % Stop the timer.

% Construct the results.
results{end+1} = struct( ...
    'name','fourier-motzkin', ...
    'time',t1 + t2, ...
    'lowerBounds',reshape(permute(cat(3,l1,l2),[1 3 2]),[n 2*bSz]), ...
    'upperBounds',reshape(permute(cat(3,u1,u2),[1 3 2]),[n 2*bSz]) ...
);
fprintf(' done\n');

% 3. Iterative Optimization of Support Function Dual Variables ------------

% Specify number of iterations to test.
numIterations = [10 100 200 500 1000];
for i=1:length(numIterations)
    % Specify the name.
    name = sprintf('dual-iter (%d)',numIterations(i));
    fprintf('Computing bounds ("%s") ...',name);
    % Specify the options.
    options.nn.conzonotope_bounding_method = 'dual-iter';
    options.nn.conzonotope_bound_step_size = 1;
    options.nn.conzonotope_bound_max_iter = numIterations(i);
    % Set default values.
    options = nnHelper.validateNNoptions(options);
    
    timerval = tic; % Start measuring time.
    % Compute the bounds.
    [l1,u1,~,~] = conZonotope.approximateBoundsWithGPU(cZsIneq1,1,options);
    t1 = toc(timerval); % Stop the timer.
    
    timerval = tic; % Start measuring time.
    % Compute the bounds.
    [l2,u2,~,~] = conZonotope.approximateBoundsWithGPU(cZsIneq2,1,options);
    t2 = toc(timerval); % Stop the timer.
    
    % Construct the results.
    results{end+1} = struct( ...
        'name',name, ...
        'time',t1 + t2, ...
        'lowerBounds',reshape(permute(cat(3,l1,l2),[1 3 2]),[n 2*bSz]), ...
        'upperBounds',reshape(permute(cat(3,u1,u2),[1 3 2]),[n 2*bSz]) ...
    );
    fprintf(' done\n');
end

% Finally, show the results. ----------------------------------------------

aux_printStats(results);

% Set the output.
res = true;

end


% Auxiliary functions -----------------------------------------------------

% Coversion Functions. ----------------------------------------------------

function cZs_ = aux_castConZonotope(cZs,x)
    % Convert a batch of constrained zonotope to the same type as x.
    cZs_.c = cast(cZs.c,'like',x);
    cZs_.G = cast(cZs.G,'like',x);
    if isfield(cZs,'dr')
        cZs_.dr = cast(cZs.dr,'like',x);
    end
    cZs_.A = cast(cZs.A,'like',x);
    cZs_.b = cast(cZs.b,'like',x);
end

function [c,G,dr,A,b] = aux_extractConZonotope(cZs,idx)
    % Extract parameters of the constraint zonotope.
    c = cZs.c;
    G = cZs.G;
    if isfield(cZs,'dr')
        dr = cZs.dr;
    else
        dr = 0;
    end
    A = cZs.A;
    b = cZs.b;
    % Extract the correct batch indices.
    if ~isempty(idx)
        c = c(:,idx);
        G = G(:,:,idx);
        if numel(dr) > 1
            dr = dr(:,idx);
        end
        A = A(:,:,idx);
        b = b(:,idx);
    end
end

function cZsEq = aux_2ConZonoWithEqConst(cZsIneq)
    % Extract parameters of the constraint zonotope.
    [c,G,dr,A,b] = aux_extractConZonotope(cZsIneq,[]);

    % We convert the inequality constraints to equality constraints by 
    % adding a slack variable.

    % Obtain number of dimensions, generators, and batch size.
    [n,~,bSz] = size(G);
    % Obtain number of constraints.
    [p,~] = size(A);

    cZsEq.c = c;
    % Add the radius to the generators.
    if any(dr ~= 0,'all')
        G = cat(2,G,permute(dr,[1 3 2]).*eye(n));
        A = cat(2,A,zeros([p n bSz]));
    end
    cZsEq.dr = 0;
    % Add a slack variable.
    cZsEq.G = cat(2,G,zeros([n p bSz]));
    % Compute scale for the slack variable.
    s = 1/2*(sum(abs(A),2) + permute(b,[1 3 2]));
    cZsEq.A = cat(2,A,eye(p).*s);
    % Compensate for the slack variable.
    cZsEq.b = b - s(:,:);
end

function cZsIneq = aux_2ConZonoWithIneqConst(cZsEq)
    % Extract parameters of the constraint zonotope.
    [c,G,dr,A,b] = aux_extractConZonotope(cZsEq,[]);

    % We convert the equality constraints to inequality constraints.

    % Obtain number of dimensions, generators, and batch size.
    [n,~,bSz] = size(G);
    % Obtain number of constraints.
    [p,~] = size(A);

    cZsIneq.c = c;
    % Add the radius to the generators.
    if any(dr ~= 0,'all')
        G = cat(2,G,permute(dr,[1 3 2]).*eye(n));
        A = cat(2,A,zeros([p n bSz]));
    end
    cZsIneq.dr = 0;
    % The generators stay the same.
    cZsIneq.G = G;
    % We duplicate the constraints to represent an equality constraint 
    % with two inequalities, i.e., 
    %   A*\beta = b <--> A*\beta <= b && A*\beta >= b.
    cZsIneq.A = [A; -A];
    cZsIneq.b = [b; -b];
end

% Results Table. ----------------------------------------------------------

function results = aux_printStats(results)
    % Compute the ratio of the radii.
    rs = [];
    for i=1:length(results)
        % Compute and append the radii for the i-th method.
        rs = cat(3,rs,1/2*(results{i}.upperBounds - results{i}.lowerBounds));
    end
    % Set NaN values (indicate empty sets) to really small radii.
    rs(isnan(rs)) = 1e-3;
    % Compute the ratios.
    qs = max(rs./permute(rs,[1 2 4 3]),[],4);
    % Set the ratios in the results.
    for i=1:length(results)
        % Compute and append the radii for the i-th method.
        results{i}.radiiRatio = reshape(qs(:,:,i),[],1);
    end

    % Setup the table
    table = CORAtable('double', ...
        {'Method','Radii Ratio','Time [ms]'}, ...
        {'s','sum{%.2f & %.2f}','.2f'}, ...
        'ColumnWidths',[20 20 10] ...
    );
    table.printHeader();

    % Print the stats for each method.
    for i=1:length(results)
        % Print the stats.
        table.printContentRow({ ...
            results{i}.name, ...
            results{i}.radiiRatio, ...
            results{i}.time*1000 ...
        });
    end
    % Finish table.
    table.printFooter();
end

% ------------------------------ END OF CODE ------------------------------
