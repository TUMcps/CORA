function [completed,res,tTotal] = example_neuralNet_reach_07_VCAS(varargin)
% example_neuralNet_reach_07_VCAS - example of reachability analysis
%    for an aircraft collision avoidance system: Vertical-CAS, where we
%    also choose the worst acceleration; note that this example has
%    discrete time steps
%
% Syntax:
%    completed = example_neuralNet_reach_07_VCAS()
%
% Inputs:
%    type - type of acceleration: 'worst' or 'middle'
%    hdot0 - e.g. one of [-19.5, -22.5, -25.5, -28.5]
%
% Outputs:
%    completed - true/false
%    res - verification result
%    tTotal - total time
%
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP21 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."
%       EPiC Series in Computing 80 (2021): 90-119.

% Authors:       Tobias Ladner
% Written:       02-June-2022
% Last update:   30-March-2022 (TL, ARCH'23 revisions)
% Last revision: 20-August-2022 (added verification)

% ------------------------------ BEGIN CODE -------------------------------

[type,hdot0] = setDefaultValues({'worst',-19.5}, varargin);
possibleTypes = {'worst', 'middle'};
inputArgsCheck({ ...
    {type, 'str', possibleTypes}; ...
    {hdot0, 'att', 'numeric', 'scalar'}; ...
})

fprintf("BENCHMARK: Vertical Collision Avoidance System (VCAS: %s)\n", type)
fprintf("dot{h}_0(0) = %.1f\n", hdot0)

% Parameters ----------------------------------------------------------

R0 = interval([-133; hdot0; 25], [-129; hdot0; 25]);
tFinal = 10;

% System Dynamics -----------------------------------------------------

% parameter
dt = 1;
g = 32.2;

% system dynamics
dx = @(x, u) [ ...
    - x(2, :) * dt - 0.5 * u(2, :) * dt^2; ...
    u(2, :) * dt; ...
    ones(1, size(x, 2)) * -dt; ...
    ];

% Neural Network ------------------------------------------------------

% load neural network controllers
% 9x [3, 20, 20, 20, 20, 20, 9]
networks = cell(9, 1);
for i = 1:length(networks)
    % load
    network = neuralNetwork.readNNetNetwork(...
        sprintf('VertCAS_noResp_pra0%d_v9_20HU_200.nnet', i));

    % normalize input
    offset = -[0; 0; 20]; % Means
    scale = diag(1./[16000; 5000; 40]); % Ranges

    % prepend to network
    network = neuralNetwork([; ...
        {nnLinearLayer(eye(3), offset, "Means")}; ...
        {nnLinearLayer(scale, zeros(3, 1), "Ranges")}; ...
        network.layers; ...
        ]);

    networks{i} = network;
end

accs = {; ...
    -g / 8, 0, g / 8; ... % COC
    -g / 3, -7 * g / 24, -g / 4; ... % DNC
    g / 4, 7 * g / 24, g / 3; ... % DND
    -g / 3, -7 * g / 24, -g / 4; ... % DES1500
    g / 4, 7 * g / 24, g / 3; ... % CL1500
    -g / 3, -g / 3, -g / 3; ... % SDES1500
    g / 3, g / 3, g / 3; ... % SCL1500
    -g / 3, -g / 3, -g / 3; ... % SDES2500
    g / 3, g / 3, g / 3; ... % SCL2500
    };

complCL = {; ...
    interval(-Inf, -Inf); ... % COC (always choose from accs)
    interval(-Inf, 0); ... % DNC
    interval(0, Inf); ... % DND
    interval(-Inf, -1500); ... % DES1500
    interval(1500, Inf); ... % CL1500
    interval(-Inf, -1500); ... % SDES1500
    interval(1500, Inf); ... % SCL1500
    interval(-Inf, -2500); ... % SDES2500
    interval(2500, Inf); ... % SCL2500
    };

% choose acc
if strcmp(type, 'worst')
    accs = accs(:, 1)';
elseif strcmp(type, 'middle')
    accs = accs(:, 2)';
else
    throw(CORAerror("CORA:wrongValue", "first", possibleTypes))
end
adv_init = 1;

% Specification -------------------------------------------------------

% A1 = [0 1 0; 0 -1 0]; b1 = [0.5; 0.5];
% A2 = [eye(3); -eye(3)]; b2 = ones(6,1);
% set = polytope([blkdiag(A1,A2),zeros(8,6)],[b1;b2]);
unsafeSet = [; ...
    -100, 100; ... % h
    -1000, 1000; ... % h^.
    25 - tFinal, 25; ... % tau
    ];

unsafeSet = interval(unsafeSet(:, 1), unsafeSet(:, 2));
spec = specification(unsafeSet, 'unsafeSet');

% Simulation ----------------------------------------------------------

numSims = 10;
adv = ones(10, 1) * adv_init;

tic
simRes = [];
time = 1:tFinal;
for j = 1:numSims
    x = R0.randPoint(1);
    xs = x';

    for i = time
        % evaluate network
        network = networks{adv(j)};
        logits = network.evaluate(x, struct);

        % determine output
        [~, adv_next] = max(logits);

        % check next adv with compliance
        if complCL{adv_next}.contains(x(2))
            acc = 0;
        else
            acc = accs{adv_next};
        end

        u_j = [adv_next; acc];
        adv(j) = adv_next;

        % update state
        x = x + dx(x, u_j);

        % store
        xs = [xs; x'];
    end
    simRes = [simRes; simResult({xs}, {[0 time]'})];
end
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);

% Check Violation -----------------------------------------------------

tic
isVio = ~check(spec, simRes);
tVio = toc;
disp(['Time to check violation in simulations: ', num2str(tVio)]);

if isVio
    res = 'VIOLATED';
    tComp = 0;
    tVeri = 0;
    R = {};
else

    % Reachability Analysis -------------------------------------------

    R0 = polyZonotope(R0);
    R = {R0};

    tic
    Ri = R0;
    adv_cur = adv_init;
    for i = 1:tFinal
        evParams = struct();
        evParams.poly_method = "singh";

        % evaluate network
        network = networks{adv_cur};
        logits = network.evaluate(Ri, evParams);

        % determine output
        logits = interval(logits);
        adv_next = -1;
        for k = 1:9
            k_ = [1:k - 1, k + 1:9];
            if all(supremum(logits(k_)) <= infimum(logits(k)))
                adv_next = k;
                break;
            end
        end
        if adv_next == -1
            disp("Unable to determine maximum!")
            break;
        end

        % check next adv with compliance
        if complCL{adv_next}.contains(interval(project(Ri, 2)))
            acc = 0;
        else
            acc = accs{adv_next};
        end

        u_i = [adv_next; acc];
        adv(j) = adv_next;

        % update state
        Ri = [1, -dt, 0; 0, 1, 0; 0, 0, 1] * Ri;
        Ri = Ri + [-0.5 * u_i(2) * dt^2; u_i(2) * dt; -1];

        R = [R, {Ri}];
        adv_cur = adv_next;
    end
    tComp = toc;
    disp(['Time to compute reachable set: ', num2str(tComp)]);

    tic
    Rend = R{end};
    Rend = interval(Rend);

    isVeri = ~isIntersecting(Rend, unsafeSet);
    tVeri = toc;
    disp(['Time to check verification: ', num2str(tVeri)]);

    if isVeri && i == tFinal
        res = 'VERIFIED';
    else
        res = 'UNKNOWN';
    end
end
tTotal = tSim+tVio+tComp+tVeri;
disp(['Total Time: ', num2str(tTotal)]);
disp(['Result: ' res])

% Visualization -------------------------------------------------------

disp("Plotting..");
figure; hold on; box on;

% title
title(sprintf('$\\dot{h}_0(0) = %.1f$', hdot0), 'interpreter', 'latex')

% plot specification
M = [0, 0, -1; 1, 0, 0];
plot(M*spec.set, [1 2], 'FaceColor', CORAcolor("CORA:unsafe"),'DisplayName','NMAC');

% plot reachable set
NVpairs = {'FaceColor',CORAcolor("CORA:reachSet"),'DisplayName','Reachable set'};
for i = 1:length(R)
    plot(M*R{i}+interval([0; 0], [1; 0]),[1,2],NVpairs{:});
    if i==1
        NVpairs = [NVpairs,{'HandleVisibility','off'}];
    end
end

% plot simulation
NVpairs = {'.k','DisplayName','Simulations'};
for i=1:length(simRes)
    scatter(-25:-15, simRes(i).x{1}(:, 1)',NVpairs{:});
    if i==1
        NVpairs = [NVpairs,{'HandleVisibility','off'}];
    end
end

% labels, limits, and legend
xlabel('Time -\tau');
ylabel('h (Relative Position of Intruder to Ownship)');
xlim([-25, -25 + tFinal]);
legend()

% example completed -------------------------------------------------------

completed = true;

% handling for ARCH competition
if nargout < 2
    clear res;
end
if nargout < 3
    clear tTotal;
end

end

% ------------------------------ END OF CODE ------------------------------
