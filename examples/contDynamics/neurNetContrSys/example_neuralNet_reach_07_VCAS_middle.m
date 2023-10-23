function completed = example_neuralNet_reach_07_VCAS_middle()
% example_neuralNet_reach_07_VCAS_middle - example of reachability analysis
%    for an aircraft collision avoidance system: Vertical-CAS, where we
%    also choose the middle acceleration; note that this example has
%    discrete time steps.
%
% Syntax:
%    completed = example_neuralNet_reach_07_VCAS_middle()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP21 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."
%       EPiC Series in Computing 80 (2021): 90-119.

% Authors:       Tobias Ladner
% Written:       02-June-2022
% Last update:   ---
% Last revision: 20-August-2022 (added verification)

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Vertical Collision Avoidance System (VCAS_middle)")

for hdot0 = [-19.5, -22.5, -25.5, -28.5]

    fprintf("dot{h}_0(0) = %.1f\n", hdot0)

    % Parameters ----------------------------------------------------------

    R0 = interval([-133; hdot0; 25], [-129; hdot0; 25]);
    tFinal = 10;

    % System Dynamics -----------------------------------------------------

    % parameters
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
        network = neuralNetwork.readNNetNetwork(sprintf('VertCAS_noResp_pra0%d_v9_20HU_200.nnet', i));

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

    % always choose middle acc
    accs = accs(:, 2)';
    adv_init = 1;


    % Specification -------------------------------------------------------

    % A1 = [0 1 0; 0 -1 0]; b1 = [0.5; 0.5];
    % A2 = [eye(3); -eye(3)]; b2 = ones(6,1);
    % set = polytope([blkdiag(A1,A2),zeros(8,6)],[b1;b2]);
    unsafeSet = [; ...
        -100, 100; ... % h
        -inf, -inf; ... % h^.
        -25, -25 + tFinal; ... % tau
        ];

    unsafeSet = interval(unsafeSet(:, 1), unsafeSet(:, 2));
    spec = specification(unsafeSet, 'unsafeSet');


    % Simulation ----------------------------------------------------------

    numSims = 10;
    x = R0.randPoint(numSims);
    simRes = {x};
    adv = ones(10, 1) * adv_init;

    tic
    for i = 1:tFinal
        X = [];
        for j = 1:numSims
            x_j = x(:, j);

            % evaluate network
            network = networks{adv(j)};
            logits = network.evaluate(x_j, struct);

            % determine output
            [~, adv_next] = max(logits);

            % check next adv with compliance
            if complCL{adv_next}.contains(x_j(2))
                acc = 0;
            else
                acc = accs{adv_next};
            end

            u_j = [adv_next; acc];
            adv(j) = adv_next;

            % update state
            x_j = x_j + dx(x_j, u_j);
            x(:, j) = x_j;
            X = [X, x_j];
        end
        simRes = [simRes; X];
    end
    tSim = toc;
    disp(['Time to compute random simulations: ', num2str(tSim)]);


    % Check Violation -----------------------------------------------------

    tic
    isVio = false;
    for i = 1:length(simRes)
        x = simRes{i};
        for j = 1:length(unsafeSet) - 1
            isVio = isVio || any( ...
                (infimum(unsafeSet(j)) <= x(j, :)) & ...
                (x(j, :) <= supremum(unsafeSet(j))));
            if isVio
                break;
            end
        end
    end
    tVio = toc;
    disp(['Time to check violation in simulations: ', num2str(tVio)]);

    if isVio
        disp("Result: VIOLATED")
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
            evParams.bound_approx = true;
            evParams.polynomial_approx = "lin";

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
            disp('Result: VERIFIED');
        else
            disp('Result: UNKNOWN');
        end
    end
    disp(['Total Time: ', num2str(tSim+tVio+tComp+tVeri)]);


    % Visualization -------------------------------------------------------

    disp("Plotting..");
    figure; hold on; box on;

    % title
    title(sprintf('$\\dot{h}_0(0) = %.1f$', hdot0), 'interpreter', 'latex')

    % plot specification
    us = plot(unsafeSet, [3, 1], 'FaceColor', [0.8, 0, 0]);
    alpha(us, .5)

    % plot reachable set
    rs = plot([], []);
    for i = 1:length(R)
        rs = plot(interval([0, 0, -1; 1, 0, 0]*R{i})+interval([0; 0], [1; 0]),...
            [1,2],'FaceColor',[.8, .8, .8],'EdgeColor',[.8, .8, .8]);
    end

    for i = 1:length(simRes)
        t = simRes{i}(3, :);
        h = simRes{i}(1, :);
        sims = scatter(-t, h, '.k');
    end

    % labels, limits, and legend
    xlabel('Time -\tau');
    ylabel('h (Relative Position of Intruder to Ownship)');
    xlim([-25, -25 + tFinal]);
    legend([us, rs, sims], "NMAC", "Reachable Set", "Simulations");

end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
