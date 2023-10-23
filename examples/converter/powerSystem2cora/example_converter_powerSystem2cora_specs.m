function completed = example_converter_powerSystem2cora_specs()
% example_converter_powerSystem2cora_specs - example for creating 
% specifications of the power system benchmarks
%
% Syntax:
%    completed = example_converter_powerSystem2cora_specs()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280

% Authors:       Matthias Althoff
% Written:       14-April-2022
% Last update:   11-September-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% this example only works if MATPOWER is installed (i.e., MATPOWER files are in the MATLAB path)
if logical(exist('test_matpower','file'))

    % list of MATPOWER models
    models_MATPOWER = {'case4gs', 'case5', 'case6ww', 'case9', 'case14', ...
        'case30', 'case57', ...
        'case10ba', 'case12da', 'case15da', 'case15nbr', ...
        'case17me', 'case18nbr', ...
        'case22', 'case28da', 'case33bw', 'case33mg', 'case34sa', ...
        'case38si', 'case51ga', 'case51he', 'case69', ...
        'case74ds', 'case85', 'case94pi', 'case118zh', 'case136ma', ...
        'case141'};

    % set path
    path = [CORAROOT filesep 'models' filesep 'powerSystems'];
    path_CORAmodel = [CORAROOT filesep 'models' filesep 'powerSystemsConverted'];
    path_spec = [CORAROOT filesep 'models' filesep 'powerSystems' filesep 'specs'];

    %% create power system models
    for i = 1:length(models_MATPOWER)
        %% obtain relevant model information
        % get model name
        model = models_MATPOWER{i}
        % load created power system description
        load([path filesep ['MATPOWER_',model]]);
        % convert to CORA
        if ~isfile([path_CORAmodel filesep ['MATPOWER_',model,'_model.mat']])
            powerSystem2cora(['MATPOWER_',model]);       
        end
        % load CORA model
        load([path_CORAmodel filesep ['MATPOWER_',model,'_model']], ['MATPOWER_',model,'_model']);
        % rename system to sys
        eval(['sys = ',['MATPOWER_',model,'_model;']]);
        % obtain parameters
        nrOfGenerators = length(bus.generator);
        nrOfLoads = length(bus.load);
        nrOfBuses = nrOfGenerators + nrOfLoads;
        nrOfInputs = sys.nrOfInputs;

        %% Data for spec MA-2014-x-1 -------------------------------------------
        % input
        u0_load = zeros(nrOfBuses,1);
        if nrOfGenerators >= 2
            u0_gen = [2.6; 0.4; zeros(nrOfGenerators-2,1)];
        else
            u0_gen = 2.6;
        end
        params.u = [u0_gen; u0_load];
        % guesses for dynamic and algebraic state variables
        x0guess = [0*ones(nrOfGenerators,1); 377*ones(nrOfGenerators,1); 1*ones(nrOfGenerators,1)];
        y0guess = [ones(nrOfBuses, 1); zeros(nrOfBuses-1, 1)];
        % does Jacobian exist for computing the steady state?
        try
            load(['jacobian_',['MATPOWER_',model]]);
        catch
            % compute Jacobian
            options.tensorOrder = 1;
            derivatives(sys, options);
        end
        % compute steady state solution
        [x0,params.y0guess] = steadyState(sys, x0guess, y0guess, params.u);
        % create initial set
        R0 = zonotope([x0, diag([5e-3*ones(nrOfGenerators,1); ...
         1e-1*ones(nrOfGenerators,1); 1e-3*ones(nrOfGenerators,1)])]);
        % create input set (zero input)
        U = zonotope(zeros(nrOfInputs,1));
        % time of fault occurance
        t_o = 0.1;
        % time of fault clearance
        t_c = 0.13;
        % name of spec
        name = ['MA2014_',num2str(nrOfBuses),'_1'];
        % save spec MA-2014-x-1
        eval([name,'.R0 = R0;']);
        eval([name,'.U = U;']);
        eval([name,'.t_o = t_o;']);
        eval([name,'.t_c = t_c;']);
        eval([name,'.t_f = [];']);
        eval([name,'.S = [];']);
        eval([name,'.X = [];']);

        eval(['save(''',path_spec,filesep,name,''',''',name,''');']);


        % Data for spec MA-2014-x-2 -------------------------------------------
        % set of acceptable states
        X = zonotope([x0, diag([pi*ones(nrOfGenerators,1); ...
         0.5*pi*ones(nrOfGenerators,1); 0.03*ones(nrOfGenerators,1)])]);

        % input set
        if nrOfBuses >= 14
            nonzeroInd = [13,14];
        else
            nonzeroInd = [nrOfBuses-1, nrOfBuses];
        end
        nonzeroInd = nonzeroInd + nrOfGenerators; % P_g^d, see [1]
        U = [zeros(nonzeroInd(1)-1,2); eye(2); zeros(nrOfInputs-nonzeroInd(2),2)]*zonotope(interval([0.04;0.04],[0.06;0.06]));

        % time horizon
        t_f = 5;

        % name of spec
        name = ['MA_2014_',num2str(nrOfBuses),'_2'];

        % save spec MA-2014-x-2
        eval([name,'.R0 = R0;']);
        eval([name,'.U_base = U;']);
        eval([name,'.U = sym(''t'')/t_f*sym(''U_base'');']);
        eval([name,'.t_o = [];']);
        eval([name,'.t_c = [];']);
        eval([name,'.t_f = t_f;']);
        eval([name,'.S = [];']);
        eval([name,'.X = X;']);

        eval(['save(''',path_spec,filesep,name,''',''',name,''');']);
    end

    % Spec for region of attraction MA2022_1_1
    x0 = [asin(1/5); 0];
    MA2022_1_1.R0 = zonotope([x0, 0.7*eye(2)]);
    MA2022_1_1.U = zonotope(0);
    MA2022_1_1.t_o = [];
    MA2022_1_1.t_c = [];
    MA2022_1_1.t_f = [];
    MA2022_1_1.S = ellipsoid([0.0201 -0.0249; -0.0249 4.9999],x0);
    MA2022_1_1.X = [];

    save([path_spec,filesep,'MA_2022_1_1'],'MA_2022_1_1');
else
    % write that MATPOWER is required
    disp('To run this example, please make sure that MATPOWER is installed.');
end
    
% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
