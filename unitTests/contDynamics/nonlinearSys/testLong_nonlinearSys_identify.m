function res = testLong_nonlinearSys_identify
% testLong_nonlinearSys_identify - unit test for system identification
%
% Syntax:
%    res = testLong_nonlinearSys_identify
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       13-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % define algorithms to test
    alg = {'dmd','opt','sparse','nested'}; 

    % define template functions to test
    basis = {'quad','cubic','trig','all'};

    % loop over all possible 
    for i = 1:length(alg)
        %for j = 1:length(basis)

            % set algorithm options
            options = [];
            options.alg = alg{i};

            if ~strcmp(options.alg,'nested')
                %options.basis = basis{i};
                options.phi = @(x,u) x;
            end

            % generate random test data
            dt = rand(); n = randi([2,3]); m = randi([1,3]); 
            N1 = randi([3,5]); N2 = randi([3,5]); s = randi([1,10]);

            t1 = (0:N1-1)*dt; t2 = (0:N2-1)*dt; t = {t1;t2};
            x1 = s*(rand(n,N1)-0.5); x2 = s*(rand(n,N2)-0.5); x = {x1;x2};
            u1 = s*(rand(m,N1)-0.5); u2 = s*(rand(m,N2)-0.5); u = {u1;u2};

            % Test 1: single trajectory without input
            sys = nonlinearSys.identify(x1,t1,options);
            assert(isa(sys,'nonlinearSys'));

            % Test 2: single trajectory with input
            sys = nonlinearSys.identify(x1,t1,u1,options);
            assert(isa(sys,'nonlinearSys'));

            % Test 3: multiple trajectories without input
            sys = nonlinearSys.identify(x,t,options);
            assert(isa(sys,'nonlinearSys'));

            % Test 4: multiple trajectories with input
            sys = nonlinearSys.identify(x,t,u,options);
            assert(isa(sys,'nonlinearSys'));

%             if strcmp(options.alg,'nested')
%                 break;
%             end
        %end
    end

    % test successfull if it runs through without errors
    res = true;

% ------------------------------ END OF CODE ------------------------------
