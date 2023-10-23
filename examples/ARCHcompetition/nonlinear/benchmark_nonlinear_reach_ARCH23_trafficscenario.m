function text = benchmark_nonlinear_reach_ARCH23_trafficscenario
% benchmark_nonlinear_reach_ARCH23_trafficscenario - example of nonlinear
%    reachability analysis.
%
% Syntax:
%    benchmark_nonlinear_reach_ARCH23_trafficscenario()
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%

% Authors:       Mark Wetzlinger
% Written:       06-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% precompute symbolic derivatives required for the computation of the

% occupancy set
tay = aux_taylorOccupancySet();

% loop over all traffic scenarios (= verification problems)
path = fileparts(which(mfilename));
filename = 'BEL_Putte-4_2_T-1_controls.csv';

[x0,u,K,t] = aux_readSolutionFile(filename);

% compute reachable set
clock = tic;
R = aux_compReachSet(K,x0,u,t);

% compute occupancy set
O = aux_compOccupancySet(R,tay);

tComp = toc(clock);

% check if input constraint are satisfied
res = aux_checkInputConstraints(R,u,K);

% export occupancy set
if res
    name = strrep(filename,'controls','occupancies');
    aux_exportOccupancySet(name,O);
else
    list{1} = strrep(filename,'_controls.csv','');
end

% display results
disp("");
if res
    disp('Input constraints satisfied!');
else
    disp('Input constraints violated!');
end
disp("Computation time: " + tComp);
disp("");

text = ['TRAF22,',num2str(res),',',num2str(tComp),', , '];

end


% Auxiliary functions -----------------------------------------------------

function R = aux_compReachSet(K,x0,u,t)
% compute the reachable set of the controlled car

    % number of reachability steps
    N = 20;

    % set of disturbances
    W = interval([-0.02; -0.3],[0.02; 0.3]);
    
    % set of measurement errors         
    V = interval([-0.0004; -0.0004; -0.006; -0.002; -0.002], ...
                 [0.0004; 0.0004; 0.006; 0.002; 0.002]);
    
    % parameter for reachability problem
    params.R0 = cartProd(zonotope(x0) + zonotope(V),zonotope(x0));
    params.U = cartProd(zonotope(W),zonotope(V));
    
    % reachability settings
    options.alg = 'lin';
    options.tensorOrder = 2;
    options.zonotopeOrder = 20;
    options.taylorTerms = 10;
    options.intermediateTerms = 4;
    
    options.lagrangeRem.simplify = 'optimize';
    
    % system dynamics
    sys = nonlinParamSys('trafficscenario_ARCH23',@aux_controlledCar,10,7,12);

    % reachability analysis
    R = [];
    
    for i = 1:size(u,2)
        
        % update parameters
        params.tStart = t(i);
        params.tFinal = t(i+1);
        options.timeStep = (params.tFinal-params.tStart)/N;
        params.paramInt = [reshape(K{i},[10,1]);u(:,i)];
        
        % compute reachable set
        Rtemp = reach(sys, params, options);
        R = add(R,Rtemp);
        
        % update initial set
        params.R0 = Rtemp.timePoint.set{end};
    end
end

function dx = aux_controlledCar(x,w,p)
% differential equation for the closed-loop system

    % extended state: x_ext = [x; x_ref]
    x_ref = x(6:10);
    x = x(1:5);

    % get feedback matrix K and reference input
    K = reshape(p(1:10),[2,5]);
    u_ref = p(11:12);

    % control law
    u = u_ref + K * (x - x_ref + w(3:end));

    % differential equation of controlled system
    dx1 = aux_car(x,w(1:2),u);
    dx2 = aux_car(x_ref,zeros(2,1),u_ref);

    dx = [dx1; dx2];
end

function dx = aux_car(x,w,u)
% differential equation for the open-loop system

    % wheelbase parameter
    l_wb = 2.578;

    % system dynamics
    dx(1, 1) = u(1) + w(1);
    dx(2, 1) = (x(3) / l_wb) * tan(x(1));
    dx(3, 1) = u(2) + w(2);
    dx(4, 1) = x(3) * cos(x(2));
    dx(5, 1) = x(3) * sin(x(2));
end

function [x0,u,K,t] = aux_readSolutionFile(path)
% read-in initial point, the control inputs, the feedback matrices, and the 
% time vector from the solution file
    
    % read initial state x0
    fid = fopen(path);
    data = fgetl(fid);
    temp = split(data);
    x0 = zeros(6,1);
    for i = 1:length(temp)
       x0(i) = str2double(strrep(temp{i},';','')); 
    end
    t0 = x0(1);
    x0 = x0(2:end);
    fclose(fid);

    % read control inputs u and time vector t
    data = readtable(path);
    t = data{:,1};
    t = [t0;t];
    u = data{:,2:3};
    u = u';

    % read feedback matrices
    K = cell(length(t)-1,1);
    
    for i = 1:length(K)
       temp = data{i,4:end};
       K{i} = [temp(1:5);temp(6:10)];
    end
end

function O = aux_compOccupancySet(R,tay)
% compute the occupancy set by adding the car dimensions to the reachable
% set

    % car length and width in [m]
    l = 4.508;
    w = 1.610;
    l_wb = 2.578;
    
    vehicle = zonotope(interval([-l/2+l_wb/2;-w/2],[l/2+l_wb/2;w/2])); 
    
    % loop over all reachable sets
    O = {};
    
    for j = 1:length(R)
        for i = 1:length(R(j).timeInterval.set)

            % define initial set
            set = project(R(j).timeInterval.set{i},[4 5 2]);
            set = reduce(set,'girard',8);
            
            X = cartProd(set,vehicle);

            % evaluate derivatives at linearization point
            p = center(X);

            [f,A,Q,T] = aux_evalDerivatives(X,p,tay);

            % compute Largrange remainder
            rem = aux_lagrangeRemainder(X,p,T);

            % compute over-approximating zonotope
            res = f + A * (X + (-p)) + 0.5*quadMap((X + (-p)),Q) + rem;

            % convert to polygon object
            if ~isa(res,'zonotope')
               res = zonotope(res); 
            end
            V = vertices(res);

            w = warning();
            warning('off');

            O{end+1}.set = polygon(V(1,:)',V(2,:)');
            O{end}.time = R(j).timeInterval.time{i};

            warning(w);
        end
    end
end

function tay = aux_taylorOccupancySet()
% compute the symbolic derivatives of the Taylor expansion of the nonlinear
% function used to compute the occupancy set

    % define transformation function
    syms x1 x2 phi z1 z2
    x = [x1;x2;phi;z1;z2];
    
    f = [x1 + cos(phi)*z1 - sin(phi)*z2; x2 + cos(phi)*z2 + sin(phi)*z1];

    % function handle for the nonlinear function
    fun =  matlabFunction(f,'Vars',{x});
    
    % first-order derivative
    A = jacobian(f,x);
    Afun =  matlabFunction(A,'Vars',{x});
    
    % second order derivative
    Qfun = cell(length(f),1);
    for i = 1:length(f)
       temp = hessian(f(i),x); 
       Qfun{i} =  matlabFunction(temp,'Vars',{x});
    end
    
    % Lagrange remainder
    Tfun = cell(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            temp = hessian(A(i,j),x);
            if any(any(temp ~= 0))
                Tfun{i,j} = matlabFunction(temp,'Vars',{x});
            end
        end
    end
    
    % store function handles
    tay.fun = fun; tay.Afun = Afun; tay.Qfun = Qfun; tay.Tfun = Tfun;
end

function [f,A,Q,T] = aux_evalDerivatives(X,p,tay)
% evaluate the derivatives at the linearization point

    % interval encluore of the set
    int = interval(X);
    
    f = tay.fun(p);
    A = tay.Afun(p);
    
    Q = cell(length(f),1);
    for i = 1:length(f)
       funHan = tay.Qfun{i};
       Q{i} = funHan(p);
    end
    
    T = cell(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if ~isempty(tay.Tfun{i,j})
                funHan = tay.Tfun{i,j};
                T{i,j} = funHan(int);
            end
        end
    end
end

function rem = aux_lagrangeRemainder(X,p,T)
% comptute the lagrange remainder of the Taylor series

    % interval enclousre of the shifted initial set
    int = interval(X) - p;

    % Lagrange remainder term
    rem = interval(zeros(size(T,1),1));
    
    for i = 1:size(T,1)
        for j = 1:size(T,2)
            if ~isempty(T{i,j})
                rem(i) = rem(i) + int(j) * transpose(int) * T{i,j} * int;
            end
        end
    end
    
    % convert to zonotope
    rem = zonotope(1/6*rem);
end

function aux_exportOccupancySet(fileSol,O)
% write the occupancy set to the specified file

    % open file
    fid = fopen(fileSol,'wt');
    
    % loop over all occupancy sets
    for i = 1:length(O)
       V = O{i}.set.set.Vertices;
       fprintf(fid,'%f;',infimum(O{i}.time)); 
       for j = 1:size(V,1)-1
           fprintf(fid,'%f;',V(j,1)); 
       end
       fprintf(fid,'%f\n',V(end,1)); 
       fprintf(fid,'%f;',supremum(O{i}.time));
       for j = 1:size(V,1)-1
           fprintf(fid,'%f;',V(j,2)); 
       end
       fprintf(fid,'%f\n',V(end,2)); 
    end

    fclose(fid);
end
 
function res = aux_checkInputConstraints(R,u,K)
% check if the input constraints are satisfied

    % set of admissible control inputs
    U = interval([-0.7;-11],[0.7;11]);
    
    % set of applied control input
    U_ = [];
    
    for i = 1:length(K)
        for j = 1:length(R(i).timeInterval.set)
           set = R(i).timeInterval.set{j};
           Utemp = u(:,i) + K{i}*zonotope(set.Z(1:5,:)-set.Z(6:10,:));
           U_ = U_ | interval(Utemp);
        end
    end
    
    % check if input constraints are satisfied
    res = contains(U,U_);
end


% ------------------------------ END OF CODE ------------------------------
