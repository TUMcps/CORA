function [X_full,zTraj] = simulate_rrt(obj,options)
% simulate_rrt - simulates a DAE system using rapidly exploring random trees
% (partially deterministic sampling)
%
% Syntax:  
%    [obj,t,x,index] = simulate_rrt_det(obj,options)
%
% Inputs:
%    obj - DAE object
%    options - options struct
%
% Outputs:
%    obj - linearSys object
%    t - time vector
%    x - state vector
%    index - returns the event which has been detected
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-December-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%nrOfSamples
nrOfSamples = 10;

%initialize
X_sample_size = 2*rad(interval(options.X_sample));
normMatrix = diag(1./X_sample_size);
xCenter = options.x0;

%possible extreme inputs
U_input_mat = get(vertices(options.U),'V');
nrOfExtrInputs = length(U_input_mat(1,:));

%obtain vertices 
V = vertices(options.R0);
Vmat = get(V,'V');
X_full = Vmat;
Y_full = options.y0*ones(1,length(X_full(1,:)));

%timeSteps
timeSteps = ceil(options.tFinal/options.timeStep);

%init
xTraj = cell(nrOfSamples,timeSteps);

for iStep = 1:timeSteps
    
    iStep
    
    for iSample = 1:nrOfSamples
        %iSample
        
        %init
        x_full = cell(1,nrOfExtrInputs);
        x_end = [];     

        %sample
        x_sample = randPoint(options.X_sample) + xCenter;

        %nearest neighbor and selected state
        ind = nearestNeighbor(x_sample,X_full,xCenter,normMatrix);
        ind
        x0 = X_full(:,ind);
        y0 = Y_full(:,ind);

        %simulate model to find out best input
        for iInput = 1:nrOfExtrInputs
            %set input
            par_options.u = options.uTrans + U_input_mat(:,iInput);
            %[x0_ss,y0] = initialState(obj, x0, y0, par_options.u);
            %simulate
            z{iInput} = par_simulate(obj,par_options,options.timeStep,x0,y0);
            %convert to states
            x_end(:,iInput) = z{iInput}(end,1:obj.dim);   
            y_end(:,iInput) = z{iInput}(end,obj.dim+1:obj.dim + obj.nrOfConstraints); 
        end

        %nearest neighbor and selected state
        ind_input = nearestNeighbor(x_sample,x_end,xCenter,normMatrix);
        ind_input

        %add selected state 
        X_next(:,iSample) = x_end(:,ind_input);
        Y_next(:,iSample) = y_end(:,ind_input);

        %add trajectory
        zTraj{iSample,iStep} = z{ind_input}';
    end
    %update variables
    X_full = X_next;
    Y_full = Y_next;
end



function ind = nearestNeighbor(x_sample,X,c,normMatrix)

%normalize set of points
X_delta = X - c*ones(1,length(X(1,:)));
X_norm = normMatrix*X_delta;

%normalize sample point
x_delta = x_sample - c;
x_norm = normMatrix*x_delta;

%norm of distance
X_rel = X_norm - x_norm*ones(1,length(X_norm(1,:)));
norm_val = vnorm(X_rel,1,2); %compute 2-norm

%find index with smallest norm
[val,ind] = min(norm_val);



%simulate
function z_full = par_simulate(obj,par_options,timeStep,x0,y0)

[dummy_a,dummy_b,z_full] = simulate(obj,par_options,0,timeStep,x0,y0,[]);



%steady state solution
function [x0,y0] = initialState(obj, x0, y0, u0)

converged = 0;

while ~converged
    k = obj.dynFile(0, x0, y0, u0);
    l = obj.conFile(0, x0, y0, u0);
    [A,B,C,D,E,F] = jacobian(x0, y0, u0);

    %new steady state solution
    M = [A C; D F];
    delta_z = M\(-[k;l]);

    %update steady state solution
    delta_y = delta_z((length(x0)+1) : end);
    
    %check convergence
    if norm(delta_y)<1e-10
        converged = 1;
    end
   
    %update y0
    y0 = y0 + delta_y;
end


%------------- END OF CODE --------------
