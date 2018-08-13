function [Rguard,tmin,tmax] = guardProjection(obj,hSpace,A,B,t_hit,tmin,tmax,loc,partialIntersectionFlag,backwards,options)
% hyperplaneMap - computes the reachable set of the system within a 
% location, detects the guard set that is hit and computes the new set on
% the hyperplane and the subsequent mapping.
%
% Syntax:  
%    [TP,R,activeGuards,Rjump,Rcont] =
%    hyperplaneMap(obj,tStart,R0,options)
%
% Inputs:
%    obj - location object
%    tStart - start time
%    R0 - initial reachable set
%    options - options struct
%
% Outputs:
%    TP - time point struct; e.g. contains time vector of minimum times for reaching guard sets
%    R - cell array of reachable sets
%    activeGuards - active guards
%    Rjump - reachable set after jump according to the reset map
%    Rcont - reachable set due to continuous evolution without guard or
%    invariant inmtersection
%
% Example: 
%
% Other m-files required: initReach, reach, potOut, potInt, guardIntersect,
% reset
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      23-August-2013
% Last update:  26-August-2013
%               28-August-2013
%               06-September-2013
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------


%separating halspace
c_new = - (hSpace.c.'*A).';
d_new = hSpace.c.'*B*options.u;
h_new = halfspace(c_new, d_new);

%check intersection time with guard set and separating halfspace
contDynamics = get(loc,'contDynamics');
[tmin,tmax] = intersectionTime_oneDirection(hSpace,obj,contDynamics,options);
if backwards
    tmax_sep = [];
else
    [tmin_sep,tmax_sep] = intersectionTime_oneDirection(h_new,obj,contDynamics,options);
end

%partial intersection?
if isempty(tmax_sep) || (tmax < tmax_sep)

    %get data
    x0 = center(obj);
    u = B*options.u;

    %compute f
    dim = length(A);
    I = eye(dim);
    %Theta
    Theta = A^2*t_hit/2;
    for i=3:options.taylorTerms
        Theta = Theta + A^i*t_hit^(i-1)/factorial(i);
    end
    %Gamma
    Gamma = I;
    for i=1:(options.taylorTerms-1)
        Gamma = Gamma + A^i*t_hit^i/factorial(i+1);
    end
    f = Theta*x0 + Gamma*u;
    flow = A*x0 + f;
    f_cell{1} = flow;

    %check flow direction
    limitRespected = 0;
    flowChanged = 0;
    while ~limitRespected
        flowDir_IH = interval(hSpace.c'*(A*obj + flow)*(1/(norm(hSpace.c)*norm(A*x0 + flow))));
        flowDir_IH_abs = min(abs(flowDir_IH(1,1)), abs(flowDir_IH(1,2)));
        %flowDir_min = 1.0;
        %flowDir_min = 0.8;
        %flowDir_min = 0.4;
        flowDir_min = 0.2;
        %flowDir_min = 0;
        if (flowDir_IH(1,1)*flowDir_IH(1,2) > 0) & (flowDir_IH_abs < 0.8*flowDir_min)
    %         %normal vector of separating halfspace
    %         n_sep = (hSpace.c.'*A).';
    %         %extract component that is perpendicular to normal vector of
    %         %original halfspace
    %         n_new = n_sep - n_sep.'*hSpace.c/norm(hSpace.c)*hSpace.c;
            %strecthing factor for flow correction
            if backwards
                alpha = (flowDir_min - flowDir_IH_abs)*norm(hSpace.c)*norm(A*x0 + flow)/(hSpace.c'*hSpace.c);
            else
                alpha = -(flowDir_min - flowDir_IH_abs)*norm(hSpace.c)*norm(A*x0 + flow)/(hSpace.c'*hSpace.c);
            end
            %flow correction
            f = f + alpha*hSpace.c;
            flow = flow + alpha*hSpace.c;
            %set flow changed
            flowChanged = 1;
        else
            limitRespected = 1;
        end
    end
    f_cell{2} = flow;

    %use reduced zonotope 
    Rred = reduce(obj,'girard',options.errorOrder);
    
    %compute trajectory error
    [c_tilde, d_tilde, E_err] = trajectoryError(A,u,obj,tmin,tmax,t_hit,options.taylorTerms,hSpace.c);

    %compute terms of Taylor series for non-error part
    [k2,L2,Q2,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(hSpace, A, f, x0, Rred);
    
    %compute terms of Taylor series for error part
    c_tilde_mid = mid(c_tilde);
    d_tilde_mid = mid(d_tilde);
    c_tilde_rad = rad(c_tilde);
    d_tilde_rad = rad(d_tilde);
    [k2_mid,L2_mid,Q2_mid,lambdaQuad_inv_zono_mid,Lambda_int_mid] = mapsQuad(hSpace, A, f, x0, Rred, c_tilde_mid, d_tilde_mid);
    [k2_rad,L2_rad,Q2_rad,lambdaQuad_inv_zono_rad,Lambda_int_rad] = mapsQuad(hSpace, A, f, x0, Rred, c_tilde_rad, d_tilde_rad);
    k2_err = zonotope([k2_mid,abs(k2_rad)]);
    L2_err = interval(L2_mid - abs(L2_rad), L2_mid + abs(L2_rad));
    for iDim = 1:length(Q2_mid)
        Q2_err{iDim} = concatenate(Q2_mid{iDim},Q2_rad{iDim});
    end
    
    
    %individual reachable sets 
    dR = obj + (-x0);
    dRred = Rred + (-x0);
    
    %linear part
    Rlin = (I+L2)*dR;
    Rlin_err = L2_err*dR;

    %quadratic part
    Rquad_zono_aux = quadraticMultiplication_zono(dRred,Q2);
    Rquad_zono = 1/2*lambdaQuad_inv_zono*Rquad_zono_aux;
    
    Rquad_zono_aux_err = quadraticMultiplication_zono(dRred,Q2_err);
    Rquad_zono_err = 1/2*lambdaQuad_inv_zono*Rquad_zono_aux_err;
    
    %Rquad_zono = 1/2*lambdaQuad_inv_zono*reduce(Rquad_zono_aux,'girard',options.zonotopeOrder);

    %combine results
    Rhit_const = x0 + k2 + exactPlus(Rlin, Rquad_zono);
    %Rhit_err = k2_err + Rlin_err + Rquad_zono_err + E_err;
    Rhit_err = k2_err + Rlin_err_comb + Rquad_zono_err + E_err;

    %for test
    Zk = k2 + zonotope([zeros(dim,1), 0.01*eye(dim)]);

    %compute mapping error
    [delta_Rhit,E] = mappingError(A,f_cell,u,obj,tmin,tmax,t_hit,options.taylorTerms,Lambda_int,hSpace.c);
    %delta_Rhit2 = mappingError_split(A,f,u,obj,tmin,tmax,t_hit,options.taylorTerms,Lambda_int,hSpace.c,0.1*options.maxTimeStep);
    
        
    %compute set on hyperplane
    %Rguard = Rhit_const + delta_Rhit;
    Rguard = Rhit_const + Rhit_err;
    
    Rguard = reduce(Rguard,'girard',options.zonotopeOrder);
    
%     %test of constant flow
%     figure
%     Rred = reduce(Rhit_const,'girard',10);
%     plot(Rred,[1 2],'g');
%     plot(obj,[1 2],'g');
%     for i = 1:100
%         if i < 50
%             x0 = randPointExtreme(obj);
%         else
%             x0 = randPoint(obj);
%         end
%         flowVec = A*x0 + f;
%         xNew = x0 + tmax*flowVec;
%         plot([x0(1),xNew(1)],[x0(2),xNew(2)],'r');
%     end

    
%     %test of exact flow
%     figure
%     Rred = reduce(Rguard,'girard',10);
%     plot(Rred,[1 2],'g');
%     plot(obj,[1 2],'g');
%     
%     %get contDynamics
%     contDynamics = get(loc,'contDynamics');
%     
%     %obtain data from location
%     inv = interval(1e6*ones(length(options.x0),1)*[-1 1]);
%     tran_old = get(loc, 'transition');
%     reset = get(tran_old{1}, 'reset');
%     tran{1}=transition(hSpace,reset,0,'a','b'); 
%     %create dummy location
%     loc = location('dummyLoc',inv,tran,contDynamics);  
%     eventOptions = odeset('Events',eventFcn(loc));
%     opt = odeset(eventOptions);
%     
%     for i = 1:100
%         if i < 50
%             options.x0 = randPointExtreme(obj);
%         else
%             options.x0 = randPoint(obj);
%         end
%         
%         %simulate continuous dynamics
%         [contDynamics,t,x] = simulate(contDynamics,options,tmin,tmax,options.x0,opt);
%         
%         plot(x(:,1),x(:,2),'r');
%     end


%     %test delta_Rhit
%     figure
%     hold on
%     
%     
%     %get contDynamics
%     contDynamics = get(loc,'contDynamics');
%     
%     %simulate
%     stepsizeOptions = odeset('MaxStep',0.05*(options.tStart-options.tFinal));
%     %generate overall options
%     opt = odeset(stepsizeOptions);
%     
%     for i = 1:3
%         x0 = randPointExtreme(obj);
%         
%         %simulate continuous dynamics
%         %[contDynamics,t,x] = simulate(contDynamics,options,tmin,tmax,x0,opt);
%         [contDynamics,t,x] = simulate(contDynamics,options,0,tmax,x0,opt);
%         
%         %xNew = x0 + (tmax-tmin)*flowVec;
%         xNew = x0 + (tmax)*(A*x0 + f);
%         plot([x0(1),xNew(1)],[x0(2),xNew(2)],'g');
%         plot(x(:,1),x(:,2),'r');
%         
%         Rred = reduce(E,'girard',10);
%         Rred_err = reduce(E_err,'girard',10);
%         %Rred2 = reduce(delta_Rhit2,'girard',10);
%         Rtest = Rred + xNew;
%         Rtest_err = Rred_err + xNew;
%         %Rtest2 = Rred2 + xNew;
%         plot(Rtest);
%         plot(Rtest_err,[1 2],'k');
%         
%         %compute delta_x
%         delta_aux = 0*x0;
%         for i = 2:options.taylorTerms
%             delta_aux = delta_aux + A^i/factorial(i)*(tmax^(i-1)*x0 - t_hit^(i-1)*x0);
%         end
%         for i = 0:options.taylorTerms
%             delta_aux = delta_aux + A^i/factorial(i+1)*(tmax^(i)- t_hit^(i))*u;
%         end
%         delta_x = tmax*delta_aux;
%         
%         x_corr = xNew + delta_x;
%         plot(x_corr(1), x_corr(2), 'or');
%     end
else
    Rguard = partialGuardProjection(obj,hSpace,A,B,loc,options);
end



%function which determines when the state of x'=Ax+u crosses the hypeplane
%c^T x = d with initial state x0
function t_hit = crossingTime(hSpace,x0,A,u,tc,tmin,tmax)

t_hit = tc;
t_lower = tmin;
t_upper = tmax;
t_old = inf;

while (abs(t_hit - t_old) > 1e-6*tmax)
    
    %remember time
    t_old = t_hit;
    
    %compute state solution
    x_upper = stateSolution(x0,A,u,0.5*(t_hit + t_upper));
    x_lower = stateSolution(x0,A,u,0.5*(t_hit + t_lower));

    %check if time is too large or too small
    if (abs(hSpace.c'*x_upper - hSpace.d) < abs(hSpace.c'*x_lower - hSpace.d))
        t_lower = t_hit;
        t_hit = 0.5*(t_hit + t_upper);
    else
        t_upper = t_hit;
        t_hit = 0.5*(t_hit + t_lower);
    end

end


%computes the solution for x0 under x'=Ax+u after time t
function x = stateSolution(x0,A,u,t)

%Gamma matrix for input solution using 10 Taylor terms
Gamma = eye(length(A))*t;
for i=1:10
    Gamma = Gamma + 1/factorial(i+1)*A^i*t^(i+1);
end

x = expm(A*t)*x0 + Gamma*u;

function Rcubic = cubicIntMultiplication(dR,Cint)

%obtain interval of dR
dx = interval(dR);

%dim
dim = length(Cint);

for i = 1:dim
    int(i,1) = interval(0,0);
    for l = 1:dim
        for m = 1:dim
            for n = 1:dim
                int(i,1) = int(i,1) + Cint{i}(l,m,n)*dx(l)*dx(m)*dx(n);
            end
        end
    end
end

%Rcubic
Rcubic = zonotope(int);


function Rcubic = cubicIntMultiplication_zono(dR,C_c,C_g)

%obtain interval of dR
Zmat = get(dR,'Z');

%dim
dim = length(C_c);
gens = length(C_g(1,:));

%init Znew
Znew = zeros(dim,(gens+1)^2);

%center computation
for i = 1:dim
    for l = 1:dim
        for m = 1:dim
            for n = 1:dim
                %compute center 
                for iGen=1:length(Zmat(1,:))
                    Znew(i,iGen) = Znew(i,iGen) + C_c{i}(l,m,n)*Zmat(l,iGen)*Zmat(m,iGen)*Zmat(n,iGen);
                end
            end
        end
    end
end

%generator computation
for iMatGen = 1:gens
    for i = 1:dim
        for l = 1:dim
            for m = 1:dim
                for n = 1:dim
                    %compute center 
                    for iGen=1:length(Zmat(1,:))
                        ind = (gens+1)*iMatGen + iGen;
                        Znew(i,ind) = Znew(i,ind) + C_g{i,iMatGen}(l,m,n)*Zmat(l,iGen)*Zmat(m,iGen)*Zmat(n,iGen);
                    end
                end
            end
        end
    end
end

%Rcubic
Rcubic = zonotope(Znew);



function delta_Rhit = splittedMappingError(A,f,u,R0,tmin,tmax,tc,order,Lambda_int,n,tmaxStep)

%compute time span
timeSpan = tmax - tmin;

%obtain necessary subintervals
subIntervals = ceil(timeSpan/tmaxStep);

%obtain time intervals
delta_t = timeSpan/subIntervals;

for i = 1:subIntervals
    %obtain new lower and upper bounds
    tmin_new = tmin + (i-1)*delta_t;
    tmax_new = tmin + i*delta_t;
    
    delta_Rhit_partial = mappingError(A,f,u,R0,tmin,tmax,tc,order,Lambda_int,n);

end


function delta_Rhit = mappingError_split(A,f,u,R0,tmin,tmax,tc,order,Lambda_int,n,tmaxStep)

%obtain center and delta set
x0 = center(R0);

%compute time span
timeSpan = tmax - tmin;

%obtain necessary subintervals
subIntervals = ceil(timeSpan/tmaxStep);

%obtain time intervals
delta_t = timeSpan/subIntervals;

for iTimeInterval = 1:subIntervals
    %obtain new lower and upper bounds
    tmin_new = tmin + (iTimeInterval-1)*delta_t;
    tmax_new = tmin + iTimeInterval*delta_t;
    
    %compute powers of tmin and tmax
    tZono = timePowers(tmin_new,tmax_new,order);
 
    %auxiliary set
    Raux_1 = zonotope(0*x0);
    Raux_2 = zonotope(0*x0);
    for i=2:order
        Raux_1 = Raux_1 + A^i/factorial(i)*(tZono{i-1}*R0 + (-tc^(i-1)*x0));
        
    end
    for i=1:order
        Raux_2 = Raux_2 + A^i/factorial(i+1)*(tZono{i}+(-tc^i))*u;
    end
    
    if iTimeInterval==1
        R_union = interval(Raux_1 + Raux_2);
    else
        R_union = R_union | interval(Raux_1 + Raux_2);
    end
end
% convert to zonotope
R_union = zonotope(R_union);

%compute time powers of full time interval
tZono = timePowers(tmin,tmax,order);

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:order
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
%compute absolute value of W for numerical stability
W = abs(W);
E = interval(-W,W);

%final additional set
G = tZono{1}*(R_union) + E*R0 + E*tmax*u;

%error due to uncertain hitting time
delta_t = interval((-n')*G)/Lambda_int;
delta_t_mid = mid(delta_t);
delta_t_rad{1} = rad(delta_t);
delta_t_zono = matZonotope(delta_t_mid,delta_t_rad);

%final uncertainty
delta_Rhit = G + delta_t_zono*(A*R0 + f);



function [delta_Rhit,G] = mappingError(A,f,u,R0,tmin,tmax,tc,order,Lambda_int,n)

%compute powers of tmin and tmax
tZono = timePowers(tmin,tmax,order);

%obtain center and delta set
x0 = center(R0);
x0zono = zonotope(x0);
uzono = zonotope(u);
dR = R0 + (-x0);

%auxiliary set: new
Raux_1 = zonotope(0*x0);
Raux_2 = zonotope(0*x0);
for i=2:order
    %Raux_1 = Raux_1 + A^i/factorial(i)*(tZono{i-1}*dR + (tZono{i-1} + (-tc^(i-1)))*x0zono);
    Raux_1 = Raux_1 + A^i/factorial(i)*(tmax^(i-1)*dR + (tZono{i-1} + (-tc^(i-1)))*x0zono);
end
for i=1:order
    Raux_2 = Raux_2 + A^i/factorial(i+1)*(tZono{i} + (-tc^i))*uzono;
end

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:order
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
%compute absolute value of W for numerical stability
W = abs(W);
E = interval(-W,W);

%final additional set
G = tZono{1}*(Raux_1 + Raux_2 + zonotope(f{1}-f{2})) + E*R0 + E*tmax*uzono;

%error due to uncertain hitting time
delta_t = interval((-n')*G)/Lambda_int;
delta_t_mid = mid(delta_t);
delta_t_rad{1} = rad(delta_t);
delta_t_zono = matZonotope(delta_t_mid,delta_t_rad);

%final uncertainty
delta_Rhit = G + delta_t_zono*(A*R0 + f{2});



function [c_tilde, d_tilde, E_err] = trajectoryError(A,u,R0,tmin,tmax,tc,order,n)

%obtain dimension
dim = length(u);

%compute powers of tmin and tmax
tZono = timePowers(tmin,tmax,order);

%obtain center and delta set
x0 = center(R0);
x0zono = zonotope(x0);
uzono = zonotope(u);
dR = R0 + (-x0);

%auxiliary matrix
A_tilde = zeros(dim);
for i=2:order
    A_tilde = A_tilde + A^i/factorial(i)*tZono{i};
end

%auxiliary vector
u_tilde = zeros(dim,1);
for i=2:order
    u_tilde = u_tilde + A^i/factorial(i)*(tZono{1} * (-tc^(i-1)))*x0zono;
end
for i=1:order
    u_tilde = u_tilde + A^i/factorial(i+1)*(tZono{i+1} + tZono{1}*(-tc^i))*uzono;
end

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:order
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
%compute absolute value of W for numerical stability
W = abs(W);
E = interval(-W,W);

%consider E
gen{1} = rad(E);
A_tilde = A_tilde + matZonotope(mid(E), gen);
u_tilde = u_tilde + E*tmax*uzono;

%test
E_err = A_tilde*R0 + u_tilde;

%obtain c_tilde and d_tilde
c_tilde = A_tilde.'*n;
c_tilde = intervalMatrix(c_tilde);
c_tilde = c_tilde.int;

d_tilde = -n.'*u_tilde;
d_tilde = interval(d_tilde);




function tZono = timePowers(tmin,tmax,order)


%first order
tminPow(1) = tmin;
tmaxPow(1) = tmax;

tradPow(1) = 0.5*(tmaxPow(1) - tminPow(1));
tmidPow(1) = 0.5*(tmaxPow(1) + tminPow(1));

tRad{1} = tradPow(1);
tZono{1} = matZonotope(tmidPow(1),tRad);

%higher orders
for i=1:order
    tminPow(i+1) = tminPow(i)*tmin;
    tmaxPow(i+1) = tmaxPow(i)*tmax;
    
    tradPow(i+1) = 0.5*(tmaxPow(i+1) - tminPow(i+1));
    tmidPow(i+1) = 0.5*(tmaxPow(i+1) + tminPow(i+1));
    
    tRad{1} = tradPow(i+1);
    tZono{i+1} = matZonotope(tmidPow(i+1),tRad);
end




%------------- END OF CODE --------------