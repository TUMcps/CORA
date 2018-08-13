function [Rguard,Rguard_noInt,split,RbeforeInt] = guardProjection(obj,hSpace,A,B,loc,partialIntersectionFlag,options)
% guardProjection - computes the ptojection of a reachable set onto a
% halfspace
%
% Syntax:  
%    Rguard = guardProjection(obj,halfspace,A,B,t_hit,tmin,tmax,location,options)
%
% Inputs:
%    obj - zonotopeBundle
%    halfspace - halfspace object
%    A - system matrix 
%    B - input matrix 
%    t_hit - time when center trajectory hits the guard set
%    t_min - minimum time when the halfspace is hit
%    t_max - maximum time when the halfspace is hit
%    location - location object
%    options - options struct
%
% Outputs:
%    Rguard - guard intersection
%
% Example: 
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      23-August-2013
% Last update:  26-August-2013
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%init
Rguard_noInt = [];
split = 0;
RbeforeInt = [];

%partial intersection?
if ~partialIntersectionFlag
    
    %determine switching time interval
    [tmin,tmax,t_hit,t_total,RbeforeInt] = switchingTimeInterval(obj,loc,hSpace,options);
    
    %if there is an intersection
    if ~isempty(RbeforeInt)
        
        %only consider first zonotope
        RbeforeInt = RbeforeInt.Z{1};

        %get data
        x0 = center(RbeforeInt);
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

        %use reduced zonotope 
        Rred = reduce(RbeforeInt,'girard',options.errorOrder);

        % flow manipulation
        %[f, f_cell] = flowManipulation_works(A,RbeforeInt,x0,f,hSpace);
        %[f, f_cell] = flowManipulation_perpendicular(A,Rred,x0,f,hSpace);
        [f, f_cell] = flowManipulation_inDir(A,Rred,x0,f,hSpace);
        flow = f_cell{1};

        %compute terms of Taylor series for non-error part
        [k2,L2,Q2,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(hSpace, A, f, x0, Rred);

        %compute trajectory error
        [c_tilde, d_tilde, E_err] = trajectoryError(A,f_cell,u,RbeforeInt,tmin,tmax,t_hit,options.taylorTerms,hSpace.c);

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
        dR = RbeforeInt + (-x0);
        dRred = Rred + (-x0);

        %linear part
        Rlin = (I+L2)*dR;
        Rlin_err = L2_err*dR;


        %quadratic part
        Rquad_zono_aux = quadraticMultiplication_zono(dRred,Q2);
        Rquad_zono = 1/2*lambdaQuad_inv_zono*Rquad_zono_aux;

        Rquad_zono_aux_err = quadraticMultiplication_zono(dRred,Q2_err);
        Rquad_zono_err = 1/2*lambdaQuad_inv_zono*Rquad_zono_aux_err;


        %combine results
        Rhit_const = x0 + k2 + exactPlus(Rlin, Rquad_zono);
        Rhit_err = k2_err + Rlin_err + Rquad_zono_err + E_err;


        %compute set on hyperplane
        %Rtmp = Rhit_const + delta_Rhit;
        Rtmp = Rhit_const + Rhit_err;

        %reduce result
        Rtmp = reduce(Rtmp,'girard',options.zonotopeOrder);

        %project result
        Rtmp = project(hSpace, Rtmp);
        
        %save result
        Rguard{1} = Rtmp;


    %     interval(Rquad_zono)
    %     interval(delta_Rhit)
    %     interval(Rhit_err)

        %error index
        errorInd = max(2*rad(interval(Rhit_err) + interval(Rquad_zono)));
        %errorInd = max(2*rad(interval(Rhit_err)));
        
        if errorInd > options.maxProjectionError
        %if errorInd > 0.1
        %if errorInd > 0.05
            split = 1
        end

    else
        Rguard_noInt = [];
        Rguard{1} = [];
        split = [];
        flow = [];
    end
else
    [Rguard_noInt, Rguard] = partialGuardProjection(obj,hSpace,A,B,loc,options);
end


function [c_tilde, d_tilde, E_err] = trajectoryError(A,f,u,R0,tmin,tmax,tc,order,n)

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
u_tilde = u_tilde + E*tmax*uzono + tZono{1}*zonotope(f{1}-f{2}); %error due to flow correction

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


function [tmin,tmax,t_hit,t_total,RbeforeInt,x_hit] = switchingTimeInterval(obj,loc_old,h,options)

%get intersection times
%[tmin,tmax,t_total,RbeforeInt] = intersectionTime(obj,loc_old,h,options);
contDynamics = get(loc_old,'contDynamics');
[tmin,tmax,t_total,RbeforeInt] = intersectionTime_oneDirection(h,obj,contDynamics,options);

if ~isempty(RbeforeInt)

    %get contDynamics
    contDynamics = get(loc_old,'contDynamics');

    %compute t_hit
    %obtain data from location
    aux = 1e6*ones(length(options.x0),1);
    inv = interval(-aux,aux);
    tran_old = get(loc_old, 'transition');
    reset = get(tran_old{1}, 'reset');
    tran{1}=transition(h,reset,0,'a','b'); 
    %create dummy location
    loc = location('dummyLoc',0,inv,tran,contDynamics);  
    eventOptions = odeset('Events',eventFcn(loc));
    opt = odeset(eventOptions);
    %simulate continuous dynamics
    [contDynamics,t,x] = simulate(contDynamics,options,tmin,tmax,center(RbeforeInt),opt);

    t_hit = t(end) - tmin;
    x_hit = x(end,:).';
else
    t_hit = [];
    x_hit = [];
end


function [f, f_cell] = flowManipulation_inDir(A,R0,x0,f,hSpace)

%obtain flow in center
flow = A*x0 + f;
f_cell{1} = flow;

%check flow direction
limitRespected = 0;

while ~limitRespected
    flowDir_IH = interval(hSpace.c'*(A*R0 + f)*(1/(norm(hSpace.c)*norm(A*x0 + f))));
    
    if flowDir_IH(1,1)*flowDir_IH(1,2) < 0
        flowDir_IH_abs = 0;
    else
        [flowDir_IH_abs] = min([abs(flowDir_IH(1,1)), abs(flowDir_IH(1,2))]);
    end
    
    %flowDir_min = 0.4;
    %flowDir_min = 0.2;
    flowDir_min = 0.1;
    %flowDir_min = 0;
    if flowDir_IH_abs < 0.8*flowDir_min

        %strecthing factor for flow correction

        direction = sign(flowDir_IH(1,1) + flowDir_IH(1,2)); %is center positive or negative?
        alpha = direction*(flowDir_min - flowDir_IH_abs)*norm(hSpace.c)*norm(A*x0 + f)/(hSpace.c'*hSpace.c);

        %flow correction
        f = f + alpha*hSpace.c;
        flow = flow + alpha*hSpace.c;
    else
        limitRespected = 1;
    end
end
f_cell{2} = flow;


%------------- END OF CODE --------------