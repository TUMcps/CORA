function [Rguard,Rguard_noInt,split,RbeforeInt] = guardProjection(obj,hSpace,A,B,loc,partialIntersectionFlag,options)
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
%               10-September-2013
%               11-September-2013
%               16-June-2016
%               25-July-2016 (intervalhull replaced by interval)
%               22-August-2016
% Last revision:---

%------------- BEGIN CODE --------------


% %separating halspace
% c_new = - (hSpace.c.'*A).';
% d_new = hSpace.c.'*B*options.u;
% h_new = halfspace(c_new, d_new);


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
        for i=1:options.taylorTerms
            Gamma = Gamma + A^i*t_hit^i/factorial(i+1);
        end
        f = Theta*x0 + Gamma*u;

        %use reduced zonotope 
        Rred = reduce(RbeforeInt,'girard',options.errorOrder);
        
        [f, f_cell] = flowManipulation_inDir(A,Rred,x0,f,hSpace);
        
        %compute trajectory error
        linSys = linearSys('abstractionErrorSys',A,1);
        [E_err, c_tilde, d_tilde] = trajectoryError(linSys,f,u,RbeforeInt,tmin,tmax,t_hit,options.taylorTerms,hSpace.c);
        %[c_tilde, d_tilde, E_err] = trajectoryError(A,f_cell,u,RbeforeInt,tmin,tmax,t_hit,options.taylorTerms,hSpace.c);

        %compute terms of Taylor series for non-error part
        [k2,L2,Q2,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(hSpace, A, f, x0, Rred);
        %[k2,L2,Q2,lambdaQuad_inv_zono,Lambda_int] = mapsQuad_unc(hSpace, A, f, x0, Rred,E_err);



        %compute terms of Taylor series for error part
        c_tilde_mid = mid(c_tilde);
        d_tilde_mid = mid(d_tilde);
        c_tilde_rad = rad(c_tilde);
        d_tilde_rad = rad(d_tilde);
        [k2_mid,L2_mid,Q2_mid,lambdaQuad_inv_zono_mid,Lambda_int_mid] = mapsQuad(hSpace, A, f, x0, Rred, c_tilde_mid, d_tilde_mid); 
        [k2_rad,L2_rad,Q2_rad,lambdaQuad_inv_zono_rad,Lambda_int_rad] = mapsQuad(hSpace, A, f, x0, Rred, c_tilde_rad, d_tilde_rad); % <-- needs to be corrected; Rred has to be bounded by a box
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


        %combine results (option 3)
        Rhit_const = x0 + k2 + exactPlus(Rlin, Rquad_zono);
        Rtmp = Rhit_const + E_err;
        
%         %option 1
%         n = hSpace.c;
%         I_den = interval(n'*(A*Rred+f));
%         Rhit_err = (A*Rred+f)*(-1)*(n'*E_err)*(1/I_den) + E_err;
        
%         %option 2 (requires above computatins that are commented)
%         Rhit_err = k2_err + Rlin_err + Rquad_zono_err + E_err; 


    %     %for test
    %     Zk = k2 + zonotope([zeros(dim,1), 0.01*eye(dim)]);

        %compute mapping error
        [delta_Rhit,E] = mappingError(A,f_cell,u,RbeforeInt,tmin,tmax,t_hit,options.taylorTerms,Lambda_int,hSpace.c);
        %delta_Rhit2 = mappingError_split(A,f,u,RbeforeInt,tmin,tmax,t_hit,options.taylorTerms,Lambda_int,hSpace.c,0.1*options.timeStep);


        %compute set on hyperplane
        Rtmp = Rhit_const + delta_Rhit;
        %Rtmp = Rhit_const + Rhit_err;

        %reduce result
        Rtmp = reduce(Rtmp,'girard',options.zonotopeOrder);

        %project result
        Rtmp = project(hSpace, Rtmp);
        
        %save result
        Rguard{1} = Rtmp;

        
%         %projected errors
%         Rhit_err_proj = project(hSpace, Rhit_err);
%         Rquad_zono_proj = project(hSpace, Rquad_zono);

%         %error index
%         eLen = 2*rad(interval(Rhit_err_proj) + interval(Rquad_zono_proj));
%         errorInd = max(eLen(1:2));
%         %errorInd = max(2*rad(interval(Rhit_err)));
%         
%         if errorInd > options.maxProjectionError
%         %if errorInd > 0.1
%         %if errorInd > 0.05
%             split = 1
%         end

    %     %test of constant flow
    %     figure
    %     Rred = reduce(Rhit_const,'girard',10);
    %     plot(Rred,[1 2],'g');
    %     plot(RbeforeInt,[1 2],'g');
    %     for i = 1:100
    %         if i < 50
    %             x0 = randPointExtreme(RbeforeInt);
    %         else
    %             x0 = randPoint(RbeforeInt);
    %         end
    %         flowVec = A*x0 + f;
    %         xNew = x0 + tmax*flowVec;
    %         plot([x0(1),xNew(1)],[x0(2),xNew(2)],'r');
    %     end


    %     %test of exact flow
    %     figure
    %     Rred = reduce(Rguard,'girard',10);
    %     plot(Rred,[1 2],'g');
    %     plot(RbeforeInt,[1 2],'g');
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
    %             options.x0 = randPointExtreme(RbeforeInt);
    %         else
    %             options.x0 = randPoint(RbeforeInt);
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
    %         x0 = randPointExtreme(RbeforeInt);
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
        Rguard_noInt = [];
        Rguard{1} = [];
        split = [];
        flow = [];
    end
else
    [Rguard_noInt, Rguard] = partialGuardProjection(obj,hSpace,A,B,loc,options);
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
    inv = interval(-aux, aux);
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
    
    if infimum(flowDir_IH)*supremum(flowDir_IH) < 0
        flowDir_IH_abs = 0;
    else
        [flowDir_IH_abs] = min([abs(infimum(flowDir_IH)), abs(supremum(flowDir_IH))]);
    end
    
    %flowDir_min = 0.4;
    flowDir_min = 0.2;
    %flowDir_min = 0.1;
    %flowDir_min = 0;
    if flowDir_IH_abs < 0.8*flowDir_min

        %strecthing factor for flow correction

        direction = sign(infimum(flowDir_IH) + supremum(flowDir_IH)); %is center positive or negative?
        alpha = direction*(flowDir_min - flowDir_IH_abs)*norm(hSpace.c)*norm(A*x0 + f)/(hSpace.c'*hSpace.c);

        %flow correction
        f = f + alpha*hSpace.c;
        flow = flow + alpha*hSpace.c;
    else
        limitRespected = 1;
    end
end
f_cell{2} = flow;

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