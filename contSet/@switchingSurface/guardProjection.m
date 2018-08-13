function Rguard = guardProjection(obj,A,B,t_hit,tmin,tmax,R0,options)
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
% Written:      08-August-2011
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%get data
x0 = center(R0);
u = B*options.u;

% tic
% t_hit_test = crossingTime(obj,x0,A,u,tc,tmin,tmax);
% toc

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
Rred = reduce(R0,'girard',options.errorOrder);

%compute terms of Taylor series
%[k,L,Q,Cint] = maps(obj, A, f, x0, R0);
%[k3,L3,Q3,C_c,C_g] = mapsCubic(obj, A, f, x0, Rred);
[k2,L2,Q2,lambdaQuad_zono,Lambda_int] = mapsQuad(obj, A, f, x0, Rred);
%[k3,L3,Q3] = mapsQuad_old(obj, A, f, x0, Rred);

% %obtain matrix zonotope
% for i=1:length(C_g(1,:))
%     C_gen{i} = C_g(:,i);
% end
% Czono = matZonotope(C_c, C_gen);

%individual reachable sets 
dR = R0 + (-x0);
dRred = Rred + (-x0);
Rlin = L2*dR;
%Rquad = 1/2*quadraticMultiplication(dRred,Q3);

Rquad_zono_aux = quadraticMultiplication_zono(dRred,Q2);
Rquad_zono = 1/2*lambdaQuad_zono*reduce(Rquad_zono_aux,'girard',options.zonotopeOrder);
%Rquad_zono_old = 1/2*quadraticMultiplication_zono(dRred,Q3);

% Rcubic_int = 1/6*cubicIntMultiplication(dR,Cint);
% Rcubic_old = 1/6*cubicIntMultiplication_zono(dRred,C_c,C_g);
% 
% Rcubic = 1/6*tensorMultiplication_zono(dRred,Czono,options);

% %test:
% dIH = zonotope(interval(dR));
% [c_int, deltaVec] = cubicMultiplication_interval(dIH,Cint);

% %combine Rlin and Rquad
% Zlin = get(Rlin,'Z');
% Zquad = get(Rquad,'Z');
% 
% gens = length(Zlin(1,:));
% Zquad_cut = Zquad(:,1:gens);
% 
% Zcomb = Zlin + Zquad_cut;
% Zquad_new = Zquad(:,(gens+1):end);
% 
% Rcomb = zonotope(Zcomb);
% Rquad_new = zonotope([0*x0, Zquad_new]);
% 
% %combine results
% Rhit_const = k + Rcomb + Rquad_new + Rcubic;

%combine results
%Rhit_const = k + Rlin + Rquad + Rcubic;
Rhit_const = k2 + Rlin + Rquad_zono;

%for test
Zk = k2 + zonotope([zeros(dim,1), 0.01*eye(dim)]);

%compute mapping error
delta_Rhit = mappingError(A,f,u,R0,tmin,tmax,t_hit,options.taylorTerms,Lambda_int,obj.c);

%compute set on hyperplane
Rguard = Rhit_const + delta_Rhit;



%function which determines when the state of x'=Ax+u crosses the hypeplane
%c^T x = d with initial state x0
function t_hit = crossingTime(obj,x0,A,u,tc,tmin,tmax)

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
    if (abs(obj.c'*x_upper - obj.d) < abs(obj.c'*x_lower - obj.d))
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


function delta_Rhit = mappingError(A,f,u,R0,tmin,tmax,tc,order,Lambda_int,n)


%compute powers of tmin and tmax-------------------------------------------
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
%--------------------------------------------------------------------------

%obtain center and delta set
x0 = center(R0);
dR = R0 + (-x0);

% %auxiliary set: old
% %Raux = A*dR;
% Raux_old = zonotope(0*x0);
% for i=2:order
%     Radd_1 = (tZono{i-1} + (-tc^(i-1)))*x0;
%     Radd_2 = tZono{i-1}*dR;
%     Raux_old = Raux_old + A^i/factorial(i)*(Radd_1 + Radd_2);
% end

%auxiliary set: new
eAt_unc_1 = 0*A;
eAt_unc_2 = 0*A;
for i=2:order
    eAt_unc_1 = eAt_unc_1 + A^(i-1)/factorial(i)*(tZono{i-1} + (-tc^(i-1)));
    eAt_unc_2 = eAt_unc_2 + A^i/factorial(i)*tZono{i-1};
end
Raux = eAt_unc_1*(A*x0 + u) + eAt_unc_2*dR;

% %auxiliary set: test
% eAt_unc_1 = A;
% eAt_unc_2 = A;
% for i=2:order
%     eAt_unc_1 = eAt_unc_1 + A^i/factorial(i)*(tZono{i-1} + (-tc^(i-1)));
%     eAt_unc_2 = eAt_unc_2 + A^i/factorial(i)*tZono{i-1};
% end
% Raux_test = eAt_unc_1*x0 + eAt_unc_2*dR;

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:(order-1)
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
Mstar = M;
%one more...
i = i+1;
Apower_abs = Apower_abs*A_abs;
M = M + Apower_abs*tmax^i/factorial(i);
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
Wstar = expm(A_abs*tmax) - Mstar;
%compute absolute value of W for numerical stability
W = abs(W);
Wstar = abs(Wstar);
E = interval(-W,W);
Estar = interval(-Wstar,Wstar);

%final additional set
%delta_Rhit = tZono{1}*Raux + E*R0 + Estar*interval(tmin,tmax)*u;
G = tZono{1}*Raux + E*R0 + Estar*interval(tmin,tmax)*u;

%error due to uncertain hitting time
delta_t = interval((-n')*G)/Lambda_int;
delta_t_mid = mid(delta_t);
delta_t_rad{1} = rad(delta_t);
delta_t_zono = matZonotope(delta_t_mid,delta_t_rad);

%final uncertainty
delta_Rhit = G + delta_t_zono*(A*R0 + f);


%------------- END OF CODE --------------