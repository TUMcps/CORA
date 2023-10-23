function text = benchmark_hybrid_reach_ARCH23_lotkaVolterra()
% benchmark_hybrid_reach_ARCH23_lotkaVolterra - hybrid Lotka-Volterra
%    benchmark from the 2023 ARCH competition
%
% Syntax:
%    text = benchmark_hybrid_reach_ARCH23_lotkaVolterra()
%
% Inputs:
%    -
%
% Outputs:
%    text - true/false
% 

% Authors:       Niklas Kochdumper
% Written:       12-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

tFinal = 3.64;
params.tFinal = tFinal;
params.startLoc = 1;

e = 0.012;
I = interval([1.3-e;1;0],[1.3+e;1;0]);
params.R0 = zonotope(I);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.005; 
options.taylorTerms = 3;
options.zonotopeOrder = 20;

options.alg = 'lin';
options.tensorOrder = 3;

options.errorOrder = 5;
options.intermediateOrder = 20;

options.guardIntersect = 'levelSet';


% System Dynamics ---------------------------------------------------------

HA = lotkaVolterra();


% Reachability Analysis ---------------------------------------------------

% start timer
t1 = tic;

% figure; hold on; box on;
% R1:	tangential crossing in upper arc (-> R1(2): 'inside')
%       full crossing in lower arc (-> R1(3): 'inside')
%       both R1(2) and R1(3) stop when they come back to 'outside'
R1 = reach(HA,params,options);

% continuation of R1(2) by R2, R3, and R4
% R2:   succeeds R1(2) with continuous dynamics until fully in 'outside'
tshift = 2;
[R2,tStop] = aux_reachRem(HA,R1(2),params,options,tshift,'sup');

% R3:   switch back to hybrid using end of R2 as start set,
%       full crossing in lower arc (-> R3(2): 'inside')
%       stops when again fully in 'outside'
params.R0 = R2.timePoint.set{end};
params.tFinal = params.tFinal - tStop;
R3 = reach(HA,params,options);

% R4:   succeeds R3 with continuous dynamics until end of time horizon
% params.tFinal = params.tFinal - tFinal;
params.tFinal = tFinal;
[R4,~] = aux_reachRem(HA,R3(2),params,options,0,'inf');

% continuation of R1(3) by R5
% R5:   succeeds R3 with continuous dynamics until end of time horizon
[R5,~] = aux_reachRem(HA,R1(3),params,options,0,'inf');

% measure full computation time
tComp = toc(t1);
disp(['Computation time: ',num2str(tComp)]);


% Verification ------------------------------------------------------------

% hyperplane at tFinal
hp = conHyperplane([0 0 1],params.tFinal);

% get interval enclosures of final reachable sets
Rfinal1 = R4.timePoint.set{end};
I1 = interval( hp & conZonotope(zonotope(Rfinal1)) );
Rfinal2 = R5.timePoint.set{end};
I2 = interval( hp & conZonotope(zonotope(Rfinal2)) );

I = I1 | I2;

% area (x-y) of final box
area = volume(project(I,1:2));
disp(['Area of final box: ',num2str(area)]);

% Visualization -----------------------------------------------------------

% show by plot:
% 1. one final set: two guard crossings
% 2. one final set: four guard crossings
% 3. no odd number of crossings

figure; hold on; box on;

% plot reachable set
plot(R1(1),[1,2],'FaceColor',colorblind('b'),'Unify',true);
plot(R1(2),[1,2],'FaceColor',colorblind('r'),'Unify',true);
plot(R1(3),[1,2],'FaceColor',colorblind('gray'),'Unify',true);
plot(R2,[1,2],'FaceColor',colorblind('y'),'Unify',true);
plot(R3(1),[1,2],'FaceColor',colorblind('r'),'Unify',true);
plot(R3(2),[1,2],'FaceColor',colorblind('y'),'Unify',true);
plot(R4,[1,2],'FaceColor',colorblind('r'),'Unify',true);
plot(R5,[1,2],'FaceColor',colorblind('b'),'Unify',true);

% plot guard set
syms x y
vars = [x;y];
eq = (x-1)^2 + (y-1)^2 - 0.161^2;

ls = levelSet(eq,vars,'==');
spec = specification(ls, 'unsafeSet');

plot(spec,[1,2]);
xlim([0.6,1.4]);
ylim([0.6,1.4]);
xlabel('x');
ylabel('y');

res = true;
% text for .csv
text = ['LOVO21, ,',num2str(res),',',num2str(tComp),',',num2str(area),', '];

end


% Auxiliary functions -----------------------------------------------------

function [R,tFinal] = aux_reachRem(HA,R,params,options,tshift,flow)

    % get time-interval reachable set
    R = R.timeInterval.set;
    
    % get guard set
    locs = HA.location;
    trans = locs(1).transition;
    guard = trans(1).guard;
    
    fun = @(x) (x(1)-1)^2 + (x(2)-1)^2 - 0.161^2;
    
    % compute interval enclosure
    int = [];
    w = [1;1;1/1000];
    
    for i = 1:length(R)
        if aux_isIntersectingPCA(fun,R{i})
            list = splitLongestGen(diag(w)*R{i});
            for j = 1:length(list)
                if aux_isIntersectingPCA(fun,list{j})
                    int = int | interval(diag(1./w)*list{j});
                end
            end
        end
    end
    
    % intersect with region x < 1 or x > 1 (depending on direction of flow)
    lb = infimum(int);
    ub = supremum(int);
    if strcmp(flow,'inf')
        lb(1) = 1;
    elseif strcmp(flow,'sup')
        ub(1) = 1;
    end
    
    int = interval(lb,ub);
    
    % contract interval to tightly enclose guard set
    int = tightenDomain(guard,int,'forwardBackward');

    % compute intersection
    Rjump = guard & polyZonotope(int);
    
    % compute reachable set
    sys = locs(1).contDynamics;
    
    params.tFinal = params.tFinal - tshift - infimum(int(3));
    params.R0 = Rjump;
    params = rmfield(params,'startLoc');
    options = rmfield(options,'guardIntersect');
    
    R = reach(sys,params,options);
    tFinal = params.tFinal + infimum(int(3));
end

function res = aux_isIntersectingPCA(fun,pZ)
% check if a polynomial zonotope intersects a level set using Principal
% Component Analysis

    % convert to taylor model
    zono = zonotope(pZ);
    zono = reduce(zono,'pca',1);
    tay = taylm(zono);
    
    % range bounding
    I = interval(fun(tay));
    
    res = contains(I,0);
end

% ------------------------------ END OF CODE ------------------------------
