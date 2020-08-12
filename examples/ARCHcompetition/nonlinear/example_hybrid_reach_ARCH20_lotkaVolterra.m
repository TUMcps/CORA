function res = example_hybrid_reach_ARCH20_lotkaVolterra()
% example_hybrid_reach_ARCH20_lotkaVolterra - hybrid Lotka-Volterra
%                                             benchmark from the 2020 ARCH 
%                                             competition
%
% Syntax:  
%    res = example_hybrid_reach_ARCH20_lotkaVolterra()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
% 
% Author:       Niklas Kochdumper
% Written:      12-May-2020
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


% Parameter  --------------------------------------------------------------

params.tFinal = 3.64;
params.startLoc = 1;

e = 0.008;
int = interval([1.3-e;1;0;0],[1.3+e;1;0;0]);
params.R0 = zonotope(int);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.005; 
options.taylorTerms = 3;
options.zonotopeOrder = 20;

options.alg = 'lin';
options.tensorOrder = 3;

options.errorOrder = 5;
options.intermediateOrder = 20;

options.guardIntersect = 'levelSet';


% Sysetem Dynamics --------------------------------------------------------

HA = lotkaVolterra();


% Reachability Analysis ---------------------------------------------------

t1 = tic;
R1 = reach(HA,params,options);
R2 = reachRem(HA,R1,params,options);
tComp = toc(t1);
disp(['Computation time: ',num2str(tComp)]);


% Verification ------------------------------------------------------------

% get interval enclosures of final reachable sets
int1 = interval(R1(1).timePoint.set{end});

Rfin = R2.timePoint.set{end};
hp = conHyperplane([0 0 0 1],params.tFinal);
Rfin = hp & conZonotope(zonotope(Rfin));
int2 = interval(Rfin);

int = int1 | int2;

% volume of final box
vol = volume(project(int,1:2));
disp(['Volume final box: ',num2str(vol)]);

% final value of the cut-value
cut = supremum(int(3));
disp(['Cnt value final box: ',num2str(cut)]);

% check if verified
res = 0;

if cut <= 0.2+eps
   res = 1; 
end

disp(['Speficiations verified: ',num2str(res)]);



% Visualization -----------------------------------------------------------

figure; hold on; box on

% plot reachable set
plot(R1(1),[1,2],'b','Filled',true,'EdgeColor','none');
plot(R1(2),[1,2],'g','Filled',true,'EdgeColor','none');
plot(R2,[1,2],'g','Filled',true,'EdgeColor','none');

% plot guard set
syms x y
vars = [x;y];
eq = (x-1)^2 + (y-1)^2 - 0.15^2;

ls = levelSet(eq,vars,'==');

plot(ls,[1,2],'r');
xlim([0.6,1.4]);
ylim([0.6,1.4]);
xlabel('x');
ylabel('y');

res = 1;

end


% Auxiliary Functions -----------------------------------------------------

function R = reachRem(HA,R,params,options)

    % get time interval reachable set
    R = R(2).timeInterval.set;
    
    % get guard set
    locs = HA.location;
    trans = locs{1}.transition;
    guard = trans{1}.guard;
    
    fun = @(x) (x(1)-1)^2 + (x(2)-1)^2 - 0.15^2;
    
    % compute interval enclosure
    int = [];
    w = [1;1;1;1/1000];
    
    for i = 1:length(R)
        if isIntersectingPCA(fun,R{i})
            list = splitLongestGen(diag(w)*R{i});
            for j = 1:length(list)
                if isIntersectingPCA(fun,list{j})
                    int = int | interval(diag(1./w)*list{j});
                end
            end
        end
    end
    
    % intersect with region x > 1
    infi = infimum(int);
    sup = supremum(int);
    infi(1) = 1;
    
    int = interval(infi,sup);
    
    % contract interval to tightly enclose guard set
    int = tightenDomain(guard,int,'forwardBackward');

    % compute intersection
    Rjump = guard & polyZonotope(int);
    
    % compute reachable set
    sys = locs{1}.contDynamics;
    
    params.tFinal = params.tFinal - infimum(int(4));
    params.R0 = Rjump;
    params = rmfield(params,'startLoc');
    options = rmfield(options,'guardIntersect');
    
    R = reach(sys,params,options);
end

function res =  isIntersectingPCA(fun,pZ)
% check if a polynomial zonotope intersects a level set using Prinzipal
% Component Analysis

    % convert to taylor model
    zono = zonotope(pZ);
    zono = reduce(zono,'pca',1);
    tay = taylm(zono);
    
    % range bounding
    int = interval(fun(tay));
    
    if in(int,0)
        res = 1;
    else 
        res = 0;
    end
end
