function completed = example_reach_2D()
% updated: 13-September-2017, AElG
% As realized in the ACC paper
% Estimating the region of attraction via forward reacahble sets

% Number of state variables of the system
dim = 2;


% Set options -------------------------------------------------------------
options.tStart   = 0;   % Start time
options.timeStep = 0.05;   % Time step size for reachable set computation
options.tFinal   = 10000; % Final time


% Initial set for the input variables
options.uTrans = 0;

options.U    = zonotope([0,0]);
options.originContained = 0;

options.zonotopeOrder = 2000; % Zonotope order
options.polytopeOrder = 10;  % Polytope order
options.taylorTerms   = 4;


% Linearization options
options.tensorOrder          = 2;
options.advancedLinErrorComp = 0;
options.intermediateOrder    = 3;
options.reductionTechnique   = 'girard';

% Avoids splitting
options.recur = 0;


options.errorOrder           = 0;
options.oldError             = zeros(dim,1);
options.maxError             = 1.0e+02 *[...
                                         0.001;...
                                         0.001];
options.reductionInterval    = inf;

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
Mod_2D = nonlinearSys(2,1,@Model_2D,options);
figure
%--------------------------------------------------------------------------
        tic
        % Enlarge the set of the equliburim point rapidly
        options.R0 = zonotope(interval([-0.001;-0.001],[0.001;0.001]));
        R0 = options.R0;
        iDomain    = 1;
        iCounter   = 0;
        R_Domain{iDomain} = interval(options.R0);
        
        while iCounter ~= iDomain
            iCounter = iDomain;
            Enlarged_R0 = enlarge(R0,1.5);            
            options.oldError = zeros(2,1);
            options.R0       = Enlarged_R0;
            [~,Rlast,options] = contReach(Mod_2D,R_Domain,options);
            if ~isempty(Rlast{1}.set)
                R_Domain{iDomain} = interval(options.R0);
                hold on
                plot(Enlarged_R0,[1 2]);
                iDomain           = iDomain + 1;
            end
            R0 = Enlarged_R0;
        end
        % Define the search domain
        c0          = [0;0];
        G0(1,1)     = 1.6;
        G0(2,2)     = 1.5;
        Z0          = zonotope([c0,G0]);
        Interval_R0 = interval(Z0);
        options.R0  = Z0;
        
         % Define the initial domain (equ. point)
        R_Domain{1} = R_Domain{iDomain-1};
        I_Set = R_Domain{1};
        iDomain     = 1;
        
        % Load system and plot the search domain
        plot(I_Set)

    
        % Start computation
        for k = 1:5
            IH_Cells = segmentInterval(partition(Interval_R0(:,:),[2^k;2^k]));
           if k == 1
               color = [.2 .2 .2];
           elseif k == 2
               color = [.4 .4 .4];
           elseif k == 3
               color = [.6 .6 .6];
           elseif k == 4
               color = [.7 .7 .7];
           else
               color = [.8 .8 .8];
           end
            for i = 1:length(IH_Cells)
                options.oldError = zeros(2,1);
                options.R0 = zonotope(IH_Cells{i});
                
                for l=1:length(R_Domain)
                    bool = IH_Cells{i} <= R_Domain{l};
                    if bool
                        disp('Been There before !');
                        break
                    end
                end
%
                if ~bool
                    [~,Rlast,~] = contReach(Mod_2D,R_Domain,options);
                    if ~isempty(Rlast{1,1}.set)
                        % Computed reachable set
                        R_Domain{iDomain} = IH_Cells{i};
                        i
                        plotFilled(R_Domain{iDomain},[1 2],color)
                        iDomain           = iDomain + 1;
                    end
                end
            end
        end
        toc
        
%example completed
completed = 1;
        

function [Rcont_x,Rlast, options] = contReach(sys,R_Domain,options)
I_Set = interval([-1e-4 1e-4;-1e-4 1e-4]);
options.oldError = zeros(2,1);
%initialize reachable set computations
[Rnext_x, options] = initReach(sys, options.R0, options);

%while final time is not reached
t=options.tStart + options.timeStep;
iSet=1;

%save reachable set in cell structure
Rcont_x{iSet}   = Rnext_x.ti{1};

while t<options.tFinal
    %increment time and set counter
    if isempty(Rnext_x.ti{1})
        %t
        disp('Initial Reachable set explosion !');
        Rlast = Rnext_x.tp;
        break
    else
        t = t+options.timeStep;
        iSet = iSet+1; 
        options.t = t;
        %compute next reachable set
        %[Rnext_x, options] = reach(sys, Rnext_x, options);
        [Rnext_x,options]=post(sys,Rnext_x,options);
    
        if isempty(Rnext_x.tp{1}.set)
            t
        	disp('Reachable set explosion !');
            Rlast = Rnext_x.tp;
            break
        else
            Rcont_x{iSet} = reduce(Rnext_x.ti{1}, 'girard', 20);
        end
    end
    
    for i=1:length(R_Domain)
        bool = interval(Rcont_x{iSet}) <= R_Domain{i};
        if bool == 1
            t
            disp('Reachable set inside domain !');
            break
        end
    end
    
    if bool == 1
        Rlast = Rnext_x.tp;
        break
    end
    
    bool = interval(Rcont_x{iSet-1}) <= interval(Rcont_x{iSet});
    if bool == 1
        Rlast = Rnext_x.tp;
        t
        disp('Invariant Set !');
        break
    end
        
    Rlast = Rnext_x.ti;
    bool = interval(Rcont_x{iSet})<= I_Set;
    if bool == 1
        t
        disp('Reachable set inside Target !');
        break
    end
end
