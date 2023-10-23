function res = test_nonlinearSysDT_observe_01_tank()
% test_nonlinearSysDT_observe_01_tank - unit_test_function for guaranteed
%    state estimation of nonlinear discrete-time tank system.
%
%    Checks the state estmation of the nonlinearSysDT class for a tank
%    example from [1]; It is checked whether the enclosing interval of the 
%    final observed set is close to an interval provided by a previous 
%    solution that has been saved
%
% Syntax:
%    res = test_nonlinearSysDT_observe_01_tank
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Reference:
%    [1] M. Althoff "Guaranteed State Estimation in CORA 2021", ARCH 2021

% Authors:       Matthias Althoff
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define tank model
dt = 0.5;
n = 6; m = 3; r = 3;
fun = @(x,u) tank6EqDT_inflow4(x,u,dt);
out_fun = @(x,u)[0,1,0,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0]*x(1:n);
tank = nonlinearSysDT('tank',fun,dt,n,m,out_fun,r);

% parameters
params.tFinal = 10;
params.R0 = zonotope(20*ones(6,1),4*eye(6));
params.U = zonotope(zeros(3,1));
params.V = zonotope(zeros(3,1),0.2*eye(3));
params.W = zonotope(zeros(6,1),0.001*eye(6));
params.u = [0.1000    0.1000    0.1000    0.1025    0.1050    0.1075    0.1100    0.1125    0.1150    0.1175    0.1200    0.1225    0.1250    0.1275    0.1300    0.1325    0.1350    0.1375    0.1400    0.1425;
    0.1000    0.1000    0.1000    0.1025    0.1050    0.1075    0.1100    0.1125    0.1150    0.1175    0.1200    0.1225    0.1250    0.1275    0.1300    0.1325    0.1350    0.1375    0.1400    0.1425;
    0.1000    0.1000    0.1000    0.1050    0.1100    0.1150    0.1200    0.1250    0.1300    0.1350    0.1400    0.1450    0.1500    0.1550    0.1600    0.1650    0.1700    0.1750    0.1800    0.1850];
params.y = [18.587107801413037  19.284508642234258  19.804172917266072  19.984844940204784  20.278953550132702  20.483743290472002  20.674406205375988  20.564287471887262  20.524616208625844  20.390858909530351  20.579505841851283  20.499347774002867  20.334643962630423 20.272233256206864  20.167743326619661  20.169739624168546  20.104461580569449  20.060900527281170  19.937554706454481  19.834875810069285;
  20.225990245627816  20.158180661305583  20.171464843998834  20.327437617339143  20.153618747112297  20.380068720276380  20.376688700844845  20.299247293957542  20.243419253755508  20.277306209637565  20.290071755637200  20.468834597540557  20.270783553258500 20.223394465701581  20.467148457713510  20.201977670029528  20.431607252968249  20.481135248015988  20.487405103697650  20.410091165821306;
  19.112251238495777  18.973762663726877  18.948170824104380  18.992192949983256  19.001881958518979  19.150117548502443  18.975701530881086  19.088283489795671  18.991871883700380  19.066458349897363  18.943805866972866  19.214278885933556  18.957978956507624 19.084269923675965  19.183361721875762  19.081041506305120  19.118244206333170  19.284683679863136  19.156503106336739  19.215720898675404];

% settings
options.timeStep = dt;
options.zonotopeOrder = 100;
options.reductionTechnique = 'pca';
options.tensorOrder = 2;
options.errorOrder = 1;
% options.points = 1;
% options.p_conf = 0.9990;

% Set of evaluated estimators
Estimator = {
    'FRad-A' 
    'FRad-B' 
    'FRad-C' 
    };

% init
resPartial = [];

% set accuracy
accuracy = 1e-6;

%% perform evaluation
% loop over estimators
for iEst = 1:length(Estimator)

    % set algorithm
    estName = Estimator{iEst};
    options.alg = estName;
    
    % if constrained zonotopes or ellipsoids are used
    paramsNew = params;
    if any(strcmp(estName,{'ESO-A','ESO-B','ESO-C','ESO-D','ESO-E'}))
        paramsNew.R0 = ellipsoid(params.R0);
        paramsNew.W = ellipsoid(params.W);
        paramsNew.V = ellipsoid(params.V);
    elseif any(strcmp(estName,{'CZN-A','CZN-B'}))
        paramsNew.R0 = conZonotope(params.R0);
    end
    
    % evaluate observer
    estSet = observe(tank,paramsNew,options);
    
    % enclose last estimated set by interval
    IH = interval(estSet.timePoint.set{end});
    
    % obtain enclosing intervals
    if strcmp(options.alg,'FRad-A')
        if params.tFinal == 10
            IH_saved = interval( ...
                [19.2976721922255834; 19.6139972122578357; 14.2446680848811198; 20.0892697254588484; 18.9581928439432339; 16.1726262416728765], ...
                [19.7876666310276192; 20.1035308207714820; 30.1688768049905462; 20.7224260892464685; 19.3788828254191650; 23.7131744879695638]);
        else
            IH_saved = interval( ...
                [17.6813648886554518; 18.2156906209731488; 21.2088894208727972; 22.3622857916818099; 22.4614208115655352; 19.2343954709766152], ...
                [18.3586147868548224; 18.9069663981664888; 23.8031380887066000; 23.2005173305999541; 23.1757768194797045; 22.9771198272602817]);
        end
    elseif strcmp(options.alg,'FRad-B')
        if params.tFinal == 10
            IH_saved = interval( ...
                [19.2712596426652674; 19.5913693494689234; 14.2056931305107454; 20.0980777032700715; 18.9587395711021074; 16.1529188853400001], ...
                [19.7900964406908493; 20.1029076589013762; 30.4543302458481833; 20.7268550247347747; 19.3789984819565717; 23.7268155612348224]);
        else
            IH_saved = interval( ...
                [17.6928498777203806; 18.2262865296234970; 21.2428263743206962; 22.3924573528288349; 22.4929791654838702; 19.2318792831541572], ...
                [18.3491195855406168; 18.8984370397519434; 23.7623779649680849; 23.1693183620486991; 23.1497672369610150; 22.9805998157471549]);
        end
    elseif strcmp(options.alg,'FRad-C')
        if params.tFinal == 10
            IH_saved = interval( ...
                [19.2882829662340818; 19.6049109089969669; 14.1375551339295349; 20.0683709752421144; 18.9552007160255158; 16.1715832824749270], ...
                [19.8034985374044119; 20.1189759828253436; 30.1814503670473009; 20.7408289024531278; 19.3772251039683709; 23.7132313182795755]);
        else
            IH_saved = interval( ...
                [17.6996609730186378; 18.2366127949558070; 21.3886795671867809; 22.3511769355149319; 22.5106975175051005; 19.2412346771689293], ...
                [18.3377371228686741; 18.8833803746827193; 23.6252998964590724; 23.2109899193938354; 23.1329874737624905; 22.9707010569532990]);
        end
    end

    %check if slightly bloated versions enclose each other
    resPartial(end+1) = isequal(IH,IH_saved,1+accuracy);
end

% final result
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
