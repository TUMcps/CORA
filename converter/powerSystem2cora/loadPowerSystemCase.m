function loadPowerSystemCase(casefilename, filename, varargin)
% loadPowerSystemCase - loads a MATPOWER model and saves it as a power system 
% description (psd)
%
%    structure of the power system description:
%    .name - name of power system
%    .Y - admittance matrix
%    .Pd - active power demand
%    .Qd - reactive power demand
%    .VM - voltage magnitude
%    .genParam - generator parameters
%    .bus.output - vector of output buses
%    .bus.input - vector of input buses
%    .bus.generator - vector of generator buses
%    .bus.load - vector of load buses
%    .bus.fault - index of faulty bus
%    .bus.slack - index of slack bus
%
% Syntax:
%    loadPowerSystemCase(casefilename, filename)
%    loadPowerSystemCase(casefilename, filename, type)
%
% Inputs:
%    casefilename - name of the MATPOWER case (specified as string)
%    filename - name of the file under which the power system description
%    should be saved
%    type - model type
%
% Outputs:
%    -
%
% Reference:
%    [1] R. D. Zimmerman, C. E. Murillo-S ́anchez. Matpower User’s Manual, 
%        Version 7.1. 2020. [Online]. Available: 
%        https://matpower.org/docs/MATPOWER-manual-7.1.pdf
%        doi: 10.5281/zenodo.4074122
%    [2] Y. C. Chen and A. D. Dom ́ınguez-Garc ́ıa. Assessing the impact of 
%        wind variability on power system small-signal reachability. In 
%        Proc. of the International Conference on System Sciences, 
%        pages 1–8, 2011.
%
% Example:
%    [x,y] = loadPowerSystemCase('case14')

%
% Other m-files required: requires installation of MATPOWER
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       05-May-2022
% Last update:   01-June-2022 (MA, PSSE models can be loaded)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% load case
if nargin == 3
    type = varargin{1};
else
    type = 'MATPOWER';
end

% MATPOWER case loaded
if strcmp(type,'MATPOWER')
    mpc = loadcase(casefilename);
% PSSE case loaded
elseif strcmp(type,'PSSE')
    mpc = psse2mpc(casefilename);
end

% set name
name = filename;

% create bus admittance matrix
Y = makeYbus(mpc); 

%% fetch bus data (see of Table B-1 on p. 139 in [1])
% active power demand
Pd = mpc.bus(:,3)/mpc.baseMVA; 
% reactive power demand
Qd = mpc.bus(:,4)/mpc.baseMVA; 
% voltage magnitude
VM = mpc.bus(:,8); 
% find load buses (entry "1" in second column of mpc.bus)
bus.load = find(mpc.bus(:,2)==1)';
% find generator buses (entry "2" in second column of mpc.bus);
% currently, we also consider the alack bus as a generator bus (entry "3" 
% in second column of mpc.bus)
bus.generator = [find(mpc.bus(:,2)==3); find(mpc.bus(:,2)==2)]';
% find slack bus (entry "3" in second column of mpc.bus)
bus.slack = find(mpc.bus(:,2)==3)';
% input buses (MATPOWER does not specify input buses)
bus.input = [];
% output buses (MATPOWER does not specify output buses)
bus.output = [];
% faulty buses (MATPOWER does not specify faulty buses)
bus.fault = [];

%% generator parameters
% MATPOWER does not specify the dynamics of generator; for this reason we
% set default parameters for one generator from [2] and copy these 
% parameters to all generators.

% damping coefficient [s/rad]
D = 0.04;
% rotational inertia [MJ/Hz2]
M = 1/(15*pi);
% proportional gain of the governor [-]
R_d = 6*pi;
% time constant of the governor [s]
T_sv = 1;
% inverse of the absolute value of the generator admittance 1/|Y_{g,i}|
Y_m = 5;
% phase angles of the generator admittance [rad]
Psi_g = -pi/2;
    
% frequency
genParam.omega_s = 120*pi;

% copy parameters for each generator
nrOfGenerators = length(bus.generator);
genParam.D = D*ones(1,nrOfGenerators);
genParam.M = M*ones(1,nrOfGenerators);
genParam.R_d = R_d*ones(1,nrOfGenerators);
genParam.T_sv = T_sv*ones(1,nrOfGenerators);
genParam.Y_m = Y_m*ones(1,nrOfGenerators);
genParam.Psi_g = Psi_g*ones(1,nrOfGenerators);

%% save power system description
% set path
path = [CORAROOT filesep 'models' filesep 'powerSystems'];
if ~isfolder(path)
    mkdir([CORAROOT filesep 'models'],'powerSystems');
end
% save 
save([path filesep filename],'bus','genParam','name','Pd','Qd','VM','Y');
% remove and add path so that file can be found
warOrig = warning;
warning('off','all');
rmpath(path);
warning(warOrig);
addpath(path);

% ------------------------------ END OF CODE ------------------------------
