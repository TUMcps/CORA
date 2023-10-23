function res = test_intervalMatrix_expm
% test_intervalMatrix_expm - unit test function of computing the
%    exponential matrix
% 
% Syntax:
%    res = test_intervalMatrix_expm
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       12-November-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create interval matrix
% center
matrixCenter = [ ...
0.1030783856804509, -0.9101278715321481, 0.6053092778292544, -0.7608595779453735, -0.2653784107456929, 0.8833290280791213, -0.5867093811050932, -0.3920964615009785, -0.9330031960700655, 0.2224896347355838; ...
-0.3358623384378412, -0.6606199258835557, 0.3453023581373951, 0.3373340161562055, 0.1728220584061801, 0.0571996064779656, -0.3622320218412414, -0.6610674364129314, 0.8292507232278290, -0.6361200702133358; ...
0.7258992470204615, 0.2226582101541923, 0.4569764282907689, -0.8375957540168222, 0.9482367162906631, 0.5867543962581900, 0.8899927549727960, -0.0438826364149949, -0.0382136940774347, 0.4472851749922520; ...
-0.3046165455476422, 0.9005982283328562, -0.6360513364531990, 0.4124682114261358, -0.8080476776070364, -0.5728913883303508, -0.5916024282199834, -0.5937651381699391, 0.6597700848648822, -0.0758241326363369; ...
0.7096807070659510, 0.4998074647484749, 0.9713758502353100, 0.9376564546258963, 0.2457649012001917, 0.8577861212555120, 0.5407432738095963, 0.4696139508867379, 0.2178340806678498, 0.6784418775414429; ...
-0.7440856538867173, -0.5786954854131867, -0.7372557698893745, 0.6013975816933068, -0.7329906793875653, -0.6534173712962281, 0.8794343222311480, -0.8762545300911007, -0.9999624970157370, -0.8504824917103764; ...
0.9084850705727758, -0.8989100454650638, 0.9920886519204546, 0.9198930027855210, 0.1288885207530313, 0.3430266672898845, 0.0650526293251279, 0.9152025750401234, -0.2713856840944358, -0.2328442840941238; ...
-0.5748584958454852, 0.7413066606288192, -0.3728964105823804, -0.4680697157913984, 0.1567639105156624, 0.0153859735195778, 0.3410331006983629, 0.1066855250868497, 0.3867543050330005, 0.1620966903388332; ...
0.8838491015893175, 0.9233098306708436, -0.9715289785483805, -0.6531267294643079, 0.5304519869096429, -0.5966071048661272, 0.8169855310441498, -0.5101349852246480, 0.5764394251020426, 0.5059990141202664; ...
-0.2456725263694224, -0.9851975999597997, -0.1965669055496231, 0.4875990931423375, 0.7099594112562833, -0.2178327228088248, -0.2044138432022624, 0.0436619915851515, -0.5099717870717175, 0.4996095777778413];

% delta
matrixDelta = [ ...
0.3751207667608223, -0.2349336009150442, -0.7114373156136358, -0.0580297417046427, -0.6875042580158681, -0.6888486679279389, 0.2898714179410671, -0.9432795972096502, -0.0744770228233731, -0.3963969063670916; ...
0.9331880749809554, 0.1355905388370493, 0.4838847359375154, -0.3342096240982033, 0.2741229250281747, 0.2297590596236823, 0.8191398259446245, 0.8837717870885617, -0.5039016392307831, 0.8107585470622514; ...
-0.7455565709123813, -0.1245963657239815, 0.8256547839760635, 0.7972918918651206, -0.4663638009902955, 0.1544071007075494, 0.4389437518529968, -0.2634678486265958, 0.4327498507690415, -0.9828520365303051; ...
-0.2769199850635151, -0.8996598433302581, 0.7604815768228936, 0.7645414388534642, -0.0569280774274030, 0.5458701694220485, -0.0310482287289708, -0.6618155490120483, -0.7304297863386380, -0.5005923658155174; ...
0.9920766192928832, 0.5215662166008941, -0.2929182386406333, 0.7079880837383690, -0.7735794896756001, -0.3082773955514337, 0.0621397979540532, -0.5505736368355505, -0.8410189658240854, -0.8198939683656246; ...
0.3860881845498279, 0.2542717222750639, 0.3572289387234904, 0.8730971642753758, 0.5564179007277741, 0.4076247420165973, 0.2434429301452767, 0.9344061714082892, 0.8477000970613400, 0.1171450667758223; ...
0.4057697372477231, 0.2462132765349827, 0.9914145978660249, -0.9380325051045260, 0.4244674376721900, -0.5639726778588379, -0.6144158462040752, -0.1228217217497798, -0.9105674333380869, 0.8831378867707247; ...
-0.6990635930836684, 0.2036635255878141, -0.4756884053950827, -0.4454674039544488, 0.5530889077827099, 0.2735185993178468, 0.5563054687550690, -0.2159551627905312, -0.4648200795383441, -0.4668460112190254; ...
-0.8296006612129700, -0.0170496509199638, -0.7850525937103894, 0.6328021784933340, 0.9314641203938441, -0.7505545889636518, -0.4079744248915713, 0.2879091025995013, -0.6669359444094343, 0.5006304092185188; ...
0.9214237728662329, -0.6681325456884544, -0.8830816306080265, -0.7804084142981342, 0.5600326808340590, -0.9978908648229325, 0.3673553292739666, 0.7261916885374695, 0.1550107415038009, -0.4856084402282732];

% instantiate interval matrix
timeStep = 0.001;
M_int = timeStep*intervalMatrix(matrixCenter, matrixDelta);

% compute exponential matrix
eM_int = expm(M_int,1,6);

% compute exponential matrices of random vertices
eM = cell(10,1);
for i=1:10
    randMatrix = randPoint(M_int);
    eM{i} = expm(randMatrix{1});
end

% check enclosure
res = false(10,1);
for i=1:length(eM)
    res(i) = all(interval(eM{i}) <= eM_int.int);
end

%result of all enclosure checks
res = all(res);

% save results if test failed
if ~res
     file_name = strcat('test_intervalMatrix_expm_', ...
                             datestr(now,'mm-dd-yyyy_HH-MM'));
                  
     file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', file_name);
     save(file_path, 'eM');
end

% ------------------------------ END OF CODE ------------------------------
