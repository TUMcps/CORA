function run_arch23_aff_category()
% run_arch23_aff_category - runs all benchmarks from the linear category
%     from the 2023 ARCH competition
%
% Syntax:
%    text = run_arch23_aff_category()
%
% Inputs:
%    -
%
% Outputs:
%    -

% Authors:       Mark Wetzlinger
% Written:       23-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% below only for CodeOcean capsule...
% addpath(genpath('../code'));

% initialize .csv file
% outputdir = '../results/'; % codeocean
outputdir = [CORAROOT filesep 'examples' filesep 'ARCHcompetition' filesep 'linear']; % git
fid = fopen([outputdir 'results.csv'],'w');

% header
text = 'benchmark,instance,result,time';
fprintf(fid,'%s\n',text);

% helper dashes
longdashes = '-------------------------------------------------------------------';
shortdashes = '------------------';

% HEAT 3D Benchmark -------------------------------------------------------

disp(['Heat 3D ' longdashes]);

disp(' ');
id = 'HEAT01';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_heat3D_HEAT01();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close;

disp(' ');
id = 'HEAT02';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_heat3D_HEAT02();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close;

disp(' ');
id = 'HEAT03';
disp([id shortdashes]);
text = example_linear_reach_ARCH23_heat3D_HEAT03();
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close;


% Beam Benchmark ----------------------------------------------------------

disp(' ');
disp(['Beam ' longdashes]);

disp(' ');
id = 'CBC01';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_beam_CBC01();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'CBC02';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_beam_CBC02();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'CBC03';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_beam_CBC03();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'CBF01';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_beam_CBF01();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'CBF02';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_beam_CBF02();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'CBF03';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_beam_CBF03();
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);


% ISS Benchmark -----------------------------------------------------------

disp(' ');
disp(['ISS ' longdashes]);

disp(' ');
id = 'ISSC01_ISS02';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_iss_ISSC01_ISS02;
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'ISSC01_ISU02';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_iss_ISSC01_ISU02;
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'ISSF01_ISS01';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_iss_ISSF01_ISS01;
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'ISSF01_ISU01';
disp([id shortdashes]);
text = example_linear_verifyFast_ARCH23_iss_ISSF01_ISU01;
fprintf(fid,'%s\n',text);
% saveas(gcf, [outputdir filesep id '.png']);
close(gcf);


% Platoon -----------------------------------------------------------------

disp(' ');
disp(['Platoon ' longdashes]);

disp(' ');
id = 'PLAA01-BND42';
disp([id shortdashes]);
text = example_linearParam_reach_ARCH23_platoon_PLAA01_BND42;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'PLAA01-BND50';
disp([id shortdashes]);
text = example_linearParam_reach_ARCH23_platoon_PLAA01_BND50;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'PLAD01-BND30';
disp([id shortdashes]);
text = example_linearParam_reach_ARCH23_platoon_PLAD01_BND30;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'PLAD01-BND42';
disp([id shortdashes]);
text = example_linearParam_reach_ARCH23_platoon_PLAD01_BND42;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'PLAN01';
disp([id shortdashes]);
text = example_linearParam_reach_ARCH23_platoon_PLAN01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

% Spacecraft Benchmark ----------------------------------------------------

disp(' ');
disp(['Spacecraft ' longdashes]);

disp(' ');
id = 'SRA01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA02';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA02;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA03';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA03;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA04';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA04;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA05';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA05;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA06';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA06;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA07';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA07;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRA08';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRA08;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRNA01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRNA01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRU01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRU01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'SRU02';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_rendezvous_SRU02;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

% Powertrain Benchmark ----------------------------------------------------

disp(' ');
disp(['Powertrain ' longdashes]);

disp(' ');
id = 'DTN01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_powerTrain_DTN01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'DTN02';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_powerTrain_DTN02;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'DTN03';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_powerTrain_DTN03;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'DTN04';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_powerTrain_DTN04;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'DTN05';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_powerTrain_DTN05;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'DTN06';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_powerTrain_DTN06;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

% Gearbox Benchmark -------------------------------------------------------

disp(' ');
disp(['Gearbox ' longdashes]);

disp(' ');
id = 'GRBX01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_gearbox_GRBX01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'GRBX02';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_gearbox_GRBX02;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

% Brake Benchmark ---------------------------------------------------------

disp(' ');
disp(['Brake ' longdashes]);

disp(' ');
id = 'BRKDC01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_brake_BRKDC01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'BRKNC01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_brake_BRKNC01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);

disp(' ');
id = 'BRKNP01';
disp([id shortdashes]);
text = example_hybrid_reach_ARCH23_brake_BRKNP01;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep id '.png']);
close(gcf);


% End ---------------------------------------------------------------------

% Close .csv file
fclose(fid);

end

% ------------------------------ END OF CODE ------------------------------
