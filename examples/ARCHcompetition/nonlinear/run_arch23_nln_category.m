function run_arch23_nln_category
% run_arch23_nln_category - runs all benchmarks from the nonlinear category
%     from the 2023 ARCH competition
%
% Syntax:
%    text = run_arch23_nln_category()
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

% only for CodeOcean capsule
% addpath(genpath('../code'));

% outputdir = '../results/'; % codeocean
outputdir = [CORAROOT filesep 'examples' filesep 'ARCHcompetition' filesep 'nonlinear' filesep]; % git
fid = fopen([outputdir 'results.csv'],'w');

% header
text = 'benchmark,instance,result,time,accuracy,timesteps';
fprintf(fid,'%s\n',text);

longdashes = '-------------------------------------------------------------------';


disp(['Robertson ' longdashes]);

text = example_nonlinear_reach_ARCH23_Robertson;
fprintf(fid,'%s\n',text{1});
fprintf(fid,'%s\n',text{2});
fprintf(fid,'%s\n',text{3});
saveas(gcf, [outputdir filesep 'Robertson.png']);
close(gcf);

disp(' ')
disp(' ')


disp(['Van der Pol ' longdashes]);
disp(' ');

text = example_nonlinear_reach_ARCH23_vanDerPol;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep 'vanderPol.png']);
close(gcf);

disp(' ')
disp(' ')


disp(['Laub Loomis ' longdashes]);
disp(' ');

text = example_nonlinear_reach_ARCH23_laubLoomis;
fprintf(fid,'%s\n',text{1});
fprintf(fid,'%s\n',text{2});
fprintf(fid,'%s\n',text{3});
saveas(gcf, [outputdir filesep 'laubloomis.png']);
close(gcf);

disp(' ')
disp(' ')


disp(['Lotka Volterra ' longdashes]);
disp(' ');

text = example_hybrid_reach_ARCH23_lotkaVolterra;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep 'lotkaVolterra.png']);
close(gcf);

disp(' ')
disp(' ')


disp(['Spacecraft ' longdashes]);
disp(' ');

text = example_hybrid_reach_ARCH23_spacecraft;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep 'spacecraft.png']);
close(gcf);

disp(' ')
disp(' ')


disp(['Traffic Scenario ' longdashes]);
disp(' ');

text = example_nonlinear_reach_ARCH23_trafficscenario;
fprintf(fid,'%s\n',text);
saveas(gcf, [outputdir filesep 'trafficscenario.png']);
close(gcf);

disp(' ')
disp(' ')

% Close .csv file
fclose(fid);

end

% ------------------------------ END OF CODE ------------------------------
