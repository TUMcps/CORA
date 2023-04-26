% load common parameters
load IEEE30

% set path
path = [CORAROOT filesep 'models' filesep 'powerSystems'];

% full system
name = 'IEEE30';
bus.output = [];
bus.input = [];
bus.generator = [1, 2, 5, 8, 11, 13];
bus.load = [3, 4, 6, 7, 9, 10, 12, 14, 15, 16:30];
bus.fault = [];
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% subsystem 1
name = 'IEEE30_sub1';
bus.output = [2, 4, 17, 20, 24, 25];
bus.input = [5, 6, 10, 22, 27];
bus.generator = [1, 2, 13];
bus.load = [3, 4, 12, 14, 15, 16:20, 23:26];
bus.fault = [];
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% subsystem 2
name = 'IEEE30_sub2';
bus.output = [5, 6, 10, 22, 27];
bus.input = [2, 4, 17, 20, 24, 25];
bus.generator = [5, 8, 11];
bus.load = [6, 7, 9, 10, 21, 22, 27:30];
bus.fault = [];
bus.slack = [];

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% four parts - subsystem 1
name = 'IEEE30_fourParts_sub1';
bus.output = [2, 4, 12];
bus.input = [5, 6, 14, 15, 16];
bus.generator = [1, 2, 13];
bus.load = [3, 4, 12];
bus.fault = [];
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% four parts - subsystem 2
name = 'IEEE30_fourParts_sub2';
bus.output = [14, 15, 16, 17, 20];
bus.input = [10, 12, 23];
bus.generator = [];
bus.load = 14:20;
bus.fault = [];
bus.slack = [];

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% four parts - subsystem 3
name = 'IEEE30_fourParts_sub3';
bus.output = [5, 6, 10, 22];
bus.input = [2, 4, 17, 20, 24];
bus.generator = [5, 8, 11];
bus.load = [6, 7, 9, 10, 21, 22];
bus.fault = [];
bus.slack = [];

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% four parts - subsystem 4
name = 'IEEE30_fourParts_sub4';
bus.output = [23, 24, 28];
bus.input = [8, 9, 15, 22];
bus.generator = [];
bus.load = 23:30;
bus.fault = [];
bus.slack = [];

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% faulty system
name = 'IEEE30_fault';
bus.output = [];
bus.input = [];
bus.generator = [1, 2, 5, 8, 11, 13];
bus.load = [3, 4, 6, 7, 9, 10, 12, 14, 15, 16:30];
bus.fault = 1;
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% four parts - faulty subsystem 1
name = 'IEEE30_fault_fourParts_sub1';
bus.output = [2, 4, 12];
bus.input = [5, 6, 14, 15, 16];
bus.generator = [1, 2, 13];
bus.load = [3, 4, 12];
bus.fault = 1;
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');