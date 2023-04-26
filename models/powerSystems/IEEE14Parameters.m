% load common parameters
load IEEE14

% set path
path = [CORAROOT filesep 'models' filesep 'powerSystems'];

% full system
name = 'IEEE14';
bus.output = [];
bus.input = [];
bus.generator = [1, 2, 3, 6, 8];
bus.load = [4, 5, 7, 9:14];
bus.fault = [];
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');


% subsystem 1
name = 'IEEE14_sub1';
bus.output = [2, 5, 6, 13];
bus.input = [3, 4, 11, 14];
bus.generator = [1, 2, 6];
bus.load = [5, 12, 13];
bus.fault = [];
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% subsystem 2
name = 'IEEE14_sub2';
bus.output = [3, 4, 11, 14];
bus.input = [2, 5, 6, 13];
bus.generator = [3, 8];
bus.load = [4, 7, 9, 10, 11, 14];
bus.fault = [];
bus.slack = [];

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% faulty system
name = 'IEEE14_fault';
bus.output = [];
bus.input = [];
bus.generator = [1, 2, 3, 6, 8];
bus.load = [4, 5, 7, 9:14];
bus.fault = 1;
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');

% faulty subsystem 1
name = 'IEEE14_fault_sub1';
bus.output = [2, 5, 6, 13];
bus.input = [3, 4, 11, 14];
bus.generator = [1, 2, 6];
bus.load = [5, 12, 13];
bus.fault = 1;
bus.slack = 1;

save([path filesep name],'name','Y','Pd','Qd','VM','genParam','bus');
