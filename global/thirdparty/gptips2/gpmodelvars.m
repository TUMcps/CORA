function hvec=gpmodelvars(gp,ID)
%GPMODELVARS Display the frequency of input variables present in the specified model.
%
%   GPMODELVARS(GP,ID) displays input frequency for the model with numeric
%   identifier ID in the GPTIPS datastructure GP.
%
%   GPMODELVARS(GP,'best') does the same for the 'best' model of the run
%   (as evaluated on the training data).
%
%   GPMODELVARS(GP,'valbest') does the same for the model that performed
%   best on the validation data (if this data exists).
%
%   GPMODELVARS(GP,'testbest') does the same for the model that performed
%   best on the test data (if this data exists).
%
%   GPMODELVARS(GP,GPMODEL) operates on the GPMODEL struct representing a
%   multigene regression model, i.e. the struct returned by the functions
%   GPMODEL2STRUCT or GENES2GPMODEL.
%
%   HITVEC = GPMODELVARS(GP,'valbest') returns a frequency vector and
%   suppresses graphical output.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPPOPVARS

hitvec = [];

if nargin < 2
    disp('Usage is GPMODELVARS(GP,ID) where ID is the numeric identifier of the selected individual.');
    disp('or GPMODELVARS(GP,''BEST'') to run the best individual of the run');
    disp('or GPMODELVARS(GP,''VALBEST'') to run the best validation individual of the run.');
    disp('or GPMODELVARS(GP,''TESTBEST'') to run the best test individual of the final population.');
    return;
end

if nargout < 1
    graph = true;
else
    graph = false;
end

if isfield(gp.userdata,'name') && ~isempty(gp.userdata.name)
    setname = ['Data: ' gp.userdata.name];
else
    setname = '';
end

%get a specified individual from the population and scan it for input
%variables
if isnumeric(ID)
    
    if ID < 1 || ID > numel(gp.pop)
        disp('Invalid value of supplied numerical identifier ID.');
        return;
    end
    
    model = gp.pop{ID};
    title_str = ['Input frequency in individual with ID: ' int2str(ID)];
    
elseif ischar(ID) && strcmpi(ID,'best')
    
    model = gp.results.best.individual;
    title_str = 'Input frequency in best individual.';
    
elseif ischar(ID) && strcmpi(ID,'valbest')
    
    % check that validation data is present
    if ~isfield(gp.results,'valbest')
        error('No validation data was found. Try GPMODELVARS(GP,''BEST'') instead.');
    end
    
    model = gp.results.valbest.individual;
    title_str = 'Input frequency in best validation individual.';
    
elseif ischar(ID) && strcmpi(ID,'testbest')
    
    % check that test data is present
    if ~isfield(gp.results,'testbest')
        error('No test data was found. Try GPMODELVARS(GP,''BEST'') instead.');
    end
    
    model = gp.results.testbest.individual;
    title_str = 'Input frequency in best test individual.';
    
    %for a supplied cell array of encoded tree expressions
elseif iscell(ID)
    model = ID;
    title_str = ('Input frequency in user model.');

elseif isstruct(ID) && isfield(ID,'genes')
     model = ID.genes.geneStrs;
     title_str = ('Input frequency in user model.');
     
else
    error('Illegal argument');
end

%perform scan
numx = gp.nodes.inputs.num_inp;
hitvec = scangenes(model,numx);

% plot results as barchart
if graph
    h = figure;
    set(h,'name','GPTIPS 2 single model input frequency','numbertitle','off');
    a = gca;
    bar(a,hitvec);shading faceted;
    
    if ~verLessThan('matlab','8.4') %R2014b
        a.Children.FaceColor = [0 0.45 0.74];
        a.Children.BaseLine.Visible = 'off';
    else
        b = get(a,'Children');
        set(b,'FaceColor',[0 0.45 0.74],'ShowBaseLine','off');
    end
    axis tight;
    xlabel('Input');
    ylabel('Input frequency');
    title({setname, title_str},'Fontweight','bold');
    grid on;
end

hitvec = hitvec';

if nargout > 0
    hvec = hitvec;
end