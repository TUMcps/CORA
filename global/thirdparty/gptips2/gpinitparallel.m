function gp = gpinitparallel(gp)
%GPINITPARALLEL Initialise the Parallel Computing Toolbox.
%
%   GP = GPINITPARALLEL(GP) initialises the Parallel Computing Toolbox if
%   it has been enabled and the Parallel Computing Toolbox license is
%   present.
%
%   IMPORTANT:
%
%   There is a known JVM bug in versions of MATLAB (all platforms) prior to
%   version R2013b (6.3). This causes a failure of the Parallel Computing
%   Toolbox in most cases.
%
%   There is a fix/workaround for this here:
%
%   http://www.mathworks.com/support/bugreports/919688
%
%   Please apply this fix if you are using a version prior to R2013b.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPINIT, GPTOOLBOXCHECK

gp.runcontrol.parallel.ok = false;
parToolbox = gp.info.toolbox.parallel;

if parToolbox && (gp.runcontrol.parallel.auto || gp.runcontrol.parallel.enable)
    mps = getMaxPoolSize;
end

%in auto mode - try to start parallel mode with auto settings but exit and
%proceed in standard mode if no license
if gp.runcontrol.parallel.auto
    if parToolbox
        gp.runcontrol.parallel.enable = true;
        gp.runcontrol.parallel.autosize = true;
        gp.runcontrol.parallel.numWorkers = mps;
    else
        gp.runcontrol.parallel.enable = false;
        return;
    end
end

if gp.runcontrol.parallel.enable
    
    if ~parToolbox
        error('You do not have a license for the Parallel Computing Toolbox or it is not installed. Please set gp.runcontrol.parallel.enable = false.');
    end
    
    ps = getCurrentPoolSize;
    
    %check if pool is already open and of right size
    if (ps > 0) && ( (gp.runcontrol.parallel.autosize && (ps == mps) ) || ...
            (~gp.runcontrol.parallel.autosize && (ps == gp.runcontrol.parallel.numWorkers) ) )
        gp.runcontrol.parallel.ok = true;
        gp.runcontrol.parallel.numWorkers = ps;
        
        if ~gp.runcontrol.quiet
            disp(' ');
            disp(['GPTIPS: proceeding in parallel mode with ' num2str(gp.runcontrol.parallel.numWorkers) ' workers.']);
            disp(' ');
        end
        
    else %otherwise stop then start correct size pool
        
        if ps > 0
            if ~gp.runcontrol.quiet
                disp(' ');
                disp('GPTIPS: stopping existing pool for parallel computation.');
                disp(' ');
            end
            
            if verLessThan('matlab','8.3.0')
                matlabpool close;
            else
                delete(gcp('nocreate'));
            end
        end
        
        try %start new pool
            
            if ~gp.runcontrol.quiet
                disp(' ');
                disp('GPTIPS: attempting to start new pool for parallel computation.');
                disp(' ');
            end
            
            if verLessThan('matlab','8.3.0')
                matlabpool close force local;
                eval(['matlabpool ' num2str(gp.runcontrol.parallel.numWorkers)]);
            else
                delete(gcp('nocreate'));
                
                if gp.runcontrol.parallel.autosize
                    a = parpool;
                    gp.runcontrol.parallel.numWorkers = a.NumWorkers;
                else
                    parpool(gp.runcontrol.parallel.numWorkers);
                end
                
            end
            
            gp.runcontrol.parallel.ok = true;
            
            if ~gp.runcontrol.quiet
                disp(' ');
                disp(['GPTIPS: proceeding in parallel mode with ' num2str(gp.runcontrol.parallel.numWorkers) ' workers.']);
                disp(' ');
            end
            
        catch
            gp.runcontrol.parallel.ok = false;
            
            if ~gp.runcontrol.quiet
                warning('There was a problem starting GPTIPS in parallel computation mode.');
                warning('To run GPTIPS in non-parallel mode set gp.runcontrol.parallel.enable = false in your config file.');
                error(['GPTIPS: could not start pool for parallel computation because: ' lasterr]);
            end
            
        end
    end
    
else
    gp.runcontrol.parallel.ok = false;
end

function x = getCurrentPoolSize
%GETCURRENTPOOLSIZE Get the current parallel pool size. Workaround
%function: calls either matlabpool('size') or pool_size hack depending on
%MATLAB version

if verLessThan('matlab', '7.7.0')
    x = pool_size;
elseif verLessThan('matlab','8.3.0')
    x = matlabpool('size');
else
    pool = gcp('nocreate');
    
    if isempty(pool)
        x = 0;
    else
        if pool.Connected
            x = pool.NumWorkers;
        else
            x = 0;
        end
    end
    
end

function x = getMaxPoolSize
%getMaxPoolSize - try to get the max number of workers in a pool for the
%default profile on the current machine

if verLessThan('matlab','8.3.0')
    c = findResource('scheduler', 'configuration', defaultParallelConfig);
    x = c.ClusterSize;
else
    c = parcluster();
    x = c.NumWorkers;
end

function x = pool_size
%POOL_SIZE - Mathworks hack to return size of current MATLABPOOL. Used for
%versions prior to MATLAB 7.7 (R2008b)
%
% See:
% http://www.mathworks.co.uk/support/solutions/en/data/1-5UDHQP/index.html

session = com.mathworks.toolbox.distcomp.pmode.SessionFactory.getCurrentSession;
if ~isempty( session ) && session.isSessionRunning() && session.isPoolManagerSession()
    client = distcomp.getInteractiveObject();
    if strcmp( client.CurrentInteractiveType, 'matlabpool' )
        x = session.getLabs().getNumLabs();
    else
        x = 0;
    end
else
    x = 0;
end