function [time,data] = usbaqc(varargin)
%% Acquire Data Using NI Devices
%
% % This example shows how to acquire data from a National Instruments 
% % device available to MATLAB(R) from the command line using the Session 
% % based interface.
%
% % EXAMPLE CODE:
% % set params
% params.background        = 0;    % [1] background, [0] foreground
% params.debug             = 0;    % [1] print extra info, [0] do not
% params.plot              = 1;    % [1] plot traces], [0] no not
% params.ScannerPin        = 13;   % to what pin is scanner connected
% params.ECGPin            = 6;    % to what 
% params.Rate              = 2000; % [Hz] sampling rate
% params.DurationInSeconds = 17;   % [sec] sampling duration
% [time,data] = usbaqc(params);
%
% [time,data] = usbaqc(); % Use default params

    global data
    global time

    if nargin == 1 && isstruct(varargin{1})
        % If one argument, and it is a struct, use it as params
        params = varargin{1};
    elseif mod(nargin,2) == 1
        % If even number of arguments, error
        error('Argments should be a single struct or name-value pairs')
    elseif nargin == 0
        % If no arguments, use defaults
        params.background        = 0;
        params.debug             = 0;
        params.plot              = 0;
        params.ScannerPin        = 13;
        params.ECGPin            = 6;
        params.Rate              = 125000;
        params.DurationInSeconds = 15;
    else
        % Else there are an even number of pairs
        % Loop through name-value pairs, and set params
        for argi=1:2:nargin
            if ~ischar(varargin{argi}), error('name then value'), end
            disp(['params.' varargin{argi} ' = ' num2str(varargin{argi+1})])
        end
    end
    
    %% Discover Analog Input Devices
    % To discover a device that supports analog input subsystems, click the
    % name of the device in the list in the Command window, or access the
    % device in the array returned by |daq.getDevices| command.
    devices = daq.getDevices;
    if isfieldtrue(params,'debug')
        disp(devices)
    end
    
    %% Create a Session and Add an Analog Input Channel
    % Create a session, and use the |addAnalogInputChannel| function to add two
    % analog input channels from this device to the session. 
    s = daq.createSession('ni');
    s.addAnalogInputChannel('Dev1', params.ScannerPin, 'Voltage');
    s.addAnalogInputChannel('Dev1', params.ECGPin, 'Voltage');
    if isfieldtrue(params,'debug')
        disp(s)
    end
    
    %% Set Session Rate
    % By default the session is configured for 1000 scans/second. 
    % Change the scan rate to acquire at 8000 scans / second.
    s.Rate = params.Rate;

    %% Acquire data for a Specified Duration  
    % By default the session runs for a duration of one second. Configure the
    % session to run for a specified duration using the |DurationInSeconds|
    % property.  
    %
    % Note, that Rate = DurationInSeconds/NumberOfScans. Changing the duration
    % changes the number of scans in the session accordingly. The acquisition
    % will now run for two seconds, acquiring 16,000 scans.
    s.DurationInSeconds = params.DurationInSeconds;

    %% Acquire a Single Scan
    % Use the |inputSingleScan| function to acquire a single scan. The result
    % is an array of size M-by-1, where M corresponds to the number of input
    % channels.
    if isfieldtrue(params,'debug')
        disp(s.inputSingleScan)
    end

    %% Start an acquisition
    if isfieldtrue(params,'background')
        close all;clf
        lh = s.addlistener('DataAvailable',@plotData);
        s.startBackground();
 
        % TODO: figure out how to use this properly...
        % Do something
        counter = 0;
        fprintf('%s ','Seconds of acq.:');
        while(~s.IsDone)
            pause(1)
            counter = counter + 1;
            fprintf('%d ',counter);
        end
        fprintf('\n ');

        % Delete the listener
        delete(lh);
    else
        % Start the Session in Foreground
        % You can acquire multiple scans using the |startForeground| function. 
        % This blocks MATLAB execution until all the data is acquired. The acquired
        % data is returned in TIME-DATA pairs. TIME is a M-by-1 matrix, where M is
        % the number of scans. DATA is a M-by-N matrix where M is the number of
            % scans and N is the number of analog input channels in the session.
        % We will now start a foreground operation.
        [data,time] = s.startForeground;
    end
    
    %% After acquisition of scans, Plot the Acquired Data
    if isfieldtrue(params,'plot')
        % This session has now acquired the scans. 
        % Plot the acquired data for each channel.
        plot(time,data);       % plot global data
        xlabel('Time (secs)');
        ylabel('Voltage')
    end
end

function plotData(src,event)
%% Listener function for appending new background data to globals
    persistent tempData tempTime
    global data time
    if(isempty(tempData)), tempData = []; end
    if(isempty(tempTime)), tempTime = []; end

    %plot(event.TimeStamps, event.Data)
    
    tempData = [tempData;event.Data];
    tempTime = [tempTime;event.TimeStamps];
    
    plot(tempTime, tempData)
    
    data = tempData;
    time = tempTime;
end

function [state] = isfieldtrue(structval,fieldName)
%% function to test if a field exists AND is true
%
% test.hi = 1; test.lo = 0;
% state1=isfieldtrue(test,'hi')
% state2=isfieldtrue(test,'lo')
% state3=isfieldtrue(test,'no')
    state = isfield(structval,fieldName) && eval(['structval.' fieldName]);
end
