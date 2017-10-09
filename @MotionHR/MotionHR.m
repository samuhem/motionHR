classdef MotionHR
    
    % MOTIONHR Collection of functions for processing Heart rate and IMU data. 
    %
    % Use with raw data URI to raw data folder (e.g., D:\univ\data\HR\raw\01_ABC):
    %  m = MotionHR;
    %  m = m.parseFromRaw('D:\univ\data\HR\raw\01_ABC'); % Replace URI folder here
    %
    %   Author: S.Hemminki
    %   All rights reserved.
    %
    %   Date: 09/10/2017
    %   Ver: 0.9
    %
    
    properties(Constant)
       
        DataIndex =                                  ...
            struct('ts',            1,               ...
                   'accX',          2,               ...
                   'accY',          3,               ...
                   'accZ',          4,               ...
                   'gyroX',         5,               ...
                   'gyroY',         6,               ...
                   'gyroZ',         7,               ...
                                                     ...
                   'accel',         2:4,             ...
                   'gyro',          5:7,             ...
                   'accelWithTime', 1:4,             ...
                   'gyroWithTime',  [1, 5:7]); 
               
        HRConstants =                                ...
            struct('ts',            1,               ...
                   'HR_BPM',        2);
               
        RESAMPLE_RATE = 50;
        
    end
    
    properties
        
        % Offsets for HR monitors
        polarOffset  = -1;
        fitbitOffset = -1;
        
        % Data
        MotionData;
        HRData;
        
        % Motion data 
        resamplingInterval;
        
        % Gravity Estimation
        initG;
        
    end
    
    methods
        
        % Empty constructor
        function obj = MotionHR()
        end
        
        % Parse from RawFolder
        function obj = parseFromRaw(obj, url)
           
            files = dir( [url, '\*.csv'] );
            
            for i=1:length(files)
               
                thisFile = files(i).name;
                
                if contains(thisFile, 'acc')
                    % Not required, since acc comes with gyro
                end
                
                if contains(thisFile, 'gyr')
                    
                    accgyro = dlmread(thisFile, ';');
                    
                    % Gyro is is degrees/sec, should be rad/sec
                    accgyro(:, obj.DataIndex.gyro) = deg2rad(accgyro(:, obj.DataIndex.gyro));
                    
                    % Move timestamps from ms --> sec
                    accgyro(:, obj.DataIndex.ts) = accgyro(:, obj.DataIndex.ts) / 1e3;
                    
                    % Resample motion sensors at constant 50hz
                    data = accgyro(:, [obj.DataIndex.accel, obj.DataIndex.gyro]);
                    time = accgyro(:, obj.DataIndex.ts);
                    [rs_data, rs_ts] = resample(data, time, obj.RESAMPLE_RATE);
                    cleanData = [rs_ts, rs_data];
                    
                    % Insert into MotionHR
                    obj = obj.addMotionData(cleanData);
                    
                end
                
                if contains(thisFile, 'hr_band')
                    obj = obj.addHRData(dlmread(thisFile, ';'), 'band2');
                end
                
                if contains(thisFile, 'hr_fitbit')
                    obj = obj.parseFitbitOffset(thisFile);
                    obj = obj.addHRData(readtable(thisFile), 'fitbit');
                end
                
                if contains(thisFile, 'hr_polar')
                    obj = obj.parsePolarOffset(thisFile);
                    obj = obj.addHRData(readtable(thisFile), 'polar');
                end
                
                if contains(thisFile, 'rr_band2')
                    % Not used for now
                end
                
            end
            
        end
        
        % Add motiondata
        function obj = addMotionData(obj, input)
            
            median_dtime = median(diff(input(:, obj.DataIndex.ts)));
            
            if median_dtime > 1
                disp('Warning: Timestamp not in seconds!');
            end
            
            % Add timezone offset
            input(:, obj.DataIndex.ts) = input(:, obj.DataIndex.ts)+3600*3;
            
            obj.MotionData = input;
        end
        
        % Add HR rata from provider
        function obj = addHRData(obj, data, provider)
            
            switch(provider)
                case('polar')
                    parsedData = obj.parsePolar(data);
                    
                case('fitbit')
                    parsedData = obj.parseFitbit(data);

                case('band2')
                    parsedData(:, obj.HRConstants.ts)     = data(:, obj.HRConstants.ts) ./ 1e3 + 3600 * 3; % Add timezone offset
                    parsedData(:, obj.HRConstants.HR_BPM) = data(:, obj.HRConstants.HR_BPM);
            end
 
            obj.HRData.(provider) = parsedData;
        end

    end
        
    methods (Access = private)
        
        % Parse offset time for Polar input table
        function obj = parsePolarOffset(obj, input)
            
            % Load data from the file
            fid = fopen(input);
            % nLines = countLines(fid);
            % lines = cell(nLines,1);
            line = 'first';
            i=1;
            while ischar(line) && i < 3
                line = fgetl(fid);
                if ischar(line)
                    headerLines{i} = line;
                end
                i = i+1;
            end
            fclose(fid);
            
             % This contains the start time
            strSplits = strsplit(headerLines{2}, ',');
            date  = strSplits{3};
            time = strSplits{4};
            
            obj.polarOffset = posixtime(datetime([date, ' ', time], 'InputFormat', 'dd-MM-yyyy HH:mm:ss'));
            
        end
        
        % Parse offset time from fileName for Fitbit
        function obj = parseFitbitOffset(obj, input)
            
            % Parse from input filename
            strSplits = strsplit(input,'_');
            year      = strSplits{1}(1:4);
            month     = strSplits{1}(5:6);
            day       = strSplits{1}(7:8);
            
            obj.fitbitOffset = posixtime(datetime(['', year, '-', month, '-', day, ' 00:00:00']));
            
        end
        
        % Parse HR data from Polar format
        function data = parsePolar(obj, input)
        
            ts = input.Time;
            
            for i=1:length(ts)
                tsSplits = strsplit(ts{i},':');
                
                tsInSeconds = str2num(tsSplits{1})*3600   + ...
                              str2num(tsSplits{2})*60     + ...
                              str2num(tsSplits{3}); 
                             
                tsUTC(i) = tsInSeconds + obj.polarOffset;
                
            end

            data(:,obj.HRConstants.ts) = tsUTC;
            data(:,obj.HRConstants.HR_BPM) = input.HR_bpm_;
            
        end
        
        % Parse HR data from Fitbit format
        function data = parseFitbit(obj, input)
            
            ts = input.Time;
            
            for i=1:length(ts)
                tsSplits = strsplit(ts{i},':');
                
                tsInSeconds = str2num(tsSplits{1})*3600   + ...
                              str2num(tsSplits{2})*60     + ...
                              str2num(tsSplits{3}); 
                             
                tsUTC(i) = tsInSeconds + obj.fitbitOffset;
                
            end
            
            data(:,obj.HRConstants.ts) = tsUTC;
            data(:,obj.HRConstants.HR_BPM) = input.HeartRate_beats_min_;
            
        end
        
    end
    
end

