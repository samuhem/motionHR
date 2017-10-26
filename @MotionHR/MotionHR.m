classdef MotionHR
    
    % MOTIONHR Collection of functions for processing Heart rate and IMU data. 
    %
    % Use with raw data URI to raw data folder.
    % For example to target data in folder 'd:\univ\data\HR\raw\01_ABC':
    %  m = MotionHR;
    %  m = m.parseFromRaw('d:\univ\data\HR\raw\01_ABC');
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
        
        segmentsUri = 'D:\Data\Univ\HR\segments_timestamp.csv';
        
    end
    
    properties
        
        % Offsets for HR monitors
        polarOffset  = -1;
        fitbitOffset = -1;
        
        % Data
        motionData;
        hrData;
        
        % Motion data 
        resamplingInterval;
        
        % Experiment
        id;
        testPattern;
        activities;
        
    end
    
    methods
        
        % =================================================================
        % CONSTRUCTOR
        % =================================================================
        
        % Empty constructor
        function obj = MotionHR()
        end
        
        % =================================================================
        % PARSING (PUBLIC)
        % =================================================================
        
        % Parse from RawFolder
        function obj = parseFromFolder(obj, uri)
           
            files = dir( [uri, '\*.csv'] );
            [obj.id, obj.testPattern] = obj.parseIdAndPatternFromFolderName(uri);
            
            for i=1:length(files)
               
                thisFile = files(i).name;
                absolutePath = [uri, '/', thisFile];
                
                if contains(thisFile, 'acc')
                    % Not required, since acc comes with gyro
                end
                
                if contains(thisFile, 'gyr')
                    
                    accgyro = dlmread(absolutePath, ';');
                    
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
                    obj = obj.addHRData(dlmread(absolutePath, ';'), 'band2');
                end
                
                if contains(thisFile, 'hr_fitbit')
                    obj = obj.parseFitbitOffset(absolutePath);
                    obj = obj.addHRData(readtable(absolutePath), 'fitbit');
                end
                
                if contains(thisFile, 'hr_polar')
                    obj = obj.parsePolarOffset(absolutePath);
                    obj = obj.addHRData(readtable(absolutePath), 'polar');
                end
                
                if contains(thisFile, 'rr_band2')
                    % Not used for now
                end
                
            end
            
        end
        
        % Segment activities
        function obj = segmentActivities(obj, uri)
            
            rawSegments  = dlmread(uri, ',', 1, 0);
            userSegments = rawSegments(rawSegments(:,1) == obj.id, :);
            n = size(userSegments, 1);
            
            for i=1:n
                a.id(i)       = userSegments(i,2);
                a.ts_start(i) = userSegments(i,3);
                a.ts_end(i)   = userSegments(i,4);
            end
            
            obj.activities = a;
        end
        
        % Add motiondata
        function obj = addMotionData(obj, input)
            
            median_dtime = median(diff(input(:, obj.DataIndex.ts)));
            
            if median_dtime > 1
                disp('Warning: Timestamp not in seconds!');
            end
            
            % Add timezone offset
            input(:, obj.DataIndex.ts) = input(:, obj.DataIndex.ts)+3600*3;
            
            obj.motionData = input;
        end
        
        % Add HRrate from provider
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
 
            obj.hrData.(provider) = parsedData;
        end

        % =================================================================
        % ANALYSIS (PUBLIC)
        % =================================================================
        
        % Analyse MotionData
        function obj = analyseMotion(obj)
           
            t   = obj.motionData(:, obj.DataIndex.ts);
            acc = obj.motionData(:, obj.DataIndex.accel);
            gyr = obj.motionData(:, obj.DataIndex.gyro);

            % Perform gravity estimation
            grav = obj.getGravityWithMahony(t, acc, gyr);
            
            % Rotate measurements
            racc = obj.rotateData(grav, acc);
            rgyr = obj.rotateData(grav, gyr);
            
            % Get Representations
            [~, accPca]     = pca(acc);                       % a_pca                
            accMagn         = getMagnitude(acc);              % a_magn
            [~, laccPca]    = pca(acc-grav);                  % al_pca1-3
            laccMagn        = getMagnitude(acc-grav);         % al_magn
            rlacc_v         = racc(:,3);                      % ar_v
            [~,rlacc_h_Pca] = pca(racc(:,1:2));               % ar_h_pca1-2
            rlacc_h_Magn    = getMagnitude(racc(:,1:2));      % ar_h_magn
            [~,totG]        = pca(grav);                      % totgrav
            [~,gyrPca]      = pca(gyr);                       % g_pca1-3
            gyrMagn         = getMagnitude(gyr);              % g_magn
            rgyr_v          = rgyr(:,3);                      % gr_v
            [~,rgyr_h_Pca]  = pca(rgyr(:,1:2));               % gr_h_pca1-2
            rgyr_h_Magn     = getMagnitude(rgyr(:,1:2));      % gr_h_magn

        end
        
    end
        
    methods (Access = private)
        
        % =================================================================
        % PARSING (PRIVATE)
        % =================================================================
        
        % Parse UserID and test order from url
        function [id,pattern] = parseIdAndPatternFromFolderName(obj, uri) 
            uriSplits = strsplit(uri,'\');
            dataFolder = uriSplits{end};
            folderSplits = strsplit(dataFolder, '_');
            id = str2num(folderSplits{1});
            pattern = folderSplits{2};
        end  
        
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
            
            % Get filename from absolutePath
            uriSplits = strsplit(input, '/');
            filename = uriSplits{2};
            
            % Parse from input filename
            strSplits = strsplit(filename,'_');
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
        
        % =================================================================
        % ANALYSIS (PRIVATE)
        % =================================================================
        
        % Get data for activity
        function activityData = getDataForActivity(obj, activityID, data)
            
            i_activity   = find(obj.activities.id == activityID, 1, 'first');
            t_start      = obj.activities.ts_start(i_activity);
            t_end        = obj.activities.ts_end(i_activity);
            
            t            = data(:, obj.DataIndex.ts);
            activityData = data(t>=t_start & t<=t_end, :);
            
        end
        
        % Gravity estimation with Mahony 
        % Part of project "On Attitude Estimation with Smartphones" 
        % URI: http://tyrex.inria.fr/mobile/benchmarks-attitude
        function grav = getGravityWithMahony(obj, t, acc, gyr)
            
            % Create dummy values for missing Magnetometer
            context.magnetic.vector = [0 0 0];
            context.magnetic.declination = 0;
            dummyMagn = zeros(length(acc),3);
            
            % Set reference gravity to 1
            context.gravity.vector = [0 0 1];
            
            % Utilise framework from the referenced project
            q = generateAttitude(t, acc, gyr, dummyMagn, 'mahony', context, 'ENU');
            
            % Map quaterns into gravity vector
            grav = obj.quat2grav(q);
            
        end
  
        % Rotate data to horizontal plane from known euler angles
        function rotated = rotateData(obj, gravity, data)
            
            Gx = gravity(:,1); 
            Gy = gravity(:,2); 
            Gz = gravity(:,3);
            
            Ax = data(:,1); 
            Ay = data(:,2); 
            Az = data(:,3);
            
            rotated = zeros(size(data));
            
            for i=1:length(Ax)
                
                % Calculate roll sin and cos
                Sin_r = Gx(i)/sqrt( Gx(i)^2 + Gz(i)^2 ); % sin of Phi
                Cos_r = Gz(i)/sqrt( Gx(i)^2 + Gz(i)^2 ); % cos of Phi
                
                % Derotate by roll
                Ax_r = Ax(i)*Cos_r - Az(i)*Sin_r;
                Az_r = Ax(i)*Sin_r + Az(i)*Cos_r;
                Gz_r = Gx(i)*Sin_r + Gz(i)*Cos_r;
                
                % Calculate pitch sin and cos
                Sin_p = -Gy(i)/sqrt( Gz_r^2 + Gy(i)^2 ); % sin of Theta
                Cos_p = Gz_r/sqrt( Gz_r^2 + Gy(i)^2 );   % cos of Theta
                
                % Derotate by pitch
                Ay_rp = Ay(i)*Cos_p + Az_r*Sin_p;
                Az_rp = -Ay(i)*Sin_p + Az_r*Cos_p;
                
                rotated(i,:) = [Ax_r, Ay_rp, Az_rp];
            end 
        end
        
        % Get gravity vector from quaternion
        function g = quat2grav(obj, q)
            
            n = size(q, 1);
            g = zeros(n, 3);
            
            for i = 1:n
                g(i,1) = 2 * (q(i,2) * q(i,4) - q(i,1) * q(i,3));
                g(i,2) = 2 * (q(i,1) * q(i,2) + q(i,3) * q(i,4));
                g(i,3) = q(i,1) * q(i,1) - q(i,2) * q(i,2) - q(i,3) * q(i,3) + q(i,4) * q(i,4);
            end
            
        end
  
        
    end
    
end

