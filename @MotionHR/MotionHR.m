classdef MotionHR
    
    % MOTIONHR Collection of functions for processing Heart rate and IMU data. 
    %
    % Use with raw data URI to raw data folder.
    % For example to target data in folder 'd:\univ\data\HR\raw\01_ABC':
    %  m = MotionHR;
    %  m = m.parseFromFolder('d:\univ\data\HR\raw\01_ABC');
    %  m = m.segmentActivities('d:\Data\Univ\HR\segments_timestamp.csv');
    %
    %   Uses attitude estimation from project "On Attitude Estimation with Smartphones" 
    %   URL: http://tyrex.inria.fr/mobile/benchmarks-attitude
    %
    %   Author: S.Hemminki
    %   All rights reserved.
    %
    %   Date: 26/10/2017
    %   Ver: 0.91
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
                a.ts_start(i) = userSegments(i,3) + 3600 * 3; % Add timezone offset
                a.ts_end(i)   = userSegments(i,4) + 3600 * 3; 
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
        
        % Main pipeline for the analysis
        function obj = analysisPipeline(obj)
            
            % Dataframe length in seconds
            frameLength = 1;
            
            % Overlap in percentage
            overlap     = 50;
            
            % Loop through activities
            for i=obj.activities.id
                
                % Get activity data
                motion   = obj.getDataForActivity(i, obj.motionData);
                polarHR  = obj.getDataForActivity(i, obj.hrData.polar);
                fitbitHR = obj.getDataForActivity(i, obj.hrData.fitbit);
                bandHR   = obj.getDataForActivity(i, obj.hrData.band2);
                
                % Get motionData representations
                t   = motion(:, obj.DataIndex.ts);
                acc = motion(:, obj.DataIndex.accel);
                gyr = motion(:, obj.DataIndex.gyro);
                rep = obj.getMotionRepresentations(t, acc, gyr);
                
                % Extract features
                motionFeatures = obj.computeMotionFeatures(t, rep, frameLength, overlap);
                
                % Analyse HR vs. motion correlations
                
                
                % Get HR error for fitbit and band, assuming Polar is ground truth
                
                
                % Train a model to predict HR Error 
                
                
            end
            
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

        % Get motion representatios
        function representations = getMotionRepresentations(obj, t, acc, gyr)

            % Perform gravity estimation
            grav = obj.getGravityWithMahony(t, acc, gyr);
            
            % Rotate measurements
            racc  = obj.rotateData(grav, acc);
            rlacc = obj.rotateData(grav, acc - grav);
            rgyr  = obj.rotateData(grav, gyr);
            
            % Get Representations
            lacc            = acc-grav;                       % a_lacc
            [~, accPca]     = pca(acc);                       % a_pca                
            accMagn         = getMagnitude(acc);              % a_magn
            [~, laccPca]    = pca(lacc);                      % al_pca1-3
            laccMagn        = getMagnitude(lacc);             % al_magn
            rlacc_v         = rlacc(:,3);                     % ar_v
            [~,rlacc_h_Pca] = pca(rlacc(:,1:2));              % ar_h_pca1-2
            rlacc_h_Magn    = getMagnitude(rlacc(:,1:2));     % ar_h_magn
            [~,totG]        = pca(grav);                      % totgrav
            [~,gyrPca]      = pca(gyr);                       % g_pca1-3
            gyrMagn         = getMagnitude(gyr);              % g_magn
            rgyr_v          = rgyr(:,3);                      % gr_v
            [~,rgyr_h_Pca]  = pca(rgyr(:,1:2));               % gr_h_pca1-2
            rgyr_h_Magn     = getMagnitude(rgyr(:,1:2));      % gr_h_magn

            % Wrap as structure
            representations.a_pca1      = accPca(:,1);
            representations.a_magn      = accMagn;
            representations.al_pca1     = laccPca(:,1);
            representations.al_magn     = laccMagn;
            representations.alr_v       = rlacc_v;
            representations.alr_h_pca1  = rlacc_h_Pca(:,1);
            representations.alr_h_magn  = rlacc_h_Magn;
            representations.totgrav     = totG(:,1);
            representations.g_pca1      = gyrPca(:,1);
            representations.g_magn      = gyrMagn;
            representations.gr_v        = rgyr_v;
            representations.gr_h_pca1   = rgyr_h_Pca(:,1);
            representations.gr_h_magn   = rgyr_h_Magn;
            
        end

        % Gravity estimation with Mahony 
        % Part of project "On Attitude Estimation with Smartphones" 
        % URL: http://tyrex.inria.fr/mobile/benchmarks-attitude
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
          
        % Compute features from the motion representations
        function features = computeMotionFeatures(obj, t, rep, frameLen, overlap)
            
            repPointers = fieldnames(rep); 
            
            for i=1:length(repPointers)
                
                rp = repPointers{i};
                repData = [t,rep.(rp)];
                [frames, t_frames] = obj.splitToFramesByTime(repData, frameLen, overlap);
 
                for j=1:length(frames)
                    frame = frames{j};
                    tdf(j,:) = obj.getTimedomainFeatures(frame);
                    fdf(j,:) = obj.getFreqDomainFeatures(frame, obj.RESAMPLE_RATE);
                end
                
                features.(rp).td = tdf;
                features.(rp).fd = fdf;
                features.(rp).t  = t_frames;
            end
        end
        
        % Split data to frames of specified duration with overlap(in percentage)
        function [frames, t_frames] = splitToFramesByTime(obj, input, frameLen, overlap) 
            t       = input(:, obj.DataIndex.ts);
            data    = input(:, 2:end);
            n       = length(t);
            i_start = 1;
            i       = 1;
            
            while i_start < n
                t_start     = t(i_start);
                i_end       = find(t > (t_start + frameLen), 1, 'first');
                % Omit the last incomplete frame
                if isempty(i_end)
                    break;
                end
                frames{i}   = data(i_start:i_end,:);
                t_frames{i} = t(i_start:i_end);
                % Apply overlap
                step        = length((i_start:i_end)) * (100-overlap)/100; 
                i_start     = floor(i_start + step + 1);
                i           = i+1;
            end     
        end
        
        % Calculate basic timedomain features
        function features = getTimedomainFeatures(obj, input)
            features = [                       ...
                mean(input),                  ...
                median(input),                ...
                std(input),                   ...
                var(input),                   ...
                min(input),                   ...
                max(input),                   ...
                range(input),                 ...
                iqr(input),                   ...
                valueCross(input, 'zero'),    ...
                valueCross(input, 'mean')];  
        end
        
        % Return a set of frequency domain features
        function features = getFreqDomainFeatures(obj, input, samplingRate)
            
            % Zero-padding requires removing mean
            input = input-mean(input);
            
            % Get zero-padded FFT of length 2^nextpow2(input)
            inputL = size(input,1);
            nextpw2 = nextpow2(inputL);
            F = fft(input,2^nextpw2);
            
            % Magnitude of FFT
            F_magn = sqrt( real(F).^2 + imag(F).^2 ) / (inputL);
            fftLen = length(F_magn);
            
            % Remove DC component
            F_magn(1) = [];
            
            % Center point
            center = floor(length(F_magn)/2);
            
            % Resample F_magn at x2 rate
            resampled = zeros(1,center*2-1);
            resampled(1:2:length(resampled)) = F_magn(1:center);
            for i=2:2:length(resampled)
                resampled(i) = (resampled(i-1)+resampled(i+1)) /2;
            end

            % Bin to Hz
            samplesPerBin = round(inputL / (samplingRate)); % 2* because of 0.5hz bins, 1* because of resampling)
            binCount = floor(length(resampled) / samplesPerBin);
            bins = reshape(resampled(1:binCount*samplesPerBin),[samplesPerBin,binCount]);
            
            if samplesPerBin > 1
                energyPer_05Hz = mean(bins);
            else
                energyPer_05Hz = bins;
            end
            
            % energyPer_05Hz = winfeats(F_magn(2:round(fftLen/2)),binsPer05Hz,'trapz');
            Ftot_1to5hz = sum(energyPer_05Hz(1:10));
            Ftot = sum(energyPer_05Hz);
            
            % Relation of frequency components
            Frel_1to10 = energyPer_05Hz(1:10)./Ftot_1to5hz;
            
            % Relation of 1-5Hz to total energy
            Frel_1to10toAll = energyPer_05Hz(1:10)./Ftot;
            
            % Sum of >5s frequencies
            Ftot_5plusHz = sum(energyPer_05Hz(11:end));
            
            % Dominant frequency
            [~,Fmax] = max(energyPer_05Hz);
            
            % Dominant frequency relation to all
            FmaxToAll = energyPer_05Hz(Fmax) / Ftot;
            
            % 50% Energy point
            c = cumsum(energyPer_05Hz);
            halfEnergyIdx = find(c>(Ftot/2),1,'first');
            
            % Spectral entropy
            Q = energyPer_05Hz./Ftot;
            H = Q.*log2(1./Q);
            spectralEntropy = sum(H)/log2(binCount);
            
            % Return result vector
            features = [
                energyPer_05Hz(1:10),   ...
                Frel_1to10,             ...
                Frel_1to10toAll,        ...
                Ftot_1to5hz,            ...
                Ftot_5plusHz,           ...
                Fmax,                   ...
                FmaxToAll,              ...
                halfEnergyIdx,          ...
                Ftot,                   ...
                spectralEntropy];
        end
        
    end
    
end

