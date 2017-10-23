% read acceleration and gyration data in
acc = dlmread('./datasets/test/20171010_02_acc_band2.csv',';');
gyr = dlmread('./datasets/test/20171010_02_gyr_band2.csv',';');

timestamp = acc(:, 1);
acc = acc(:, [2 3 4])*9.81;     % convert from g unit to m/s^2
gyr = gyr(:, [5 6 7]); 
mag = zeros(size(acc));
linear_acc = zeros(size(acc));
angle = zeros(size(acc));

% construct context info
location = struct('latitude',60.1699, 'longitude', 24.9384, 'altitude', 50);    % geo information of Helsinki
date = struct('year', 2017, 'month', 10, 'day', 10);
coordinateSystem = 'ned';
algorithm = 'Mahony';

context = createContextFromLocationAndDate(location, date, coordinateSystem);

attitude = generateAttitude(timestamp, acc, gyr, mag, algorithm, context, coordinateSystem);

% rotate acc vector by quarternions and convert quaternion to angle
for i = 1:length(acc)
    q = attitude(i, :);
    v = acc(i, :);
    linear_acc(i, :) = quatrotate(q, v);
%     linear_acc(i, :) = qvrot(q, v);
    angle(i, :) = toEulerAngle(q);
end
