function angle=toEulerAngle(q)
	% roll (x-axis rotation)
	sinr = +2.0 * (q(1) * q(2) + q(3) * q(4));
	cosr = +1.0 - 2.0 * (q(2) * q(2) + q(3) * q(3));
	roll = atan2(sinr, cosr);

	% pitch (y-axis rotation)
	sinp = +2.0 * (q(1) * q(3) - q(4) * q(2));
        if abs(sinp) >= 1
            pitch = 1.57*sign(sinp); % use 90 degrees if out of range
        else
            pitch = asin(sinp);
        end

	% yaw (z-axis rotation)
	siny = +2.0 * (q(1) * q(4) + q(2) * q(3));
	cosy = +1.0 - 2.0 * (q(3) * q(3) + q(4) * q(4));  
	yaw = atan2(siny, cosy);
    angle = [pitch*57.2 roll*57.2 yaw*57.2];
end