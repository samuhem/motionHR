% function ZR=rotz(angle)
% % Returns rotation matrix for rotation about Z-axis
% %
% % Usage: ZR=rotz(angle)
% %
% % angle....Rotation angle (radians)
% %
% % e.g. ZR=rotz(pi/2) % +ve Rotation of 90 Deg around the Z-axis
% %
% % Positive rotation is defined as clockwise as looking from the
% % end of the Z-axis towards the origin.
% 
% ZR=[ cos(degtorad(angle))   sin(degtorad(angle)) 0    
%     -sin(degtorad(angle))   cos(degtorad(angle)) 0
%          0             0     1];
