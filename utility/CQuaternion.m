classdef CQuaternion < handle
    %QUATERNION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
    end
    
    methods (Static)
        %% Set by Quat array
        function Q = SetByArray( Quats )
            Q = Quats / norm(Quats);
        end
        
        %% Display quat, fixed angles, angle & axis
        function disp(Q)
            fixed_angles = CQuaternion.QuatToFixedAngles(Q);
            fixed_angles = rad2deg(fixed_angles);
            [angle, axis] = CQuaternion.GetAngleAxis(Q);
            angle = rad2deg(angle);
            fprintf('Quat: %f, %f, %f, %f\n', Q(1), Q(2), Q(3), Q(4));
            fprintf('Fixed Angles: %f, %f, %f\n', fixed_angles(1), fixed_angles(2), fixed_angles(3));
            fprintf('Angle: %f, Axis: %f, %f, %f\n', angle, axis(1), axis(2), axis(3));
          
        end
    
        %% CQuaternion Multiplication
        function Q = Multiply( Q1, Q2 )
            a = Q1(1);
            u = Q1(2:4);
            t = Q2(1);
            v = Q2(2:4);
            
            Q(1) = a*t - dot(u, v);
            Q(2:4) = a*v + t*u  + cross(u, v);
            
            Q = Q / norm( Q );
        end
        
        %% setting quaternion by fixed angles / XYZ in radian
        function Q = FixedAnglesToQuat( angles )
            R = angles(1);
            P = angles(2);
            Y = angles(3);
            
            cr = cos( 0.5 * R );
            sr = sin( 0.5 * R );
            cp = cos( 0.5 * P );
            sp = sin( 0.5 * P );
            cy = cos( 0.5 * Y );
            sy = sin( 0.5 * Y );
            
            Q(1) = cy*cr*cp + sy*sr*sp;
            Q(2) = cy*sr*cp - sy*cr*sp;
            Q(3) = cy*cr*sp + sy*sr*cp;
            Q(4) = sy*cr*cp - cy*sr*sp;
            
            Q = Q / norm(Q);
        end
        
        %% setting fixed angle roll pitch yaw / XYZ in radians by quaternion
        function FixedAngles = QuatToFixedAngles( Q )
            temp = Q(1)*Q(3) - Q(2)*Q(4);
            if abs(temp) >= 0.4999
                if temp < 0
                    FixedAngles(1) = 0;
                    FixedAngles(2) = -0.5 * pi;
                    FixedAngles(3) = 2 * atan2( Q(2), Q(1) );
                    
                    if FixedAngles(3) > pi
                        FixedAngles(1) = FixedAngles(3) - pi;
                        FixedAngles(3) = pi;
                    elseif FixedAngles(3) < -pi
                        FixedAngles(1) = FixedAngles(3) + pi;
                        FixedAngles(3) = -pi;
                    end
                else
                    FixedAngles(1) = 0;
                    FixedAngles(2) = 0.5 * pi;
                    FixedAngles(3) = -2 * atan2( Q(2), Q(1) );
                    
                    if FixedAngles(3) > pi
                        FixedAngles(1) = FixedAngles(3) + pi;
                        FixedAngles(3) = pi;
                    elseif FixedAngles(3) < -pi
                        FixedAngles(1) = FixedAngles(3) - pi;
                        FixedAngles(3) = -pi;
                    end
                end
            else
                FixedAngles = zeros(1, 3);
                FixedAngles(1) = atan2( 2*( Q(1)*Q(2) + Q(3)*Q(4) ), 1 - 2*( Q(2)^2 + Q(3)^2 ) );
                FixedAngles(2) = asin( 2*( Q(1)*Q(3) - Q(2)*Q(4) ) );
                FixedAngles(3) = atan2( 2*( Q(1)*Q(4) + Q(2)*Q(3) ), 1 - 2*( Q(3)^2 + Q(4)^2 ) );
            end
        end
        
        %% Quaternion conjugate
        function Q = Conj( Q1 )
            Q = Q1 .* [1 -1 -1 -1];
        end
        
        %% Quat difference
        function QDiff = Diff( Q1, Q2 )
            Q1 = CQuaternion.Conj( Q1 );
            QDiff = CQuaternion.Multiply( Q2, Q1 );
        end
        
        %% Convert Quat to Rotation matrix
        function R = Quat2Rot( Q )
            w = Q(1);
            x = Q(2);
            y = Q(3);
            z = Q(4);
            
            R = [ 1-2*y^2-2*z^2   2*(x*y-z*w)     2*(x*z+y*w);
                2*(x*y+z*w)     1-2*x^2-2*z^2   2*(y*z-x*w);
                2*(x*z-y*w)     2*(y*z+x*w)     1-2*x^2-2*y^2;];
        end
        
        %% Convert rotation matrix to quat
        function Q = Rot2Quat( R )
            Q(1) = sqrt( 1 + R(1,1) + R(2,2) + R(3,3) ) / 2;
            Q(2) = ( R(3,2) - R(2,3) ) / ( 4 * Q(1) );
            Q(3) = ( R(1,3) - R(3,1) ) / ( 4 * Q(1) );
            Q(4) = ( R(2,1) - R(1,2) ) / ( 4 * Q(1) );
        end
        
        
        %% Roll pitch yaw transforming to rotation matrix
        function Rot = FixedAngleToRot( angles )
            angles = rad2deg(angles);
            Rot = rotz(angles(3)) * roty(angles(2)) * rotx(angles(1));
        end
        
        %% rotation matrix with self-defined order and fixed or euler method
        %order: [1 pi/3; 2 pi/6;...] n*2 matrix; each row indicates a rotation
        %step; first col = axis; 2nd col = angle;
        %fixedEuler: 0: Fixed angle, 1: Euler angle
        %20161206 Modified: can input syms
        function outM = rotationMod(order,fixedEuler)
            axisCode = 'xyz';
            outM = eye(3);
            for i = 1:size(order,1)
                axis = axisCode(order(i,1));
                angle = order(i,2);
                switch axis
                    case 'x'
                        rotatM = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle);];
                    case 'y'
                        rotatM = [cos(angle) 0 sin(angle);0 1 0;-sin(angle) 0 cos(angle);];
                    case 'z'
                        rotatM = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1;];
                end
                if fixedEuler == 0
                    outM = rotatM*outM; %Fixed
                else
                    outM = outM*rotatM; %Euler
                end
            end
        end
        
        %% rotaion matrix transforming to roll pitch yaw
        function fixed_angles = Rot2FixedAngles( Rot )
            Quat = CQuaternion.Rot2Quat( Rot );
            fixed_angles = CQuaternion.Quat2RPY( Quat ) * 180 / pi;
        end
        
        %% difference between two RPY angles
        % by inputing start and end RPY angles, return the difference between them
        % which can be used when constructing an INCMOVL command
        function angles_diff = RotDiff( a1, a2 )
            Q1 = CQuaternion.FixedAnglesToQuat( a1( 1 ), a1( 2 ), a1( 3 ) );
            Q2 = CQuaternion.FixedAnglesToQuat( a2( 1 ), a2( 2 ), a2( 3 ) );
            
            QDiff = CQuaternion.Diff( Q1, Q2 );
            angles_diff = CQuaternion.QuatToFixedAngles( QDiff );
        end
        
        %% Slerp interpolation of two quats
        function q = SlerpQuat( q0, q1, h )
            Angle = acos( dot( q0, q1 ) );
            q = ( q0 * sin( ( 1 - h ) * Angle ) + q1 * sin( h * Angle ) ) / sin( Angle );    
        end
        
        %% Set by angle and axis
        function q = SetByAngleAndAxis( AngleR, Axis )
           q = [0 0 0 0];
           Axis = Axis / norm(Axis);
           sinHalfAngle = sin( 0.5 * AngleR );
           q(1) = cos( 0.5 * AngleR );
           q(2:end) = Axis' * sinHalfAngle;
        end
        
        %% Get axis from quaternion
        function axis = GetAxis( q )
           angle = acos(q(1)) * 2;
           sinHalfAngle = sin(0.5 * angle);
           axis = q(2:end) / sinHalfAngle;
        end
        
        %% Get angle
        function [angle, axis] = GetAngleAxis( q )
           angle = acos(q(1)) * 2;
           axis = CQuaternion.GetAxis( q );
        end
        
        %% Get Quaternion interpolation 
        function qOutput = GetIntrp( q, r )
           angle = acos(q(1)) * 2;
           if angle == 0 
               qOutput = [1 0 0 0];
               return;
           end
           sinHalfAngle = sin(0.5 * angle);
           axis = q(2:end) / sinHalfAngle;
           intrpAngle = angle * r;
           qOutput = CQuaternion.SetByAngleAndAxis( intrpAngle, axis );
        end
    end
end

