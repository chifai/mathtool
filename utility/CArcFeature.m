classdef CArcFeature < handle
    properties (SetAccess = private)
        m_TotalQuat;
        m_SwingQuat;
        m_TwistQuat;
        m_TwistQuat2;
        m_ArcAngle;
        m_Radius;
        m_AngleSpan;
        m_FrameZ;
        m_Center = [0, 0];
    end
    
    methods
        function obj = CArcFeature(MidPosition, EndPosition, EulerDeg1, EulerDeg2)
            % mid and end position are at 2D plane, starting position is
            % (0,0)
            x2 = MidPosition(1);
            y2 = MidPosition(2);
            x3 = EndPosition(1);
            y3 = EndPosition(2);
            det = x2*y3 - x3*y2;
            if det == 0
                error( 'Not a valid circle' );
            end
            obj.m_Center(1) = ( (x2^2+y2^2) * y3 - (x3^2+y3^2) * y2 ) / ( 2 * det );
            obj.m_Center(2) = ( (x2^2+y2^2) * -x3 + (x3^2+y3^2) * x2 ) / ( 2 * det );
            
            % use vector to calculate radius, arc angle etc
            v1 = [MidPosition 0];
            v2 = [EndPosition 0] - [MidPosition 0];
            v3 = [-EndPosition 0];
            
            obj.m_ArcAngle = acos( dot(-v1, v2) / norm(v1) / norm(v2) );
            obj.m_Radius = norm(v3) / (2 * sin(obj.m_ArcAngle));
            obj.m_AngleSpan = 2 * ( pi - obj.m_ArcAngle );
            
            % obtain local Z frame
            xAxis = v1' / norm(v1);
            zAxis = cross(xAxis, v2' / norm(v2));
            zAxis = zAxis / norm(zAxis);
            yAxis = cross(zAxis, xAxis);
            obj.m_FrameZ = [xAxis yAxis zAxis];
            
            % Total Quat
            Q1 = CQuaternion.RPY2Quat(EulerDeg1(1), EulerDeg1(2), EulerDeg1(3));
            Q2 = CQuaternion.RPY2Quat(EulerDeg2(1), EulerDeg2(2), EulerDeg2(3));
            obj.m_TotalQuat = CQuaternion.Diff(Q1, Q2);
            
            % Swing Quat
            obj.m_SwingQuat = CQuaternion.SetByAngleAndAxis( obj.m_AngleSpan, zAxis );
            
            % Check Total and Swing Quat axis should be at same direction
            % i.e. dot product of two axes >= 0
            TotalAxis = CQuaternion.GetAxis( obj.m_TotalQuat );
            SwingAxis = CQuaternion.GetAxis( obj.m_SwingQuat );
            
            if dot( TotalAxis, SwingAxis ) < 0
                obj.m_TotalQuat = -1 * obj.m_TotalQuat;
            end
            
            % Twist Quat
            swingConj = CQuaternion.Conj( obj.m_SwingQuat );
            obj.m_TwistQuat = CQuaternion.Multiply( swingConj, obj.m_TotalQuat );
            
            % Alternative TwistQuat2
            TwistAxis = CQuaternion.GetAxis( obj.m_TwistQuat );
            TwistAxis2 = -TwistAxis;
            angle2 = 2 * pi - acos(obj.m_TwistQuat(1)) * 2;
            obj.m_TwistQuat2 = CQuaternion.SetByAngleAndAxis( angle2, TwistAxis2 );
            
            displayFeature(obj);
        end
        
        function displayFeature(obj)
           fprintf( 'Center: [%.2f, %.2f]\n', obj.m_Center(1), obj.m_Center(2) );
           fprintf( 'ArcAngle: %f\n', obj.m_ArcAngle*180/pi );
           fprintf( 'Radius: %f\n', obj.m_Radius );
           fprintf( 'AngleSpan: %f\n', obj.m_AngleSpan*180/pi );
           fprintf( 'FrameZ:\n' );
           disp( obj.m_FrameZ );
           disp( 'Total Quat:' );
           disp( obj.m_TotalQuat );
           disp( 'Swing Quat:' );
           disp( obj.m_SwingQuat );
           disp( 'Twist Quat:' );
           disp( obj.m_TwistQuat );
           disp( 'Twist Quat2:' );
           disp( obj.m_TwistQuat2 );
        end
    end
    
    methods (Static)
        
    end
end
