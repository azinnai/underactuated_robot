function [q1DD, q2DD, constraintsViolated] = computeAcceleration(q, qD, active_joints, desiredTaskAcceleration, constraintsViolated)

    n_joints = size(active_joints,1);
    n_joints_unactive = sum(active_joints(:) ==0);


    B = BFunc(q(2),q(3),q(4),q(5));
    C = CFunc(q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));
    h = hFunc(q(1),q(2),q(3),q(4),q(5));
    
    B11 = B(1:n_joints_unactive,1:n_joints_unactive);
    B12 = B(1:n_joints_unactive,n_joints_unactive+1:n_joints);


    C1 = C(1:n_joints_unactive);


    h1 = h(1:n_joints_unactive);


    J = JFunc(q(1),q(2),q(3),q(4),q(5));
    
    J1 = J(:,1:n_joints_unactive);
    J2 = J(:,n_joints_unactive+1:n_joints);  

    Jdot = JdotFunc(q(1),q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));
    
    JbarCurrent = J2 - J1*(inv(B11)*B12);
    try 
        JbarPinvCurrent = pinv(JbarCurrent);
    catch
        JbarPinvCurrent = zeros(size(JbarCurrent'));
        constraintsViolated = true;
    end
    
    q2DD = JbarPinvCurrent * (desiredTaskAcceleration - Jdot*qD + J1*inv(B11)*(C1 + h1));
    
    q1DD = -inv(B11)*(B12*q2DD + C1 + h1);
    

    

end