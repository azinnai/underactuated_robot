function newNode = checkConstraintsRK(state, desiredTaskAcceleration, deltaT, tauLimit, jointLimit, active_joints)

    constraintsViolated = false;
    
    n_joints = size(active_joints,1);
    n_joints_unactive = sum(active_joints(:) ==0);
    q = state(1:n_joints)';
    qD = state(n_joints+1 : n_joints*2)';
    
    B = BFunc(q(2),q(3),q(4),q(5));
    C = CFunc(q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));
    h = hFunc(q(1),q(2),q(3),q(4),q(5));
    
    B21 = B(n_joints_unactive+1:n_joints,1:n_joints_unactive);
    B22 = B(n_joints_unactive+1:n_joints,n_joints_unactive+1:n_joints);
    
    C2 = C(n_joints_unactive+1:n_joints);
    
    h2 = h(n_joints_unactive+1:n_joints);
    
    qk1 = q;
    qDk1 = qD;
    
    q = reorderingMatrix(q, active_joints);
    qD = reorderingMatrix(qD, active_joints);

    [q1DDk1, q2DDk1, constraintsViolated] = computeAcceleration(qk1, qDk1, active_joints, desiredTaskAcceleration,constraintsViolated);
    if (constraintsViolated == false)
        qDDk1 = [q1DDk1;q2DDk1];

        qDk1 = reorderingMatrix(qDk1, active_joints);

        qDk2 = qD + 0.5*qDDk1*deltaT;
        qk2 = q + 0.5*qDk1*deltaT;

        qk2 = INVorder(qk2, active_joints);
        qDk2 = INVorder(qDk2, active_joints);

        [q1DDk2, q2DDk2, constraintsViolated] = computeAcceleration(qk2, qDk2, active_joints, desiredTaskAcceleration,constraintsViolated);
        if (constraintsViolated == false)
            qDDk2 = [q1DDk2;q2DDk2];


            qDk2 = reorderingMatrix(qDk2, active_joints);

            qDk3 = qD + 0.5*qDDk2*deltaT;
            qk3 = q + 0.5*qDk2*deltaT;

            qk3 = INVorder(qk3, active_joints);
            qDk3 = INVorder(qDk3, active_joints);

            [q1DDk3, q2DDk3, constraintsViolated] = computeAcceleration(qk3, qDk3, active_joints, desiredTaskAcceleration,constraintsViolated);
            if (constraintsViolated == false)
                qDDk3 = [q1DDk3;q2DDk3];


                qDk3 = reorderingMatrix(qDk3, active_joints);

                qDk4 = qD + qDDk3*deltaT;
                qk4 = q + qDk3*deltaT;

                qk4 = INVorder(qk4, active_joints);
                qDk4 = INVorder(qDk4, active_joints);

                [q1DDk4, q2DDk4, constraintsViolated] = computeAcceleration(qk4, qDk4, active_joints, desiredTaskAcceleration,constraintsViolated);
                if (constraintsViolated == false)
                    qDk4 = reorderingMatrix(qDk4, active_joints);
                    
                    
                    q1DDRK = (1/6) * (q1DDk1 + q1DDk2 +q1DDk3 + q1DDk4);
                    q2DDRK = (1/6) * (q2DDk1 + q2DDk2 +q2DDk3 + q2DDk4);
                 

                    %Once I get the RK accelerations, I have to check if
                    %is admissibile according to the torque constraints.
                    
                    tau = B21*q1DDRK + B22*q2DDRK + C2 + h2;
                    for i=1:size(tau,1)
                       if (tau(i) > tauLimit)
                            constraintsViolated = true;
                        elseif (tau(i) < - tauLimit)
                            constraintsViolated = true;
                       end
                    end
                    
                    qDRK = (1/6) * (qDk1 + 2*qDk2 + 2*qDk3 + qDk4);

                    nodeQD = qD + (deltaT)*[q1DDRK;q2DDRK];
                    nodeQ = q + (deltaT)*qDRK;


                    for i=2:size(nodeQ,1)
                        if ((nodeQ(i) < pi + jointLimit) && (nodeQ(i) > pi -jointLimit))
                            constraintsViolated = true;
                        end
                    end
                    
                end
            end
        end
    end


    if (constraintsViolated == false)
        newNode = [INVorder(nodeQ,active_joints)', INVorder(nodeQD,active_joints)'];

    else
            newNode = 9999;
    end



end