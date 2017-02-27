function newNode = checkConstraints(state, desiredTaskAcceleration, Knull, deltaTPlanning, deltaT, tauLimit, jointLimit, active_joints)

    constraintViolated = false;
    n_joints = size(active_joints,1);
    n_joints_unactive = sum(active_joints(:) ==0);
    n_joints_active = n_joints - n_joints_unactive;

    q = state(1:n_joints)';
    qD = state(n_joints+1 : n_joints*2)';
    
    epsilon = 0.000000008;


    for i = 0 : deltaT : deltaTPlanning


        B = BFunc(q(2),q(3),q(4),q(5));
        C = CFunc(q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));
        h = hFunc(q(1),q(2),q(3),q(4),q(5));

        B11 = B(1:n_joints_unactive,1:n_joints_unactive);
        B12 = B(1:n_joints_unactive,n_joints_unactive+1:n_joints);
        B21 = B(n_joints_unactive+1:n_joints,1:n_joints_unactive);
        B22 = B(n_joints_unactive+1:n_joints,n_joints_unactive+1:n_joints);

        C1 = C(1:n_joints_unactive);
        C2 = C(n_joints_unactive+1:n_joints);

        h1 = h(1:n_joints_unactive);
        h2 = h(n_joints_unactive+1:n_joints);

        J = JFunc(q(1),q(2),q(3),q(4),q(5));

        J1 = J(:,1:n_joints_unactive);
        J2 = J(:,n_joints_unactive+1:n_joints);  

        Jdot = JdotFunc(q(1),q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));

        %Riordino qua, altrimenti avrei problemi ad assegnare parametri alle
        %funzioni
        q = reorderingMatrix(q, active_joints);
        qD = reorderingMatrix(qD, active_joints);

        JbarCurrent = J2 - J1*(inv(B11)*B12);
        
        try
            JbarPinvCurrent = pinv(JbarCurrent);
        catch
            JbarPinvCurrent = zeros(size(JbarCurrent'));
            constraintViolated = true;
        end    

        
        if (Knull > 0)
            nullSpaceAcceleration = zeros(n_joints_active,1);
            for p=1:n_joints_active

                qPlus = q;
                qPlus(p + n_joints_unactive) = qPlus(p + n_joints_unactive) + epsilon;
                qMinus = q;
                qMinus(p + n_joints_unactive) = qMinus(p + n_joints_unactive) - epsilon;

                qPlus = INVorder(qPlus, active_joints);
                qMinus = INVorder(qMinus, active_joints);

                JbarPlus = JbarFunc(qPlus(1),qPlus(2),qPlus(3),qPlus(4),qPlus(5));
                JbarMinus = JbarFunc(qMinus(1),qMinus(2),qMinus(3),qMinus(4),qMinus(5));


                HPlus = sqrt(det(JbarPlus*JbarPlus'));
                HMinus = sqrt(det(JbarMinus*JbarMinus'));

                gradHp = (HPlus - HMinus)/(2*epsilon);

                nullSpaceAcceleration(p) = gradHp;
            end
            
            projectedGradient = Knull*(eye(n_joints_active) - JbarPinvCurrent*JbarCurrent)*nullSpaceAcceleration;
            
        else
            projectedGradient = zeros(n_joints_active,1);
        end
        
        
        
        q2DDCurrent = JbarPinvCurrent * (desiredTaskAcceleration - Jdot*qD + J1*inv(B11)*(C1 + h1)) + projectedGradient;

        q1DDCurrent = -inv(B11)*(B12*q2DDCurrent + C1 + h1);

        tauCurrent = B21*q1DDCurrent + B22*q2DDCurrent + C2 + h2;


        for k=1:size(tauCurrent,1)
            if (tauCurrent(k) > tauLimit)
                constraintViolated = true;
                disp('tau')
                tauCurrent(k)
            elseif (tauCurrent(k) < - tauLimit)
                constraintViolated = true;
                disp('tau')
                tauCurrent(k)
            end
        end


            qD = qD + [q1DDCurrent; q2DDCurrent]*deltaT;

            q = mod(q + qD*deltaT + (1/2)*[q1DDCurrent; q2DDCurrent]*deltaT^2, 2*pi);

            for j=1:size(q,1)
                if ((q(j) < pi + jointLimit) && (q(j) > pi -jointLimit) && (active_joints(j) == 1))
                    constraintViolated = true;
                    disp('joint limit')
                    q(j)
                end
            end
           
            
            q = INVorder(q,active_joints);
            qD = INVorder(qD,active_joints);
            
            if(constraintViolated) 
                break;
            end
            
    end
    
    if (constraintViolated == false)
        newNode = [INVorder(q,active_joints)', INVorder(qD,active_joints)'];
    else
        newNode = 9999;
    end
    
end