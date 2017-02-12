function newNode = checkConstraints(state, desiredTaskAcceleration, deltaT, tauLimit, jointLimit, active_joints)

constraintViolated = false;


n_joints = size(active_joints,1);
n_joints_unactive = sum(active_joints(:) ==0);
q = state(1:n_joints)';
qD = state(n_joints+1 : n_joints*2)';


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
    
    
    JbarPinvCurrent = pinv(JbarCurrent);
    
    q2DDCurrent = JbarPinvCurrent * (desiredTaskAcceleration - Jdot*qD + J1*inv(B11)*(C1 + h1));
    
    q1DDCurrent = -inv(B11)*(B12*q2DDCurrent + C1 + h1);
    
    tauCurrent = B21*q1DDCurrent + B22*q2DDCurrent + C2 + h2;

    
    for i=1:size(tauCurrent,1)
        if (tauCurrent(i) > tauLimit)
            constraintViolated = true;
        elseif (tauCurrent(i) < - tauLimit)
            constraintViolated = true;
        end
    end

    
    if (constraintViolated == false)
        
        nodeQD = qD + [q1DDCurrent; q2DDCurrent]*deltaT;
    
        nodeQ = mod(q + qD*deltaT + (1/2)*[q1DDCurrent; q2DDCurrent]*deltaT^2, 2*pi);
    
        %%%CHECK NODE SATURATION
    
        if (constraintViolated == false)
            newNode = [INVorder(nodeQ,active_joints)', INVorder(nodeQD,active_joints)'];
        else
            newNode = 9999;
        end

    else
        newNode = 9999;
    end



end