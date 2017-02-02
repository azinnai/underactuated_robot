function C = Cmatrix(B)

syms q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 real
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9 dq10 real


joint_velocities = [dq1; dq2; dq3; dq4; dq5; dq6; dq7; dq8; dq9; dq10];

joints = [q1 q2 q3 q4 q5 q6 q7 q8 q9 q10];


N_joints = size(B,1);


joint_velocities = joint_velocities(1:N_joints);

joints = joints(1:N_joints);

C = sym('C',[N_joints, 1]);

for i=1:N_joints
    fprintf('JOINT %d ::::::::::::::::::::::\n', i);
    jacob = jacobian(B(:,i),joints)
    specific_Derivative = diff(B,joints(i))
   K = (1/2)* (jacobian(B(:,i),joints) +  jacobian(B(:,i),joints)' - diff(B,joints(i)));
   C(i) = joint_velocities'* K * joint_velocities; 
  %pause;
end

disp('C MATRIX ::::::::::::::::::::::')
C = simplify(C);
end