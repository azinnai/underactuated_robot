function C = Cmatrix(B)

syms q1 q2 q3 q4 q5 q1D q2D q3D q4D q5D real;


joint_velocities = [q1D;q2D;q3D;q4D;q5D];

joints = [q1; q2; q3; q4; q5];


N_joints = size(B,1);


joint_velocities = joint_velocities(1:N_joints);

joints = joints(1:N_joints);

C = sym('C',[N_joints, 1]);

for i=1:N_joints
   K = (1/2)* (jacobian(B(:,i),joints) +  jacobian(B(:,i),joints)' - diff(B,joints(i)));
   C(i) =joint_velocities' * K * joint_velocities; 
end

C = simplify(C);
end