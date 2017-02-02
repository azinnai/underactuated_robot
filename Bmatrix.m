function B = Bmatrix(T,N_joints)

syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9 dq10 real

joint_velocities = [dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9 dq10];

joint_velocities = joint_velocities(1:N_joints);

B2 = sym('B2',N_joints);
Bt = sym('Bt',N_joints);



for i = 1:N_joints
   for  k = 1: N_joints
      if k == i 
          B2(k,k) = diff(T,joint_velocities(i),2);
      
      elseif i < k 
          Bt(i,k) = diff(T,joint_velocities(i));
      
      elseif i > k
          B2(i,k) = diff(Bt(k,i),joint_velocities(i));
          B2(k,i) = B2(i,k); 
      end
      
      
             
   end

end



B = simplify(B2);

%B = collect(B2,2);

end


