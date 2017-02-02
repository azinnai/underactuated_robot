function B = Bmatrix(T,N_joints)

syms q1D q2D q3D q4D q5D real;

joint_velocities = [q1D;q2D;q3D;q4D;q5D];


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

B = simplify(collect(B2,5));

end


