
function [U,T] = computeEnergies(node)

    m = 0.2;
    l = 0.2;
    I = 0.05;
    lc = 0.1; 
   
    q1 = node(1);
    q2 = node(2);
    q3 = node(3);
    q4 = node(4);
    q5 = node(5);
    q1D = node(6);
    q2D = node(7);
    q3D = node(8);
    q4D = node(9);
    q5D = node(10);

    omega = [q1D; q1D+q2D; q1D+q2D+q3D; q1D+q2D+q3D+q4D; q1D+q2D+q3D+q4D+q5D];
    
    
    %Compact definitions
    s1 = sin(q1);
    s12 = sin(q1+q2);
    s123 = sin(q1+q2+q3);
    s1234 = sin(q1+q2+q3+q4);
    s12345 = sin(q1+q2+q3+q4+q5); 

    c1 = cos(q1);
    c12 = cos(q1+q2);
    c123 = cos(q1+q2+q3);
    c1234 = cos(q1+q2+q3+q4);
    c12345 = cos(q1+q2+q3+q4+q5); 

    %robot link COM
    pc1 = [lc*c1;lc*s1];
    pc2 = [l*c1+lc*c12;l*s1+lc*s12];
    pc3 = [l*c1+l*c12+lc*c123;l*s1+l*s12+lc*s123];
    pc4 = [l*c1+l*c12+l*c123+lc*c1234;l*s1+l*s12+l*s123+lc*s1234];
    pc5 = [l*c1+l*c12+l*c123+l*c1234+lc*c12345;l*s1+l*s12+l*s123+l*s1234+lc*s12345];

    
    Jc1 = Jc1Func(q1);
    Jc2 = Jc2Func(q1,q2);
    Jc3 = Jc3Func(q1,q2,q3);
    Jc4 = Jc4Func(q1,q2,q3,q4);
    Jc5 = Jc5Func(q1,q2,q3,q4,q5);
    
    

    vc1 = Jc1*[q1D;q2D;q3D;q4D;q5D];
    vc2 = Jc2*[q1D;q2D;q3D;q4D;q5D];
    vc3 = Jc3*[q1D;q2D;q3D;q4D;q5D];
    vc4 = Jc4*[q1D;q2D;q3D;q4D;q5D];
    vc5 = Jc5*[q1D;q2D;q3D;q4D;q5D];


    g = [0;9.81];
    U = m*g'*(pc1 + pc2 + pc3 + pc4 + pc5);  


    Tv = (1/2)*m*(vc1'*vc1 + vc2'*vc2 + vc3'*vc3 + vc4'*vc4 + vc5'*vc5);
    Tomega = (1/2)*I*(omega'*omega);

    T = Tv + Tomega;
    
    
end