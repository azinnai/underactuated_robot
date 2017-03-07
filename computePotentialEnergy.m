
function U = computePotentialEnergy(node)

%    m = 0.2;
%    l = 0.2;
%    I = 0.05;
%    lc = 0.1; 
   
%    q1 = node(1);
%    q2 = node(2);
%    q3 = node(3);
%    q4 = node(4);
%    q5 = node(5);



    
    
    %Compact definitions
%    s1 = sin(q1);
%    s12 = sin(q1+q2);
%    s123 = sin(q1+q2+q3);
%    s1234 = sin(q1+q2+q3+q4);
%    s12345 = sin(q1+q2+q3+q4+q5); 

%    c1 = cos(q1);
%    c12 = cos(q1+q2);
%    c123 = cos(q1+q2+q3);
%    c1234 = cos(q1+q2+q3+q4);
%    c12345 = cos(q1+q2+q3+q4+q5); 

    %robot link COM
%    pc1 = [lc*c1;lc*s1];
%    pc2 = [l*c1+lc*c12;l*s1+lc*s12];
%    pc3 = [l*c1+l*c12+lc*c123;l*s1+l*s12+lc*s123];
%    pc4 = [l*c1+l*c12+l*c123+lc*c1234;l*s1+l*s12+l*s123+lc*s1234];
%    pc5 = [l*c1+l*c12+l*c123+l*c1234+lc*c12345;l*s1+l*s12+l*s123+l*s1234+lc*s12345];


%    g = [0;9.81];
%    U = m*g'*(pc1 + pc2 + pc3 + pc4 + pc5);  

taskValue = taskFunc(node(1),node(2),node(3),node(4),node(5));

U = taskValue(2)*sin(taskValue(1));
    
    
end