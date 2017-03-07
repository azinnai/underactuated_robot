function qRand = randomConfig(jointLimits)



%a = pi - jointLimits;
%b = -a;

%qRand = a  + (b-a) * rand(5,1);

qRand = zeros(2,1);

a1 = 0;
b1 = 2*pi;
qRand(1) = a1  + (b1-a1) * rand(1);

a2 = 0;
b2 = 0.5;

qRand(2) = a2  + (b2-a2) * rand(1);



end