function difference = boxMinus(a,b)
    
if (size(a,1)~=size(b,1))
    disp('ERROR IN BOXMINUS, DIFFERENT SIZES ')
    size(a)
    size(b)
    difference = -1;
else
    difference = zeros(size(a));
    for i = 1 : size(a,1)
        currentAngleMat = [cos(a(i)), - sin(a(i)); sin(a(i)), cos(a(i))];
        referenceAngleMat = [cos(b(i)), - sin(b(i)); sin(b(i)), cos(b(i))];
        errorAngleMat = (currentAngleMat)\referenceAngleMat;
        difference(i) = atan2(errorAngleMat(2,1),errorAngleMat(1,1));
    end
    
end
    
    
end