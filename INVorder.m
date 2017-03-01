function M = INVorder(T,v)
    if (size(T,1) ~= size(v,1)) 
        error('The matrix and the vector are not compatible!');
        M=0;
    else
        M = zeros(size(T));

        m = sum(v(:) == 0);
        r = zeros(size(v));
        j=1;
        k=m+1;
        for i=1:1:size(v,1)
            if (v(i) == 0)
                r(i) = j;
                j = j+1;
            else 
                r(i) = k;
                k = k + 1;
            end
        end

        for i=1:size(v,1)
           M(i,:) = T(r(i),:);
        end
        
       
    end
end

        

