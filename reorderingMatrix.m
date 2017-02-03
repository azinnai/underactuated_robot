function M = reorderingMatrix(T,v)
    
    if (size(T,1) ~= size(v,1)) 
        error('The matrix and the vector are not compatible!');
        M=0;
    else
        if isa(T,'sym')
            M = sym('M',size(T));
        else
        M = zeros(size(T));
        end
        m = sum(v(:) == 0);

        j=1;
        k=m+1;
        for i=1:1:size(v,1)
            if (v(i) == 0)
                M(j,:) = T(i,:);
                j = j+1;
            else 
                M(k,:) = T(i,:);
                k = k + 1;
            end
        end
    end
end

        

