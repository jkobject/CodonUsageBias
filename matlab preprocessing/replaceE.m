function Sr = replaceE(N,ct)

r = randi([1,length(ct)],1,N);

switch(length(ct))
    case 2
        j=1;
        for j=1:N
            if r(j)==1
                SR(j)=ct(1);
                
            else
                SR(j)=ct(2);
            end
        end
        
    case 3
        j=1;
        for j=1:N
            if r(j)==1
                SR(j)=ct(1);
                
            elseif r(j)==2
                SR(j)=ct(2);
                
            else SR(j)=ct(3);
            end
        end
        
    case 4
        j=1;
        for j=1:N
            if r(j)==1
                SR(j)=ct(1);
                
            elseif r(j)==2
                SR(j)=ct(2);
                
            elseif r(j)==3
                SR(j)=ct(3);
            else SR(j)=ct(4);
            end
        end
        
    case 6
        j=1;
        for j=1:N
            if r(j)==1
                SR(j)=ct(1);
                
            elseif r(j)==2
                SR(j)=ct(2);
                
            elseif r(j)==3
                SR(j)=ct(3);
                
            elseif r(j)==4
                SR(j)=ct(4);
                
            elseif r(j)==5
                SR(j)=ct(5);
                
            else SR(j)=ct(6);
            end
        end
end

Sr=cellstr(SR');
end

