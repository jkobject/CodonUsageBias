function Sr = replaceT(N,ct,cf)

r = rand([1,N]);

switch(length(ct))
    
    case 2
        for j=1:N
            if r(j)<cf(1)
                SR(j)=ct(1);
            else SR(j)=ct(2);
            end
        end
        
    case 3
        for j=1:N
            if r(j)<cf(1)
                SR(j)=ct(1);
                
            elseif r(j)>=cf(1) && r(j)<(cf(1)+cf(2))
                SR(j)=ct(2);
            else SR(j)=ct(3);
            end
        end
        
    case 4
        for j=1:N
            if r(j)<cf(1)
                SR(j)=ct(1);
                
            elseif r(j)>=cf(1) && r(j)<(cf(1)+cf(2))
                SR(j)=ct(2);
                
            elseif r(j)>=(cf(1)+cf(2)) && r(j)<(cf(1)+cf(2)+cf(3))
                SR(j)=ct(3);
            else SR(j)=ct(4);
            end
        end
        
    case 6 
        for j=1:N
            if r(j)<cf(1)
                SR(j)=ct(1);
                
            elseif r(j)>=cf(1) && r(j)<(cf(1)+cf(2))
                SR(j)=ct(2);
                
            elseif r(j)>=(cf(1)+cf(2)) && r(j)<(cf(1)+cf(2)+cf(3))
                SR(j)=ct(3);
            
            elseif r(j)>=(cf(1)+cf(2)+cf(3)) && r(j)<(cf(1)+cf(2)+cf(3)+cf(4))
                SR(j)=ct(4);
            
            elseif r(j)>=(j)<(cf(1)+cf(2)+cf(3)+cf(4)) && r(j)<(cf(1)+cf(2)+cf(3)+cf(4)+cf(5))
                SR(j)=ct(5);
                
            else SR(j)=ct(6);
            end
        end
end

Sr=cellstr(SR');

end


