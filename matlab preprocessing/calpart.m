function mnvect = calpart(sum,level) %% i is to keep the bit for mnvect,j is the original level amount

bitmin=ceil(sum/level);

if level>1  %% if it's not the last bit, preserve the first bit in mnvect

    for n=(sum-1):-1:bitmin
        
        temp(1)=sum-1;%% the first bit is between [],exclude same combinations
 
        temp(2)=sum-bit1; %% calculate the second bit
        

    %     minB=min(mnvect(1:i));

        if (temp(2)<=1 && length(find(diff(mnvect(1:i))<=0))==(i-1))
            if i<j
            temp(i+1:j)=0;
            end
        end
        mnvect=temp;
        disp(mnvect);%%try to use j to adjust i as the index of mnvect   
        

        if (bit2>1 && length(find(diff(mnvect(1:i))<=0))==(i-1))
            if i<j
            temp(i+1:j)=0;
            end
        mnvect=temp;
        disp(mnvect);
            
        temp=calpart(sum-,level-1);
        if i<j
        i=i+1;
        end
        mnvect=[mnvect;[temp,repmat()];
        end
        i=1;
    end
end
end
   







    


    