function xx = histY2(n,q)  %% n: length; q: breaking points to assign probability value

for i=1:length(n)
    
l(:,1:100000)=n(i); 

xx(:,i) = arrayfun(@(x) entropyD(x,q),l); %% for a certain length draw its histogram


end
end

