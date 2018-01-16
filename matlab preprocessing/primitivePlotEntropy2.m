clc
clear

%% generate mRNA sequence with length 200-2000 at 100 intervals
l = 200:100:2000;  % length matrix
nl = length(l);    
S=zeros(1,nl);     % entropy matrix for individual length

for i=1:nl
    B = char('a' + floor(rand(1,l(i))./0.25));
    S(i) = entropy2(B);  
    clear B
end

S
%% plot relationship between entropy and length
plot(l,S)

