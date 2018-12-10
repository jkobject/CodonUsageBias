function s = entropy2(B)
%% initialize

L = size(B,2);           %% B: aimed mRNA sequence matrix  
                          % L: length of aimed mRNA sequence
Q = zeros(1,L);          %% coresponding probability matrix for aimed mRNA sequece

s = 0;                   %% entropy of aimed mRNA 



% Distribution Pa Pb Pc Pd 
% Pa=binornd(counta,0.1,[1,counta]);
% Pb=rand(1,countb);
% Pc=gamrnd(5,10,[1 countb]);
% Pd=chi2rnd(6,[1,countd]);

%%main loop

a = find (B =='a');
La = length (a);
Q(a) = binornd(La,0.1,[1,La]);

b = find (B =='b');
Lb = length (b);
Q(b) = rand(1,Lb);

c = find (B =='c');  
Lc = length (c);
Q(c) = gamrnd(5,10,[1,Lc]);

d = find (B=='d');
Ld = length (d);
Q(d) = chi2rnd(6,[1,Ld]);

s = sum(Q .* log(Q));
    
end
