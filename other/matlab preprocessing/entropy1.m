function S =entropy1(B,P)
%% initialize
A = ['a','b','c','d'];   %% 4 codon types 
L = size(B,2);           %% B: aimed mRNA sequence matrix  
                          % L: length of aimed mRNA sequence
Q = zeros(1,L);          %% coresponding probability matrix for aimed mRNA sequece
S = 0;                   %% entropy of aimed mRNA 


%%main loop
for i=1:L
    for r=1:4
        if B(i) == A(r);
           Q(i) = P(r);   
            break;
        end
    end
end
    S = sum(Q .* log(Q));
end
