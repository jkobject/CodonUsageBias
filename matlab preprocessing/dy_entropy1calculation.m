clc
clear

%% generate coresponding probability to a,b,c,d
P= rand(1,4)

%% generate mRNA sequence with length 200 codons
B = char('a' + floor(rand(1,200)./0.25))

%% calculate entropy for the whole sequence
entropy1(B,P)