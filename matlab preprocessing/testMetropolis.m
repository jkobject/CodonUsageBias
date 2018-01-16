%% test Metropolis 

SpSize=10000000; %%SpSize: sample size

% [Pfin,OutX]=sampSequenceD(3,20,SpSize); %% call function of Metropolis algorithm to find configuration 
% [Pfin,OutX]=saAveEntropy6f201t300.matmpSequenceD(6,500,SpSize);

% A=[10,5,5];    %% find how many [5,5,10]?
% P=[1/3,1/3,1/3];
A=[100,100,100,100,50,50];    %% find how many [5,5,10]?
P=[1/6,1/6,1/6,1/6,1/6,1/6];



for i=1:SpSize
CompareX{i}=OutX500{i}-A;
B(i)=CompareX{i}(1,1)||CompareX{i}(1,2)||CompareX{i}(1,3)||CompareX{i}(1,4)||CompareX{i}(1,5)||CompareX{i}(1,6);   %% only if all the elements match [5,5,10] have 0 value in B
% B(i)=CompareX{i}(1,1)||CompareX{i}(1,2)||CompareX{i}(1,3);
end

CountId=find(~B);   %% find [5,5,10] indicesc
Count=length(CountId);

Proportion=Count/SpSize;  %% [5,5,10] propotion among all the samples

Pmpf=mnpdf(A,P);

Ratio=Proportion/Pmpf  %% compare with mnpdf value to test whether Metroplis work


%% when sampling 1000000, test AveEntropy6(56) and p11Sample6(101) since we know there is a cycle relationship between length and entropy values