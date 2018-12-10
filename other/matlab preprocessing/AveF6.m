ave=zeros(1,300);
for l=101:400
    filename=['partition6o',num2str(l),'p.csv'];
    pf=dlmread(filename);
    pmax=Efor(6,l);
    pref=log(pf./pmax);
    ave(1,l-100)=sum(pref)/length(pref);
end
    AveEntropy6f(1:100)=AveEntropy6(1:100);
    AveEntropy6f(101:400)=ave(1:300);
    
    save('AveEntropy6f.mat',AveEntropy6f)