function pref = getPref(NN,syno)
switch syno
    case 2
        filename=['partition2o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(2,NN);
        pref=(log(pmax./pf))/NN;
    case 3 
        filename=['partition3o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(3,NN);
        pref=(log(pmax./pf))/NN;
    case 4
        filename=['partition4o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(4,NN);
        pref=(log(pmax./pf))/NN;
    case 6
        filename=['partition6o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(6,NN);
        pref=(log(pmax./pf))/NN;
end
        
end