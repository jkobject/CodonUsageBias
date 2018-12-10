%%%%%%%find the real entropy location in histogram (selection pressure: accumulative probability up to this location)
function PlocationF = entropy4Location(syno,NN,realEntropy)

switch syno
    case 2
        filename=['partition2o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(2,NN);
        pref=(log(pmax./pf))/NN;
        [NHIST,edges]=histcounts(pref,'BinWidth',0.0001);
        lhist=length(edges);
        if edges(1)<0
            NHISTnew(1)=sum(NHIST(1:2)); %%previously delete all the negative edges;but not good
            NHISTnew(2:(lhist-2))=NHIST(3:(lhist-1));
            edgesNew(1:(lhist-2))=edges(2:(lhist-1)); %% since -0 take up big portion and should consider as 0;-0 due to accuration loss like 4/3*3-1~=0;
            %%% equal to find the midpoint corresponding to NHIST
        else
            NHISTnew=NHIST;
            edgesNew(1:(lhist-1))=edges(1:(lhist-1));
        end
        PiOriginal=pmax./exp(edgesNew.*NN);
        NHISTfinal=NHISTnew.*PiOriginal;
        CumSumNHIST=cumsum(NHISTfinal);
        Plocationid=find(CumSumNHIST<1,1,'last');
        
        if realEntropy>edgesNew(end)
            PlocationF=CumSumNHIST(Plocationid);
        else
            [~,PlocationId,~]=find(realEntropy<=edgesNew,1);
            Plocation=CumSumNHIST(PlocationId);
            
            if PlocationId==1
                PlocationF=Plocation;
            else
                if Plocation<=0.95
                    PlocationF=Plocation;
                elseif Plocation>0.95 && Plocation<1
                    PlocationF=CumSumNHIST(PlocationId-1);
                else
                    PlocationF=CumSumNHIST(Plocationid);
                end
            end
        end
        
    case 3
        filename=['partition3o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(3,NN);
        pref=(log(pmax./pf))/NN;
        [NHIST,edges]=histcounts(pref,'BinWidth',0.0001);
        lhist=length(edges);
        if edges(1)<0
            NHISTnew(1)=sum(NHIST(1:2));   %%%%% add count of negative -0 to 0
            NHISTnew(2:(lhist-2))=NHIST(3:(lhist-1));
            edgesNew(1:(lhist-2))=edges(2:(lhist-1));
        else
            NHISTnew=NHIST;
            edgesNew(1:(lhist-1))=edges(1:(lhist-1));
        end
        PiOriginal=pmax./exp(edgesNew.*NN);
        NHISTfinal=NHISTnew.*PiOriginal; %% actually midpoint should be used, but time cost and accuracy
        CumSumNHIST=cumsum(NHISTfinal);
        Plocationid=find(CumSumNHIST<1,1,'last');
        
        if realEntropy>edgesNew(end)
            PlocationF=CumSumNHIST(Plocationid);
        else
            [~,PlocationId,~]=find(realEntropy<=edgesNew,1);
            Plocation=CumSumNHIST(PlocationId);
            
            if PlocationId==1
                PlocationF=Plocation;
            else
                if Plocation<=0.95
                    PlocationF=Plocation;
                elseif Plocation>0.95 && Plocation<1
                    PlocationF=CumSumNHIST(PlocationId-1);
                else
                    PlocationF=CumSumNHIST(Plocationid);
                end
            end
        end
        
        
        
    case 4
        filename=['partition4o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(4,NN);
        pref=(log(pmax./pf))/NN;
        [NHIST,edges]=histcounts(pref,'BinWidth',0.0001);
        
        lhist=length(edges);
        if edges(1)<0
            NHISTnew(1)=sum(NHIST(1:2));
            NHISTnew(2:(lhist-2))=NHIST(3:(lhist-1));
            edgesNew(1:(lhist-2))=edges(2:(lhist-1));
        else
            NHISTnew=NHIST;
            edgesNew(1:(lhist-1))=edges(1:(lhist-1));
        end
        
        PiOriginal=pmax./exp(edgesNew.*NN);
        NHISTfinal=NHISTnew.*PiOriginal;
        CumSumNHIST=cumsum(NHISTfinal);
        Plocationid=find(CumSumNHIST<1,1,'last');
        
        if realEntropy>edgesNew(end)
            PlocationF=CumSumNHIST(Plocationid);
        else
            [~,PlocationId,~]=find(realEntropy<=edgesNew,1);
            Plocation=CumSumNHIST(PlocationId);
            
            if PlocationId==1
                PlocationF=Plocation;
            else
                if Plocation<=0.95
                    PlocationF=Plocation;
                elseif Plocation>0.95 && Plocation<1
                    PlocationF=CumSumNHIST(PlocationId-1);
                else
                    PlocationF=CumSumNHIST(Plocationid);
                end
            end
        end
        
        
    case 6
        filename=['partition6o',num2str(NN),'p.csv'];
        pf=csvread(filename);
        pmax=Efor(6,NN);
        pref=(log(pmax./pf))/NN;
        if NN>100    %%%%sampling method, needn't adjust pi value
            [NHIST,edges]=histcounts(pref,'BinWidth',0.0001,'Normalization','cdf');
            lhist=length(edges);
            if edges(1)<0          %%%%%let edgesNew and NHIST match each other
                NHISTnew(1)=NHIST(2);
                NHISTnew(2:(lhist-2))=NHIST(3:(lhist-1));
                edgesNew(1:(lhist-2))=edges(2:(lhist-1));
            else
                NHISTnew=NHIST;
                edgesNew(1:(lhist-1))=edges(1:(lhist-1));
            end
            
            CumSumNHIST=NHISTnew;
            Plocationid=find(CumSumNHIST<1,1,'last');
            
            if realEntropy>edgesNew(end)
                PlocationF=CumSumNHIST(Plocationid);
            else
                [~,PlocationId,~]=find(realEntropy<=edgesNew,1);
                Plocation=CumSumNHIST(PlocationId);
                
                if PlocationId==1
                    PlocationF=Plocation;
                else
                    if Plocation<=0.95
                        PlocationF=Plocation;
                    elseif Plocation>0.95 && Plocation<1
                        PlocationF=CumSumNHIST(PlocationId-1);
                    else
                        PlocationF=CumSumNHIST(Plocationid);
                    end
                end
            end
            
            
            
        else
            [NHIST,edges]=histcounts(pref,'BinWidth',0.0001);
            lhist=length(edges);
            if edges(1)<0
                if lhist>2
                    NHISTnew(1)=sum(NHIST(1:2));
                    NHISTnew(2:(lhist-2))=NHIST(3:(lhist-1));
                    edgesNew(1:(lhist-2))=edges(2:(lhist-1));
                else
                    NHISTnew(1)=NHIST(1);
                    edgesNew(1)=edges(2);                    
                end
                
            else
                NHISTnew=NHIST;
                edgesNew(1:(lhist-1))=edges(1:(lhist-1)); %%%choose which edge point to match NHIST? My choice is here
            end
            PiOriginal=pmax./exp(edgesNew.*NN);   %%%%consider which edge contribute more
            NHISTfinal=NHISTnew.*PiOriginal;
            CumSumNHIST=cumsum(NHISTfinal);
            Plocationid=find(CumSumNHIST<1,1,'last');
            
            if realEntropy>edgesNew(end)
                PlocationF=CumSumNHIST(Plocationid);
            else
                [~,PlocationId,~]=find(realEntropy<=edgesNew,1);
                Plocation=CumSumNHIST(PlocationId);
                
                if PlocationId==1
                    PlocationF=Plocation;
                else
                    if Plocation<=0.95
                        PlocationF=Plocation;
                    elseif Plocation>0.95 && Plocation<1
                        PlocationF=CumSumNHIST(PlocationId-1);
                    else
                        PlocationF=CumSumNHIST(Plocationid);
                    end
                end
            end
            
            
        end
        
end

end


