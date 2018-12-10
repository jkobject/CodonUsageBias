clc
clear

setSynonymousCodonTable;

pasteCodon=getCodonSequence; %% call function getCodonSequence to creat input

NumSynonymous={'4','2'};

NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)'};

NameAminoAcid={'Gly','Asp'};


for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
    G(i,:)=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA 
    D(i,:)=AspAminoAcidH(pasteCodon{1,i}); %%get Ydd for each mRNA
   
end


X={G,D};

for j=1:length(X)
% Xpr=zeros(length(X{j}),1);
% Xpr=X{j};    
% Xok=Xpr(~(isnan(Xpr))); %% for mean calculation matrix
% Xst=sort(Xok);
% xlen=abs(Xst(length(Xst)-Xst(1)));
% mean=sum(Xok)/xlen;  

figure 

h=histogram(X{j},'BinWidth',0.04); %% normalize bin size

xlim([0,1]);

xlabel([NameAminoAcid(j),'ratio of probability']);

ylabel('frequency');

title([NameSynonymous(j), 'Histogram in all mRNAs']);


[Nbin,edges]=histcounts(X{j}); %% prepare to plot midpoint,slope.


k=1;
mid=zeros(1,(length(edges)-1)); %%get (mid) midpoints vector
for k=1:(length(edges)-1)
    if log(Nbin(k))<=1;    %%discard outliers
        mid(k)=NaN;
        Nbin(k)=NaN;
    else
        mid(k)=(edges(k)+edges(k+1))/2;
    end
end


% pointPosition=[mid',Nbin'];
% 
% [hn,cn]=size(pointPosition);

% p=1;
% for p=1:hn
%     pointLabel=['(',num2str(mid(p)),',',num2str(Nbin(p)),')']; %%label points
%     text(mid(p),Nbin(p)+0.1,pointLabel);
% end


figure

plot(mid,log(Nbin),'o');

hold on;

Good = ~(isnan(mid) | isnan(log(Nbin))) ; %%discard NaN values for further plot
Ngbin=log(Nbin);
midG=mid(Good);
NgbinG=Ngbin(Good);

mean=sum(midG.*NgbinG)./sum(NgbinG);

slopeOneC = polyfit(midG,NgbinG,1); %%coefficiencies of fitted line

filin=polyval(slopeOneC,midG);

% yresid=filin-NgbinG;     %% calculate fitting error
% SSresid=sum(yresid.^2);
% SStotal=(length(NgbinG)-1)*var(NgbinG);
% rsq=1-SSresid/SStotal;

slopeOne=slopeOneC(1);

%plot(midG,filin,'g-',midG,filin+2*delta,'r:',midG,filin-2*delta,'r:');

plot(midG,filin);

legend(['slope:',num2str(slopeOne)]);

xlabel('midpoint of bin');

ylabel('log of frequency');

title([NameAminoAcid(j),'bin plot']);

mdl = fitlm((midG)',(NgbinG)');

figure
plotResiduals(mdl,'probability');
title(NameAminoAcid(j));
% 
% B=log(Nbin);
% B(isnan(B))=0; %%delete NaN and calculate mean
% Bid=find(B); %% delete 0
% % b=length(Bid)*(h.BinWidth);
% mean=sum(B(Bid))/length(Bid);

meanAll(1,j)=mean;
RsqAll(1,j)=mdl.Rsquared.Ordinary;
SlopeAll(1,j)=slopeOne;

end

tbl = table((NumSynonymous)',(SlopeAll)',(meanAll)',(RsqAll)','VariableNames',{'NumSynonymous' 'Slope' 'mean' 'Rsq'},'RowNames',(NameAminoAcid)');
    
Tbl = sortrows(tbl,{'NumSynonymous'},'ascend')
