NNsp=[NN,Y];
NNs=sortrows(NNsp);
lenNN=length(NNs(:,1));
NNs1=NNs(:,1);
NNs2=NNs(:,2);
NNnew=NNs1(0.2*lenNN:0.8*lenNN);
Ynew=NNs2(0.2*lenNN:0.8*lenNN);

figure
histogram(NNnew);
Ynewl=log(Ynew);
figure
h=histogram(Ynewl);

[Nbin,edges]=histcounts(Ynewl); %% prepare to plot midpoint,slope..

k=1;
mid=zeros(1,(length(edges)-1)); %%get (mid) midpoints vector
for k=1:(length(edges)-1)
%      if Nbin(k)==0    %%discard outliers
%         mid(k)=NaN;
%         Nbin=NaN;
%     else
        mid(k)=(edges(k)+edges(k+1))/2;
%     end
end

Good = find(~(isinf(log(Nbin)))); %%discard NaN values for further plot
midG=mid(Good);
NbinG=Nbin(Good); %% NbinG is the histogram value excluding NaN
NgbinG=log(NbinG);
figure
plot(midG,NgbinG,'o'); 
xlabel('midpoint of bin');
ylabel('log of frequency');
title('Arg, Sah6,(Pr/Pm) probability ratio with certain lengths deleted');