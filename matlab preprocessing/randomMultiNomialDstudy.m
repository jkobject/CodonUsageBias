% function [midG,NgbinG]=randomMultiNomialDstudy(m,Nl,N)
Nl=4;    %% Nl:sequence length
N=1000;
m=4;         %% category number
 %% N&m jointly decide the sketch of the graph (for each length) so sample size matters
P=zeros(1,m);
P(1:m)=1/m;
% Y=zeros(1,1000);
Yp=zeros(1,N);
Xm=zeros(1,m);
Xm(1:m)=Nl/m;
Pm=zeros(1,m);
Pm(1:m)=1/m;
ym=mnpdf(Xm,Pm);
log(ym)
for i=1:N
Xr=randi(m,[1,Nl]);
for j=1:m
X=zeros(1,m);
X(j)=length(find(Xr==j));
end
Yp(i)=mnpdf(X,P);
end

% Y=log(Yp/ym);
Y=log(Yp);
% Y=Yp;
% figure
% h=histogram(Y,1000);


[Nbin,edges]=histcounts(Y,1000); 

k=1;
mid=zeros(1,(length(edges)-1)); 

for k=1:(length(edges)-1)

        mid(k)=(edges(k)+edges(k+1))/2;

end

Good = find(~(isinf(log(Nbin)))); 
midG=mid(Good);
NbinG=Nbin(Good);
NgbinG=log(NbinG);
% figure
plot(midG,NgbinG,'o'); 
xlabel('midpoint of bin');
ylabel('frequency');
title(['sample ', 'N=',num2str(N),' m=',num2str(m),' nl=',num2str(Nl)]);
    
hold on
% end