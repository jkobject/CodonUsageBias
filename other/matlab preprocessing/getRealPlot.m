function [NHISTreal,EdgesReal]=getRealPlot(X)
Ddone1=X(~isnan(X));
Ddone2=Ddone1(Ddone1>0);
[NHISTrealp,EdgesRealp]=histcounts(Ddone2,'BinWidth',0.001,'Normalization','Probability');
EdgesReal=EdgesRealp(2:end);
NHISTreal=NHISTrealp./sum(NHISTrealp);
end