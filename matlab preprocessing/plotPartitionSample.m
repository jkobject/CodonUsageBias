function [] = plotPartitionSample (pin)
[NHIST,edges]=histcounts(pin,'BinWidth',0.00000001);
EdgesNew=edges(1:end-1);
NHISTnew=NHIST.*EdgesNew;
figure
% plot(EdgesNew,NHISTnew,'.');
histogram(NHISTnew);
end