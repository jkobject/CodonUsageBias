function savePlot

h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['figure' num2str(length(h)+1-i)], 'pdf');
end
end