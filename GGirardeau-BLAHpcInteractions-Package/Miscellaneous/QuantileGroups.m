function groups = QuantileGroups(x,quant)
% Gives a vector with groups for quantiles. x=data, quant=quantiles obtained using matlab quantile or prctile

groups=zeros(length(x),1);
groups(x<=quant(1))=1;
for i=2:length(quant)
  groups(x>quant(i-1)&x<=quant(i))=i;
end
groups(x>quant(end))=length(quant)+1;