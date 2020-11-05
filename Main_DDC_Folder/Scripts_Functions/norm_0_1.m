function [Array2]=norm_0_1(Array)


max1=max(Array);
min1=min(Array);

Array2=[];
for i=1:length(Array)
    Array2(i)=(Array(i)-min1)/(max1-min1);
end