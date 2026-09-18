function [] = rowdata_columnumber_plot(inputmat)
% input: [—————
%         ——x——
%         —————]    num*dim

hold on
rowsize=size(inputmat,1);
colsize=size(inputmat,2);

for i=1:rowsize
   %scatter(1:colsize,inputmat(i,:))
end

% errorbar
standard_deviation=[];
for i=1:colsize
    standard_deviation=[standard_deviation,std(inputmat(:,i))];
end
%sem= std/sqrt(n) using sem!
errorbar(mean(inputmat,1),standard_deviation./sqrt(rowsize))

hold off
end