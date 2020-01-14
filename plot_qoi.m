% HOW TO USE THIS SCRIPT:
% 1. put a break-point before Line 49 of sizeMeanStd.m
% 2. run main.m
% 3. once at the break-point, run the following code

subplot(2,1,1)
xx = (0:numel(sample)-1)*dx-40;
plot(xx,sample)
hold on
plot(xx,0*xx+cmin(i),'--','linewidth',2)
plot(xx,0*xx+cmax(i),'--','linewidth',2)
hold off
%legend({'example field','avg local min','avg local max'})
subplot(2,1,2)
histogram(S,7)
%legend('particle size distribution')