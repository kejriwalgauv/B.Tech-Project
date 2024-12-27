X = [peak maxtomin];
opts = statset('Display','final');
[idx,C] = kmeans(X,3,...
'Replicates',100,'Options',opts);
figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
'Location','NW')
title 'Max to Min ratio vs Peak power' ;
ylabel 'Max to Min power ratio';
xlabel 'Peak power consumption (kWhh)';
%colorbar
hold off
