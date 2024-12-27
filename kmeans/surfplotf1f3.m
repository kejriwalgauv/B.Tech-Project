for i = 1:count1
    user1(i) = i;
    i=i+1;
end;
for i = 1:count2
    user2(i) = i;
    i=i+1;
end;
for i = 1:count3
    user3(i) = i;
    i=i+1;
end;
%------------------surface plot
mymap = [0 0 0.7
         0 0.7 0.7
         0 0.5 0
         1 1 0       
         1 0 0];
%surf(time,user3,clus3','FaceAlpha',1,'EdgeColor','k','FaceColor','interp','FaceLighting','gouraud')
%surf(time,user2,clus2','FaceAlpha',1,'EdgeColor','k','FaceColor','interp','FaceLighting','gouraud')
surf(time,user1,clus1','FaceAlpha',1,'EdgeColor','k','FaceColor','interp','FaceLighting','gouraud')
colormap(mymap)
colorbar;

