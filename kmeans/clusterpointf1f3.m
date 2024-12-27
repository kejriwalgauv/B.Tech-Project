for i=1:186
    if idx(i,1) ==1
        region1(i,:) = X(i,:);
    elseif idx(i,1) ==2
        region2(i,:) = X(i,:);
    else
        region3(i,:) = X(i,:);
    end;
    i=i+1;
end;
sz1 = size(region1);
sz2 = size(region2);
sz3 = size(region3);
%-----------cluster1
a1=0;
a2=0;
count1 = 0;
for r=1:sz1(1,1)
    a1 = a1 + region1(r,1);
    a2 = a2 + region1(r,2);
    if region1(r,1) ~=0
        count1 = count1+1;
        user1(count1,:) = region1(r,:);
    else
        count1 = count1;
    end;
    
    r=r+1;
end;
region1avg(1,1) = a1/count1;
region1avg(1,2) = a2/count1;
%--------------------- cluster2
a3=0;
a4=0;
count2 = 0;
for r=1:sz2(1,1)
    a3 = a3 + region2(r,1);
    a4 = a4 + region2(r,2);
    if region2(r,1) ~=0
        count2 = count2+1;
        user2(count2,:) = region2(r,:);
    else
        count2 = count2;
    end;
    r=r+1;
end;
region2avg(1,1) = a3/count2;
region2avg(1,2) = a4/count2;
%-------------------------- cluster3
a5=0;
a6=0;
count3 = 0;
for r=1:sz3(1,1)
    a5 = a5 + region3(r,1);
    a6 = a6 + region3(r,2);
    if region3(r,1) ~=0
        count3 = count3+1;
        user3(count3,:) = region3(r,:);
    else
        count3 = count3;
    end;
    r=r+1;
end;
region3avg(1,1) = a5/count3;
region3avg(1,2) = a6/count3;

%----------centroids of clusters
C1 = C(1,:);
C2 = C(2,:);
C3 = C(3,:);
%----------distance between centroids
S1 = [C1 ; C2];
d1 = pdist(S1 ,'euclidean');
S2 = [C2 ; C3];
d2 = pdist(S2 ,'euclidean');
S3 = [C1 ; C3];
d3 = pdist(S3 ,'euclidean');

%-------spread of points within clusters
s1 = std(user1);
s2 = std(user2);
s3 = std(user3);

%-------------- AREA OF CLUSTERS
A1 =pi *s1(1,1)*s1(1,2);
A2 =pi *s2(1,1)*s2(1,2);
A3 =pi *s3(1,1)*s3(1,2);
%----------PI
p1 = d1/(A1+A2);
p2= d2/(A2+A3);
p3 = d3/(A3+A1);

P = p1+p2+p3; 

