sz1 = size(region1Copy);
sz2 = size(region2Copy);
sz3 = size(region3Copy);
%-----------cluster1
count1 = 0;
for r=1:sz1(1,1)
    if region1Copy(r,1) ~=0
        count1 = count1+1;
        clus1(:,count1) = totalperinst(:,r);
        userid1(count1,:) = r;
    else
        count1 = count1;
    end;
    
    r=r+1;
end;
%--------------------- cluster2
count2 = 0;
for r=1:sz2(1,1)
    if region2Copy(r,1) ~=0
        count2 = count2+1;
        clus2(:,count2) = totalperinst(:,r);
        userid2(count2,:) = r;
    else
        count2 = count2;
    end;
    r=r+1;
end;
%-------------------------- cluster3
count3 = 0;
for r=1:sz3(1,1)
    if region3Copy(r,1) ~=0
        count3 = count3+1;
        clus3(:,count3) = totalperinst(:,r);
        userid3(count3,:) = r;
    else
        count3 = count3;
    end;
    r=r+1;
end;

