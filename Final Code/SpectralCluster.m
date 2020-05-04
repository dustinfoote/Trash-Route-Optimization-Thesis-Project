function [node_clusters,nodes_in_cluster,cluster_len]=SpectralCluster(c,allnodes,input_nodes)
% INPUT: c is 1 if clustering DLs. c is a 5x1 vector of closest DL (tracked by node # to
% the cluster means

%DATA=load('Final Distance Matrix.mat');
%DATA=DATA.DATA;
%DATACopy=DATA;
%[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C263');
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet3','Range','P2:Q154');

%compsData=load('DistComps - Compactors x DLs.mat');
%compsData=compsData.DistComps;
%[LongComps, LatComps] = readvars('Landfill_Compactors.xlsx','Sheet','Sheet1','Range','X67:Y86');

%DATA=[DATA,[10000*ones(1,length(compsData(:,1)));compsData'];[10000*ones(length(compsData(:,1)),1),compsData],10000*ones(length(compsData(:,1)))-10000*eye(length(compsData(:,1)))];

if length(c)==1
    X=[Longitude(1:133), Latitude(1:133),linspace(2,134,133)'] ;%[Longitude, Latitude,linspace(2,154,153)'] ;%load('Project Data.txt');
else
    X=[Longitude(134:153), Latitude(134:153),linspace(135,154,20)'] ;
end
%Comp_data=load('Project Data Compacters.txt');
%[a b]=size(BinsandComps_loc);
%divider=0.15;

% Just for 3 cluster version
X=X(input_nodes,:);



% BinsandComps_loc(:,1)=(BinsandComps_loc(:,1)+111.935).*100;%Clean up the data
% BinsandComps_loc(:,2)=(BinsandComps_loc(:,2)-33.418).*1000;

%X=[BinsandComps_loc(:,2)./5,BinsandComps_loc(:,1).*1.75];
N=length(X);

%N=238;
%X=[Bins_data(:,2),7.8.*Bins_data(:,1)];




%one=ones(5,2);
%r=12.*(rand(5,2)-0.5.*ones(5,2));
if length(c)==1
    %r=2.*(rand(5,2)-0.5.*ones(5,2));%randperm(153,5);%
    r=2.*(rand(3,2)-0.5.*ones(3,2));
else
    %rloc=[];
    r=[];
    for i=1:5
        rloc=find(allnodes==c(i));
        %rloc(i)=rloc(i)-1;
        r(i,1)=Longitude(rloc-1);
        r(i,2)=Latitude(rloc-1);
    end
end
%r=X(ran,:);
figure
scatter(X(:,1),X(:,2),'b*');
hold on
scatter(r(:,1),r(:,2),'r*')
hold off
len=length(r);
errx=1;erry=1;
%g=1;
if length(c)==1
while errx>.005 && erry>.005
%g=g+1;

eucdist=[];
%dist=[];
for i=1:len
    eucdist(:,i)=((r(i,1).*ones(N,1)-X(:,1)).^2+(r(i,2).*ones(N,1)-X(:,2)).^2).^.5;
%     for j=2:134
%         if j==ran(j)
%             % do nothing
%         else
%             dist(j,i)=DATA(allnodes(j),ran(i));
%         end
%     end
end
[~,B]= min(eucdist');
B=B';
i1=0;
i2=0;
i3=0;
i4=0;
i5=0;
Group1=[];
Group2=[];
Group3=[];
Group4=[];
Group5=[];
for k=1:N
    if B(k)==1
        i1=i1+1;
        Group1(i1,:)=X(k,:);
    elseif B(k)==2
        i2=i2+1;
        Group2(i2,:)=X(k,:);
    elseif B(k)==3
        i3=i3+1;
        Group3(i3,:)=X(k,:);
    elseif B(k)==4
        i4=i4+1;
        Group4(i4,:)=X(k,:);
    elseif B(k)==5
        i5=i5+1;
        Group5(i5,:)=X(k,:);        
    end
end





len1=length(Group1);
len2=length(Group2);
len3=length(Group3);
len4=length(Group4);
len5=length(Group5);

rold=r;

r(1,:)=[sum(Group1(:,1)),sum(Group1(:,2))]./len1;
r(2,:)=[sum(Group2(:,1)),sum(Group2(:,2))]./len2;
r(3,:)=[sum(Group3(:,1)),sum(Group3(:,2))]./len3;
%r(4,:)=[sum(Group4(:,1)),sum(Group4(:,2))]./len4;
%r(5,:)=[sum(Group5(:,1)),sum(Group5(:,2))]./len5;

errx=norm(abs(r(:,1)-rold(:,1)));
erry=norm(abs(r(:,2)-rold(:,2)));
end
figure
scatter(Group1(:,1),Group1(:,2),'k*')
hold on
scatter(Group2(:,1),Group2(:,2),'b*')
scatter(Group3(:,1),Group3(:,2),'y*')
%scatter(Group4(:,1),Group4(:,2),'g*')
%scatter(Group5(:,1),Group5(:,2),'m*')
hold off
else
   % g=g+1;

eucdist=[];
%dist=[];
for i=1:len
    eucdist(:,i)=((r(i,1).*ones(N,1)-X(:,1)).^2+(r(i,2).*ones(N,1)-X(:,2)).^2).^.5;
%     for j=2:134
%         if j==ran(j)
%             % do nothing
%         else
%             dist(j,i)=DATA(allnodes(j),ran(i));
%         end
%     end
end
[~,B]= min(eucdist');
B=B';
i1=0;
i2=0;
i3=0;
i4=0;
i5=0;
Group1=[];
Group2=[];
Group3=[];
Group4=[];
Group5=[];
for k=1:N
    if B(k)==1
        i1=i1+1;
        Group1(i1,:)=X(k,:);
    elseif B(k)==2
        i2=i2+1;
        Group2(i2,:)=X(k,:);
    elseif B(k)==3
        i3=i3+1;
        Group3(i3,:)=X(k,:);
    elseif B(k)==4
        i4=i4+1;
        Group4(i4,:)=X(k,:);
    elseif B(k)==5
        i5=i5+1;
        Group5(i5,:)=X(k,:);        
    end
end



figure
scatter(Group1(:,1),Group1(:,2),'k*')
hold on
scatter(Group2(:,1),Group2(:,2),'b*')
scatter(Group3(:,1),Group3(:,2),'y*')
%scatter(Group4(:,1),Group4(:,2),'g*')
%scatter(Group5(:,1),Group5(:,2),'m*')
hold off

len1=length(Group1(:,1));
len2=length(Group2(:,1));
len3=length(Group3(:,1));
len4=length(Group4(:,1));
len5=length(Group5(:,1));

end
    
    %%
cluster_len=[len1,len2,len3]%,len4,len5];
if length(c)==1
[~,B2]=min(eucdist);
node_clusters=[];
for k=1:3%5
    node_clusters(k,:)=X(B2(k),:);
end
end
nodes_in_cluster=zeros(max([len1,len2,len3]),3);%,len4,len5]),5);
nodes_in_cluster(1:len1,1)=Group1(:,3);
nodes_in_cluster(1:len2,2)=Group2(:,3);
nodes_in_cluster(1:len3,3)=Group3(:,3);
%nodes_in_cluster(1:len4,4)=Group4(:,3);
%nodes_in_cluster(1:len5,5)=Group5(:,3);

figure;
for i=1:4%6
    if i==1
        scatter(Group1(:,1),Group1(:,2),'k*');
        hold on;
    elseif i==2
        scatter(Group2(:,1),Group2(:,2),'b*');
    elseif i==3
        scatter(Group3(:,1),Group3(:,2),'y*');
%     elseif i==4
%         scatter(Group4(:,1),Group4(:,2),'g*');
%     elseif i==5
%         scatter(Group5(:,1),Group5(:,2),'m*');
    else
        if length(c)==1
        scatter(node_clusters(:,1),node_clusters(:,2),'r*');
        node_clusters=node_clusters(:,3);
        else
            node_clusters=1;
        end
        hold off;
    end
end    
    

end


% k=5;
% idx=spectralcluster(X,k);                %(X,48,'SimilarityGraph','knn','Radius',5)
% figure
% gscatter(X(:,1),X(:,2),idx)

