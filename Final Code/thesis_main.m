%function []=thesis_main()
% INPUTS: 
clear;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data points (without the repeats)
allnodes=load('DLandComp_nodes_without_Duplicates.mat');
allnodes=allnodes.DLandComp_nodes_without_Duplicates;


%% calculate clusters (5 zones)
[node_clusters,nodes_in_cluster,cluster_len]=SpectralCluster(1,1);
% outputs 5 clusters along with nodes that are in each cluster

%% use intlinprog without comps on clusters 
cluster_order=myintlinprog(node_clusters,allnodes);
%%
[~,nodes_in_cluster_comp,cluster_len_comp]=SpectralCluster(cluster_order(2:6),allnodes);

%%
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C263');
node_clusters=[node_clusters,linspace(1,5,5)'];
neworder=[];
for j=1:5
    spot=find(cluster_order(2:6)==allnodes(node_clusters(j,1)));
    neworder=[neworder;j,spot]; %left column is node_cluster group, right column is cluster_order group
end
size=max(cluster_len)+max(cluster_len_comp);
nodes_comp_dl=zeros(size,5);
for k=1:5
    nodes_comp_dl(1:cluster_len(k)+cluster_len_comp(neworder(k,2)),k)=[nodes_in_cluster(1:cluster_len(k),k);nodes_in_cluster_comp(1:cluster_len_comp(neworder(k,2)),neworder(k,2))];  
end

figure
hold on;
datanodes=allnodes(nodes_comp_dl(1:cluster_len(1)+cluster_len_comp(neworder(1,2)),1));
scatter(Longitude(datanodes),Latitude(datanodes),'k*')
datanodes=allnodes(nodes_comp_dl(1:cluster_len(2)+cluster_len_comp(neworder(2,2)),2));
scatter(Longitude(datanodes),Latitude(datanodes),'b*')
datanodes=allnodes(nodes_comp_dl(1:cluster_len(3)+cluster_len_comp(neworder(3,2)),3));
scatter(Longitude(datanodes),Latitude(datanodes),'y*')
datanodes=allnodes(nodes_comp_dl(1:cluster_len(4)+cluster_len_comp(neworder(4,2)),4));
scatter(Longitude(datanodes),Latitude(datanodes),'g*')
datanodes=allnodes(nodes_comp_dl(1:cluster_len(5)+cluster_len_comp(neworder(5,2)),5));
scatter(Longitude(datanodes),Latitude(datanodes),'m*')
hold off

%%
% use intlinprog with comps on all nodes in each cluster
%for i=1:5
% x=load('xtest.mat');
% x=x.x;

    i=5;
    input_nodes=nodes_comp_dl(1:cluster_len(i)+cluster_len_comp(neworder(i,2)),i);%nodes_in_cluster(1:cluster_len(i),i);
    %input_nodes=input_nodes(17:22);
    x=myintlinprog_compactor(allnodes,cluster_len,input_nodes);
%end
%%
j=2;
i=3;
    input_nodes1=nodes_comp_dl(1:cluster_len(i)+cluster_len_comp(neworder(i,2)),i);%nodes_in_cluster(1:cluster_len(i),i);
    input_nodes2=nodes_comp_dl(1:cluster_len(i)+cluster_len_comp(neworder(i,2)),i);%nodes_in_cluster(1:cluster_len(i),i);
    input_nodes=[input_nodes1;input_nodes2];
    x=myintlinprog_compactor(allnodes,cluster_len,input_nodes);

%%
nDL=10;
inputnodes=linspace(1,nDL-1,nDL-1);
 allnodes=[randperm(238,nDL-3)+ones(1,nDL-3),249,259]';% [1,139,214,73,88,103,232,249];%% row vector, must always start with node 1
     x=myintlinprog_compactor(allnodes,cluster_len,inputnodes); 
%%
zone5=allnodes(nodes_comp_dl(1:cluster_len(5)+cluster_len_comp(neworder(5,2)),5));
AllZoneCoords(1:length(zone1),1:2)=[Longitude(zone1),Latitude(zone1)];




%end

