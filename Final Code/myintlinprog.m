function [cluster_order]=myintlinprog(node_clusters,allnodes)
%clear;clc;
m=1; %number of salemen/vehicles
DATA=load('Final Distance Matrix.mat');
DATA=DATA.DATA;
DATACopy=DATA;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C263');
compsData=load('DistComps - Compactors x DLs.mat');
compsData=compsData.DistComps;
%[LongComps, LatComps] = readvars('Landfill_Compactors.xlsx','Sheet','Sheet1','Range','X67:Y86');

DATA=[DATA,[10000*ones(1,length(compsData(:,1)));compsData'];[10000*ones(length(compsData(:,1)),1),compsData],10000*ones(length(compsData(:,1)))-10000*eye(length(compsData(:,1)))];

%datanodes=linspace(1,239,239);
%LINE BELOW IS NEW
 %datanodes=[1,randperm(238,5)+ones(1,5)]; % row vector, must always start with node 1
% %DATA=DATA([1,5,76,108,223],[1,5,76,108,223]);
datanodesALL=load('DLandComp_nodes_without_Duplicates.mat');
datanodesALL=datanodesALL.DLandComp_nodes_without_Duplicates;


%datanodes=linspace(1,259,259);
datanodes=[1;allnodes(node_clusters)];%datanodes(1:134);
Long=Longitude(datanodes);
Lat=Latitude(datanodes);
%x=load('x2.mat');
%x=x.x2;

nStops = length(datanodes); % you can use any number, but the problem size scales as N^2
nDLs=nStops-1;


%%% Calculate Distances Between Points
idxs1 = nchoosek(1:259,2);
idxs1=[idxs1;idxs1(:,2),idxs1(:,1)]; 
idxs=[];
id=[];
for p=1:nStops
    for w=1:nStops
if w~=p
        id=find(idxs1(:,1)==datanodes(p) & idxs1(:,2)==datanodes(w));
idxs=[idxs;idxs1(id,:)];
%id=find(idxs1(:,2)==datanodes(p) & idxs1(:,1)==datanodes(w)); 
%idxs=[idxs;idxs1(id,:)];
end
    end
end


nCombs=length(idxs);

dist=zeros(nCombs,1);

for i=1:nCombs
    dist(i)=DATA(idxs(i,1),idxs(i,2));
end

%%
f=[dist;zeros(nDLs,1)]; % concatenate f with zeros at the end to account for u at the end of x vector (length #nodes)
%Aeq=spalloc(2*nStops,nCombs+nDLs,125000);
Aeq=zeros(2*nStops,nCombs+nDLs);
beq=zeros(2*nStops,1);


%% Constraint 2
Aeq(1,find(idxs(:,1)==1))=ones(1,nDLs);
beq(1)=m;
%% Constraint 3
Aeq(2,find(idxs(:,2)==1))=ones(1,nDLs);
beq(2)=m;
%% Constraint 4

for j=2:nStops
    rowidxs=find(idxs(:,2)==datanodes(j));
    Aeq(j+1,rowidxs)=ones(1,nDLs);
    beq(j+1)=1;
end

%% Constraint 5
cnt=nStops+2;
for i=2:nStops
    rowidxs=find(idxs(:,1)==datanodes(i));
    Aeq(cnt,rowidxs)=ones(1,nDLs);
    beq(cnt)=1;
    cnt=cnt+1;
end

%% Constraint 6, 7
%A=spalloc(2*nDLs+nDLs^2,nCombs+nDLs,252000);
A=zeros(2*nDLs+nDLs^2,nCombs+nDLs);
b=zeros(2*nDLs+nDLs^2,1);
L=7;
K=3;
cnter=nStops;
for i=2:nStops
    x1i=find(idxs(:,1)==1 & idxs(:,2)==datanodes(i));
    xi1=find(idxs(:,2)==1 & idxs(:,1)==datanodes(i));
    ui=nCombs+i-1;
    A(i-1,[x1i,xi1,ui])=[L-2,-1,1];
    b(i-1)=L-1;
    A(cnter,[x1i,xi1,ui])=[-1,K-2,-1];
    b(cnter)=-2;
    cnter=cnter+1;
end

%% Constraint 9
counter=2*nStops-1;
for i=2:nStops
    for j=2:nStops
        if i==j
            % do nothing
        else
            xij=find(idxs(:,1)==datanodes(i) & idxs(:,2)==datanodes(j));
            xji=find(idxs(:,2)==datanodes(i) & idxs(:,1)==datanodes(j));
            ui=nCombs+i-1;
            uj=nCombs+j-1;
            A(counter,[xij,xji,ui,uj])=[L,L-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(counter)=L-1;
            counter=counter+1;
        end
    end
end
%% Constraint 8
    xstart=find(idxs(:,1)==1);
    xend=find(idxs(:,2)==1);
    for i=1:nDLs

    A(counter,[xstart(i),xend(i)])=[1 1];
    b(counter)=1;
    counter=counter+1;
end


%%
int=ones(length(f),1);
intcon=find(int);

lb=[zeros(nCombs,1);ones(nDLs,1)];
ub=ones(nCombs,1);
% %options = optimoptions('intlinprog','MaxTime',172800);%,'BranchRule','maxpscost','MaxNodes',10000);%,'MaxFeasiblePoints',10000,'Heuristics','diving','NodeSelection','minobj');
x=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

%%
segments = find(x(1:nCombs)<1.05 & x(1:nCombs)>.95); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
truetripscopy=truetrips;
xplot=zeros(1,2);
yplot=zeros(1,2);
%%
% figure;hold on;
% plot(Longitude,Latitude,'r*')
% plot(Long,Lat,'g*')
% for i=1:length(truetripscopy)
%     xplot=[Latitude(truetripscopy(i,1)),Latitude(truetripscopy(i,2))];
%     yplot=[Longitude(truetripscopy(i,1)),Longitude(truetripscopy(i,2))];
%     plot(yplot,xplot,'b')
% end
% 
% hold off




%%
figure;hold on;
%plot(Longitude,Latitude,'r*')

% REMOVE NEXT 2 LINES TO RUN FULL SCRIPT
%truetrips=truetripscopy;
%L=45;
plot(Long,Lat,'r*')
route_long_lat=zeros(L+2,4*m);
for h=1:m
q=1;
yplot(2)=0;
while yplot(2)~=Longitude(1)
    if q==1
        i=find(truetrips(:,1)==1);
        i=i(1);
    else
        i=find(truetrips(:,1)==nextDL);
        i=i(1);
    end
    
    xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
    route_long_lat(q,4*h-2)=xplot(1);
    yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
    route_long_lat(q,4*h-3)=yplot(1);
    route_long_lat(q,[4*h-1,4*h])=truetrips(i,:);
    if h==1
        plot(yplot,xplot,'b')    
    elseif h==2
        plot(yplot,xplot,'c')
    elseif h==3
        plot(yplot,xplot,'g')
    elseif h==4
        plot(yplot,xplot,'m')  
    elseif h==5
        plot(yplot,xplot,'k')  
    end
    nextDL=truetrips(i,2);
    truetrips(i,:)=[];
    q=q+1;
end
route_long_lat(q,4*h-2)=xplot(2);
route_long_lat(q,4*h-3)=yplot(2);
%route_long_lat(q,[4*h-1,4*h])=[nextDL,1];
end

hold off

% calculate order of stopping at clusters in vector output
show=find(route_long_lat(:,3));
cluster_order= [route_long_lat(show,3);1] ;
end
%%
%vecg=[13,234,19,210,15,11,235];
% xcopy=x;
% truetrips=truetripscopy;
% total=0;
% vecg=route_long_lat(1:21,11:12);%[14,15,16,17,18];
% leng=length(vecg);
% for g=1:leng
%     call1=find(idxs(:,1)==vecg(g,1) & idxs(:,2)==vecg(g,2));
%     call2=find(idxs(:,1)==vecg(g,2) & idxs(:,2)==vecg(g,1));
%     if x(call1)~=0
%         x(call1)=0;
%         total=total+1;
%         nodenum=find(datanodes==vecg(g,1));
%         x(nCombs+nodenum-1)=0;
%         
%     end
%     if x(call2)~=0
%         x(call2)=0;
%         nodenum=find(datanodes==vecg(g,1));
%         x(nCombs+nodenum-1)=0;
%         
%         total=total+1;
%     end
% end
% vecb=route_long_lat(1:40,3:4);%[14,15,16,17,18];
% lenb=length(vecb);
% for g=1:lenb
%     call1=find(idxs(:,1)==vecb(g,1) & idxs(:,2)==vecb(g,2));
%     call2=find(idxs(:,1)==vecb(g,2) & idxs(:,2)==vecb(g,1));
%     if x(call1)~=0
%         x(call1)=0;
%         total=total+1;
%         nodenum=find(datanodes==vecb(g,1));
%         x(nCombs+nodenum-1)=0;
%         
%     end
%     if x(call2)~=0
%         x(call2)=0;
%         total=total+1;
%         nodenum=find(datanodes==vecb(g,1));
%         x(nCombs+nodenum-1)=0;
%        
%     end
% end
