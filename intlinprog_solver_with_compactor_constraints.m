clear;clc;
m=4; %number of salemen/vehicles
H=3; %number of DL stops b4 visiting compactor
L=15;
K=2;
fullscale=0; % if problem is full scale, put 1, otherwise 0

DATA=load('Final Distance Matrix.mat');
DATA=DATA.DATA;
%DATACopy=DATA;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C263');
compsData=load('DistComps - Compactors x DLs.mat');
compsData=compsData.DistComps;
[LongComps, LatComps] = readvars('Landfill_Compactors.xlsx','Sheet','Sheet1','Range','X67:Y86');

DATA=[DATA,[10000*ones(1,length(compsData(:,1)));compsData'];[10000*ones(length(compsData(:,1)),1),compsData],10000*ones(length(compsData(:,1)))-10000*eye(length(compsData(:,1)))];
%LINE BELOW IS NEW
if fullscale==0
    datanodes=[1,randperm(238,20)+ones(1,20),248,254,251,242]; % row vector, must always start with node 1
    %datanodes=[1,56,100,13,231,88,99,200,66,44,133,156,180]; % row vector, must always start with node 1
    Long=Longitude(datanodes);
    Lat=Latitude(datanodes);
    count=1;
    activecomps=[];
    for k=length(datanodes)
        if datanodes(k)>239
            activecomps(count)=datanodes(k);
            count=count+1;
        end
    end
    lenactivecomps=length(activecomps);
else
    for j=1:length(DATA)
        datanodes=[idxscomps,j];
    end
    activecomps=linspace(1,20,1);
    lenactivecomps=length(activecomps);
end
%DATA=DATA([1,5,76,108,223],[1,5,76,108,223]);


nStops = length(datanodes); % you can use any number, but the problem size scales as N^2
nDLs=nStops-1;


%%% Calculate Distances Between Points
idxs1 = nchoosek(1:259,2);
idxs1=[idxs1;idxs1(:,2),idxs1(:,1)]; 
if fullscale==0
idxs=[];
id=[];
for p=1:nStops
    for w=1:nStops
        if p~=w
id=find(idxs1(:,1)==datanodes(p) & idxs1(:,2)==datanodes(w));
idxs=[idxs;idxs1(id,:)];
        end
%id=find(idxs1(:,2)==datanodes(p) & idxs1(:,1)==datanodes(w)); 
%idxs=[idxs;idxs1(id,:)];
    end
end
else
    idxs=idxs1;
end


nCombs=length(idxs);

dist=zeros(nCombs,1);

for i=1:nCombs
    dist(i)=DATA(idxs(i,1),idxs(i,2));
end

f=[dist;zeros(2*nDLs,1)]; % concatenate f with zeros at the end to account for u at the end of x vector (length #nodes)

Aeq=zeros(nDLs+2+lenactivecomps,nCombs+2*nDLs);
beq=zeros(nDLs+2+lenactivecomps,1);
A=zeros(5*nDLs+nDLs^2,nCombs+2*nDLs);
b=zeros(5*nDLs+nDLs^2,1);


%% Constraint 2
%Aeq(1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
Aeq(1,find(idxs(:,1)==1))=ones(1,nDLs);
beq(1)=m;
%% Constraint 3
%Aeq(2,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
Aeq(2,find(idxs(:,2)==1))=ones(1,nDLs);
beq(2)=m;

%% Constraint 3.1, 4

for j=2:nStops
    rowidxs=find(idxs(:,2)==datanodes(j));
    %A(j-1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    A(j-1,rowidxs)=-1*ones(1,nDLs);
    b(j-1)=-1;
    
    rowidxs=find(idxs(:,2)==datanodes(j));
    rowidxs2=find(idxs(:,1)==datanodes(j));
    
    %Aeq(j+1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    Aeq(j+1,rowidxs)=ones(1,nDLs);
    Aeq(j+1,rowidxs2)=-1*ones(1,nDLs);
    beq(j+1)=0;    
end

%% Constraint 5
cnt=nStops;
for i=2:nStops
    rowidxs=find(idxs(:,1)==datanodes(i));
    %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    A(cnt,rowidxs)=-1*ones(1,nDLs);
    b(cnt)=-1;
    cnt=cnt+1;
end

%% Constraint 6, 7, 6.1
%cnter=nStops;
for i=2:nStops
    x1i=find(idxs(:,1)==1 & idxs(:,2)==datanodes(i));
    xi1=find(idxs(:,2)==1 & idxs(:,1)==datanodes(i));
    ui=nCombs+i-1;
    %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    A(cnt,[x1i,xi1,ui])=[L-2,-1,1];
    b(cnt)=L-1;
    cnt=cnt+1;
    %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    A(cnt,[x1i,xi1,ui])=[-1,K-2,-1];
    b(cnt)=-2;
    cnt=cnt+1;
    si=nCombs+nDLs+i-1;
    %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    A(cnt,[x1i,xi1,si])=[H-2,-1,1];
    b(cnt)=H-1;
    cnt=cnt+1;
end


%% Constraints 9, 9.1,9.2
counter=2+nStops;
for i=2:nStops
    for j=2:nStops
        if i==j
            % do nothing
        else
            xij=find(idxs(:,1)==datanodes(i) & idxs(:,2)==datanodes(j));
            xji=find(idxs(:,2)==datanodes(i) & idxs(:,1)==datanodes(j));
            ui=nCombs+i-1;
            uj=nCombs+j-1;
            si=nCombs+nDLs+i-1;
            sj=nCombs+nDLs+j-1; 
            %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
            A(cnt,[xij,xji,ui,uj])=[L,L-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=L-1;
            cnt=cnt+1;
            %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
            A(cnt,[xij,xji,si,sj])=[H,H-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=H-1;
            cnt=cnt+1;
        end

    end
    if DATA(1,i)==10000
         %Aeq(cnter,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
         Aeq(cnter,si)=1;
         beq(cnter)=0;
         cnter=cnter+1;
    end
end
%% Constraint 8
    xstart=find(idxs(:,1)==1);
    xend=find(idxs(:,2)==1);
    for i=1:nDLs

    A(cnt,[xstart(i),xend(i)])=[1 1];
    b(cnt)=1;
    counter=cnt+1;
end


%%
int=ones(length(f),1);
intcon=find(int);

lb=[zeros(nCombs,1);ones(nDLs,1);zeros(nDLs,1)];
ub=ones(nCombs,1);

x=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

%%
segments = find(x(1:nCombs)<1.05 & x(1:nCombs)>.95); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
truetripscopy=truetrips;

%%
figure;hold on;
plot(Longitude,Latitude,'r*')
plot(Long,Lat,'g*')
for i=1:length(truetripscopy)
    xplot=[Latitude(truetripscopy(i,1)),Latitude(truetripscopy(i,2))];
    yplot=[Longitude(truetripscopy(i,1)),Longitude(truetripscopy(i,2))];
    plot(yplot,xplot,'b')
end

hold off


%%
% figure;hold on;
% if fullscale==0
%     plot(Long,Lat,'g*')
% end
% 
% plot(Longitude,Latitude,'r*')
% 
% route_long_lat=zeros(L+2,2*m);
% for h=1:m
% q=1;
% yplot(2)=0;
% while yplot(2)~=Longitude(1)
%     if q==1
%         i=find(truetrips(:,1)==1);
%         i=i(1);
%     else
%         i=find(truetrips(:,1)==nextDL);
%         i=i(1);
%     end
%     
%     xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
%     route_long_lat(q,2*h)=xplot(1);
%     yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
%     route_long_lat(q,2*h-1)=yplot(1);
%     if h==1
%         plot(yplot,xplot,'b')    
%     elseif h==2
%         plot(yplot,xplot,'c')
%     elseif h==3
%         plot(yplot,xplot,'g')
%     elseif h==4
%         plot(yplot,xplot,'m')  
%     elseif h==5
%         plot(yplot,xplot,'k')  
%     end
%     nextDL=truetrips(i,2);
%     truetrips(i,:)=[];
%     q=q+1;
% end
% route_long_lat(q,2*h)=xplot(2);
% route_long_lat(q,2*h-1)=yplot(2);
% end
% 
% hold off