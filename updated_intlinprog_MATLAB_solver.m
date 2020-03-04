clear;clc;
m=4; %number of salemen/vehicles
DATA=load('Final Distance Matrix.mat');
DATA=DATA.DATA;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C243');

%LINE BELOW IS NEW
DATA=DATA(1:17,1:17);
Longitude=Longitude(1:17);
Latitude=Latitude(1:17);

nStops = length(DATA); % you can use any number, but the problem size scales as N^2
nDLs=nStops-1;


%%% Calculate Distances Between Points
idxs = nchoosek(1:nStops,2);
idxs=[idxs;idxs(:,2),idxs(:,1)]; 
nCombs=length(idxs);

dist=zeros(nCombs,1);

for i=1:nCombs
    dist(i)=DATA(idxs(i,1),idxs(i,2));
end

f=[dist;zeros(nDLs,1)]; % concatenate f with zeros at the end to account for u at the end of x vector (length #nodes)

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
    rowidxs=find(idxs(:,2)==j);
    Aeq(j+1,rowidxs)=ones(1,nDLs);
    beq(j+1)=1;
end

%% Constraint 5
cnt=nStops+2;
for i=2:nStops
    rowidxs=find(idxs(:,1)==i);
    Aeq(cnt,rowidxs)=ones(1,nDLs);
    beq(cnt)=1;
    cnt=cnt+1;
end

%% Constraint 6, 7
A=zeros(2*nDLs+nDLs^2,nCombs+nDLs);
b=zeros(2*nDLs+nDLs^2,1);
L=6;
K=2;
cnter=nStops;
for i=2:nStops
    x1i=find(idxs(:,1)==1 & idxs(:,2)==i);
    xi1=find(idxs(:,2)==1 & idxs(:,1)==i);
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
            xij=find(idxs(:,1)==i & idxs(:,2)==j);
            xji=find(idxs(:,2)==i & idxs(:,1)==j);
            ui=nCombs+i-1;
            uj=nCombs+j-1;
            A(counter,[xij,xji,ui,uj])=[L,L-2,1,-2];
            b(counter)=L-1;
            counter=counter+1;
        end
    end
end
%% Constraint 8
for i=1:nDLs
    xstart=find(idxs(:,1)==1);
    xend=find(idxs(:,2)==1);
    A(counter,[xstart(i),xend(i)])=[1 1];
    b(counter)=1;
    counter=counter+1;
end


%%
int=ones(length(f),1);
intcon=find(int);

lb=[zeros(nCombs,1);ones(nDLs,1)];
ub=ones(nCombs,1);

x=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

%%
segments = find(x(1:nCombs)); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
xplot=zeros(length(truetrips),1);
yplot=zeros(length(truetrips),1);
figure;hold on;
for i=1:length(truetrips)
    xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
    yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
    plot(yplot,xplot,'b')
end

plot(Longitude,Latitude,'r*')
hold off