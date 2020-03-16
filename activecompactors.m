function [y]=activecompactors(x,u)
%%
compsData=load('DistComps - Compactors x DLs.mat');
compsData=compsData.compsData;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C243');
[LongComps, LatComps] = readvars('Landfill_Compactors.xlsx','Sheet','Sheet1','Range','X67:Y86');

%LINE BELOW IS NEW
DATA=DATA(1:20,1:20);
%Longitude=Longitude(1:17);
%Latitude=Latitude(1:17);

%nStops = length(DATA); % you can use any number, but the problem size scales as N^2
%nDLs=nStops-1;


%%% Calculate Distances Between Points
%idxs = nchoosek(1:nStops,2);
%idxs=[idxs;idxs(:,2),idxs(:,1)]; 

%%

lentour=length(u);
n=ceil(lentour/20); 



DLcomps=[];

nComps=length(compsData);
idxscomps=[];
for j=1:nComps
    idxscomps=[idxscomps;j];
end


for k=1:n
    if k<n
    idxu=find(u==20*k);
    DLcomps=[DLcomps;idxu];
    else
        idxu=find(u(end));
        DLcomps=[DLcomps;idxu];
    end
end
idxsCOMPS=[];
onescomps=ones(nComps,1);
for q=1:length(DLcomps)
    idxsCOMPS=[idxsCOMPS;DLcomps(q)*onescomps,idxscomps];
end

nCombs=length(idxsCOMPS);

dist=zeros(nCombs,1);

for i=1:nCombs
    dist(i)=compsData(idxs(i,1),idxs(i,2));
end

f=[dist]; % only dist data for one direction (col 1 idxs to col 2)


cnt=1;
for k=1:n
    if k<n
        idxu=find(u==20*k);
        s=length(idxu);
        for i=1:s
            Aeq(cnt,:)=zeros(1,nCombs);
            %beq(cnt)=zeros(1,1);
            idxx=find(idxsCOMPS(:,1)==idxu(i));
            Aeq(cnt,idxx)=ones(length(idxx),length(idxx));
            %beq(cnt)=1;
            cnt=cnt+1;
        end
    else 
        %copy above ^^ but replace idxu with idxu=find(u(end));
        idxu=find(u(end));
        s=length(idxu);
        for i=1:s
            idxx=find(idxsCOMPS(:,1)==idxu(i));
            Aeq(cnt,idxx)=ones(length(idxx),length(idxx));
            %beq(cnt)=1;
            cnt=cnt+1;
        end
    end
end
beq=ones(cnt-1,1);

%%
int=ones(length(f),1);
intcon=find(int);

lb=[zeros(nCombs,1)];
ub=ones(nCombs,1);

y=intlinprog(f,intcon,[],[],Aeq,beq,lb,ub);

%%
segments = find(x(1:nCombs)); % Get indices of lines on optimal path
segmentsy = find(y(1:nCombs));

truetrips=idxs(segments',:);
truetripsy=idxsCOMPS(segmentsy',:);

xplot=zeros(length(truetrips),2);
yplot=zeros(length(truetrips),2);

xplot2=zeros(length(truetripsy),2);
yplot2=zeros(length(truetripsy),2);

figure;hold on;
for i=1:length(truetrips)
    xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
    yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
    plot(yplot,xplot,'b')
end
for i=1:length(truetripsy)
    xplot2=[Latitude(truetripsy(i,1)),Latitude(truetripsy(i,2))];
    yplot2=[Longitude(truetripsy(i,1)),Longitude(truetripsy(i,2))];
    plot(yplot2,xplot2,'g')
end

plot(Longitude,Latitude,'r*')
plot(LongComps,LatComps,'k*')
hold off


end
        