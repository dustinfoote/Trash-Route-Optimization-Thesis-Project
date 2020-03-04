
m=1; %number of salemen/vehicles
DATA=load('Final Distance Matrix.mat');
DATA=DATA.DATA;
DATA=DATA(1:20,1:20);
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C243');
nStops = length(DATA); % you can use any number, but the problem size scales as N^2

%%% Calculate Distances Between Points
idxs = nchoosek(1:nStops,2);

dist=zeros(length(idxs),1);
for i=1:length(idxs)
    dist(i)=DATA(idxs(i,1),idxs(i,2));
end

f=dist;
Aeq=zeros(19,length(idxs));
beq=2*ones(19,1);
%beq(1)=2*m;

for stops = 1:nStops
        whichIdxs = (idxs == stops);
        whichIdxs = any(whichIdxs,2); 
        %constrtrips(stops) = sum(trips(whichIdxs)) == 2*m; % m vehicles leave and return to the first "depot" node
        whichIdxs=find(whichIdxs);
        Aeq(stops,whichIdxs)=ones(1,length(whichIdxs));
end

x=inlinprog(f,incon,[],[],Aeq,beq,0,1);
%%
segments = find(x); % Get indices of lines on optimal path

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