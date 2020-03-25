function [yy]=updatedactivecompactors(m,x,route_long_lat,rate,priorCombs,idxs)
%%
compsData=load('DistComps - Compactors x DLs.mat');
compsData=compsData.DistComps;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C243');
[LongComps, LatComps] = readvars('Landfill_Compactors.xlsx','Sheet','Sheet1','Range','X67:Y86');
%compsData=compsData(:,truetripscopy(m:end,1)');


%%
truetripsy=[];

%for v=1:m
v=3; 
p=any(route_long_lat(:,3*v),2);
    p=sum(p)-2;
Aeq=[];   
beq=[];

DLcomps=[];

nComps=length(compsData(:,1));
idxscomps=[];


for j=1:nComps
    idxscomps=[idxscomps;j];
end
for k=1:floor(p/rate)%+1

    idxu=route_long_lat(rate*k+1,3*v);
    DLcomps=[DLcomps;idxu];

end
   % u=idxscomps(1:p);%define for each trip
cnt=1;
%solution=[];
%for r=1:length(DLcomps)
    
idxsCOMPS=[];
idxsCOMPS2=[];
onescomps=ones(nComps,1);
for q=1:length(DLcomps)  %length(compsData)
    idxsCOMPS=[idxsCOMPS;DLcomps(q)*onescomps,idxscomps];
    idxx10=find(route_long_lat(:,3*v)==DLcomps(q));
	idxsCOMPS2=[idxsCOMPS2;route_long_lat(idxx10+1,3*v)*onescomps,idxscomps];
end
idxsCOMPS=idxsCOMPS-[ones(length(idxsCOMPS),1),zeros(length(idxsCOMPS),1)];
idxsCOMPS2=idxsCOMPS2-[ones(length(idxsCOMPS2),1),zeros(length(idxsCOMPS2),1)];
nCombs=length(idxsCOMPS);

dist=zeros(nCombs,1);
dist2=zeros(nCombs,1);
for i=1:nCombs
    dist(i)=compsData(idxsCOMPS(i,2),idxsCOMPS(i,1));
    dist2(i)=compsData(idxsCOMPS2(i,2),idxsCOMPS2(i,1));
end
idxsCOMPS=idxsCOMPS+[ones(length(idxsCOMPS),1),zeros(length(idxsCOMPS),1)];
idxsCOMPS2=idxsCOMPS2+[ones(length(idxsCOMPS2),1),zeros(length(idxsCOMPS2),1)];
f=dist+dist2; % only dist data for one direction (col 1 idxs to col 2)

for w=1:length(DLcomps)
    
[val,loc]=min(f(20*w-19:20*w));
truetripsy([2*cnt-1:2*cnt],[1:2])=[idxsCOMPS(loc+20*(w-1),:);idxsCOMPS2(loc+20*(w-1),:)];
cnt=cnt+1;
  %solution=[dist,DL # (from 239 list), Compactor]
%solution(r,:)=[val,idxsCOMPS(loc,:),idxsCOMPS2(loc,:)];
end
%end
% yy=solution;
% cnt=1;
% for k=1:floor(p/rate)
% 
%  %% Sum of Full DL and subsequent DL active DL-Comp trips is 2       
%             Aeq(cnt,:)=zeros(1,nCombs);
%             %beq(cnt)=zeros(1,1);
%             idxx=find(idxsCOMPS(:,1)==DLcomps(k));
%             Aeq(cnt,idxx)=ones(1,length(idxx));
%             beq(cnt)=1;
%             
%             %%
%             cnt=cnt+1;
%             Aeq(cnt,:)=zeros(1,nCombs);
%             %beq(cnt)=zeros(1,1);
%             idxx2=find(route_long_lat(:,3*v)==DLcomps(k));
%             idxx2=find(idxsCOMPS(:,1)==route_long_lat(idxx2+1,3*v));
%             Aeq(cnt,idxx2)=ones(1,length(idxx2));
%             %beq(cnt)=2;
%             cnt=cnt+1;
%             bcnt=cnt;
% %% Each compactor i paired with an active DL (both full and subsequent) combination must sum to 2
% for i=1:20
%             Aeq(cnt,:)=zeros(1,nCombs);
%             idxx3=find(idxsCOMPS([idxx,idxx2],2)==i);
%             %idxx=find(idxsCOMPS(:,1)==DLcomps(k) & idxsCOMPS(:,2)
%             A(cnt,idxx3)=ones(1,2);
%             cnt=cnt+1;
% end
% 
% end
% beq=ones(cnt-1,1);
% beq(bcnt:end)=2*beq(bcnt:end);
% %% Constraint each compactor visited must be stopped at and left from at least once (summing to two)
% % cnt2=1;
% % A=[];
% % for i=1:nComps
% %     A(cnt2,:)=zeros(1,nCombs);
% %     idxx=find(idxsCOMPS(:,2)==i);
% %     A(cnt2,idxx)=ones(1,length(idxx));
% %     cnt2=cnt2+1;
% % end
% % b=2*ones(cnt2-1,1);
% 
% 
% %%
% int=ones(length(f),1);
% intcon=find(int);
% 
% lb=zeros(nCombs,1);
% ub=ones(nCombs,1);
% 
% y=intlinprog(f,intcon,[],[],Aeq,beq,lb,ub);
% %if v==1
%     yy=zeros(length(nCombs),1);
% %end
% yy=y+yy;
% %end

%%
segments = find(x(1:priorCombs)<1.05 & x(1:priorCombs)>.95); % Get indices of lines on optimal path
%segmentsy = find(yy(1:nCombs)<1.05 & yy(1:nCombs)>.95);

truetrips=idxs(segments',:);
%truetripsy=idxsCOMPS(segmentsy',:);

figure;hold on;
for i=1:length(truetrips)
    xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
    yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
    plot(yplot,xplot,'b')
end
for i=1:length(truetripsy)
    xplot2=[Latitude(truetripsy(i,1)),LatComps(truetripsy(i,2))];
    yplot2=[Longitude(truetripsy(i,1)),LongComps(truetripsy(i,2))];
    plot(yplot2,xplot2,'g')
end

plot(Longitude,Latitude,'r*')
plot(LongComps,LatComps,'k*')
hold off


end