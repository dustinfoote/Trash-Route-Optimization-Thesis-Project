function [x]=myintlinprog_compactor(allnodes,cluster_len,nodes_in_cluster)

m=1; %number of salemen/vehicles
H=9; %number of DL stops b4 visiting compactor
L=150;
K=5;
fullscale=0; % if problem is full scale, put 1, otherwise 0

DATA=load('Final Distance Matrix.mat');
DATA=DATA.DATA;
%DATACopy=DATA;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C263');
compsData=load('DistComps - Compactors x DLs.mat');
compsData=compsData.DistComps;
[LongComps, LatComps] = readvars('Landfill_Compactors.xlsx','Sheet','Sheet1','Range','X67:Y86');

DATA=[DATA,[10000*ones(1,length(compsData(:,1)));compsData'];[10000*ones(length(compsData(:,1)),1),compsData],10000*ones(length(compsData(:,1)))-10000*eye(length(compsData(:,1)))];
%DATA=[DATA,[2000*rand(1,length(compsData(:,1)));compsData'];[2000*rand(length(compsData(:,1)),1),compsData],2000*rand(length(compsData(:,1)))];


%for y=1:5
%LINE BELOW IS NEW
if fullscale==0
    %datanodes=[1,randperm(238,6)+ones(1,6),249];% [1,139,214,73,88,103,232,249];%% row vector, must always start with node 1
    %datanodes=[1,56,100,13,231,88,99,200,66,44,133,156,180]; % row vector, must always start with node 1
    datanodes=[1;allnodes(nodes_in_cluster)];
    Long=Longitude(datanodes);
    Lat=Latitude(datanodes);
    count=1;
    activecomps=[];
    for k=1:length(datanodes)
        if datanodes(k)>239
            activecomps(count)=datanodes(k);
            compLong(count)=Longitude(datanodes(k));
            compLat(count)=Latitude(datanodes(k));
            count=count+1;
            
        end
    end
    lenactivecomps=length(activecomps);
else
    count=1;
    datanodes=linspace(1,259,259);%[idxscomps,j];
    for j=1:length(DATA)
        if datanodes(j)>239
            activecomps(count)=datanodes(j);
            compLong(count)=Longitude(datanodes(j));
            compLat(count)=Latitude(datanodes(j));
            count=count+1;
        end
    
    lenactivecomps=length(activecomps);
    end
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
%nDLs-(lenactivecomps-floor((nDLs-lenactivecomps)/H))
f=[dist;zeros(2*nDLs,1)];%zeros(2*nDLs+2*floor((nDLs-lenactivecomps)/H)+2,1)]; % concatenate f with zeros at the end to account for u at the end of x vector (length #nodes)
%HERE (line above)
Aeq=zeros(2*nStops+2*lenactivecomps,nCombs+2*nDLs);%nCombs+2*nDLs+2*floor((nDLs-lenactivecomps)/H)+2); %HERE
beq=zeros(2*nStops+2*lenactivecomps,1);
A=zeros(2*nDLs^2+nDLs+lenactivecomps,nCombs+2*nDLs);%nCombs+2*nDLs+2*floor((nDLs-lenactivecomps)/H)+2);   %HERE
b=zeros(2*nDLs^2+nDLs+lenactivecomps,1);%(5*nDLs+nDLs^2,1);


%% Constraint 2 , 2.1
%Aeq(1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
Aeq(1,find(idxs(:,1)==1))=ones(1,nDLs);
beq(1)=m;

% Constraints below prevent subtours but require number of travelers to =
% number of compactors!
count=2;
% for i=1:lenactivecomps
%     rowidx=find(idxs(:,1)==activecomps(i) & idxs(:,2)==1);
%     Aeq(count,rowidx)=1;
%     beq(count)=1;
%     count=count+1;
% end
%% Constraint 3
%Aeq(2,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
Aeq(count,find(idxs(:,2)==1))=ones(1,nDLs);
beq(count)=m;
count=count+1;
%% Constraint 3.1, 4

cnt=1;
for j=2:nStops-lenactivecomps
    rowidxs=find(idxs(:,2)==datanodes(j));% & idxs(:,1)<240);
%     for k=1:length(rowidxs)
%         rowi=idxs(rowidxs(k),1);
%         if rowi>239
%             rowidxs(k)=0;
%         end
%     end
%     remove=find(rowidxs==0);
%     rowidxs(remove)=[];
%     %A(j-1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
%     Aeq(count,rowidxs)=ones(1,length(rowidxs));
%     beq(count)=1;
%     count=count+1;
    A(cnt,rowidxs)=-1*ones(1,nDLs);
    b(cnt)=-1;
    cnt=cnt+1;
    
%     rowidxs=find(idxs(:,2)==datanodes(j));
%     rowidxs2=find(idxs(:,1)==datanodes(j));
%     
%     %Aeq(j+1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
%     Aeq(count,rowidxs)=ones(1,nDLs);
%     Aeq(count,rowidxs2)=-1*ones(1,nDLs);
%     beq(count)=0; 
%     count=count+1;
end

for j=nStops-lenactivecomps+1:nStops
%     rowidxs=find(idxs(:,2)==datanodes(j));
%     %A(j-1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
%     A(j-1,rowidxs)=-1*ones(1,nDLs);
%     b(j-1)=-1;

% visit a compactor one or less times
%     rowidxs=find(idxs(:,2)==datanodes(j));
%     %A(j-1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
%     A(cnt,rowidxs)=ones(1,nDLs);
%     b(cnt)=1;
%     cnt=cnt+1;
%     
    rowidxs=find(idxs(:,2)==datanodes(j));
    rowidxs2=find(idxs(:,1)==datanodes(j));
    
    %Aeq(j+1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
    Aeq(count,rowidxs)=ones(1,nDLs);
    Aeq(count,rowidxs2)=-1*ones(1,nDLs);
    beq(count)=0; 
    count=count+1;
end
%% Constraint 5
%cnt=nStops;
 for i=2:nStops-lenactivecomps
     rowidxs=find(idxs(:,1)==datanodes(i));% & idxs(:,1)<240);
%     for k=1:length(rowidxs)
%         rowi=idxs(rowidxs(k),2);
%         if rowi>239
%                rowidxs(k)=0;
%         end
%     end
%     remove=find(rowidxs==0);
%     rowidxs(remove)=[];
    %A(j-1,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
%     Aeq(count,rowidxs)=ones(1,length(rowidxs));
%     beq(count)=1;
%     count=count+1;
    A(cnt,rowidxs)=-1*ones(1,nDLs);
    b(cnt)=-1;
    cnt=cnt+1;
end

%% Constraint 6, 7, 6.1
%cnt=1;
for i=2:nStops%-lenactivecomps
    x1i=find(idxs(:,1)==1 & idxs(:,2)==datanodes(i));
    xi1=find(idxs(:,2)==1 & idxs(:,1)==datanodes(i));
    ui=nCombs+i-1;
%             iloc=find(datanodes==datanodes(i));
%             if floor((iloc-1)/H)==0
%                 ci=0;
%             else
%                 ci=1;
%             end
%             ui=nCombs+iloc-1+floor((iloc-1)/H)+1*ci;
%             si=nCombs+nDLs+iloc-1+floor((iloc-1)/H)+1*ci;%floor((iloc-1-lenact)/H);
%             
%     
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
    if i<=nStops-lenactivecomps
    A(cnt,[x1i,xi1,si])=[H-2,-1,1];
    b(cnt)=H-1;
    cnt=cnt+1;
    else
    A(cnt,[x1i,xi1,si])=[H-1,H-1,1];
    b(cnt)=H-1;
    cnt=cnt+1;   
    end
end


%% Constraints 9, 9.1,9.2
%counter=2+nStops;
for i=2:nStops%-lenactivecomps
    for j=2:nStops%-lenactivecomps
        if i==j
            % do nothing
        else
            %if datanodes(i)<239 && datanodes(j)<239
            xij=find(idxs(:,1)==datanodes(i) & idxs(:,2)==datanodes(j));
            xji=find(idxs(:,2)==datanodes(i) & idxs(:,1)==datanodes(j));

%             iloc=find(datanodes==datanodes(i));
%             jloc=find(datanodes==datanodes(j));
%             if floor((iloc-1)/H)==0
%                 ci=0;
%             else
%                 ci=1;
%             end
%             if floor((jloc-1)/H)==0
%                 cj=0;
%             else
%                 cj=1;
%             end
            ui=nCombs+i-1;%iloc-1+floor((iloc-1)/H)+1*ci;
            uj=nCombs+j-1;%jloc-1+floor((jloc-1)/H)+1*cj;
            si=nCombs+nDLs+i-1;%iloc-1+floor((iloc-1)/H)+1*ci;%floor((iloc-1-lenact)/H);
            sj=nCombs+nDLs+j-1;%jloc-1+floor((jloc-1)/H)+1*cj; 
            %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
            
            %A(cnt,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
            if i<=nStops-lenactivecomps && j<=nStops-lenactivecomps
            A(cnt,[xij,xji,si,sj])=[H,H-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=H-1;
            cnt=cnt+1;            
            

            A(cnt,[xij,xji,ui,uj])=[L,L-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=L-1;
            cnt=cnt+1;
            elseif i<=nStops-lenactivecomps && j>nStops-lenactivecomps
            A(cnt,[xij,xji,si,sj])=[-1,H-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=H;%H-1;
            cnt=cnt+1; 
            A(cnt,[xij,xji,ui,uj])=[L,L-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=L-1;
            cnt=cnt+1;
            elseif j<=nStops-lenactivecomps && i>nStops-lenactivecomps
            A(cnt,[xij,xji,si,sj])=[H,2*H-1,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=H-1;
            cnt=cnt+1;
            A(cnt,[xij,xji,ui,uj])=[L,L-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=L-1;
            cnt=cnt+1;
            elseif j>nStops-lenactivecomps && i>nStops-lenactivecomps
            A(cnt,[xij,xji,si,sj])=[H-1,H-1,1,-1]; % On March 11 changed from -2 to -1 for u_j
            b(cnt)=H-1;
            cnt=cnt+1; 
            end
            %end
% COMBINATION OF BELOW 2 CONSTRAINTS WORKS!
% %BELOW 167-169 WORKS
%             A(cnt,[xij,xji,si,sj])=[-1,H-2,1,-1]; % On March 11 changed from -2 to -1 for u_j
%             b(cnt)=H-1;
%             cnt=cnt+1;
% %BELOW 171-173 WORKS
%             A(cnt,[xij,xji,si,sj])=[H,2*H-1,1,-1]; % On March 11 changed from -2 to -1 for u_j
%             b(cnt)=2*H;
%             cnt=cnt+1;
        end

    end
%     if datanodes(i)>239
%          %Aeq(cnter,nCombs+2*nDLs)=zeros(1,nCombs+2*nDLs);
%          Aeq(count,si)=1;
%          beq(count)=0;
%          count=count+1;
%     end
end

%% Constraint 8 , 8.1
%     xstart=find(idxs(:,1)==1);
%     xend=find(idxs(:,2)==1);
% for i=1:nDLs
% 
%     A(cnt,[xstart(i),xend(i)])=[1 1];
%     b(cnt)=1;
%     cnt=cnt+1;
%     %counter=cnt+1;
% end

%% xij+xji<=1
% for i=2:nStops-lenactivecomps
%     for j=2:nStops-lenactivecomps
%         if j~= i
%     xstart=find(idxs(:,1)==datanodes(i) & idxs(:,2)==datanodes(j));
%     xend=find(idxs(:,1)==datanodes(j) & idxs(:,2)==datanodes(i));
%     A(cnt,[xstart,xend])=[1 1];
%     b(cnt)=1;
%     cnt=cnt+1;
%         end
%     %counter=cnt+1;
%     end
% end

%%
int=ones(length(f),1);
intcon=find(int);

lb=[zeros(nCombs,1);ones(nDLs-lenactivecomps,1);zeros(lenactivecomps,1);ones(nDLs-lenactivecomps,1);zeros(lenactivecomps,1)];%ones(nDLs+floor((nDLs-lenactivecomps)/H)+1,1);zeros(nDLs+floor((nDLs-lenactivecomps)/H)+1,1)];
ub=[ones(nCombs,1);(nDLs)*ones(nDLs,1);H*ones(nDLs-lenactivecomps,1);zeros(lenactivecomps,1)];%nStops*ones(nDLs+floor((nDLs-lenactivecomps)/H)+1,1);H*ones(nDLs+floor((nDLs-lenactivecomps)/H)+1,1)];
% for c=1:floor((nDLs-lenactivecomps)/H)
%     ub(nCombs+nDLs+c*H+c)=0;
% end
% ub(end)=0;
options=optimoptions('intlinprog','MaxTime',172800);
[x,fval,exitflag,output]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,[],options);

%%
segments = find(x(1:nCombs)<1.05 & x(1:nCombs)>.95); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
truetripscopy=truetrips;

%%
figure;
plot(Longitude,Latitude,'r*');
hold on;
plot(Long,Lat,'g*');
plot(compLong,compLat,'k*');
for i=1:length(truetripscopy)
    xplot=[Latitude(truetripscopy(i,1)),Latitude(truetripscopy(i,2))];
    yplot=[Longitude(truetripscopy(i,1)),Longitude(truetripscopy(i,2))];
    plot(yplot,xplot,'b')
end


hold off


%% MATLAB TSP SUBTOUR ELIMINATION

truetrips=1;


while length(truetrips(:,1))>0
    x_tsp = logical(round(x));
Gsol = graph(idxs(find(x_tsp(1:nCombs)==1),1),idxs(find(x_tsp(1:nCombs)==1),2));
truetrips=table2array(Gsol.Edges);
q=1;
nextDL=0;

while nextDL~=1
    if q==1
        i=find(truetrips==1);
        i=i(1);
        if i>length(truetrips(:,1))
            i=i-length(truetrips(:,1));
        end
        if truetrips(i,1)==1
            nextDL=truetrips(i,2);
        elseif truetrips(i,2)==1
            nextDL=truetrips(i,1);
        end
    else
        i=find(truetrips==nextDL);
        if i>length(truetrips(:,1))
            i=i-length(truetrips(:,1));
        end
        %i=i(1);
        place=truetrips(i,:)';
        col=find(place==nextDL);
        if col==1
            nextDL=truetrips(i,2);
        elseif col==2
            nextDL=truetrips(i,1);
        end
    end
   
    truetrips(i,:)=[];
    q=q+1;
end        
        
b(cnt)=0;
for j=1:length(truetrips(:,1))
    rowidx=find(idxs(:,1)==truetrips(j,1) & idxs(:,2)==truetrips(j,2));
    rowidx2=find(idxs(:,2)==truetrips(j,1) & idxs(:,1)==truetrips(j,2));
    A(cnt,[rowidx,rowidx2])=[1,1];
    b(cnt)=b(cnt)+1;
end

b(cnt)=b(cnt)-1;
cnt=cnt+1;
[x,fval,exitflag,output]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,x,options);


x_tsp = logical(round(x));
Gsol = graph(idxs(find(x_tsp(1:nCombs)==1),1),idxs(find(x_tsp(1:nCombs)==1),2));
truetrips=table2array(Gsol.Edges);

q=1;
nextDL=0;

while nextDL~=1
    if q==1
        i=find(truetrips==1);
        i=i(1);
        if i>length(truetrips(:,1))
            i=i-length(truetrips(:,1));
        end
        if truetrips(i,1)==1
            nextDL=truetrips(i,2);
        elseif truetrips(i,2)==1
            nextDL=truetrips(i,1);
        end
    else
        i=find(truetrips==nextDL);
        if i>length(truetrips(:,1))
            i=i-length(truetrips(:,1));
        end
        %i=i(1);
        place=truetrips(i,:)';
        col=find(place==nextDL);
        if col==1
            nextDL=truetrips(i,2);
        elseif col==2
            nextDL=truetrips(i,1);
        end
    end
   
    truetrips(i,:)=[];
    q=q+1;
end    
end
%%
segments = find(x(1:nCombs)<1.05 & x(1:nCombs)>.95); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
truetripscopy=truetrips;

%%
figure;
plot(Longitude,Latitude,'r*');
hold on;
plot(Long,Lat,'g*');
plot(compLong,compLat,'k*');
for i=1:length(truetripscopy)
    xplot=[Latitude(truetripscopy(i,1)),Latitude(truetripscopy(i,2))];
    yplot=[Longitude(truetripscopy(i,1)),Longitude(truetripscopy(i,2))];
    plot(yplot,xplot,'b')
end


hold off



% 
% 
% tourIdxs = conncomp(Gsol.Edges(1:15,1:2));
% numtours = max(tourIdxs); % number of subtours
% fprintf('# of subtours: %d\n',numtours);
% 
% while numtours > ceil(nStops/H) % Repeat until there is just one subtour
%     % Add the subtour constraints
%     b = [b;zeros(numtours,1)]; % allocate b
%     A = [A;zeros(numtours,size(A,2))]; % A guess at how many nonzeros to allocate
%     for ii = 2:numtours
%         rowIdx = size(A,1) + 1; % Counter for indexing
%         subTourIdx = find(tourIdxs == ii); % Extract the current subtour
% %         The next lines find all of the variables associated with the
% %         particular subtour, then add an inequality constraint to prohibit
% %         that subtour and all subtours that use those stops.
%         variations = nchoosek(1:length(subTourIdx),2);
%         for jj = 2:length(variations)
%             whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
%                        (sum(idxs==subTourIdx(variations(jj,2)),2));
%             A(rowIdx,whichVar) = 1;
%         end
%         b(rowIdx) = length(subTourIdx) - 1; % One less trip than subtour stops
%     end
% 
%     % Try to optimize again
%     [x_tsp,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,[],options);
%     x_tsp = logical(round(x_tsp));
%     Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));
%     
%     % Visualize result
%     hGraph.LineStyle = 'none'; % Remove the previous highlighted path
%     highlight(hGraph,Gsol,'LineStyle','-')
%     drawnow
%     
%     % How many subtours this time?
%     tourIdxs = conncomp(Gsol);
%     numtours = max(tourIdxs); % number of subtours
%     fprintf('# of subtours: %d\n',numtours)
% end
% 
% 


%end
end
%%
% figure;hold on;
% if fullscale==0
%plot(Long,Lat,'g*')
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