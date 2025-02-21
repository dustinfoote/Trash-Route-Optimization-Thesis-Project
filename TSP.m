%% mTSP

%% Load Data and initialize

m=5; %number of salemen/vehicles
DATA=load('Final Distance Matrix.mat');
DATA=DATA.DATA;
DistComps=load('DistComps - Compactors x DLs.mat');
DistComps=DistComps.DistComps;
[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong_Distance Matrix with Compactors +Depot.xlsx','Sheet','Sheet2','Range','B5:C243');
nStops = length(DATA); % you can use any number, but the problem size scales as N^2

%%% Calculate Distances Between Points
idxs = nchoosek(1:nStops,2);

dist=zeros(length(idxs),1);
for i=1:length(idxs)
    dist(i)=DATA(idxs(i,1),idxs(i,2));
end

lendist = length(dist);

%%% Create Variables and Problem
tsp = optimproblem;
trips = optimvar('trips',lendist,1,'Type','integer','LowerBound',0,'UpperBound',1);

tsp.Objective = dist'*trips;

%%% Equality Constraints
constrips = sum(trips) == nStops+m-1; % restricts the total number of trips needed for m round trip subtours
tsp.Constraints.constrips = constrips;

constrtrips = optimconstr(nStops,1);
for stops = 1:nStops
    if stops==1 
        whichIdxs = (idxs == stops);
        whichIdxs = any(whichIdxs,2); 
        constrtrips(stops) = sum(trips(whichIdxs)) == 2*m; % m vehicles leave and return to the first "depot" node
        
    else
        whichIdxs = (idxs == stops);
        whichIdxs = any(whichIdxs,2); 
        constrtrips(stops) = sum(trips(whichIdxs)) == 2; % one trip to and from each DL node
    end
end
tsp.Constraints.constrtrips = constrtrips;

%constrcomp=
% %ATTEMPT TO FIX OPT VAR PROBLEM
% % routes=[];
% % aux = optimexpr(lendist); % Allocate aux
% % for idx=1:lendist
% % aux(idx)=trips(idx)*idx;
% % if aux(idx)~=0
% %     routes(j,:)=[idxs(idx,1),idxs(idx,2)];
% %     j=j+1;
% % end
% % end
% % routescons=sum(aux)==sum(trips);
% % tsp.Constraints.consroutes = routescons;
% 
% 
% 
% %ATTEMPT TO FIX OPT VAR PROBLEM
% % idx = 1:(nPeriods-1);
% % w(idx,:) = y(idx+1,:,'Low') - y(idx,:,'Low') + y(idx+1,:,'High') - y(idx,:,'High');
% % w(nPeriods,:) = y(1,:,'Low') - y(nPeriods,:,'Low') + y(1,:,'High') - y(nPeriods,:,'High');
% % switchcons = w - z <= 0;
% 
% %ATTEMPT TO FIX OPT VAR PROBLEM
% % routes=[];
% % rout=zeros(lendist,1);
% % for pairs=1:lendist
% %     rout(pairs)=trips(pairs)*pairs;
% % end
% % for pairs=1:lendist
% %     if rout(pairs)*pairs~=0
% %     routes(j,:)=[idxs(pairs,1),idxs(pairs,2)];
% %     j=j+1;
% %     end
% % end


%% Inequality constraints

% routes=optimconstr(length(DATA)+m-1,2);%zeros(length(sum(trips)),2); 
% % trips=optimconstr(length(dist),1);
% j=1;
% for pairs=1:length(trips)
%     if trips(pairs)==1                              % ERROR - cannot use opt var trips
%         routes(j,:)=[idxs(pairs,1),idxs(pairs,2)]; % Stores active trip node pairs but not in final route sequence
%         j=j+1;
%     end
% end
% 
% origroutes=routes; 
% totaltrips=length(routes);   % number of active node pairs for all 5 routes     
% for numpaths=1:m
%     
%     count=1;
%     routeX=[];
%     routeX(count,:)=[0,0];
%     while routeX(count,2)~=1 || routeX(1,1)==0 % create each vehicle's individual path vector: route##
% 
%         if count==1
%             [row,col,val]=find(routes==1);
%             routeX(1,:)=routes(row(1),:);
%         else
%             [row,col,val]=find(routes==routeX(count-1,2));
%             routeX(count,:)=routes(row(1),:);
%             if routeX(count,1)~=routeX(count-1,2)
%                 routeX(count,:)=[routeX(count,2),routeX(count,1)];
%             end
% 
%         end
%         count=count+1;
% 
%     end
%     
%  
%     if numpaths==1
%         route1=routeX;
%         elseif numpaths==2
%         route2=routeX;
%         elseif numpaths==3
%         route3=routeX;
%         elseif numpaths==4
%         route4=routeX;
%         elseif numpaths==5
%         route5=routeX;
%     end
%    
%    [row,col,val]=find(routes==1);
%    routes(row(1),:)=[];              % removes the first trip from depot used in routeX from routes
%    [row,col,val]=find(routes==routeX(length(routeX),1));
%    if routes(row(1),1)==1            % removes the last trip back to depot used in routeX from routes
%         routes(row(1),:)=[];
%    else
%         routes(row(2),:)=[]; 
%    end
%     
% end
% %%
% numpickups=[length(route1),length(route2),length(route3),length(route4),length(route5)];
% conspickupsmin=min(numpickups); 
% conspickupsmax=max(numpickups);
% 
% tsp.Constraints.conspickupsmin = conspickupsmin>=33;% sets K, minimum pickups per vehicle
% tsp.Constraints.conspickupsmax = conspickupsmax<=63;% sets L, maximum pickups per vehicle
% 
% %%
% routesvector=[route1;route2;route3;route4;route5];
% u=zeros(totaltrips,1);
% u(1)=1;
% cnt2=1; cnt3=1; cnt4=1; cnt5=1;
% for order=2:totaltrips     % determines u vector of elements=pickup # in subtour, row=node of DL to pickup at
%     if 1<order<length(route1)+1
%         uorder=routesvector(order,1);
%         u(uorder)=order;
%     elseif length(route1)<order<length(route1)+length(route2)+1
%         uorder=routesvector(order,1);
%         u(uorder)=cnt2;
%         cnt2=cnt2+1;
%     elseif length(route2)<order<length(route1)+length(route2)+length(route3)+1
%         uorder=routesvector(order,1);
%         u(uorder)=cnt3;
%         cnt3=cnt3+1;
%     elseif length(route3)<order<length(route1)+length(route2)+length(route3)+length(route4)+1
%         uorder=routesvector(order,1);
%         u(uorder)=cnt4;
%         cnt4=cnt4+1;
%     elseif length(route4)<order<length(route1)+length(route2)+length(route3)+length(route4)+length(route5)+1
%         uorder=routesvector(order,1);
%         u(uorder)=cnt5;
%         cnt5=cnt5+1;
%     end
% end
% 
% 
% L=63; %conspickupsmax;
% K=33; %conspickupsmin;
% 
% 
% firsttrips=[route1(1,:);route2(1,:);route3(1,:);route4(1,:);route5(1,:)];
% lasttrips=[route1(length(route1),:);route2(length(route2),:);route3(length(route3),:);route4(length(route4),:);route5(length(route5),:)];
% 
% consmax=optimconstr(nstops-1,1);
% consmin=optimconstr(nstops-1,1);
% for maxpts=2:nstops                       % determines if trip between nodes 1 and j or j and 1 is active
%     [xonei,col,val]=find(firsttrips==maxpts);
%     [xione,col,val]=find(lasttrips==maxpts);
%     if 0<xonei<6 %240
%         xonei=1;
%     else
%         xonei=0;
%     end
%     if 0<xione<6  %240
%         xione=1;
%     else
%         xione=0;
%     end
%     consmax(maxpts-1)=u(maxpts)+(L-2)*xonei-xione;
%     consmin(maxpts-1)=u(maxpts)+xonei+(2-K)*xione;
% end
% tsp.Constraints.consmax=consmax<=L-1;
% tsp.Constraints.consmin=consmin>=2;
% 
% %%
% conSECs=optimconstr(nstops-1,nstops);
% 
% for ii=2:nstops % determines if trip between i and j exists for i>=2, 1=<j=<nstops
%     for jj=1:nstops
%         if ii~=jj
%             xij=find(routesvector(:,1)==ii & routesvector(:,2)==jj);
%             xji=find(routesvector(:,2)==ii & routesvector(:,1)==jj);
%             if isempty(xij)==1
%                 xij=0;
%             else
%                 xij=1;
%             end
%             if isempty(xji)==1
%                 xji=0;
%             else
%                 xji=1;
%             end
%             conSECs(ii-1,jj)=u(ii)-u(jj)+L*xij+(L-2)*xji;
%         end
%     end
% end
%  
% 
% tsp.Constraints.conSECs=conSECs<=L-1;




%% Solve the Initial Problem
opts = optimoptions('intlinprog','Display','off');
tspsol = solve(tsp,'options',opts)




%% Visualize the Solution

segments = find(tspsol.trips); % Get indices of lines on optimal path

truetrips=idxs(segments',:);
% xplot=zeros(length(truetrips),1);
% yplot=zeros(length(truetrips),1);
figure;hold on;
for i=1:length(truetrips)
    xplot=[Latitude(truetrips(i,1)),Latitude(truetrips(i,2))];
    yplot=[Longitude(truetrips(i,1)),Longitude(truetrips(i,2))];
    plot(xplot,yplot,'b')
end
%%
hold on
plot(Latitude,Longitude,'r*')

lh=zeros(nStops,1);
lh = updateSalesmanPlot(lh,tspsol.trips,idxs,Latitude,Longitude); %CHANGE ORDER OF LAT AND LONG HERE AND IN PLOT!!!
title('Solution with Subtours');

tours = detectSubtours(tspsol.trips,idxs);
numtours = length(tours); % number of subtours
fprintf('# of subtours: %d\n',numtours);

% Index of added constraints for subtours
k = 1;
while numtours > 5 % repeat until there is just one subtour
    % Add the subtour constraints
    for ii = 1:numtours
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        a = false(length(idxs),1);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            a = a | whichVar;
        end
        tsp.Constraints.(sprintf('subtourconstr%i',k)) = sum(trips(a)) <= length(subTourIdx)-1;
        k = k + 1;
    end
    
        % Try to optimize again
    [tspsol,fval,exitflag,output] = solve(tsp,'options',opts);

    % Visualize result
    %lh = updateSalesmanPlot(lh,tspsol.trips,idxs,stopsLon,stopsLat);

    % How many subtours this time?
    tours = detectSubtours(tspsol.trips,idxs);
    numtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',numtours);
end

title('Solution with Subtours Eliminated');
hold off

disp(output.absolutegap)



