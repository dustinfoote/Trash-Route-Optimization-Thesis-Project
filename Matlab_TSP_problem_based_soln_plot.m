
segments = find(tspsol.trips); % Get indices of lines on optimal path

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