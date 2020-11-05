%This is a little function to grab distances for the calculation of the
%real pairwise distance distribution.
function [DistanceMatrix]=Make_Distance_Matrix_Individual(Localizations, DistanceMatrix, frame, B)

[~, Inds]=sort(frame);
frame=frame(Inds);
Localizations=Localizations(Inds,:);
for i=1:length(Localizations)
    for jj=i:length(frame)
        if abs(frame(i)-frame(jj))<B+1 && (frame(jj)-frame(i))>0
            Datapoint=(Localizations(i,:)-Localizations(jj,:));
            %Changing this so it works in 3 dim
            
            DistanceMatrix{abs(frame(i)-frame(jj))}(end+1)=(sum(Datapoint.^2))^.5;
        end
    end
end


