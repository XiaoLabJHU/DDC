%This is going to be a script that thresholds the localizations in time and
%space.

function [Thresh_Loc, Thresh_frame]=Threshold(LocalizationsFinal, frame, Thresh_Time, Thresh_dist)

Thresh_Loc=[];
[mmms, index]=sort(frame);
Thresh_frame=mmms;

    Thresh_Loc=LocalizationsFinal(index,:);
temps2=0;
temps=length(Thresh_Loc);
%while temps-temps2~=0
temps=length(Thresh_Loc);
for i=length(Thresh_frame):-1:2
   % i
    elim=0;
    for ii=i-1:-1:1
        Datapoint=(Thresh_Loc(i,:)-Thresh_Loc(ii,:));
        Distance=(Datapoint(1)^2+Datapoint(2)^2)^.5;
        Time=abs(Thresh_frame(i)-Thresh_frame(ii));
        
        if Distance<=Thresh_dist && Time<=Thresh_Time && Time>0
            elim=1;
            break
        end
    end
    
    if elim==1
        %{
        hold on
        subplot(2,1,1)
        plot(Thresh_Loc(i,1),Thresh_Loc(i,2),'r.','MarkerSize',10)
        drawnow
        %}
        Thresh_Loc(i,:)=[];
        Thresh_frame(i)=[];
    end
    %{
     final_cord=Thresh_Loc;
final_Frame_Information=Thresh_frame;
subplot(2,1,2)
scatter(final_cord(:,1),final_cord(:,2),10, (final_Frame_Information),'filled')
drawnow
%}
%end
temps2=length(Thresh_Loc);
end