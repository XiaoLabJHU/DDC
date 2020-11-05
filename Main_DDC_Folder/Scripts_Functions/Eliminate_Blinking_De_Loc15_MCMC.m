%This is going to be a script that thresholds the localizations in time and
%space.

function [Loc, Probs]=Eliminate_Blinking_De_Loc15_MCMC(LocalizationsFinal, frame, Resolution, A, Deviation_in_Probability, Traj_num)

%The three different arrays below are used to determine the blinks that
%belong to each molecule.

Linked_Loc=zeros(1,length(LocalizationsFinal));
Loc_is_Blink=zeros(1,length(LocalizationsFinal));
Trajectory=zeros(1,length(LocalizationsFinal));
Probs=zeros(1,length(LocalizationsFinal));
%These should all ready be in the proper order.
Loc=LocalizationsFinal;
Frame=frame;

counts=0;
%The frame differences to investigate

for i=1:round(A)
    %hey=1;
    %If there is no more that can be linked get out to save time!
    for ii=1:length(Frame)
        
        %Get up to the current investigation spot
        investigate=0;
        
        while Frame(ii+investigate)-Frame(ii)<i
            
            
            investigate=investigate+1;
            if ii+investigate>length(Frame)
                break
            end
        end
        
        if ii+investigate>length(Frame)
            break
        end
        
        while Frame(ii+investigate)-Frame(ii)==i
            counts=counts+1;
            
            %             if counts>length(Thresh_dist)
            %                 Thresh_dist(counts)=1;
            %             end
            
            %Make sure the localization
            if Linked_Loc(ii)~=1 && Loc_is_Blink(ii+investigate)~=1
                
                %Datapoint=(Loc(ii,:)-Loc(ii+investigate,:));
                
                %I am changing this so that it works in three dim
                % Distance=(sum(Datapoint.^2))^.5;
               % Distance=(Datapoint(1)^2+Datapoint(2)^2)^.5;
                
                if Trajectory(ii)~=0
                    iit=find(Trajectory==Trajectory(ii));
                else
                    iit=ii;
                end
                
                if  Trajectory(ii+investigate)~=0
                    iit2=find(Trajectory==Trajectory(ii+investigate));
                else
                    iit2=ii+investigate;
                end
                
                %Datapoint= (Loc(iit,:))-Loc(iit2,:);
                
                Distance = sqrt( bsxfun(@plus,sum(Loc(iit,:).^2,2),sum(Loc(iit2,:).^2,2)') - 2*(Loc(iit,:)*Loc(iit2,:)') );
                %Distance=pdist2(Loc(iit,:),Loc(iit2,:));
                %Framer=pdist2(Frame(iit)',Frame(iit2)');
                AAAA=Frame(iit);
                AAAA=[AAAA(:), AAAA(:)*0];
                BBBB=Frame(iit2);
                BBBB=[BBBB(:), BBBB(:)*0];
                
                Framer=sqrt( bsxfun(@plus,sum(AAAA.^2,2),sum(BBBB.^2,2)') - 2*((AAAA*BBBB')) );
                
                %Distance=(Datapoint(1).^2+Datapoint(2).^2).^.5;
                %Framer=Frame(iit2)-Frame(iit);
                
                DDt=floor(Distance(:)/Resolution)+1;
                DDt(DDt>size(Deviation_in_Probability,2))=size(Deviation_in_Probability,2);
                
                if max(Framer)<A %&& min(Deviation_in_Probability(Framer,floor(Distance/Resolution)+1))>0
                    idx = sub2ind(size(Deviation_in_Probability), [Framer(:)], [DDt]);
                   
                    if mean(Deviation_in_Probability(idx))>.5%(counts)% && ggo>0
                        
                        %This loc is now considered a blink
                        Loc_is_Blink(ii+investigate)=1;
                        
                        Probs(ii+investigate)=mean(Deviation_in_Probability(idx));
                        %The one is considered to be linked
                        Linked_Loc(ii)=1;
                        
                        %Lets get a better estimate of the mean locaiton of the
                        %true molecule
                        % Loc(ii+investigate,:)=Loc(ii,:).*(Means_Adjust(ii)/(Means_Adjust(ii)+1))+Loc(ii+investigate,:)./(Means_Adjust(ii)+1);
                        
                        if Trajectory(ii)==0
                            Trajectory(ii)=Traj_num(ii);
                        end
                        
                        if Trajectory(ii+investigate)~=0
                            inds_of_traj=(Trajectory==Trajectory(ii+investigate));
                        else
                            inds_of_traj=ii+investigate;
                        end
                        
                        Trajectory(inds_of_traj)=Trajectory(ii);
                        
                        % Loc(ii,:)=Loc(ii+investigate,:);
                        
                        %Means_Adjust(ii)=Means_Adjust(ii+investigate);
                        %                     scatter(Loc(:,1),Loc(:,2),10,log(Frame),'filled')
                        %                     axis equal
                        %                     drawnow
                        %                     pause(.1)
                    elseif mean(Deviation_in_Probability(idx))>Probs(ii+investigate)
                        Probs(ii+investigate)=mean(Deviation_in_Probability(idx)); 
                    end
                    
                end
            end
            
            investigate=investigate+1;
            
            if ii+investigate>length(Frame)
                break
            end
            
        end
        
    end
    
    %Time to eliminate the locs that are blinks and have been linked
    %You may want to include this step or not, depending on how many
    %localizaitons you have and their blinking kinetics.
    
    %     Elim=Linked_Loc+Loc_is_Blink;
    %     Linked_Loc(Elim>1)=[];
    %     Loc_is_Blink(Elim>1)=[];
    %     Loc(Elim>1,:)=[];
    %     Frame(Elim>1)=[];
    %     Means_Adjust(Elim>1)=[];
    %     density(Elim>1)=[];
end
%Trajectoryt=Trajectory;
% scatter(Loc(:,1),Loc(:,2),10,Trajectory,'filled')
% %scatter(Loc(:,1),Loc(:,2),10,Frame,'filled')
% axis equal
% drawnow
% colormap jet
% caxis([0 length(Loc)])
% pause(.01)

Loc=ones(1,length((LocalizationsFinal)));
Loc(Loc_is_Blink>0)=0;


