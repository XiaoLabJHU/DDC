%This is a function that links the localizations into trajectories and is
%at the heart of DDC. This is ALG 1 of supporting material!!!!!!!!!!!!

function [Loc, Frame, Trajectory]=Eliminate_Blinking_De_Loc15(LocalizationsFinal, frame, Constants, Resolution, A, Deviation_in_Probability, Traj_num)

%The three different arrays below are used to determine the blinks that
%belong to each molecule.

Linked_Loc=zeros(1,length(LocalizationsFinal));
Loc_is_Blink=zeros(1,length(LocalizationsFinal));
%I have updated this CHB-2020
Trajectory= Traj_num;%zeros(1,length(LocalizationsFinal));

%These should all ready be in the proper order.
Loc=LocalizationsFinal;
Frame=frame;

counts=0;
%The frame differences to investigate

for i=1:round(A)

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
            
            %Make sure the localization is not allready considered a blink
            if Linked_Loc(ii)~=1 && Loc_is_Blink(ii+investigate)~=1
                
                %Datapoint=(Loc(ii,:)-Loc(ii+investigate,:));
                %Distance=(Datapoint(1)^2+Datapoint(2)^2)^.5;
                
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
                
                AAAA=Frame(iit);
                AAAA=[AAAA(:), AAAA(:)*0];
                BBBB=Frame(iit2);
                BBBB=[BBBB(:), BBBB(:)*0];
                
                Framer=sqrt( bsxfun(@plus,sum(AAAA.^2,2),sum(BBBB.^2,2)') - 2*((AAAA*BBBB')) );
                
                DDt=floor(Distance(:)/Resolution)+1;
                DDt(DDt>size(Deviation_in_Probability,2))=size(Deviation_in_Probability,2);
                
                if max(Framer)<A %&& min(Deviation_in_Probability(Framer,floor(Distance/Resolution)+1))>0
                    
                    %Make sure you use the ind of the matrix M
                    idx = sub2ind(size(Deviation_in_Probability), Framer(:), DDt);
                    
                    if mean(Deviation_in_Probability(idx))/abs(1+Constants(ii+investigate))>.5%(counts)% && ggo>0
                        
                        %This loc is now considered a blink
                        Loc_is_Blink(ii+investigate)=1;
                        
                        %The one is considered to be linked
                        Linked_Loc(ii)=1;
                        
                        if Trajectory(ii)==0
                            Trajectory(ii)=Traj_num(ii);
                        end
                        
                        if Trajectory(ii+investigate)~=0
                            inds_of_traj=(Trajectory==Trajectory(ii+investigate));
                        else
                            inds_of_traj=ii+investigate;
                        end
                        
                        Trajectory(inds_of_traj)=Trajectory(ii);
                        
                    end
                end
            end
            
            investigate=investigate+1;
            
            if ii+investigate>length(Frame)
                break
            end
            
        end
        
    end
    
   
end

%If anyone cares to look at the trajectories, please feel free
%Trajectoryt=Trajectory;
% figure(1)
%  scatter(Loc(:,1),Loc(:,2),10,Trajectory,'filled')
%  drawnow
%   colormap jet
%  figure(2)
% scatter(Loc(:,1),Loc(:,2),10,Frame,'filled')
%  axis equal
%  drawnow
%  colormap jet
% % caxis([0 length(Loc)])
% % pause(.01)
% pause

Loc(Loc_is_Blink>0,:)=[];
Frame(Loc_is_Blink>0)=[];


