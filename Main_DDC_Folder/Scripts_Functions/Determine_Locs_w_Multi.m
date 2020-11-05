%This is going to be a little script that will find localizaitons that are
%on the same frame and within 7 Resollution of each other.

Frames_w_Multi=[];
for nnnn2=1:length(Frame_Information)
    
    ffs=Frame_Information(nnnn2);
    inds42=find(Frame_Information==ffs);
    
    if length(inds42)>1 && sum(Frames_w_Multi==ffs)==0
        %Now we need to determine whether any of these localizaitons are
        %within a  resonable distance of each other.
        Go_w_Multi=0;
        for iosnei=inds42
            for iinoeubg=inds42
                if iinoeubg~=iosnei
                    Datapoint=(LocalizationsFinal(iosnei,:)-LocalizationsFinal(iinoeubg,:));
                    Distance=(Datapoint(1)^2+Datapoint(2)^2)^.5;
                    if Distance<250
                        Go_w_Multi=1;
                        break
                    end
                end
            end
        end
        
        if Go_w_Multi==1
            Frames_w_Multi(end+1)=ffs;
        end
        
    end
end

