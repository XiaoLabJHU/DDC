%This is going to be a little script that will find localizaitons that are
%on the same frame and within 7 Resollution of each other.

Frames_w_Multi=[];
for nnnn2=length(Frame_Information):-1:1
    
    ffs=Frame_Information(nnnn2);
    inds42=find(Frame_Information==ffs);
    
    if length(inds42)>1 && isempty(Frames_w_Multi==ffs)
        %Now we need to determine whether any of these localizaitons are
        %within a  resonable distance of each other.
        
        Go_w_Multi=0;
        for iosnei=inds42(:)'
            inds42=find(Frame_Information==ffs);
            if length(inds42)>1
            for iinoeubg=inds42(:)'
%                 iosnei
%                 iinoeubg
%                 pause(.1)
                if iinoeubg~=iosnei
                    Datapoint=(LocalizationsFinal(iosnei,:)-LocalizationsFinal(iinoeubg,:));
                    Distance=(Datapoint(1)^2+Datapoint(2)^2)^.5;
                    
                    if Distance<Resolution*5
                     
                        LocalizationsFinal(iinoeubg,:)=[];
                        Frame_Information(iinoeubg)=[];
                        iinoeubg
                    end
                end
            end
            end
        end
        
        
    end
end

