%Here we will dertermine the "relative" density for each localization

DistanceDens=2*Resolution;
density=zeros(length(LocalizationsFinal),1);
for isonepp3=1:length(LocalizationsFinal)
    
    for onerrrr=1:length(LocalizationsFinal)
    
    Datapoint=(LocalizationsFinal(isonepp3,:)-LocalizationsFinal(onerrrr,:));
    Distance=(sum(Datapoint.^2))^.5;
    
    if Distance<=DistanceDens && abs(Frame_Information(isonepp3)-Frame_Information(onerrrr))>A
        density(isonepp3)=density(isonepp3)+1;
    end
    end
    
end

density=(density-min(density))./(max(density)-min(density));