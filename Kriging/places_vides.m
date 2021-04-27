function [posvide] = places_vides(cpos,A)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
posvide=[];
for i=1:size(A,1)
    for j=1:size(A,2)
        if isempty(find(ismember(cpos,[i,j],'rows')))
            posvide=[posvide;i,j];
        end
    end
end

end

