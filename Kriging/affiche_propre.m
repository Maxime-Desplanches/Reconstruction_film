function [] = affiche_propre(type,A,cpos,spos,color)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        figure,
    
        imagesc(A)
        hold on
        for i=1:size(cpos,1)
            plot(cpos(i,2),cpos(i,1),color(i),'Marker','*')
        end
        plot(spos(:,2),spos(:,1),'ro')
        hold off
        colorbar
        type=['Krigeage avec variogramme ' type];
        title(type)
end

