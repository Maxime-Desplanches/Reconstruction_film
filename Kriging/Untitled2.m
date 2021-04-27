clear all
close all

p=1;
for i=2:40
    p=p*(1-i/365);
end

[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
figure,
surf(X,Y,Z)

%% verif proba
    cpos=zeros(nbc,2);
    % cval=rand(nbc,1);
    chance=1;

    for i=2:nbc
        chance=chance*(1-(i-1)/nbpoint);
    end
    
    %% affichage en surface
        [X,Y] = meshgrid(1:dimp(2)/pasg,1:dimp(1)/pasg);
    
%     figure,
%     surf(X,Y,A);
%     colormap winter

%% sparse

A=sparse(2,5,18);