clear all
close all

%% génération de la grille

nbc=20;
dimp=[32 32];
pasg=1;
A=zeros(dimp(1)/pasg,dimp(2)/pasg);
nbpoint=size(A,1)*size(A,2);
dossier='indsources_32x32.mat';
type='exp';
position='connu';
k=zeros(3,1);
Cexp=0;
Cgauss=0;
Csph=0;
n=50;
erreur=zeros(3);
for kjq=1:n
    %% création source
    
    nbs=3;
    [spos]=creation_position('random',nbs,A);
    
    %% création capteur
    
    [cpos]=creation_position(position,nbc,A,dossier);
    
    %% évaluation valeur au capteur
    
    cval = creation_valeur(cpos,spos,nbc,nbs,type);
    
    %% lisibilité de la carte
    
%     cval=cval/mean(cval);
    [~,indice]=sort(cval);
    color =[''];
    
    for i=1:nbc
        if i<floor(nbc/5)
            color(indice(i))='c';
        elseif i<floor(2*nbc/5)
            color(indice(i))='g';
        elseif i<floor(3*nbc/5)
            color(indice(i))='m';
        elseif i<floor(4*nbc/5)
            color(indice(i))='b';
        else
            color(indice(i))='k';
        end
    end
    % figure,
    % hold on
    % for i=1:nbc
    %    plot(cpos(i,1),cpos(i,2),color(i),'Marker','*')
    %
    %     axis equal
    % end
    % plot(spos(:,1),spos(:,2),'ro')
    % axis equal
    % hold off
    
    for i=1:nbc
        A(cpos(i,1),cpos(i,2))=cval(i);
    end
    
    
    %% Variogram exp
    % figure,
    % [v_struct] =variogram(cpos,cval,'plotit',true,'nrbins',20);
    [v_struct] =variogram(cpos,cval,'nrbins',20);
    posvide=places_vides(cpos,A);
    
    %% kriging exp
    
    str_model='exponential';
    [~,~,~,vfit] = synthetic_variogram(v_struct.distance,v_struct.val,1,var(cval),[],'model',str_model);
    posvide=places_vides(cpos,A);
    
    [zi,~] = kriging(vfit,cpos(:,1),cpos(:,2),cval,posvide(:,1),posvide(:,2),18);
    
    for i=1:nbpoint-nbc
        A(posvide(i,1),posvide(i,2))=zi(i);
    end
    
%      affiche_propre('exponentiel',A,cpos,spos,color)
    Aexp=normalisation(A,'normal');
    %%  gauss
%     str_model='gaussian';
%     [~,~,~,vfit] = variogramfit(v_struct.distance,v_struct.val,1,var(cval),[],'model',str_model);
%     
%     [zi,~] = kriging(vfit,cpos(:,1),cpos(:,2),cval,posvide(:,1),posvide(:,2),18);
%     
%     for i=1:nbpoint-nbc
%         A(posvide(i,1),posvide(i,2))=zi(i);
%     end
% %     affiche_propre('gaussien',A,cpos,spos,color)
%     Agauss=normalisation(A,'normal');
%     
%     %%  sphérique
%     
%     str_model='spherical';
%     [~,~,~,vfit] = variogramfit(v_struct.distance,v_struct.val,1,var(cval),[],'model',str_model);
%     
%     [zi,~] = kriging(vfit,cpos(:,1),cpos(:,2),cval,posvide(:,1),posvide(:,2),18);
%     
%     for i=1:nbpoint-nbc
%         A(posvide(i,1),posvide(i,2))=zi(i);
%     end
% %      affiche_propre('spherique',A,cpos,spos,color)
%     Asph=normalisation(A,'normal');
        
    %% vraie carte
    
    B=vrai_carte(spos,nbs,type,A);
    B=B/mean(mean(B));
%     affiche_propre('la vrai',B,cpos,spos,color)
    
    %% comparaison
%     [D] = diff_derive(A,B)

         Dg=abs(Agauss-B); 
%          affiche_propre('différence',Dg,cpos,spos,color)
         erreur(1,2)=erreur(1,2)+min(Dg(:));
         erreur(2,2)=erreur(2,2)+mean(Dg(:));
         erreur(3,2)=erreur(3,2)+max(Dg(:));
         Dg=abs(Aexp-B);
         erreur(1,1)=erreur(1,1)+min(Dg(:));
         erreur(2,1)=erreur(2,1)+mean(Dg(:));
         erreur(3,1)=erreur(3,1)+max(Dg(:));
         Dg=abs(Asph-B);
         erreur(1,3)=erreur(1,3)+min(Dg(:));
         erreur(2,3)=erreur(2,3)+mean(Dg(:));
         erreur(3,3)=erreur(3,3)+max(Dg(:));


    
    %% corrélation
    
    Cexp=Cexp+mean(mean(corr(Aexp,B)));
    Cgauss=Cgauss+mean(mean(corr(Agauss,B)));
    Csph=Csph+mean(mean(corr(Asph,B)));
    %% fiabilite
    
    seuil=2/3;% on fixe un critère de détection d'une source
    
    for i=1:nbs
        if Aexp(spos(i,1),spos(i,2))>=seuil*max(Aexp(:))
            k(1)=k(1)+1;
        end
        if Agauss(spos(i,1),spos(i,2))>=seuil*max(Agauss(:))
            k(2)=k(2)+1;
        end
        if Asph(spos(i,1),spos(i,2))>=seuil*max(Asph(:))
            k(3)=k(3)+1;
        end
    end
    
end
erreur=erreur/n;
k=k/3/n;
Cexp=Cexp/n;
 Cgauss=Cgauss/n;
Csph=Csph/n;

data=[type '_et_' position];
save(data)