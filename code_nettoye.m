clear all
close all

%% Initialisation

nbc=20;
dimp=[120 60];
pasg=1;
pasmq=5;
A=zeros(dimp(1)/pasg,dimp(2)/pasg);

nbpoint=size(A,1)*size(A,2);
dossier='indsources_32x32.mat';
type='normal';
type_carte='normal';
position='grille';
n=50;
Ahsi=zeros(dimp(1)/pasmq,dimp(2)/pasmq,n);
nbpointhsi=size(Ahsi,1)*size(Ahsi,2);
Agauss=zeros(dimp(1)/pasg,dimp(2)/pasg,n/5-1);
B=zeros(dimp(1)/pasg,dimp(2)/pasg,n);
Ecart=5;
vect=zeros(1,n);
%% création source

nbs=2;
[spos]=creation_position('source_eloigne',nbs,A);
%% création capteur

[cpos,posl,posc]=creation_position(position,nbc,A,dossier);
Ahsi=zeros(size(posl,2),size(posc,2),n);
%%
for kjq=1:n
    %% évaluation valeur au capteur
    
    cval = creation_valeur_3d(cpos,spos,nbc,nbs,type,1.5*kjq);
    
    %% lisibilité de la carte
    
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
    for i=1:nbc
        A(cpos(i,1),cpos(i,2))=cval(i);
    end
    
    Ahsi(:,:,kjq)=reshape(cval,[size(Ahsi,1),size(Ahsi,2)]);
    B(:,:,kjq)=vrai_carte(spos,nbs,type_carte,A,1.5*kjq);
    switch type_carte
        case 'normal'
            vect(kjq)=1.5*kjq;
        case 'diffuse'
            vect(kjq)=1.5*kjq*10;
    end
end

%% affiche beau
% figure,
% for i=1:10
%     subplot(2,5,i)
%     imagesc(Agauss(:,:,i))
%     hold on
%     for i=1:size(cpos,1)
%         plot(cpos(i,2),cpos(i,1),color(i),'Marker','*')
%     end
%     plot(spos(:,2),spos(:,1),'ro')
%     hold off
%     colorbar
% end

%% affiche moche
% figure,
% for i=1:n
%     subplot(5,10,i)
%     imagesc(Ahsi(:,:,i))
%     hold on
%     for i=1:size(cpos,1)
%         plot(cposhsi(i,2),cposhsi(i,1),color(i),'Marker','*')
%     end
%     plot(sposhsi(:,2),sposhsi(:,1),'ro')
%     hold off
%
% end

%% on tente des trucs

PM=zeros(n,size(Agauss,3));
for i=0:9
    PM(5*i+1,i+1)=1;
end
PM=sparse(PM');
Bgauss=squeeze(B(:,:,1:5:n));
P1=zeros(size(Agauss,1),size(posl,2))';
P2=zeros(size(Agauss,2),size(posc,2))';
for i=1:size(posl,2)
    P1(i,posl(end-i+1))=1;
end
for i=1:size(posc,2)
    P2(i,posc(end-i+1))=1;
end
[Ih,Jh,K]=size(Ahsi);
H3=reshape(Ahsi,[Ih*Jh,K]);
t_rank=2; % rank of the tensor -- change for different noise levels
maxit=50; % number of CPD iterations

[A_hat,B_hat,C_hat,A_tilde,B_tilde,C_tilde]=TenRec(Bgauss,H3,maxit,t_rank,P1,P2) ;

% S1_hat1=khatri_rao(C_hat,B_hat)*A_hat';
lamda=0.001;

s_iter=10; % blind stereo iterations -- play between 2-15

[ A1,B1,C1,~ ] = STEREO( Ahsi,Bgauss,P1,P2,PM,s_iter,lamda,A_hat,B_hat,C_hat,C_tilde);

S1_hat1=khatri_rao(C1,B1)*A1';


Res=zeros(size(A,1),size(A,2),n);
[G,H,L]=size(Res);
for i=1:G
    Res(i,:,:)=reshape(S1_hat1(:,i),[H,L]);
end

figure,
for i=1:n
    subplot(5,10,i)
    imagesc(Res(:,:,i))
    hold on
    for i=1:size(cpos,1)
        plot(cpos(i,2),cpos(i,1),color(i),'Marker','*')
    end
    plot(spos(:,2),spos(:,1),'ro')
    hold off
     colorbar
end
%%
figure,
for i=1:n
    subplot(5,10,i)
    imagesc(B(:,:,i))
    hold on
    for i=1:size(cpos,1)
        plot(cpos(i,2),cpos(i,1),color(i),'Marker','*')
    end
    plot(spos(:,2),spos(:,1),'ro')
    hold off
     colorbar
end

%% corrélation

correlation=zeros(1,n);

for i=1:n
    x=Res(:,:,i);
    x=x(:);
    x_tilde=B(:,:,i);
    x_tilde=x_tilde(:);
    correlation(i)=x_tilde'*x/norm(x_tilde)/norm(x);
end
figure,
plot(vect,correlation)

x=Res(:);
x_tilde=B(:);