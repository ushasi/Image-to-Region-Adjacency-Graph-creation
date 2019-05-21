function segmented_image = segmentation(image,alpha,max_num_of_regions)

 
%addpath('Segmentation');
k = max_num_of_regions;

sz = size(image);
Hc=ones(sz(1:2)); 
Vc=Hc;

%i_ground = 0; % rank of the background for plotting, 0: the darkest; 
%k-1 the brightest; 99: nowhere

%diff=10000;
%an_energy=999999999;
iter=0;
iter_v=0;
energy_global_min=99999999;
distance = 'sqEuclidean'; % Feature space distance

% Initialization: cluster the data into k regions
%tic,
disp('Starting segmentation..');

%image = imread('agricultural00.tif');
%sz = size(image);

% im contains the input RGB image as a SINGLE array
regionSize = 20 ;
regularizer = 50 ;
image = im2single(image);

segments = vl_slic(image, regionSize, regularizer) ;

%data = ToVector(image);
% % data = reshape(image, 1, []);
% % data = data';
% % [~, c] = kmeans(data, k, 'distance', distance,'EmptyAction','drop');
% % 
% % c=c(find(isfinite(c(:,1))),:);
% % k=size(c,1);
% % k_max=k;
%kmean_time = toc,

%Dc = zeros([sz(1:2) k],'single');   

%tic
%{
while iter < 5
    iter=iter+1;
    if k == 1
        continue;
    end
    Dc = zeros(sz(1),sz(2),k);
    clear K
    
    for ci=1:k
        K=kernel_RBF(image,c(ci,:));
        Dc(:,:,ci)=1-K;
    end   
    clear Sc
    clear K
    
    %% The smoothness term
    
    Sc = alpha*(ones(k) - eye(k)); 
    gch = GraphCut('open', Dc, Sc, Vc, Hc);
    [gch, L] = GraphCut('swap',gch);
    [gch, se, de] = GraphCut('energy', gch);
    nv_energy = de + se;
    gch = GraphCut('close', gch);
 
    if (nv_energy<=energy_global_min)
        %diff = abs(energy_global_min-nv_energy)/energy_global_min;
        energy_global_min=nv_energy;
        L_global_min=L;
        k_max=k;
        
        %nv_energy;
        iter_v=0;
        % Calculate region Pl of label l
        
        if size(image, 3)==3 % Color image
        for l=0:k-1
            Pl=find(L==l);
            card=length(Pl);
            K1=kernel_RBF(image(Pl),c(l+1,1));K2=kernel_RBF(image(Pl),c(l+1,2));K3=kernel_RBF(image(Pl),c(l+1,3));
            smKI(1)=sum(image(Pl).*K1); smKI(2)=sum(image(Pl+prod(sz(1:2))).*K2); smKI(3)=sum(image(Pl+2*prod(sz(1:2))).*K3);
            smK1=sum(K1);smK2=sum(K2);smK3=sum(K3);
            if (card~=0)
                c(l+1,1)=smKI(1)/smK1;c(l+1,2)=smKI(2)/smK2;c(l+1,3)=smKI(3)/smK3;
            else
                c(l+1,1)=999999999;c(l+1,2)=999999999;c(l+1,3)=999999999;
            end
        end
        end
        
        if size(image, 1)==1 % Gray-level image
        for l=0:k-1
            Pl=find(L==l);
            card=length(Pl);
            K=kernel_RBF(image(Pl),c(l+1,1));
            smKI=sum(image(Pl).*K);
            smK=sum(K);
            if (card~=0)
                c(l+1,1)=smKI/smK;
            else
                c(l+1,1)=999999999;
            end
        end
        end
        
        
        c=c(find(c(:,1)~=999999999),:);
        c_global_min=c;
        %k_global=
        k=length(c(:,1));

    else
        iter_v=iter_v+1;
        %---------------------------------
        %       Begin updating labels
        %---------------------------------
        % Calculate region Pl of label l
        
        if size(image, 3)==3 % Color image 
        for l=0:k-1           
            Pl=find(L==l);
            card=length(Pl);
            K1=kernel_RBF(image(Pl),c(l+1,1));K2=kernel_RBF(image(Pl),c(l+1,2));K3=kernel_RBF(image(Pl),c(l+1,3));
            smKI(1)=sum(image(Pl).*K1); smKI(2)=sum(image(Pl+prod(sz(1:2))).*K2); smKI(3)=sum(image(Pl+2*prod(sz(1:2))).*K3);
            smK1=sum(K1);smK2=sum(K2);smK3=sum(K3);
            % Calculate contour Cl of region Pl
            if (card~=0)
                c(l+1,1)=smKI(1)/smK1;c(l+1,2)=smKI(2)/smK2;c(l+1,3)=smKI(3)/smK3;
            else
                c(l+1,1)=999999999;c(l+1,2)=999999999;c(l+1,3)=999999999;
                area(l+1)=999999999;
            end
        end
        end
        
        if size(image, 3)== 1 % Gray-level image 
        for l=0:k-1           
            Pl=find(L==l);
            card=length(Pl);
            K=kernel_RBF(image(Pl),c(l+1,1));
            smKI=sum(image(Pl).*K);
            smK=sum(K);
            % Calculate contour Cl of region Pl
            if (card~=0)
                c(l+1,1)=smKI/smK;
            else
                c(l+1,1)=999999999;
                area(l+1)=999999999;
            end
        end
        end
              
        c=c(find(c(:,1)~=999999999),:);
        k=length(c(:,1));
end
   
end
 %}
%L=L_global_min;
L_global_min = 0;
L = 0;
%energy_global_min;
c= min(min(segments));
%c=c_global_min;

%number_of_regions = size(c,1)
%iter;

%% Show the results

[n, k]=size(segments);
 k_max=k;
if size(image, 3)==3 % Color image 
img=zeros(sz(1),sz(2),3);
j=1;
segmented_image = segments;
%imagesc(segments); axis off; hold on; 


%{
for i=0:k_max-1
    LL=(L_global_min==i);
    is_zero=sum(sum(LL));
    if is_zero
         img(:,:,1)=img(:,:,1)+LL*c(j,1);
         img(:,:,2)=img(:,:,2)+LL*c(j,2);
         img(:,:,3)=img(:,:,3)+LL*c(j,3);
         j=j+1;
    end
    
end
end

if size(image, 3)==1 % Gray-level image 
img=zeros(sz(1),sz(2));
j=1;
%imagesc(im); axis off; hold on; colormap gray; 

for i=0:k_max-1
    LL=(L_global_min==i);
    is_zero=sum(sum(LL));
    if is_zero
         img(:,:,1)=img(:,:,1)+LL*c(j,1);
         j=j+1;
    end
    
end
end

segmented_image = rgb2gray(img);
uniq_label = unique(segmented_image(:));
for i = 1:length(uniq_label)
    segmented_image((segmented_image == uniq_label(i))) = i;
end
%}

end