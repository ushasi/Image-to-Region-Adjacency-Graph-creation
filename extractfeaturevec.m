function [ graph ] = extractfeaturevec(img,fullfilename,maximum,minimum,index)
workspace;


%fullfilename = 'agricultural00';
%img =  imread([fullfilename,'.tif']);
% x = randi(10,1); featurevec = rand(5,x);
disp(['Processing ',fullfilename,'..']);
hsv_image = rgb2hsv(img); lab_image = rgb2lab(img); image = im2double(img);

% load the filtered images
load(['Filtered Images\image2_',fullfilename,'.mat']);
load(['Filtered Images\image3_',fullfilename,'.mat']);
load(['Filtered Images\image4_',fullfilename,'.mat']);
load(['Filtered Images\image5_',fullfilename,'.mat']);
load(['Filtered Images\image6_',fullfilename,'.mat']);
load(['Filtered Images\image7_',fullfilename,'.mat']);

mini = load('minimum.mat');
minimum = mini.minimum;
maxi = load('maximum.mat');
maximum = maxi.maximum;

%minimum(1:21) = 0.3;
%maximum(1:21) = 1;
x = (maximum - minimum)/8;

% Preprocessing for LBP & LPQ features
load(['binarypatterns\lbp_',fullfilename,'.mat']);
load(['binarypatterns\rilpq_',fullfilename,'.mat']);
load(['binarypatterns\cont_',fullfilename,'.mat']);
histedges = 0:32:256;
if size(image,1)>256 || size(image,2)>256
   sz = size(image,1);
else
   sz = 256;
end
lbpimage = lbp(rgb2gray(img),'i');
image1 = zeros(256);
image1(2:255,2:255) = lbpimage;

image4 = zeros(sz);
image4(3:254,3:254) = rilpqimage;

image5 = zeros(sz);
image5(2:255,2:255) = contimage;

% Preprocessing for Tamura features
% [SbestR,thetaR] = tamura_features(img(:,:,1));
% [SbestG,thetaG] = tamura_features(img(:,:,2));
% [SbestB,thetaB] = tamura_features(img(:,:,3));
% histedges = -0.5:1:255.5; graylevels = 0:255;
% phai=0:0.001:pi;
% THRESHOLD=0.01;


%% Define the attributes of the ARG

% Define the node attributes
% GENERAL FEATURES
node_field1 = 'Area';
node_field2 = 'ConvexArea';
node_field3 = 'Perimeter';
node_field4 = 'Extent';
node_field5 = 'Solidity';
node_field6 = 'Eccentricity';
node_field7 = 'Entropy';
node_field8 = 'FilledArea';
node_field9 = 'colorHistogram';
node_field10 = 'colorMoments';
node_field11 = 'fourierdescriptors';
node_field12 = 'spectralhistogram';
node_field13 = 'LBPfeatures';
node_field14 = 'LPQfeatures';
node_field15 = 'densesiftfeatures';
node_field16 = 'EquivDiameter';
node_field17 = 'EulerNumber';
node_field18 = 'MajorAxisLength';
node_field19 = 'MinorAxisLength';
node_field20 = 'Orientation';
node_field21 = 'BoundingBox';
node_field22 = 'Centroid';

%node_field1 = 'Area';
%node_field2 = 'ConvexArea';
%node_field3 = 'FilledArea';
%node_field3 = 'EquivDiameter';
%node_field4 = 'Perimeter';
%node_field5 = 'Extent';
%node_field6 = 'Solidity';
%node_field8 = 'EulerNumber';
%node_field7 = 'Eccentricity';
%node_field8 = 'Entropy';

% COLOR FEATURES
%node_field9 = 'colorHistogram';
%node_field10 = 'colorMoments';
%node_field13 = 'colorAutoCorrelogram';

% SHAPE FEATURE
%node_field11 = 'fourierdescriptors';

% TEXTURE FEATURES
%node_field12 = 'spectralhistogram';
%node_field13 = 'LBPfeatures';
%node_field14 = 'LPQfeatures';
%node_field15 = 'densesiftfeatures';

image_properties = struct(node_field1,[],node_field2,[],node_field3,[],node_field4,[],node_field5,[],node_field6,[],node_field7,[],node_field8,[],node_field9,[],node_field10,[],node_field11,[],node_field12,[],node_field13,[],node_field14,[],node_field15,[],node_field16,[],node_field17,[],node_field18,[],node_field19,[],node_field20,[],node_field21,[]);

% Define the edge attributes
% edge_field1 = 'Centroid-Distance'; edge_field2 = 'Orientation';

max_num_of_regions=300;
alpha = 0.1;
imshow(image);
segmented_image = segmentation(image,alpha,max_num_of_regions);
%% NODE PROCESSING

labels = unique(segmented_image)';
node_properties = [];
features = [];
weights = [];
%s = strel('disk',2);    % define structuring element so that most of the regions have area more than 50 pixels

for label = labels
    
    disp(['Processing label ',num2str(label),'..']);    
    binary_image = segmented_image; % initialize new binary image for each label
    binary_image(binary_image~=label) = false; binary_image(binary_image==label) = true;
           
    % extract region properties
    CC = bwconncomp(binary_image);
    boundary = bwboundaries(binary_image,'noholes');
    region_properties = regionprops(CC,node_field5,node_field6,node_field6,node_field16,node_field17,node_field18,node_field19,node_field20,node_field21);
    region_properties1 = regionprops(CC,'Centroid','Orientation',node_field1,node_field2,node_field3,node_field4,node_field5,node_field6,node_field8);

    
    % calculate the number of disjoint components for each label
    % n = CC.NumObjects
        
    for idx = 1:CC.NumObjects
        
%         if (region_properties(idx).(node_field1)<250)
%             break
%         end
        image_properties.(node_field1) = region_properties1(idx).(node_field1);
        image_properties.(node_field2) = region_properties1(idx).(node_field2);
        image_properties.(node_field3) = region_properties1(idx).(node_field3);
        image_properties.(node_field4) = region_properties1(idx).(node_field4);
        image_properties.(node_field5) = region_properties1(idx).(node_field5);
        image_properties.(node_field6) = region_properties1(idx).(node_field6);
        image_properties.(node_field8) = region_properties1(idx).(node_field8);
        image_properties.(node_field16) = region_properties(idx).(node_field16);
        image_properties.(node_field17) = region_properties(idx).(node_field17);
        image_properties.(node_field18) = region_properties(idx).(node_field18);
        image_properties.(node_field19) = region_properties(idx).(node_field19);
        image_properties.(node_field20) = region_properties(idx).(node_field20);
        image_properties.(node_field21) = region_properties(idx).(node_field21);
        %image_properties.(node_field7) = region_properties(idx).(node_field7);
                      
        % get the intensity values of each region
        len = length(CC.PixelIdxList{idx});
        cropped_region_hsv = zeros(len,1,3);
        cropped_region_lab = zeros(len,1,3);
        cropped_region_rgb = zeros(len,1,3);
        %cropped_region = zeros(len,128);
        [rowsub,colsub] = ind2sub(size(image(:,:,1)),CC.PixelIdxList{idx});
        
        imgd = single(rgb2gray(img));
        cropped_region = zeros(1,300);
        cropped_region1 = zeros(1,300);
        [siftImage, descriptors] = vl_sift(imgd);
        [n,k] = size(siftImage);
        %for i = 1:n
            cropped_region1 = reshape(siftImage,1,[]);
            for i = 1:300
                if (n >= 300)
                cropped_region(1,i) = cropped_region1(i);
                end
            end
            
            %cropped_region = reshape(siftImage,[1,300]);
            %cropped_region = detectSURFFeatures(imgd);
        %end
        
        for i = 1:len
            cropped_region_hsv(i,1,:) = hsv_image(rowsub(i),colsub(i),:);
            cropped_region_rgb(i,1,:) = image(rowsub(i),colsub(i),:);
            cropped_region_lab(i,1,:) = lab_image(rowsub(i),colsub(i),:);
            %[x,y] = centroid(rowsub(i),colsub(i),:)
            %[siftImage, descriptors] = vl_sift(rowsub(i),colsub(i),:);
            
        end
        
        % calculate entropy of regions
        entropy_vec = [entropy(cropped_region_rgb(:,:,1)),entropy(cropped_region_rgb(:,:,2)),entropy(cropped_region_rgb(:,:,3))];
        image_properties.(node_field7) = entropy_vec;
        
        % add the color histogram attribute
        image_properties.(node_field9) = hsvHistogram(cropped_region_hsv);
           
        % add the color moments (spectral) attribute
        R = cropped_region_lab(:,:,1);
        G = cropped_region_lab(:,:,2);
        B = cropped_region_lab(:,:,3);
        image_properties.(node_field10) = [mean(R(:)),std(R(:)),skewness(R(:)),mean(G(:)),std(G(:)),skewness(G(:)),mean(B(:)),std(B(:)),skewness(B(:))];
        
         
        % add the shape feature
        b = boundary{idx};
        t = (b(:, 1) - mean(b(:,1))) + 1i*(b(:, 2) - mean(b(:, 2)));   % Convert coordinates to complex numbers
        z = fft(abs(t));    % Compute the descriptors
        % Normalize the descriptors
        if abs(z(1)) ~= 0
            z = z/abs(z(1));
        end
        fourier_descriptors = z;
        abs_fourier_descriptors = abs(fourier_descriptors);
        angle_fourier_descriptors = angle(fourier_descriptors);
        image_properties.(node_field11) = [mean(abs_fourier_descriptors), std(abs_fourier_descriptors), mean(angle_fourier_descriptors), std(angle_fourier_descriptors)];
               
        % add spectral histogram feature
        imagehist = spectralhistogram(image,log_image1,log_image2,gaborfilter1,gaborfilter2,gaborfilter3,gaborfilter4,minimum,x,CC.PixelIdxList{idx});
        %imagehist(9:24) = image_properties.(node_field9);
        %imagehist = imagehist(9:end);
        image_properties.(node_field12) = imagehist;
        
        % add LBP and LPQ features
        lbphist = histcounts(image1(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        rilpqhist = histcounts(image4(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        conthist = histcounts(image5(CC.PixelIdxList{idx}),8,'Normalization','probability');
        image_properties.(node_field13) = lbphist;
        image_properties.(node_field14) = [rilpqhist,conthist];
        
        %{ 
        % add the color autocorrelogram feature
        %distance_vector = 1;
        %image_properties.(node_field13) = (color_auto_correlogram(cropped_region_rgb,distance_vector))';      
        
        % add the Tamura features        
        [FcoarsenessR, FcontrastR, FdirectionR] = tamura_calc(SbestR, thetaR, CC.PixelIdxList{idx}, img(:,:,1), len, histedges, graylevels, phai, THRESHOLD);
        [FcoarsenessG, FcontrastG, FdirectionG] = tamura_calc(SbestG, thetaG, CC.PixelIdxList{idx}, img(:,:,2), len, histedges, graylevels, phai, THRESHOLD);
        [FcoarsenessB, FcontrastB, FdirectionB] = tamura_calc(SbestB, thetaB, CC.PixelIdxList{idx}, img(:,:,3), len, histedges, graylevels, phai, THRESHOLD);              
        image_properties.(node_field15) = [FcoarsenessR, FcoarsenessG, FcoarsenessB, FcontrastR, FcontrastG, FcontrastB, FdirectionR, FdirectionG, FdirectionB];
           
        % add the wavelet attributes
        coeff1 = dwt(double(cropped_region_rgb(:,:,1)),'coif1');
        coeff2R = dwt(coeff1,'coif1');
        coeff1 = dwt(double(cropped_region_rgb(:,:,2)),'coif1');
        coeff2G = dwt(coeff1,'coif1');
        coeff1 = dwt(double(cropped_region_rgb(:,:,2)),'coif1');
        coeff2B = dwt(coeff1,'coif1');
        image_properties.(node_field16) = [mean(coeff2R(:)), std(coeff2R(:)), mean(coeff2G(:)), std(coeff2G(:)), mean(coeff2B(:)), std(coeff2B(:))];
        %}
        %feature_vector = [image_properties.(node_field1), image_properties.(node_field2), image_properties.(node_field3), image_properties.(node_field4),image_properties.(node_field5), image_properties.(node_field6), image_properties.(node_field7),image_properties.(node_field8),image_properties.(node_field9), image_properties.(node_field10), image_properties.(node_field11), image_properties.(node_field12), image_properties.(node_field13), image_properties.(node_field14),image_properties.(node_field15), image_properties.(node_field16),image_properties.(node_field17), image_properties.(node_field18),image_properties.(node_field19), image_properties.(node_field20),image_properties.(node_field21)];
        %features = [features;feature_vector];
        % add D-SIFT descriptors
        
        %centroid = features(:,end-2:end-1);
        %orientation = features(:,end);
        %indices = knnsearch(centroid,cropped_region);
        %siftvec = hist(indices, 1:size(centroids,1))/len;
        image_properties.(node_field15) = cropped_region;
        
        
        %create feature vector
        feature_vector = [image_properties.(node_field1), image_properties.(node_field2), image_properties.(node_field3), image_properties.(node_field4),image_properties.(node_field5), image_properties.(node_field6), image_properties.(node_field7),image_properties.(node_field8),image_properties.(node_field9), image_properties.(node_field10), image_properties.(node_field11), image_properties.(node_field12), image_properties.(node_field13), image_properties.(node_field14),image_properties.(node_field15), image_properties.(node_field16),image_properties.(node_field17), image_properties.(node_field18),image_properties.(node_field19), image_properties.(node_field20),image_properties.(node_field21)];
        features = [features;feature_vector];
        weights = [weights;image_properties.(node_field1)];
        
        % create a structure array with all the region properties
         node_properties = [node_properties;image_properties];                                   
    end
    
end


%% EDGE PROCESSING

disp('Starting edge processing..');
number_of_nodes = length(node_properties);
edge_matrix1 = zeros(number_of_nodes,number_of_nodes);      %containing centroid values
edge_matrix2 = zeros(number_of_nodes,number_of_nodes);      %containing orientation values

% Create the Region Adjacency graph (RAG)
edge_info = imRAG(segmented_image);

% create the matrices of region centroids and orientation
centroid = features(:,end-2:end-1);
orientation = features(:,end);

% calculate the Euclidean distance & relative orientation between the centroids of each pair of regions
for i = 1:size(edge_info,1)
    edge_matrix1(edge_info(i,1),edge_info(i,2)) = pdist2(centroid(edge_info(i,1),:),centroid(edge_info(i,2),:));
    edge_matrix2(edge_info(i,1),edge_info(i,2)) = abs(orientation(edge_info(i,1)) - orientation(edge_info(i,2)));    
end

% Normalize the values
if number_of_nodes > 1
    edge_matrix1 = edge_matrix1/max(edge_matrix1(:));
    edge_matrix2 = edge_matrix2/max(edge_matrix2(:));
end

% Create final adjacency matrix with 80% weightage to centroid distance value and 20% weightage to orientation difference
edge_matrix = 0.8*edge_matrix1 + 0.2* edge_matrix2;

% Convert the directed graph into undirected graph
edge_matrix = edge_matrix + edge_matrix';



%% FINAL GRAPH CREATION

%graph = struct('features',features,'weights',weights,'imagefilename',filename);
graph = struct('nodes',node_properties,'edges',edge_matrix,'features',features(:,1:end-3),'imagefilename',fullfilename,'index',index);

disp('Graph created..');

end