%function [graph] = createfilteredgraph(segmented_image,filename,maximum,minimum)
function [graph] = createfilteredgraph(segmented_image,filename)

%addpath(genpath('Images'));
load(['binarypatterns\lbp_',filename,'.mat']);
load(['binarypatterns\rilbp_',filename,'.mat']);
load(['binarypatterns\lpq_',filename,'.mat']);
load(['binarypatterns\rilpq_',filename,'.mat']);
load(['binarypatterns\cont_',filename,'.mat']);
histedges = 0:16:256;

%if strcmp(filename,'harbor10') == 1
%    graph = [];return;
%end
image1 = zeros(256);
image1(2:end-1,2:end-1) = lbpimage;

image2 = zeros(256);
image2(2:end-1,2:end-1) = rilbpimage;

image3 = zeros(256);
image3(2:end-1,2:end-1) = lpqimage;

image4 = zeros(256);
image4(3:end-2,3:end-2) = rilpqimage;

image5 = zeros(256);
image5(2:end-1,2:end-1) = contimage;

%image = im2double(imread([filename,'.tif']));
%x = (maximum - minimum)/10;

%number_of_spectra = size(image,3);
%image = im2double(image);
%image = uint8(image);
%image = im2double(image);
%gray_image = rgb2gray(image);
node_field1 = 'lbphist';
node_field2 = 'rilbphist';
node_field3 = 'lpqhist';
node_field4 = 'rilpqhist';
node_field5 = 'conthist';
%node_field19 = 'Centroid'; node_field20 = 'Orientation';

%image_properties = struct('filterhistogram',[],node_field19,[],node_field20,[]);
image_properties = struct(node_field1,[],node_field2,[],node_field3,[],node_field4,[],node_field5,[]);


%% NODE PROCESSING

number_of_labels = length(unique(segmented_image));
new_image = segmented_image;  % initialize new labelled image after region grouping
current_label = 1;
node_properties = [];
features = [];

for label = 1:number_of_labels
    
    disp(['Processing label ',num2str(label),'..']);
    
    binary_image = segmented_image; % initialize new binary image for each label
    binary_image(binary_image~=label) = false; binary_image(binary_image==label) = true;
   
    % Perform morphological operations
    %{
    s = strel('disk',1);    % define structuring element % so that most of the regions have area more than 100 pixels
    temp_image = uint8(imopen(binary_image,s));   % perform opening operation
    morphed_image = uint8(imclose(temp_image,s));   % perform closing operation
    %}
    morphed_image = binary_image;
    
    % extract region properties
    CC = bwconncomp(morphed_image);
    region_properties = regionprops(CC,'Centroid','Orientation');
    % calculate the number of disjoint components for each label
    n = CC.NumObjects;
    
    for idx = 1:n
        
        %image_properties.(node_field19) = region_properties(idx).(node_field19);
        %image_properties.(node_field20) = region_properties(idx).(node_field20);
                    
        % create the new image with separate labels for separate regions
        new_image(CC.PixelIdxList{idx}) = current_label;    
        current_label = current_label + 1;
        
        %{
        %create filter histogram
        mat = image(:,:,1);
        histedges = minimum(1,1) + x(1)*[0:10];
        hist1 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = image(:,:,2);
        histedges = minimum(1,2) + x(2)*[0:10];
        hist2 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = image(:,:,3);
        histedges = minimum(1,3) + x(3)*[0:10];
        hist3 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = log_image1(:,:,1);
        histedges = minimum(1,4) + x(4)*[0:10];
        hist4 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = log_image1(:,:,2);
        histedges = minimum(1,5) + x(5)*[0:10];
        hist5 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = log_image1(:,:,3);
        histedges = minimum(1,6) + x(6)*[0:10];
        hist6 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = log_image2(:,:,1);
        histedges = minimum(1,7) + x(7)*[0:10];
        hist7 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = log_image2(:,:,2);
        histedges = minimum(1,8) + x(8)*[0:10];
        hist8 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = log_image2(:,:,3);
        histedges = minimum(1,9) + x(9)*[0:10];
        hist9 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        
        mat = gaborfilter1(:,:,1);
        histedges = minimum(1,10) + x(10)*[0:10];
        hist10 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter1(:,:,2);
        histedges = minimum(1,11) + x(11)*[0:10];
        hist11 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter1(:,:,3);
        histedges = minimum(1,12) + x(12)*[0:10];
        hist12 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter2(:,:,1);
        histedges = minimum(1,13) + x(13)*[0:10];
        hist13 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter2(:,:,2);
        histedges = minimum(1,14) + x(14)*[0:10];
        hist14 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter2(:,:,3);
        histedges = minimum(1,15) + x(15)*[0:10];
        hist15 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter3(:,:,1);
        histedges = minimum(1,16) + x(16)*[0:10];
        hist16 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter3(:,:,2);
        histedges = minimum(1,17) + x(17)*[0:10];
        hist17 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter3(:,:,3);
        histedges = minimum(1,18) + x(18)*[0:10];
        hist18 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter4(:,:,1);
        histedges = minimum(1,19) + x(19)*[0:10];
        hist19 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter4(:,:,2);
        histedges = minimum(1,20) + x(20)*[0:10];
        hist20 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        mat = gaborfilter4(:,:,3);
        histedges = minimum(1,21) + x(21)*[0:10];
        hist21 = histcounts(mat(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        
        imagehist = [hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9,hist10,hist11,hist12,hist13,hist14,hist15,hist16,hist17,hist18,hist19,hist20,hist21];
        image_properties.filterhistogram = imagehist;
        clear('hist1','hist2','hist3','hist4','hist5','hist6','hist7','hist8','hist9','hist10','hist11','hist12','hist13','hist14','hist15','hist16','hist17','hist18','hist19','hist20','hist21','histedges','mat');
        %}
        
        %create texture histograms
        lbphist = histcounts(image1(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        image_properties.(node_field1) = lbphist;
        rilbphist = histcounts(image2(CC.PixelIdxList{idx}),10,'Normalization','probability');
        image_properties.(node_field2) = rilbphist;
        lpqhist = histcounts(image3(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        image_properties.(node_field3) = lpqhist;
        rilpqhist = histcounts(image4(CC.PixelIdxList{idx}),histedges,'Normalization','probability');
        image_properties.(node_field4) = rilpqhist;
        conthist = histcounts(image5(CC.PixelIdxList{idx}),8,'Normalization','probability');
        image_properties.(node_field5) = conthist;
        
        %create feature vector
        feature_vector = [image_properties.(node_field1), image_properties.(node_field2), image_properties.(node_field3), image_properties.(node_field4), image_properties.(node_field5), region_properties(idx).Centroid, region_properties(idx).Orientation];
        
        % create a structure array with all the region properties
        node_properties = [node_properties;image_properties];
        features = [features;feature_vector];
    end
    
end


%% EDGE PROCESSING

disp('Starting edge processing..');
number_of_nodes = length(node_properties);
edge_info = imRAG(new_image);
%edge_matrix = zeros(number_of_nodes,number_of_nodes);
%{
% create the matrices of region centroids and orientation
centroid = zeros(number_of_nodes,2);
orientation = zeros(number_of_nodes,1);

for i = 1:number_of_nodes
    centroid(i,:) = node_properties(i).Centroid;
    orientation(i) = node_properties(i).Orientation;
end

% calculate the Euclidean distance & relative orientation between the centroids of each pair of regions
for i = 1:size(edge_info,1)
    edge_matrix(edge_info(i,1),edge_info(i,2)) = 0.8*(pdist([centroid(edge_info(i,1)); centroid(edge_info(i,2))])) + 0.2*(abs(orientation(edge_info(i,1)) - orientation(edge_info(i,2))));    
end
%}

edge_matrix1 = zeros(number_of_nodes,number_of_nodes);
edge_matrix2 = zeros(number_of_nodes,number_of_nodes);

% create the matrices of region centroids and orientation
centroid = features(:,end-2:end-1);
orientation = features(:,end);
% calculate the Euclidean distance & relative orientation between the centroids of each pair of regions
for i = 1:size(edge_info,1)
    edge_matrix1(edge_info(i,1),edge_info(i,2)) = pdist2(centroid(edge_info(i,1),:),centroid(edge_info(i,2),:));
    edge_matrix2(edge_info(i,1),edge_info(i,2)) = abs(orientation(edge_info(i,1)) - orientation(edge_info(i,2)));    
end

if number_of_nodes > 1
    edge_matrix1 = edge_matrix1/max(edge_matrix1(:));
    edge_matrix2 = edge_matrix2/max(edge_matrix2(:));
end
edge_matrix = 0.8*edge_matrix1 + 0.2* edge_matrix2;

% Convert the directed graph into undirected graph
edge_matrix = edge_matrix + edge_matrix';

%% FINAL GRAPH CREATION

%node_properties = rmfield(node_properties,node_field19);
%node_properties = rmfield(node_properties,node_field20);

graph = struct('nodes',node_properties,'edges',edge_matrix,'features',features(:,1:end-3),'imagefilename',filename);
%graph = struct('nodes',node_properties,'edges',edge_matrix,'imagefilename',filename);
disp('Graph created..');

%clear('image','segmented_image','new_image','binary_image','morphed_image','original_image','log_image1','log_image2','gaborResult','node_properties','edge_matrix','filename','maximum','minimum','image_properties','region_properties');

end