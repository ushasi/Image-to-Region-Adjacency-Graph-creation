%function filterimage()
clc;
workspace;
% Include the folder containing the image database
addpath(genpath('Images'));

filenames = cell(21,1);     % names of the categories in image database

filenames{1} = 'agricultural';
filenames{2} = 'airplane';
filenames{3} = 'baseballdiamond';
filenames{4} = 'beach';
filenames{5} = 'buildings';
filenames{6} = 'chaparral';
filenames{7} = 'denseresidential';
filenames{8} = 'forest';
filenames{9} = 'freeway';
filenames{10} = 'golfcourse';
filenames{11} = 'harbor';
filenames{12} = 'intersection';
filenames{13} = 'mediumresidential';
filenames{14} = 'mobilehomepark';
filenames{15} = 'overpass';
filenames{16} = 'parkinglot';
filenames{17} = 'river';
filenames{18} = 'runway';
filenames{19} = 'sparseresidential';
filenames{20} = 'storagetanks';
filenames{21} = 'tenniscourt';


minimum = 999*ones(1,21);    % to store minimum values of the filtered images
maximum = -999*ones(1,21);    % to store maximum values of the filtered images

size_of_filter = 3;    % filter size n (for n*n filter)

h1 = fspecial('log',size_of_filter,0.2);
h2 = fspecial('log',size_of_filter,1);
gaborArray = gaborFilterBank(4,8,size_of_filter,size_of_filter);


for j = 1:21
    for i = 0:99
        %disp('**START**');
        % create the filenames required
        if i<10
            fullfilename = strcat(filenames{j},'0',num2str(i));
        else
            fullfilename = strcat(filenames{j},num2str(i));
        end
        
        original_image = imread([fullfilename,'.tif']);
        original_image = imresize(original_image, [256, 256]);

        saveimagefilename2 = strcat('Filtered Images\image2_',fullfilename,'.mat');
        saveimagefilename3 = strcat('Filtered Images\image3_',fullfilename,'.mat');
        saveimagefilename4 = strcat('Filtered Images\image4_',fullfilename,'.mat');
        saveimagefilename5 = strcat('Filtered Images\image5_',fullfilename,'.mat');
        saveimagefilename6 = strcat('Filtered Images\image6_',fullfilename,'.mat');
        saveimagefilename7 = strcat('Filtered Images\image7_',fullfilename,'.mat');
        
        %disp('Saving file..');
        %gray_image = rgb2gray(im2double(original_image));
        original_image = im2double(original_image);
        sz = size(original_image);
        log_image1 = zeros(sz);
        log_image2 = zeros(sz);
        gaborfilter1 = zeros(sz);
        gaborfilter2 = zeros(sz);
        gaborfilter3 = zeros(sz);
        gaborfilter4 = zeros(sz);
        
        minimum(1:3) = 0;
        maximum(1:3) = 1;
        
        % Filtering of image by Laplacian of Gaussian with sigma 0.2
        log_image1(:,:,1) = imfilter(original_image(:,:,1),h1);
        log_image1(:,:,2)= imfilter(original_image(:,:,2),h1);
        log_image1(:,:,3) = imfilter(original_image(:,:,3),h1);
        minimum(1,4) = min(minimum(1,4),min(min(log_image1(:,:,1))));
        maximum(1,4) = max(maximum(1,4),max(max(log_image1(:,:,1))));
        minimum(1,5) = min(minimum(1,5),min(min(log_image1(:,:,2))));
        maximum(1,5) = max(maximum(1,5),max(max(log_image1(:,:,2))));
        minimum(1,6) = min(minimum(1,6),min(min(log_image1(:,:,3))));
        maximum(1,6) = max(maximum(1,6),max(max(log_image1(:,:,3))));
        
        save(saveimagefilename2,'log_image1');
        
        % Filtering of image by Laplacian of Gaussian with sigma 1
        log_image2(:,:,1) = imfilter(original_image(:,:,1),h2);
        log_image2(:,:,2)= imfilter(original_image(:,:,2),h2);
        log_image2(:,:,3) = imfilter(original_image(:,:,3),h2);
        minimum(1,7) = min(minimum(1,7),min(min(log_image2(:,:,1))));
        maximum(1,7) = max(maximum(1,7),max(max(log_image2(:,:,1))));
        minimum(1,8) = min(minimum(1,8),min(min(log_image2(:,:,2))));
        maximum(1,8) = max(maximum(1,8),max(max(log_image2(:,:,2))));
        minimum(1,9) = min(minimum(1,9),min(min(log_image2(:,:,3))));
        maximum(1,9) = max(maximum(1,9),max(max(log_image2(:,:,3))));
        
        save(saveimagefilename3,'log_image2');
        
        % Filtering of image by Gabor filters with 4 orientations
        gaborfilter1(:,:,1) = abs(imfilter(original_image(:,:,1),gaborArray{1,1}));
        gaborfilter1(:,:,2) = abs(imfilter(original_image(:,:,2),gaborArray{1,1}));
        gaborfilter1(:,:,3) = abs(imfilter(original_image(:,:,3),gaborArray{1,1}));
        gaborfilter2(:,:,1) = abs(imfilter(original_image(:,:,1),gaborArray{1,2}));
        gaborfilter2(:,:,2) = abs(imfilter(original_image(:,:,2),gaborArray{1,2}));
        gaborfilter2(:,:,3) = abs(imfilter(original_image(:,:,3),gaborArray{1,2}));
        gaborfilter3(:,:,1) = abs(imfilter(original_image(:,:,1),gaborArray{1,3}));
        gaborfilter3(:,:,2) = abs(imfilter(original_image(:,:,2),gaborArray{1,3}));
        gaborfilter3(:,:,3) = abs(imfilter(original_image(:,:,3),gaborArray{1,3}));
        gaborfilter4(:,:,1) = abs(imfilter(original_image(:,:,1),gaborArray{1,4}));
        gaborfilter4(:,:,2) = abs(imfilter(original_image(:,:,2),gaborArray{1,4}));
        gaborfilter4(:,:,3) = abs(imfilter(original_image(:,:,3),gaborArray{1,4}));
        
        save(saveimagefilename4,'gaborfilter1');
        minimum(1,10) = min(minimum(1,10),min(min(gaborfilter1(:,:,1))));
        maximum(1,10) = max(maximum(1,10),max(max(gaborfilter1(:,:,1))));
        minimum(1,11) = min(minimum(1,11),min(min(gaborfilter1(:,:,2))));
        maximum(1,11) = max(maximum(1,11),max(max(gaborfilter1(:,:,2))));
        minimum(1,12) = min(minimum(1,12),min(min(gaborfilter1(:,:,3))));
        maximum(1,12) = max(maximum(1,12),max(max(gaborfilter1(:,:,3))));
        
        save(saveimagefilename5,'gaborfilter2');
        minimum(1,13) = min(minimum(1,13),min(min(gaborfilter2(:,:,1))));
        maximum(1,13) = max(maximum(1,13),max(max(gaborfilter2(:,:,1))));
        minimum(1,14) = min(minimum(1,14),min(min(gaborfilter2(:,:,2))));
        maximum(1,14) = max(maximum(1,14),max(max(gaborfilter2(:,:,2))));
        minimum(1,15) = min(minimum(1,15),min(min(gaborfilter2(:,:,3))));
        maximum(1,15) = max(maximum(1,15),max(max(gaborfilter2(:,:,3))));
        
        save(saveimagefilename6,'gaborfilter3');
        minimum(1,16) = min(minimum(1,16),min(min(gaborfilter3(:,:,1))));
        maximum(1,16) = max(maximum(1,16),max(max(gaborfilter3(:,:,1))));
        minimum(1,17) = min(minimum(1,17),min(min(gaborfilter3(:,:,2))));
        maximum(1,17) = max(maximum(1,17),max(max(gaborfilter3(:,:,2))));
        minimum(1,18) = min(minimum(1,18),min(min(gaborfilter3(:,:,3))));
        maximum(1,18) = max(maximum(1,18),max(max(gaborfilter3(:,:,3))));
        
        save(saveimagefilename7,'gaborfilter4');
        minimum(1,19) = min(minimum(1,19),min(min(gaborfilter4(:,:,1))));
        maximum(1,19) = max(maximum(1,19),max(max(gaborfilter4(:,:,1))));
        minimum(1,20) = min(minimum(1,20),min(min(gaborfilter4(:,:,2))));
        maximum(1,20) = max(maximum(1,20),max(max(gaborfilter4(:,:,2))));
        minimum(1,21) = min(minimum(1,21),min(min(gaborfilter4(:,:,3))));
        maximum(1,21) = max(maximum(1,21),max(max(gaborfilter4(:,:,3))));
        %{
        % Filtering of image by averaging filter
        
        avg_image(:,:,1) = imfilter(original_image(:,:,1),h3);
        avg_image(:,:,2)= imfilter(original_image(:,:,2),h3);
        avg_image(:,:,3) = imfilter(original_image(:,:,3),h3);
        minimum(1,22) = min(minimum(1,22),min(min(avg_image(:,:,1))));
        maximum(1,22) = max(maximum(1,22),max(max(avg_image(:,:,1))));
        minimum(1,23) = min(minimum(1,23),min(min(avg_image(:,:,2))));
        maximum(1,23) = max(maximum(1,23),max(max(avg_image(:,:,2))));
        minimum(1,24) = min(minimum(1,24),min(min(avg_image(:,:,3))));
        maximum(1,24) = max(maximum(1,24),max(max(avg_image(:,:,3))));
        
        save(saveimagefilename8,'avg_image');
        
        % Filtering of image by sobel operator (horizontal)
        
        grad_image1(:,:,1) = imfilter(original_image(:,:,1),h4);
        grad_image1(:,:,2)= imfilter(original_image(:,:,2),h4);
        grad_image1(:,:,3) = imfilter(original_image(:,:,3),h4);
        minimum(1,22) = min(minimum(1,22),min(min(grad_image1(:,:,1))));
        maximum(1,22) = max(maximum(1,22),max(max(grad_image1(:,:,1))));
        minimum(1,23) = min(minimum(1,23),min(min(grad_image1(:,:,2))));
        maximum(1,23) = max(maximum(1,23),max(max(grad_image1(:,:,2))));
        minimum(1,24) = min(minimum(1,24),min(min(grad_image1(:,:,3))));
        maximum(1,24) = max(maximum(1,24),max(max(grad_image1(:,:,3))));

        save(saveimagefilename9,'grad_image1');
        
        % Filtering of image by sobel operator (vertical)
        
        grad_image2(:,:,1) = imfilter(original_image(:,:,1),h4');
        grad_image2(:,:,2)= imfilter(original_image(:,:,2),h4');
        grad_image2(:,:,3) = imfilter(original_image(:,:,3),h4');
        minimum(1,25) = min(minimum(1,25),min(min(grad_image2(:,:,1))));
        maximum(1,25) = max(maximum(1,25),max(max(grad_image2(:,:,1))));
        minimum(1,26) = min(minimum(1,26),min(min(grad_image2(:,:,2))));
        maximum(1,26) = max(maximum(1,26),max(max(grad_image2(:,:,2))));
        minimum(1,27) = min(minimum(1,27),min(min(grad_image2(:,:,3))));
        maximum(1,27) = max(maximum(1,27),max(max(grad_image2(:,:,3))));

        save(saveimagefilename10,'grad_image2');
        %}
    end
    j
end

% save the variables maximum and minimum
save('minimum.mat','minimum');
save('maximum.mat','maximum');

%end

LPQfilters=createLPQfilters(3);
mapping=getmapping(8,'riu2');

for j = 1:21
    %filename = filenames{j};

    for i = 0:99
        
       % disp('**START**');
        % create the filenames required
        if i<10
            fullfilename = strcat(filenames{j},'0',num2str(i));
        else
            fullfilename = strcat(filenames{j},num2str(i));
        end
        
        original_image = imread([fullfilename,'.tif']);
        image = rgb2gray(original_image);
        %if isequal(size(image),[256,256]) == 0
            
            image = imresize(image,[256 256]);
            saveimagefilename1 = strcat('binarypatterns\lbp_',fullfilename,'.mat');
            saveimagefilename2 = strcat('binarypatterns\cont_',fullfilename,'.mat');
            saveimagefilename3 = strcat('binarypatterns\lpq_',fullfilename,'.mat');
            saveimagefilename4 = strcat('binarypatterns\rilpq_',fullfilename,'.mat');
            saveimagefilename5 = strcat('binarypatterns\rilbp_',fullfilename,'.mat');
            
            charOri=charOrientation(image);
            
            lbpimage = lbp(image,1,8,0,'i');
            save(saveimagefilename1,'lbpimage');
            
            rilbpimage = lbp(image,1,8,mapping,'i');
            save(saveimagefilename5,'rilbpimage');
            
            
            contimage = cont(image);
            save(saveimagefilename2,'contimage');
            
            lpqimage = lpq(image,3,1,1,'im');
            save(saveimagefilename3,'lpqimage');
            
            rilpqimage = ri_lpq(image,LPQfilters,charOri,'im');
            save(saveimagefilename4,'rilpqimage');
        %end
        
    end
    j
end


load('minimum.mat','minimum');
load('maximum.mat','maximum');

%totgra = [];
%totgra = cell(2100,1);
index = 1;
for j = 1:21  %21
    %filename = filenames{j};

    for i = 0:99 %
        
        if i<10
            %fullfilename = strcat('Images\',filenames{j},'\',filenames{j},'0',num2str(i));
            fullfilename = strcat(filenames{j},'0',num2str(i));
        else
            %fullfilename = strcat('Images\',filenames{j},'\',filenames{j},num2str(i));
             fullfilename = strcat(filenames{j},num2str(i));
        end
        
        img = imread([fullfilename,'.tif']);
        img = imresize(img, [256, 256]);
        
         disp('**START**');

        graph  = extractfeaturevec(img,fullfilename,maximum,minimum,index);
        totgra(index) = graph;
        
        
        index = index + 1;
        
    end
end

%seg = superpixels(image,5);
%bw = boundarymask(seg);
%imshow(imoverlay(image,bw,'red'),'InitialMagnification',67)
