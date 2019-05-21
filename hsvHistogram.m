function  hsvhist  =  hsvHistogram(cropped_region_hsv)
%filenames = dir(fullfile(image_folder, '*.TIF'));  % read all images with specified extention, its jpg in our case

%total_images = numel(filenames);    % count total number of photos present in that folder

 %for n = 1:total_images
 % full_name= fullfile(image_folder, filenames(n).name);         % it will specify images names with full path and extension
   % our_images = imread(full_name);                 % Read images  
    %figure (n)                           % used tat index n so old figures are not over written by new new figures
    %imshow(our_images)                  % Show all images
 

     
    image= cropped_region_hsv;
    Red = image(:,:,1);
    Green = image(:,:,2);
    Blue = image(:,:,3);
    
    %Get histValues for each channel
    yRed = imhist(Red,32);
    yGreen = imhist(Green,32);
    yBlue = imhist(Blue,32);
    hsvhist(:) = [yRed' yGreen'  yBlue'];
    %Plot them together in one plot
    %plot(x, yRed, 'Red', x, yGreen, 'Green', x, yBlue, 'Blue');
    
end