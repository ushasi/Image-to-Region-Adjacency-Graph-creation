function LPQdesc=ri_lpq(img,LPQfilters,charOri,mode)
% Funtion LPQdesc=ri_lpq(img,winSize,LPQfilters,charOri,mode) computes the rotation invariant Local Phase Quantization (LPQ) 
% descriptor for the input image img. descriptors are calculated using only valid pixels i.e. size(img)-(sqrt(2)*winSize-1).
%
% Inputs: (All empty or undefined inputs will be set to default values)
% img = N*N uint8 or double, format gray scale image to be analyzed.
% LPQfilters = 8*neigborhoodPix*numAngle double, Precomputed LPQ filters (computed using createLPQfilters.m function). 
%                                                If no filters are given they are automatically created with default parameters.
% charOri = N*N double, Characteristic orientation for each pixel (computed using charOrientation.m)
%                       If no orientation are given, they are computed using default parameters.
% mode = 1*n char, defines the desired output type. Possible choices are:
%        (default) 'nh' -> normalized histogram of LPQ codewords (1*256 double vector, for which sum(result)==1)
%                  'h'  -> un-normalized histogram of LPQ codewords (1*256 double vector)
%                  'im' -> LPQ codeword image ([size(img,1)-r,size(img,2)-r] double matrix)
%
% Output:
% LPQdesc = 1*256 double or size(img)-(sqrt(2)*winSize-1) uint8, LPQ descriptors histogram or LPQ code image (see "mode" above)
%
% Example usage:
%
% LPQfilters=createLPQfilters(9);
% img=imread('cameraman.tif');
% charOri=charOrientation(img);
% LPQhist=ri_lpq(img,LPQfilters,charOri);
% figure; bar(LPQhist);
%

% Version published in 2010 by Janne Heikkil? Esa Rahtu, and Ville Ojansivu 
% Machine Vision Group, University of Oulu, Finland

%% Defaut parameters
% LPQ filters
if nargin<2 || isempty(LPQfilters)
    LPQfilters=createLPQfilters; % Creat LPQ filters if they are not given. Use default parameters (see createLPQfilters.m)
end

% Characteristic orientation filters
if nargin<3 || isempty(charOri)
    charOri=charOrientation(img); % Compute characteristic orientations if they are not given. Use default parameters (see charOrientation.m)
end

% Output mode
if nargin<4 || isempty(mode)
    mode='nh'; % return normalized histogram as default
end

%% Check inputs
if size(img,3)~=1
    error('Only gray scale image can be used as input');
end
if sum(strcmp(mode,{'nh','h','im'}))==0
    error('mode must be nh, h, or im. See help for details.');
end

%% Initialize
img=single(img);
[imgRow,imgCol]=size(img);
numAngle=size(LPQfilters,3); % Number of different oriented LPQ filters
winSize=sqrt(size(LPQfilters,2)); % Size of the window (read from pre-computed LPQ filters) This is "enlarged" window size, which fits all rotated filters.
r=ceil((winSize-1)/2); % Get radius from window size

%% Get characteristic orientations for the valid area (i.e. pixels that have full neighborhood)
rt=charOri((r+1):(end-r),(r+1):(end-r));
rt=rt(:)';

%% Reformat image data to matrix containing neighborhoods as columns
F=zeros(winSize^2,(imgRow-winSize+1)*(imgCol-winSize+1));
ii=1;
for i=1:winSize
   for k=1:winSize       
       d=img(i:(end-winSize+i),k:(end-winSize+k));
       F(ii,:)=d(:)';
       ii=ii+1;
   end
end

%% Quantize rotations to the available filter angles (quantized angles must be in interval [0,2*pi])
ang=0:(2*pi/numAngle):(2*pi-2*pi/numAngle);
[temp,rtqind]=min(abs(rt(:)*ones(1,length(ang)+1)-ones(length(rt(:)),1)*[ang,2*pi]),[],2);
rtqind(rtqind==(length(ang)+1))=1;

%% Compute LPQ response using oriented filters 
G=zeros(8,(imgRow-winSize+1)*(imgCol-winSize+1));
for i=1:numAngle
    ii=find(rtqind==i);
    G(:,ii)=LPQfilters(:,:,i)*F(:,ii);
end

%% Quantize values and form LPQ codewords
LPQdesc=(G(1,:)>=0)+(G(2,:)>=0)*2+(G(3,:)>=0)*4+(G(4,:)>=0)*8+(G(5,:)>=0)*16+(G(6,:)>=0)*32+(G(7,:)>=0)*64+(G(8,:)>=0)*128;

%% Resize and switch format to uint8 if LPQ code image is required as output
if strcmp(mode,'im')
    LPQdesc=uint8(LPQdesc);
    LPQdesc=reshape(LPQdesc,[(imgRow-winSize+1),(imgCol-winSize+1)]);
end

%% Histogram if needed
if strcmp(mode,'nh') || strcmp(mode,'h')
    LPQdesc=hist(LPQdesc(:),0:255);
end

%% Normalize histogram if needed
if strcmp(mode,'nh')
    LPQdesc=LPQdesc/sum(LPQdesc);
end


