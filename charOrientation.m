function charOri=charOrientation(img,winSize,numAngl,useApprox)
% charOri=charOrientation(img,winSize,numAngl,useApprox) computes the characteristic orientation 
% for the input image img. Orientations are calculated only for valid pixels i.e. size(img)-(winSize-1). Borders default ot zero.
%
% Inputs: (All empty or undefined inputs will be set to default values)
% img = N*N uint8 or double, format gray scale image to be analyzed.
% winSize = 1*1 double, size of the local window. winSize must be odd number and greater or equal to 3 (default winSize=5).
% numAngl = 1*1 double, defines the number of angle samples used in orientation estimation. 
%                       numAngl must be even number and greater or equal to 2. (Has no effect if approximation is used).
% useApprox = 1*1 double, indicates whether approximation scheme is used or not. Possible values are:
%              (default)  0 -> no approximation, 
%                         1 -> use approximation
%
% Output:
% charOri = size(img) double matrix, estimated characteristic orientations. Border values (outside valid area) are set to zero.
%
% Example usage:
% img=imread('cameraman.tif');
% charOri=estimate_rotation(img);
%

% Version published in 2010 by Janne Heikkil? Esa Rahtu, and Ville Ojansivu 
% Machine Vision Group, University of Oulu, Finland


%% Defaut parameters
% Local window size
if nargin<2 || isempty(winSize)
    winSize=5; % default window size 5
end

% Number of angle samples
if nargin<3 || isempty(numAngl)
    numAngl=36; % default number of angle samples 36 (has no effect if approximation is used)
end

% Approximation scheme (approximating response using cosine)
if nargin<4 || isempty(useApprox)
    useApprox=0; % do not use approximation as default
end

%% Check inputs
if size(img,3)~=1
    error('Only gray scale image can be used as input');
end
if numAngl<2 || rem(numAngl,2)~=0
   error('Number of angle samples (numAngl) must be even number and greater or equal to 2.');
end
if winSize<3 || rem(winSize,2)~=1
   error('Window size winSize must be odd number and greater than equal to 3');
end
if sum(useApprox==[0 1])==0
    error('useApprox parameter must be set to 0->do not use approximation or 1->use approximation. See help for details.');
end

%% Initialize
img=single(img);
[imgRow,imgCol]=size(img);
r=(winSize-1)/2; % Get radius from window size

%% Get rotation estimation filters
Wi=getRotationEstimationMasks_(winSize,numAngl,useApprox);

%% Compute estimation angles (if approximation is used, only 0 and pi/2 angles are applied).
if useApprox == 1
    ang=[0;pi/2];
else
    numAngl=numAngl/2; % We only need angles up to pi because of the symmetry. 
    ang=(0:(pi/numAngl):(pi-pi/numAngl))'; % Imaginary part is antisymmetric.
end

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

%% Run characteristic orientation estimation using either approximation or sampling method.
if useApprox==1 % Approximation using cosine fitting
    % Rotation estimation
    R=(Wi*F); % imag frequency components at freq 1/n at angles ang
    tmp1=(R(2,:)*cos(ang(1))-R(1,:)*cos(ang(2)));
    tmp2=(R(2,:)*sin(ang(1))-R(1,:)*sin(ang(2)));
    
    % Get phase angle
    phi=atan(tmp1./tmp2);

    % Get characteristic orientation as peak location (equivalent to phase angle of the complex moment as below).
    charOri=mod(-phi,2*pi); 

else % Characteristic orientation using sampling method
    % Get signs of the imaginary part for each angle sample for each valid pixel position
    R=(Wi*F)>=0;

    % Calculate complex moment
    fullang=exp(1i*[ang; ang+pi]);
    tmp=sum(([R; 1-R]).*(fullang*ones(1,size(R,2))));
    
    % Get characteristic orientation as phase angle of the complex moment.
    charOri=mod(atan2(imag(tmp),real(tmp)),2*pi);

end

%% Reformating estimates
charOri=reshape(charOri,[imgRow-2*r,imgCol-2*r]);
    
%% Zero padd to mach the image size (image borders result zero orientation, since no well defined neigborhood is available)
charOri=[zeros(r,size(charOri,2)+2*r);[zeros(size(charOri,1),r),charOri,zeros(size(charOri,1),r)];zeros(r,size(charOri,2)+2*r)];





%%%%%%%%%%%%%%%%%%%%%%%%
% Additional functions %
%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute filter masks for rotation estimation
function Wi = getRotationEstimationMasks_(winSize,numAngl,useApprox)

% Default parameters
stdev = 2; % Standard deviation of the Gaussian window
%sigmaS=(winSize-1)/4; % Tuottaa eri tuloksen winSize=5 -> stdev=1;

% Initialize
r=(winSize-1)/2; % Get radius from window size
x=-r:r; % Form coordinate values

% Compute estimation angles (if approximation is used, only 0 and pi/2 angles are applied).
if useApprox == 1
    ang=[0;pi/2];
else
    numAngl=numAngl/2; % We only need angles up to pi because of the symmetry. 
    ang=(0:(pi/numAngl):(pi-pi/numAngl))'; % Imaginary part is antisymmetric.
end

% Form Gaussian window
gs=exp(-(x.^2)/(2*stdev^2));
H=gs'*gs; H=H/sum(H(:));

% Compute filter masks
xi=(1/winSize)*[cos(ang),-sin(ang)]';
Wi=zeros(size(xi,2),winSize^2);
for i=1:size(xi,2)
    M=exp(-2*pi*1i*xi(1,i)*x(:)) * exp(-2*pi*1i*xi(2,i)*(x(:).'));
    M=M.';
    M=H.*M; % Gaussian window
    Wi(i,:)=imag(M(:).');
end





