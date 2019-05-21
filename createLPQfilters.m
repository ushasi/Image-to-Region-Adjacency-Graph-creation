function LPQfilters = createLPQfilters(winSize,numAngl)
% function LPQfilters = createLPQfilters(winSize,numAngl) precomputes the LPQ filters to be applied with 
% rotation invariant LPQ descriptors (ri_LPQ.m).
%
% Inputs: (All empty or undefined inputs will be set to default values)
% winSize = 1*1 double, size of the local window. winSize must be odd number and greater or equal to 3 (default winSize=9).
% numAngl = 1*1 double, defines the number of angles for which LPQ filters are generated. 
%                       numAngl must be even number and greater or equal to 2. (default value numAngl=36).
%
% Outputs:
% LPQfilters = 8*neigborhoodPixNum*numAngl double, precomputed LPQ filters for numAngl evenly distributed orientation.
%
% Example usage:
% 
% LPQfilters=createLPQfilters(9,36);
%

% Version published in 2010 by Janne Heikkil? Esa Rahtu, and Ville Ojansivu 
% Machine Vision Group, University of Oulu, Finland


%% Default parameters
% Local window size
if nargin<1 || isempty(winSize)
    winSize=9; % default window size 9
end

% Number of angle samples
if nargin<2 || isempty(numAngl)
    numAngl=36; % default number of angle samples 36 
end

% Local frequency estimation (Frequency points used [alpha,0], [0,alpha], [alpha,alpha], and [alpha,-alpha]) 
STFTalpha=1/winSize;  % alpha in STFT approaches 

%% Check inputs
if winSize<3 || rem(winSize,2)~=1
   error('Window size winSize must be odd number and greater than equal to 3');
end
if numAngl<2 || rem(numAngl,2)~=0
   error('Number of angle samples (numAngl) must be even number and greater or equal to 2.');
end

%% Initialize
r=(winSize-1)/2; % Get radius from window size
x=-r:r; % Form spatial coordinates
rw=ceil(sqrt(2)*r); % Enlarge window size to fit rotated filters
winSizeEnl=2*rw+1; % Size of enlarged window
[xbb,ybb]=meshgrid(-rw:rw,-rw:rw); % 2-D coordinates for the enlarged window
ang=0:(2*pi/numAngl):(2*pi-2*pi/numAngl); % Initialize rotation angles. Rotation is counter clock wise in image.
LPQfilters=zeros(8,winSizeEnl^2,numAngl); % Initialize transform matrices. 
xi=STFTalpha*[1 1 1 0; -1 0 1 1]; % Frequency points used [alpha,0], [0,alpha], [alpha,alpha], and [alpha,-alpha]

%% Form filter matrices for each rotation 
for i=1:numAngl
    % Rotation matrix
    R=[cos(ang(i)), -sin(ang(i)); sin(ang(i)), cos(ang(i))]; 
    
    % Rotated 2-D coordinates
    temp=R*[xbb(:)';ybb(:)'];
    xbbi=reshape(temp(1,:),winSizeEnl,winSizeEnl);
    ybbi=reshape(temp(2,:),winSizeEnl,winSizeEnl);
    
    % Form filter for each frequency point and rotated it to the current angle (ang(i))
    cnt=1; % Initialize counter
    for k=1:size(xi,2)
        M=zeros(winSizeEnl); % Initialize filter to zero
        M(((winSizeEnl-winSize)/2+1):(end-(winSizeEnl-winSize)/2),((winSizeEnl-winSize)/2+1):(end-(winSizeEnl-winSize)/2))= exp(-2*pi*1i*xi(1,k)*x(:)) * exp(-2*pi*1i*xi(2,k)*(x(:).'));
        M=M.'; 
        Mi=interp2(xbb,ybb,M,xbbi,ybbi,'*linear',0); % Interpolate filter to new orientation
        LPQfilters(cnt,:,i)=real(Mi(:).'); % Store real part of the rotated filter
        LPQfilters(cnt+1,:,i)=imag(Mi(:).'); % Store complex part of the rotated filter
        cnt=cnt+2; % Increase counter
    end

end



