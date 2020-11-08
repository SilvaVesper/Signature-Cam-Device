 clear; clc; 
% Testing image processing
%%
Orig_Sign = imread('Test Signature.png');        % Reading the image and putting it into MatLab
%%% NOTE: The image has to be 192 x 96 pixels to be able to scale it. 
%%%       The letters also have to be relatively connected/ no vast breaks
%%%       in curvature. (I.e. having a T as two separate sticks)


Sig_Bin = imbinarize(Orig_Sign);                   % Binarizing the file (no grayscale, only black and white)
Sig_Inv = imcomplement(Sig_Bin);                   % Color inverting the binarization (black -> white and white -> black)

[B,L,N,A] = bwboundaries(Sig_Inv);
      % B = Number of traces needed ; L = Label matrix
      % N = number of objects found ; A = Adjacency matrix
      

%% Citation: MathWorks bwboundaries 
%%% https://www.mathworks.com/help/images/ref/bwboundaries.html

imshow(Sig_Inv);                        % Open the inverted plot 
hold on;                                % Hold on to plot the trace plots as well
colors = ['b' 'g' 'r' 'c' 'm' 'y'];     % Color options (blue, green, red, cyan, magenta, yellow)

for k = 1:length(B);                    % Plotting B number of traces

  boundary = B{k};                      % Store in a cell array  
  cidx = mod(k,length(colors))+1;       % Modulo to change the color
  plot(boundary(:,2), boundary(:,1),... % Plot these things on top of the NS_Inv plot
       colors(cidx),'LineWidth',2);      
end


%% X and Y trace calculations

% Initializing the x values
Ax = 0;
Bx = 0;
Cx = 0;
Dx = 0;
Ex = 0;
Fx = 0;
Gx = 0;
Hx = 0;
Jx = 0;
Kx = 0;
Lx = 0;
Mx = 0;

% Initializing the y values
Ay = 0;
By = 0;
Cy = 0;
Dy = 0;
Ey = 0;
Fy = 0;
Gy = 0;
Hy = 0;
Jy = 0;
Ky = 0;
Ly = 0;
My = 0;

n = 1;   % Initializing the counter
for n = 1:length(B)
% For however many traces the signature needs (compensates for up to 12 traces)
thcam = 100;
     if n == 1
         AA = B{1,1};  % First trace    
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Ax = AA(:,2);
         Ay = -AA(:,1);
         Ath = linspace(1,thcam,length(Ax));
     elseif n == 2
         BB = B{2,1};  % Second trace
         % Breaking those individual trace arrays into x and y components 
         % and storing it in a corresponding array
         Bx = BB(:,2);
         By = -BB(:,1);
         Bth = linspace(1,thcam,length(Bx));
     elseif n == 3
         CC = B{3,1};  % Third trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Cx = CC(:,2);
         Cy = -CC(:,1);
         Cth = linspace(1,thcam,length(Cx));
     elseif n == 4
         DD = B{4,1};  % Fourth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Dx = DD(:,2);
         Dy = -DD(:,1);
         Dth = linspace(1,thcam,length(Dx));
     elseif n == 5
         EE = B{5,1};  % Fifth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Ex = EE(:,2);
         Ey = -EE(:,1);
         Eth = linspace(1,thcam,length(Ex)); 
     elseif n == 6
         FF = B{6,1};  % Sixth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Fx = FF(:,2);
         Fy = -FF(:,1);
         Fth = linspace(1,thcam,length(Fx));
     elseif n == 7
         GG = B{7,1};  % Seventh trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Gx = GG(:,2);
         Gy = -GG(:,1);
         Gth = linspace(1,thcam,length(Gx));
     elseif n == 8
         HH = B{8,1};  % Eighth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Hx = HH(:,2);
         Hy = -HH(:,1);
         Hth = linspace(1,thcam,length(Hx));
     elseif n == 9
         JJ = B{9,1};  % Ninth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Jx = JJ(:,2);
         Jy = -JJ(:,1);
         Jth = linspace(1,thcam,length(Jx));
     elseif n == 10
         KK = B{10,1}; % Tenth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Kx = KK(:,2);
         Ky = -KK(:,1);
         Kth = linspace(1,thcam,length(Kx));
     elseif n == 11
         LL = B{11,1}; % Eleventh trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Lx = LL(:,2);
         Ly = -LL(:,1);
         Lth = linspace(1,thcam,length(Lx));
     elseif n == 12
         MM = B{12,1}; % Twelfth trace
         % Breaking those individual trace arrays into x and y components
         % and storing it in a corresponding array
         Mx = MM(:,2);
         My = -MM(:,1);
         Mth = linspace(1,thcam,length(Mx));
     else
         % Do nothing (Do you really need more than twelve???)
     end
end


%% Cam shape calculations
% NOTE: Only need the first 3 traces as those are the outer shape of the letters

PA = thcam/(length(Ax)-1);  % Cam division (in this case 100 degrees) divided
PB = thcam/(length(Bx)-1);  % by the length of the corresponding array - 1
PC = thcam/(length(Cx)-1);  % (will be used to calculate the slope between points)

%%% Converting the A, B, and C values from original array size to 100 data points
xA = 1;                     % Initializer
for xA = 1:101              % For 100 degrees of the cam (1 extra data point for later)
    Adiv = length(Ax)/101;  % Divions needed to convert Ax into 100 data pts
    Axn = nearest(Adiv*xA); 
    % Since arrays only deal with integers, finds the nearest integer to
    % the calculated decimal
    if Axn == 0
        Axnew(xA) = Ax(1);  % Prevents the first point from rounding to 0
    else
        Axnew(xA) = Ax(Axn);
    % New Ax array (total of 100 data points) and setting that n space
    % equivalent to the original Ax array at that nearest array spot
    end
    xA = xA + 1;            % Increase counter
end

yA = 1;                     % Initializer
for yA = 1:101              % For 100 degrees of the cam (1 extra data point for later)
    Adiv = length(Ay)/101;  % Divions needed to convert Ay into 100 data pts
    Ayn = nearest(Adiv*yA);
    % Since arrays only deal with integers, finds the nearest integer to
    % the calculated decimal
    if Ayn == 0
        Aynew(yA) = Ay(1);  % Prevents the first point from rounding to 0
    else
        Aynew(yA) = Ay(Ayn);
    % New Ay array (total of 100 data points) and setting that n space
    % equivalent to the original Ay array at that nearest array spot
    end
    yA = yA + 1;            % Increase counter
end

xB = 1;                     % Initializer
for xB = 1:101              % For 100 degrees of the cam (1 extra data point for later)
    Bdiv = length(Bx)/101;  % Divions needed to convert Bx into 100 data pts
    Bxn = nearest(Bdiv*xB);
    % Since arrays only deal with integers, finds the nearest integer to
    % the calculated decimal
    if Bxn == 0
        Bxnew(xB) = Bx(1);  % Prevents the first point from rounding to 0
    else
        Bxnew(xB) = Bx(Bxn);
    % New Bx array (total of 100 data points) and setting that n space
    % equivalent to the original Bx array at that nearest array spot
    end
    xB = xB + 1;            % Increase counter
end

yB = 1;                     % Initializer
for yB = 1:101              % For 100 degrees of the cam (1 extra data point for later)
    Bdiv = length(By)/101;  % Divions needed to convert By into 100 data pts
    Byn = nearest(Bdiv*yB);
    % Since arrays only deal with integers, finds the nearest integer to
    % the calculated decimal
    if Byn == 0
        Bynew(yB) = By(1);  % Prevents the first point from rounding to 0
    else
        Bynew(yB) = By(Byn);
    % New By array (total of 100 data points) and setting that n space
    % equivalent to the original By array at that nearest array spot
    end
    yB = yB + 1;
end

xC = 1;                     % Initializer
for xC = 1:101              % For 100 degrees of the cam (1 extra data point for later)
    Cdiv = length(Cx)/101;  % Divions needed to convert Cx into 100 data pts
    Cxn = nearest(Cdiv*xC);
    % Since arrays only deal with integers, finds the nearest integer to
    % the calculated decimal
    if Cxn == 0
        Cxnew(xC) = Cx(1);  % Prevents the first point from rounding to 0
    else
        Cxnew(xC) = Cx(Cxn);
    % New Cx array (total of 100 data points) and setting that n space
    % equivalent to the original Cx array at that nearest array spot
    end
    xC = xC + 1;
end

yC = 1;                     % Initializer
for yC = 1:101              % For 100 degrees of the cam (1 extra data point for later)
    Cdiv = length(Cy)/101;  % Divions needed to convert Cy into 100 data pts
    Cyn = nearest(Cdiv*yC);
    % Since arrays only deal with integers, finds the nearest integer to
    % the calculated decimal
    if Cyn == 0
        Cynew(yC) = Cy(1);  % Prevents the first point from rounding to 0
    else
        Cynew(yC) = Cy(Cyn);
    % New Cy array (total of 100 data points) and setting that n space
    % equivalent to the original Cy array at that nearest array spot
    end
    yC = yC + 1;
end


%% Calculating the Rx and Ry cam/radius 

BRad = 192; % Base radius of the cam in pixels (Cam diameter of 4")

DeltRAx = zeros(100,1); % Creating array to put the delta (change in) radius values in the X direction
DeltRAy = zeros(100,1); % Creating array to put the delta (change in) radius values in the Y direction
for n = 1:100
    DeltAx = (Axnew(n+1)- Axnew(n));  % Change in x values
    DeltAy = (Aynew(n+1)- Aynew(n));  % Change in y values
    if n == 1
        DeltRAx(n) = (DeltAx ./ PA) + BRad; % Slope in the x direction (increasing or decreasing) from the base radius 
        DeltRAy(n) = (DeltAy ./ PA) + BRad; % Slope in the y direction (increasing or decreasing) from the base radius
    else
        DeltRAx(n) = (DeltAx ./ PA) + DeltRAx(n-1); % Slope in the x direction (increasing or decreasing) from the previous radius
        DeltRAy(n) = (DeltAy ./ PA) + DeltRAy(n-1); % Slope in the y direction (increasing or decreasing) from the previous radius
    end
end

BRadx = DeltRAx(100) + 150; 
    % Alternate base radius (adding 150 pixels) to shift over the cam so the pen does not trace
    %   the second letter on top of the first. 
DeltRBx = zeros(100,1); % Creating array to put the delta (change in) radius values in the X direction
DeltRBy = zeros(100,1); % Creating array to put the delta (change in) radius values in the Y direction
for n = 1:100
    DeltBx = (Bxnew(n+1)- Bxnew(n));  % Change in x values
    DeltBy = (Bynew(n+1)- Bynew(n));  % Change in y values
    if n == 1
        DeltRBx(n) = (DeltBx ./ PB) + BRadx; % Slope in the x direction (increasing or decreasing) from the alternate base radius 
        DeltRBy(n) = (DeltBy ./ PB) + BRad;  % Slope in the y direction (increasing or decreasing) from the base radius
    else
        DeltRBx(n) = (DeltBx ./ PB) + DeltRBx(n-1); % Slope in the x direction (increasing or decreasing) from the previous radius
        DeltRBy(n) = (DeltBy ./ PB) + DeltRBy(n-1); % Slope in the y direction (increasing or decreasing) from the previous radius
    end
end

BRadx = DeltRBx(100) + 150;
    % Alternate base radius (adding 150 pixels) to shift over the cam so the pen does not trace
    %   the third letter on top of the second. 
DeltRCx = zeros(100,1); % Creating array to put the delta (change in) radius values in the X direction
DeltRCy = zeros(100,1); % Creating array to put the delta (change in) radius values in the Y direction
for n = 1:99
    DeltCx = (Cxnew(n+1)- Cxnew(n));  % Change in x values
    DeltCy = (Cynew(n+1)- Cynew(n));  % Change in y values
    if n == 1
        DeltRCx(n) = (DeltCx ./ PC) + BRadx; % Slope in the x direction (increasing or decreasing) from the alternate base radius
        DeltRCy(n) = (DeltCy ./ PC) + BRad;  % Slope in the y direction (increasing or decreasing) from the base radius
    else
        DeltRCx(n) = (DeltCx ./ PC) + DeltRCx(n-1); % Slope in the x direction (increasing or decreasing) from the previous radius
        DeltRCy(n) = (DeltCy ./ PC) + DeltRCy(n-1); % Slope in the y direction (increasing or decreasing) from the previous radius
    end
end

DeltRCx(100) = DeltRAx(1); 
DeltRCy(100) = DeltRAy(1); 
% Setting the last value of C equal to the first value of A (in the x and 
%   y direction) such that the final value of the cam (at 360*) does not  
%   spike to zero, but rather closes up the cam.

%%%% TRANSITION PTS %%%%
CTOnex = linspace(DeltRAx(100),DeltRBx(1),30); % Creating an even spacing over 30 degrees 
CTOney = linspace(DeltRAy(100),DeltRBy(1),30); %   from the end of the first letter to the
CTTwox = linspace(DeltRBx(100),DeltRCx(1),30); %   beginning of the second letter (and the
CTTwoy = linspace(DeltRBy(100),DeltRCy(1),30); %   end of the second to the start of the third)


%% Combining the pieces from above to have the cam radii in terms of theta from 0 to 360 degrees

Rx = zeros(1,360);     % Array to hold the X cam values
Rx(1:100)   = DeltRAx; % First letter 
Rx(101:130) = CTOnex;  % First transition
Rx(131:230) = DeltRBx; % Second letter
Rx(231:260) = CTTwox;  % Second transition
Rx(261:360) = DeltRCx; % Third letter

Ry = zeros(1,360);     % Array to hold the Y cam values
Ry(1:100)   = DeltRAy; % First letter 
Ry(101:130) = CTOney;  % First transition
Ry(131:230) = DeltRBy; % Second letter
Ry(231:260) = CTTwoy;  % Second transition
Ry(261:360) = DeltRCy; % Third letter

Rz = ones(1,360);      % Array to hold the Z cam values
Rz = Rz .* BRad;       % Multiplying the ones array by the base radius value
Rz(101:130) = BRad + BRad/25;   % First transition
Rz(231:260) = BRad + BRad/25;   % Second transition


%% Plots
% NOTE: If the guidelines for the image stated above were followed, the 
%       first 3 traces will be the 3 letters for the initials

figure
subplot(2,2,1), imshow(Orig_Sign); title('Original .png file')
subplot(2,2,2), imshow(Sig_Bin); title('Binarized file')
subplot(2,2,3), imshow(Sig_Inv); title('Inverted Binarization')

figure
subplot(2,2,1), plot(Ax,Ay); title('Trace 1')
subplot(2,2,2), plot(Bx,By); title('Trace 2')
subplot(2,2,3), plot(Cx,Cy); title('Trace 3')
sgtitle('First three plots')

figure
subplot(2,3,1), plot(Ath,Ax); title('Trace 1 -> Ax wrt 100 degrees')
subplot(2,3,2), plot(Bth,Bx); title('Trace 2 -> Bx wrt 100 degrees')
subplot(2,3,3), plot(Cth,Cx); title('Trace 3 -> Cx wrt 100 degrees')
subplot(2,3,4), plot(Ath,Ay); title('Trace 1 -> Ay wrt 100 degrees')
subplot(2,3,5), plot(Bth,By); title('Trace 2 -> By wrt 100 degrees')
subplot(2,3,6), plot(Cth,Cy); title('Trace 3 -> Cy wrt 100 degrees')
sgtitle('X and Y versus theta trace plots')

figure
polarplot(Rx)
title('Rx cam - Base diameter of 4"')

figure
polarplot(Ry)
title('Ry cam - Base diameter of 4"')

figure
polarplot(Rz)
title('Rz cam - Base diameter of 4"')


%%% Extra Plots %%%
% If desired, these can be uncommented to look at more traces on the
%   curvature of the initials. 
%
%   WARNING: The program will throw an error if you uncomment more plots
%   than traces exist. i.e. Signature only has 8 traces, and all are
%   uncommented, the traces 9 - 12 have no data and therefore are unable to
%   be plotted and thus shatters the program. 

% figure
% subplot(3,3,1), plot(Dx,Dy); title('Trace 4')
% subplot(3,3,2), plot(Ex,Ey); title('Trace 5')
% subplot(3,3,3), plot(Fx,Fy); title('Trace 6')
% subplot(3,3,4), plot(Gx,Gy); title('Trace 7')
% subplot(3,3,5), plot(Hx,Hy); title('Trace 8')
% subplot(3,3,6), plot(Jx,Jy); title('Trace 9')
% subplot(3,3,7), plot(Kx,Ky); title('Trace 10')
% subplot(3,3,8), plot(Lx,Ly); title('Trace 11')
% subplot(3,3,9), plot(Mx,My); title('Trace 12')
% sgtitle('Traces 4 - 12')

% figure
% subplot(3,3,1), plot(Dth,Dx); title('Trace 4 -> Dx wrt 100 degrees')
% subplot(3,3,2), plot(Eth,Ex); title('Trace 5 -> Ex wrt 360 degrees')
% subplot(3,3,3), plot(Fth,Fx); title('Trace 6 -> Fx wrt 360 degrees')
% subplot(3,3,4), plot(Gth,Gx); title('Trace 7 -> Gx wrt 360 degrees')
% subplot(3,3,5), plot(Hth,Hx); title('Trace 8 -> Hx wrt 360 degrees')
% subplot(3,3,6), plot(Jth,Jx); title('Trace 9 -> Jx wrt 360 degrees')
% subplot(3,3,7), plot(Kth,Kx); title('Trace 10 -> Kx wrt 360 degrees')
% subplot(3,3,8), plot(Lth,Lx); title('Trace 11 -> Lx wrt 360 degrees')
% subplot(3,3,9), plot(Mth,Mx); title('Trace 12 -> Mx wrt 360 degrees')
% sgtitle('Traces 4 - 12, X vals versus theta')

% figure
% subplot(3,3,1), plot(Dth,Dy); title('Trace 4 -> Dy wrt 100 degrees')
% subplot(3,3,2), plot(Eth,Ey); title('Trace 5 -> Ey wrt 360 degrees')
% subplot(3,3,3), plot(Fth,Fy); title('Trace 6 -> Fy wrt 360 degrees')
% subplot(3,3,4), plot(Gth,Gy); title('Trace 7 -> Gy wrt 360 degrees')
% subplot(3,3,5), plot(Hth,Hy); title('Trace 8 -> Hy wrt 360 degrees')
% subplot(3,3,6), plot(Jth,Jy); title('Trace 9 -> Jy wrt 360 degrees')
% subplot(3,3,7), plot(Kth,Ky); title('Trace 10 -> Ky wrt 360 degrees')
% subplot(3,3,8), plot(Lth,Ly); title('Trace 11 -> Ly wrt 360 degrees')
% subplot(3,3,9), plot(Mth,My); title('Trace 12 -> My wrt 360 degrees')
% sgtitle('Traces 4 - 12, Y vals versus theta')
