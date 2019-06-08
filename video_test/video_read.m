%% video_read.m


%% record camera
clear;clc;close all;
camList = webcamlist
cam = webcam(1)
preview(cam);
img = snapshot(cam);

% Display the frame in a figure window.
figure(111), image(img);

for idx = 1:5
    figure(idx+1), 
    img = snapshot(cam);
    image(img);
end
clear cam

%% record video
% Connect to the webcam. 
cam = webcam 
vidWriter = VideoWriter('frames.avi'); 
open(vidWriter); 
% for index = 1:20 
for index = 1:50
    % Acquire frame for processing 
    img = snapshot(cam); 
    
    % Write frame to video 
    writeVideo(vidWriter, img); 
end
close(vidWriter); 
clear cam