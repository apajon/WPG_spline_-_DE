format long
clc
clear all
close all

addpath ./input
addpath ./FEM
addpath ./iso2mesh
addpath ./geom3d/geom3d
addpath ./plane_line_intersect
addpath ./gradient

resize = [1200 1600]; %4:3
ImagesFolder='results/CMAME/Video/1-1Neoprene/';
VideoFile=strcat(ImagesFolder,'\1-1Neoprene_simulation_view2');
writerObj = VideoWriter(VideoFile);
writerObj.Quality = 75;
fps = 20; 
writerObj.FrameRate = fps;
open(writerObj);

ImagesFolder='results/CMAME/Video/1-1Neoprene/Simulation/view2';
jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
S = [jpegFiles(:).datenum]; 
[S,S] = sort(S);
jpegFilesS = jpegFiles(S);
for t = 1:length(jpegFilesS)
%for t = 1:10
    t
    Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
    %J = imresize(Frame, resize);
    writeVideo(writerObj,im2frame(Frame));
end
close(writerObj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Create Initial Page for Video                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1 = figure(1); clf;
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 40 30])
% print(fig1,'results/CMAME/Video/InitialSlidePaper','-dpng','-r100')
% textColor    = [0, 0, 0]; % [red, green, blue]
% textLocation = [591 315];       % [x y] coordinates
% textInserter = vision.TextInserter('Shape Optimization of Flexible Soles', 'Color', textColor, 'FontSize', 150, 'Location', textLocation);
% % for Humanoid Robots
% I = imread('results/CMAME/Video/InitialSlidePaper.png');
% J = step(textInserter, I);
% textColor    = [0, 0, 0]; % [red, green, blue]
% textLocation = [1195 515]; % [x y] coordinates
% textInserter = vision.TextInserter('for Humanoid Robots', 'Color', textColor, 'FontSize', 150, 'Location', textLocation);
% J1 = step(textInserter, J);
% 
% textColor    = [0, 0, 0]; % [red, green, blue]
% textLocation = [705 2515]; % [x y] coordinates
% textInserter = vision.TextInserter('G. De Magistris, S. Miossec, A. Escande, A. Kheddar', 'Color', textColor, 'FontSize', 100, 'Location', textLocation);
% J2 = step(textInserter, J1);
% imshow(J2);
% imwrite(J2,'results/CMAME/Video/InitialSlidePaper.png')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   Initial Page for the video 		           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ImagesFolder='results/CMAME/Video/';
% VideoFile=strcat(ImagesFolder,'\CMAME');
% writerObj = VideoWriter(VideoFile);
% writerObj.Quality = 75;
% fps = 20; 
% writerObj.FrameRate = fps;
% open(writerObj);
% for t = 1:60
%     t
%     Frame=imread('results/CMAME/Video/InitialSlidePaper.png');
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J))    
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       Type = 1                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Create Initial Page for Neoprene, Trajectory=1
% ImagesFolder='results/CMAME/Video/';
% for t = 1:60
%     t
%     Frame=imread('results/CMAME/Video/IniSlideNeTr1.png');
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J))    
% end 
% 
% ImagesFolder='results/CMAME/Video/Neoprene/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1Neoprene/3d view/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1Neoprene/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1000Neoprene/3d view/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1000Neoprene/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1NoDir/3d view/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1NoDir/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       Type = 2                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Create Initial Page for ButylRubber, Trajectory=1
% ImagesFolder='results/CMAME/Video/';
% for t = 1:60
%     Frame=imread('results/CMAME/Video/IniSlideBRTr1.png');
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J))    
% end 
% 
% ImagesFolder='results/CMAME/Video/ButylRubber/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1ButylRubber/3d view/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1ButylRubber/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       Type = 3                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Create Initial Page for Neoprene, Trajectory=2
% ImagesFolder='results/CMAME/Video/';
% for t = 1:60
%     Frame=imread('results/CMAME/Video/IniSlideNeTr2.png');
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J))    
% end 
% 
% ImagesFolder='results/CMAME/Video/Straight/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1Straight/3d view/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% 
% ImagesFolder='results/CMAME/Video/1-1Straight/Simulation/';
% jpegFiles = dir(strcat(ImagesFolder,'\*.png'));
% S = [jpegFiles(:).datenum]; 
% [S,S] = sort(S);
% jpegFilesS = jpegFiles(S);
% for t = 1:length(jpegFilesS)
% %for t = 1:10
%     t
%     Frame=imread(strcat(ImagesFolder,'\',jpegFilesS(t).name));
%     J = imresize(Frame, resize);
%     writeVideo(writerObj,im2frame(J));
% end
% close(writerObj);