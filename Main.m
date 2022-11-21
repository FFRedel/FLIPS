%% FLIPS
%  These files are alowed to be adjusted. However, without permission of
%  the authors, it is not allowed to publish or distrubute these files.

clc, clear all
close all 
restoredefaultpath
addpath 'Fista package'
addpath 'inputs'
%% Initialize and import 

% Import image
im.original    = imread('cameraman.tif');
patchsize      = 8;                                                             
rescale_min    = 0;                                                    
rescale_max    = 1;                                                        
im.original    = rescale( im.original,rescale_min,rescale_max);           

% Selecting patch dimensions
bigimdim1    = 200;                                                     
bigimdim2    = 200;     
var_noise    = 0.0055;
start        = 1;                       % Starting pixel
endd         = start + bigimdim1 - 1;   % Final pixel 

% Resize image and adding noise 
im.original              = imresize(im.original,[bigimdim1 bigimdim2]);                   
im.original              = im.original(start:endd,start:endd);
[imdim1, imdim2]         = size(im.original);
rng(1)                                                             
mean_noiselevel          = 0;                                                    
im.noise                 = imnoise(im.original,'gaussian',mean_noiselevel,var_noise); 
indices_store_patches_inp= reshape(1:imdim1*imdim2,[imdim1 imdim2]);        

% Extract patches
X.noise      = double(im2col(im.noise,   [patchsize patchsize],'Sliding'));            
X.no_noise   = double(im2col(im.original,[patchsize patchsize],'Sliding'));              
store_inp    = double(im2col(indices_store_patches_inp, [patchsize patchsize],'Sliding'));                         
                                           
% Dimensions 
N  = size(X.no_noise,2);                                               
n  = size(X.no_noise,1);                                            
k  = n;                                                                 

%% Selecting phi(f)

% DCT Dictionary
D = DCT(n);

% Measurement matrix 
C = eye(n);                                                               

% Linear mapping
phi = C*D;                                                                                       

%% Calculating appriopriate epsilon

epsilon = sqrt(var_noise)*sqrt(n);   
%% Initial least squares solution

H0_LL      = phi\X.noise;
H0_LL_norm = H0_LL;
for i=1:N
    H0_LL_norm(:,i) = H0_LL_norm(:,i)/norm(H0_LL_norm(:,i),1);
end

%% Chambolle-Pock

L_CP       = svd(phi);
L_CP       = L_CP(1);
tau        = 1/L_CP;
sigma      = 1/L_CP;
maxiter_CP = 50;
theta      = 0.8;

[F_CP,store_CP] = ChambollePock(X.noise,phi,epsilon,N,k,H0_LL,maxiter_CP,tau,sigma,theta);

%% C-SALSA 

maxiter_csalsa=50;

% Selecting tuned values for mu (tuned for cameraman image)
if patchsize == 64 || patchsize == 128 || patchsize == 150
    mu = 3;  
elseif patchsize == 8 || patchsize == 16 || patchsize==32
    mu = 2.5;
end

[F_csalsa,store_csalsa] = C_SALSA(X.noise,phi,epsilon,N,k,H0_LL,maxiter_csalsa,mu);

%% FLIPS (Quadratic Oracle)
% For this part, the function ProjectionOntoL1Ball can be downloaded from Author: John Duchi (jduchi@cs.berkeley.edu) https://stanford.edu/~jduchi/projects/DuchiShSiCh08.html

% Select appropriate smoothness parameter for current patch size 
% (tuned values for cameraman image)
if patchsize ==8
    betainv = 0.002;  
elseif patchsize == 16
    betainv = 0.0016; 
elseif patchsize ==32
    betainv = 0.00028; 
elseif patchsize == 64
    betainv = 0.00001; 
elseif patchsize == 128
    betainv = 0.000001; 
end
maxiter_FLIPS  = 50;                                                          
TT             = 1;       
manual_FLIPS   ='off';
theta          = 0.7;

[F_FLIPS,store_FLIPS] = FLIPS_Solver(X.noise,phi,epsilon,N,k,H0_LL_norm,maxiter_FLIPS,TT,betainv,theta,manual_FLIPS,H0_LL);

%% FISTA
% The FISTA package can be downloaded from: https://www.personal.psu.edu/thv102/
% Or directly from: https://filesender.surf.nl/?s=download&token=dc84a726-5ea5-43e9-8f4f-4d8178d466ec 

%{
lambda = 0.15; 
opts.max_iter = 75;

[F_fista, iter, store_fista]= lasso_fista(X.noise, phi, H0_LL, lambda, opts);
store_fista.maxiter = opts.max_iter;

%}
%% FLIPS (Linear Oracle)  
%{
manual_FW  ='on';
maxiter_FW = 500; 
       
[F_FW,store_FW] = Frank_Wolfe(X.noise,phi,epsilon,N,k,maxiter_FW,manual_FW,H0_LL_norm);
   
%}
%% Recover images 

% Image recovery 
Recover_input                       = X.no_noise;
Recover_input_noise                 = X.noise;
Recover_output_FLIPS                = phi*F_FLIPS;
Recover_output_csalsa               = phi*F_csalsa;
Recover_output_CP                   = phi*F_CP;
% Recover_output_FW                   = phi*F_FW;
% Recover_output_fista                = phi*F_fista;

reordered_im_input                  = patch2image(Recover_input,im.noise,store_inp);
reordered_im_input_noise            = patch2image(Recover_input_noise,im.noise,store_inp);
reordered_im_output_FLIPS           = patch2image(Recover_output_FLIPS,im.noise,store_inp);
reordered_im_output_csalsa          = patch2image(Recover_output_csalsa,im.noise,store_inp);
reordered_im_output_CP              = patch2image(Recover_output_CP,im.noise,store_inp);
% reordered_im_output_FW              = patch2image(Recover_output_FW,im.noise,store_inp);
% reordered_im_output_fista           = patch2image(Recover_output_fista,im.noise,store_inp);

%% Plot

imagee=zeros(imdim1,imdim2,2);
imagee(:,:,1) = reordered_im_output_FLIPS;
imagee(:,:,2) = reordered_im_output_csalsa;
imagee(:,:,3) = reordered_im_output_CP;
% imagee(:,:,4) = reordered_im_output_fista;
% imagee(:,:,5) = reordered_im_output_FW;

storeall(1) = store_FLIPS;
storeall(2) = store_csalsa; 
storeall(3) = store_CP; 
% storeall(4) = store_fista;
% storeall(5) = store_FW;

colorr(1,1:3,1) =[0,0.7,0.9];
colorr(1,1:3,2) =[1,0,1];
colorr(1,1:3,3) =[0.5,0.5,0.5];
colorr(1,1:3,4) =[0, 1, 0];
colorr(1,1:3,5) =[0, 0, 1];
% linestyles(1,1:2) ='--';
% linestyles(2,1:2) ='-.';
% linestyles(3,1)   =':';
linestyles(1,1) ='-';
linestyles(2,1) ='-';
linestyles(3,1) ='-';
marker(1,1:2)   ='-o';
marker(2,1:2)   ='-s';
marker(3,1:2)   ='-*';
LW = 1.5;

for i=1:3
    imageee = imagee(:,:,i);
    fig12=figure(14);
    subplot(3,3,1+3*(i-1))
    imagesc(imageee);
    colormap('gray');
    if i ==1
        title('FLIPS')
    elseif i==2 
        title('C-SALSA')
    elseif i==3
        title('Chambolle-Pock')     
    elseif i==4
        title('FISTA')
    elseif i==5
        title('FW')
    end
    subplot(3,3,2+3*(i-1))
    imagesc(reordered_im_input_noise);
    title('Noisy input')
    subplot(3,3,3+3*(i-1))
    imagesc(reordered_im_input);
    title('Input')
    drawnow;

    plotstore = storeall(i);
    fig6=figure(61);
    plot(1:size(plotstore.eta,2),(plotstore.eta)-store_FLIPS.eta(end),marker(i,:),"Color",colorr(:,:,i),"LineStyle",linestyles(i,:),'LineWidth',LW) % or loglog instead of plot
    hold on
    title('\eta vs iterations')
    grid on 
    legend('FLIPS','C-SALSA','Chambolle-Pock')
%     legend('FLIPS','C-SALSA','Chambolle-Pock','FISTA','FW')

    if i<=3
        fig4=figure(4);
        if i>=2
            hold on
        end
        plot(1:plotstore.maxiter,plotstore.gamma,marker(i,:),"Color",colorr(:,:,i),"LineStyle",linestyles(i,:),'LineWidth',LW)
        grid on
        title('Step sizes')
        ylabel('Step size')
        xlabel('Iterations')
    end
    legend('FLIPS','C-SALSA','Chambolle-Pock')
%     legend('FLIPS','C-SALSA','Chambolle-Pock','FISTA','FW')
 
    fig5=figure(5);
    normss=zeros(1,plotstore.maxiter);
    for j = 1:plotstore.maxiter
        normss(j)=norm(plotstore.h(:,j)-store_FLIPS.h(:,end),2);
    end
    plot(1:plotstore.maxiter,normss,marker(i,:),"Color",colorr(:,:,i),"LineStyle",linestyles(i,:),'LineWidth',LW)
    hold on
    grid on
    title('||f_t-f*||_2')
    ylabel('Step size')
    xlabel('Iterations')

    legend('FLIPS','C-SALSA','Chambolle-Pock')
%     legend('FLIPS','C-SALSA','Chambolle-Pock','FISTA','FW')
end    

