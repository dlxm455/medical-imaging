% Simulation of Parallel-Ray (1st Generation) Computed Tomography (CT) 
% Shijia Geng

%% 1 Read in 4 "cross-sectional" images, convert to gray scale and resize
IMG=['image1.jpg','image2.jpg','image3.jpg','image4.jpg'];
for ii=1:4
A=imread(IMG(ii*10-9:ii*10));
B=rgb2gray(A);
s=size(B);
img1=imresize(B,96/max(s));
[row, col]=size(img1);
figure(1);
subplot (3,4,ii);
imshow(A);
title(strcat('Image ',int2str(ii),' Orginal'));
subplot(3,4,ii+4);
imshow(B);
title(strcat('Image',int2str(ii),'Gray'));
subplot(3,4,ii+8);
imshow(img1);
title(strcat('Image',int2str(ii),'Resized'));

%% 2 Numerical Approach
RRA=imresize(img1,[5 4]);
RA=double(RRA);
AA=[  1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 ;
      0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 ;
      0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 ;
      0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 ;
      1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
      0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 ;
      0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 ;
      0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 ;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 ;
      0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 ;
      0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 ;
      0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 ;
      1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 ;
      0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 ;
      0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ;
      0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
      0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;
      0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 ;
      0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 ;
      0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 ;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 ];

bb2=sum(RA,2)';
bb1=sum(RA,1);
bb4=[RA(1,2)+RA(2,1),RA(1,3)+RA(2,2)+RA(3,1),RA(1,4)+RA(2,3)+RA(3,2)+RA(4,1),RA(2,4)+RA(3,3)+RA(4,2)+RA(5,1),RA(3,4)+RA(4,3)+RA(5,2),RA(4,4)+RA(5,3)];
bb3=[RA(4,1)+RA(5,2),RA(3,1)+RA(4,2)+RA(5,3),RA(2,1)+RA(3,2)+RA(4,3)+RA(5,4),RA(1,1)+RA(2,2)+RA(3,3)+RA(4,4),RA(1,2)+RA(2,3)+RA(3,4),RA(1,3)+RA(2,4)];
bb=[bb1 bb2 bb3 bb4]';

IA=pinv(AA);
XX=IA*bb;
X=reshape(XX,4,5);
X_rec=X';
DX=X_rec-RA;

figure(2);
subplot(3,4,ii);
imshow(RRA,[]);
title(strcat('Image',int2str(ii),' Resized Image'))
% subplot(3,4,ii);
% plot(bb);
% title(strcat('Image',int2str(ii),' Projection'))
subplot(3,4,ii+4);
imshow(X_rec,[]);
title(strcat('Image',int2str(ii),' Reconstruction'))
subplot(3,4,ii+8);
imshow(DX,[])
title(strcat('Image',int2str(ii),' Difference'))

%% 3 Radon Transform
p=zeros(round(sqrt(row^2+col^2))+5,180);
row_p=size(p,1);
for theta=1:180
    img1_r=imrotate(img1,theta,'loose');
    [row_r, col_r]=size(img1_r);
    p(1+round((row_p-col_r)/2):col_r+round((row_p-col_r)/2),theta)=sum(img1_r,1)';
end
figure(3);
subplot(3,4,ii);
imshow(p,[]);
title(strcat('Image',int2str(ii),' Projectioin Image'));

%% 4 Direct Backprojection
g=zeros(row,col);
    for i=1:row
        for j=1:col
            for theta=1:180
            g(i,j)=g(i,j)+p(round((j-col/2)*cos(theta/180*pi)+(i-row/2)*sin(theta/180*pi)+row_p/2),theta);
            end
        end
   end
figure(4);
subplot(3,4,ii);
imshow(g,[]);
title(strcat('Image',int2str(ii),' Direct BP Reconstruction'));

%% 5 Filtered Backprojection
f_radon=zeros(row_p,180);
ff=zeros(row_p,1);
for theta=1:180
    for i=1:row_p
        ff(i,1)=abs(i-round(row_p/2));
    end
    P=fft(p(:,theta));
    ft_piece=ff.*fftshift(P);
    ift_piece=ifft(fftshift(ft_piece));
    f_radon(:,theta)=ift_piece;
end

h=zeros(row,col);
    for i=1:row
        for j=1:col
            for theta=1:180
            h(i,j)=h(i,j)+f_radon(round((j-col/2)*cos(theta/180*pi)+(i-row/2)*sin(theta/180*pi)+row_p/2),theta);
            end
        end
    end
d=h-g;
figure(5);
% subplot(4,4,ii);
% plot(ff);
% title('Frequency Filter');
% subplot(3,4,ii);
% imshow(abs(f_radon),[]);
% title(strcat('Image', int2str(ii),'FD Filtered Projection Image'));
subplot(3,4,ii+4);
imshow(abs(h),[]);
title(strcat('Image', int2str(ii),'FD Filtered BP Reconstruction'));
subplot(3,4,ii+8);
imshow(abs(d),[]);
title(strcat('Image', int2str(ii),'FD Filtered Reconstruction Difference'));

%% 6 Convolution Backprojection
s_radon=zeros(row_p,180);
sf=fftshift(ifft(fftshift(ff)));
for theta=1:180
    s_piece=p(:,theta);
    sf_piece=conv(s_piece,sf,'same');
    s_radon(:,theta)=sf_piece;
end
k=zeros(row,col);
for i=1:row
    for j=1:col
        for theta=1:180
            k(i,j)=k(i,j)+s_radon(round((j-col/2)*cos(theta/180*pi)+(i-row/2)*sin(theta/180*pi)+row_p/2),theta);
        end
    end
end
dd=k-g;
figure(6);
% subplot(4,4,ii);
% imshow(abs(sf),[]);
% title('Spatial Filter');
% subplot(4,4,ii+4);
% imshow(abs(s_radon),[]);
% title(strcat('Image', int2str(ii),' SD Filtered Projection Image'));
subplot(3,4,ii);
imshow(abs(k),[]);
title(strcat('Image', int2str(ii),' SD Filtered BP Reconstruction'));
subplot(3,4,ii+4);
imshow(abs(dd),[])
title(strcat('Image', int2str(ii),' SD Filtered Reconstruction Difference'));

%% 7 Testing
load('proj2_10.mat');
q=Proj';
[row_q, col_q]=size(q);
row=row_q-50;
col=row_q-50;
bq1=zeros(row,col);
for i=1:row
    for j=1:col
        for k=1:col_q
            theta=(180/col_q)*k;
            bq1(i,j)=bq1(i,j)+q(round((j-col/2)*cos(theta/180*pi)+(row/2-i)*sin(theta/180*pi)+row_q/2),k);
        end
    end
end
figure(7)
subplot(3,4,1);
imshow(q,[]);
title('Projection');
subplot(3,4,2);
imshow(bq1,[]);
title ('Direct BP Reconstruction');

qf_radon=zeros(row_q,col_q);
qff=zeros(row_q,1);
for k=1:col_q
    for i=1:row_q
        qff(i,1)=abs(i-round(row_q/2));
    end
    Q=fft(q(:,k));
    qft_piece=qff.*fftshift(Q);
    qift_piece=ifft(fftshift(qft_piece));
    qf_radon(:,k)=qift_piece;
end
bq2=zeros(row,col);
    for i=1:row
        for j=1:col
            for  k=1:col_q
                 theta=(180/col_q)*k;
                 bq2(i,j)=bq2(i,j)+qf_radon(round((j-col/2)*cos(theta/180*pi)+(row/2-i)*sin(theta/180*pi)+row_q/2),k);
            end
        end
    end
    
subplot(3,4,3);    
imshow(abs(bq2),[]);
title('FD Filtered BP Reconstruction');


qs_radon=zeros(row_q,col_q);
qsf=fftshift(ifft(fftshift(qff)));
for k=1:col_q
    qs_piece=q(:,k);
    qsf_piece=conv(qs_piece,qsf,'same');
    qs_radon(:,k)=qsf_piece;
end
bq3=zeros(row,col);
    for i=1:row
        for j=1:col
            for  k=1:col_q
                 theta=(180/col_q)*k;
                 bq3(i,j)=bq3(i,j)+qs_radon(round((j-col/2)*cos(theta/180*pi)+(row/2-i)*sin(theta/180*pi)+row_q/2),k);
            end
        end
    end
subplot(3,4,4);
imshow(abs(bq3), []);
title('SD Filtered BP Reconstruction');
end
