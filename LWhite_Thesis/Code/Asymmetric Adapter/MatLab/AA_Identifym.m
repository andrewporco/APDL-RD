clc
clear
compx = 0.5;% Adjust stl in x axis through compensation
compy =0.5;% Adjust stl in y axis through compensation
compz = 0.5;% Adjust stl in z axis through compensation
gridX =(compx:39+compx);%22;
gridY =(compy:105+compy);%96;
gridZ =(compz:51+compz);%63;
% gridX = (0.5:10);
% gridY =(0.5:5.5);
% gridZ =(0:5);
Xx =[ ]; %initialize array
C = [ ];
P = [ ];
[gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,'BulkSup_Part.stl'); %Completely solid structure
[faces,vertices] = CONVERT_voxels_to_stl('AsymBulkCombined.stl',gridOUTPUT,gridCOx,gridCOy,gridCOz,'ascii'); %Voxelized structure
kk=1;
Kk = 1;
KK=1;
%obtain x, y, and z coordinates of all vertices
size = 1; %cell size
% gridOUTPUT(:,:,1)
% gridOUTPUT(3,1,1)
% gridOUTPUT(:,:,1)
% flipud(gridOUTPUT(:,:,1)')
for k = 1:length(gridZ)%Y-direction
    for j = 1:length(gridY) %Z-direction
        for i=1: length(gridX) %X-direction
            if k<length(gridZ)

                if gridOUTPUT(i,j,k)==1 && gridOUTPUT(i,j,k+1)~=0
                    X = gridCOx(i)-compx;
                    Y = gridCOy(j)-compy;
                    Z = gridCOz(k)-compz;
                    for Cc = 0:1
                        for Bb = 0:1
                            for Aa = 0:1
                                P(KK,1) =  X+Aa*size;
                                P(KK,2) =  Y+Bb*size;
                                P(KK,3) =  Z+Cc*size;
                                KK= KK+1;
                            end
                        end
                    end
            end                
            end
        end
    end
end
% [vol_handle]=VoxelPlotter(gridOUTPUT,1);
% alpha(0.4)

hold on;
% P=unique(P,'rows');
% % % pts = unique( C 'rows',);
pts = P;
% Import an STL mesh, returning a PATCH-compatible face-vertex structure
[F,V,N]=stlread('Part.stl');
V(V<1e-6)=0;
V(:,1) = V(:,1);
V(:,2) = V(:,2);%-1.5;
V(:,3) = V(:,3);%-4.5;

in = inpolyhedron(F,V,pts,'FlipNormals',false);
Part = [ ]; %initialize for part point cloud

Raw = [ ];  %point cloud of all vertices
E = [P in]; %each vertex that is the part
Pp = 1;
n = 1;
N = 1;
R = 1;
out   = permute(reshape(E',4,8,[]),[2,1,3]);%sort into cubes with 8 vertices
%% Get points of potential support structure locations
for I=1:length(out)
    if ismember(1,out(:,4,I))
        Raw(R,1) = out(1,1,I);
        Raw(R,2) = out(1,2,I);
        Raw(R,3) = out(1,3,I);
        R = R+1;
        for ii = 1:8
            Part(Pp,1) = out(ii,1,I);
            Part(Pp,2) = out(ii,2,I);
            Part(Pp,3) = out(ii,3,I);
            Pp = Pp+1;
        end       
    else
        for ii = 1:8
            t(N,1) = out(ii,1,I);
            t(N,2) = out(ii,2,I);
            t(N,3) = out(ii,3,I);
            N = N+1;
        end
    end
end
Raw = unique(Raw,'rows');
% P(Del,:) = [ ];
t = unique(t,'rows');

%% Final  support structure list %%%%%%
% Sort unit cells from point cloud
xmax = max(P(:,1));
ymax = max(P(:,2));
zmax = max(P(:,3));
S_ = [ ]; %initialize 2 mm unit cell size
Sing_ = [ ]; %initialize 1 mm unit cell size
x_c =0; %minimum x coordinate value
y_c= 0; %minimum y coordinate value
z_c = 0;
Nn = 1;

while z_c < zmax
    while y_c < ymax
        while x_c < xmax
            A = [x_c y_c z_c];
            if ismember(A,t,'rows')==1
%                 disp('In1')
                x_p=x_c;
                y_p=y_c+2;
                z_p=z_c;
                B = [x_p y_p z_p];
                if ismember(B,t,'rows')==1
%                     disp('In2')
                    x_p=x_c+2;
                    y_p=y_c+2;
                    z_p=z_c+2;
                    C = [x_p y_p z_p];
                    if ismember(C,t,'rows')==1
%                         disp('In3')
                        x_p=x_c+2;
                        y_p=y_c+2;
                        z_p=z_c+0;
                        D = [x_p y_p z_p];
                        if ismember(D,t,'rows')==1
                            x_p=x_c+2;
                            y_p=y_c;
                            z_p=z_c;
                            E = [x_p y_p z_p];
                            if ismember(E,t,'rows')==1
                                x_p=x_c+2;
                                y_p=y_c;
                                z_p=z_c+2;
                                F = [x_p y_p z_p];
                                if ismember(F,t,'rows')==1
                                    x_p=x_c;
                                    y_p=y_c;
                                    z_p=z_c+2;
                                    G = [x_p y_p z_p];
                                    if ismember(G,t,'rows')==1
                                        x_p=x_c;
                                        y_p=y_c+2;
                                        z_p=z_c+2;
                                        H = [x_p y_p z_p];
                                        if ismember(H,t,'rows')==1
                                            S_(Nn,1) = x_c;
                                            S_(Nn,2) = y_c;
                                            S_(Nn,3) = z_c;
                                            x_c = x_c+2;
                                            Nn = Nn+1;
                                        else
                                            Temp = sing(x_c, y_c, z_c,t);
                                            Sing_ = [Sing_;Temp];
                                             x_c=x_c+1;     
                                        end
                                    else
                                        Temp = sing(x_c, y_c, z_c,t);
                                        Sing_ = [Sing_;Temp];
                                         x_c=x_c+1;     
                                    end
                                else
                                    Temp = sing(x_c, y_c, z_c,t);
                                    Sing_ = [Sing_;Temp];
                                    x_c=x_c+1;     
                                end
                            else
                            Temp = sing(x_c, y_c, z_c,t);
                            Sing_ = [Sing_;Temp];
                            x_c=x_c+1;      
                        end
                        else
                        Temp = sing(x_c, y_c, z_c,t);
                        Sing_ = [Sing_;Temp];
                        x_c = x_c+1;
                    end
                    else
                        Temp = sing(x_c, y_c, z_c,t);
                        Sing_ = [Sing_;Temp];
                    x_c = x_c+1;
                    end
                else
                Temp = sing(x_c, y_c, z_c,t);
                Sing_ = [Sing_;Temp];

                x_c = x_c+1;
                end
            else
                Temp = sing(x_c, y_c, z_c,t);
                Sing_ = [Sing_;Temp];
                x_c=x_c+1;     
            end
        end
        y_c = y_c+2;
        x_c = 0;
    end
    z_c = z_c+2;
    x_c = 0;
    y_c = 0;
end
% 
% 
% pts1 = S_;
% pts2 = Sing_;
% pts4 = [S_; Sing_];
WSing_ = [ ];
Nn = 1;
for i =1:length(Sing_)
    X = Sing_(i,1);
    Y = Sing_(i,2);
    Z = Sing_(i,3);
    for j =0:1
        for k = 0:1
            for l = 0:1
                WSing_(Nn,1)=X+l;
                WSing_(Nn,2)=Y+k;
                WSing_(Nn,3)=Z+j;
                Nn=Nn+1;
            end
        end
    end
end
