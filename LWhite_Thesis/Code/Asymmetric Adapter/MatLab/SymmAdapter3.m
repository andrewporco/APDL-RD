clc
clear
compx = 0.5;% Adjust stl in x axis through compensation
compy = 7.5;% Adjust stl in y axis through compensation
compz = 24.5;% Adjust stl in z axis through compensation
gridX =(compx:12+compx);%22;
gridY =(compy:48+compy);%96;
gridZ =(compz:35+compz);%63;

Xx =[ ]; %initialize array
C = [ ];
P = [ ];
[gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,'AlignedScaledSymm_Sup_12.stl'); %Completely solid structure
[faces,vertices] = CONVERT_voxels_to_stl('SymmScaledAdapt2.stl',gridOUTPUT,gridCOx,gridCOy,gridCOz,'ascii');
kk=1;
Kk = 1;
KK=1;
%obtain x, y, and z coordinates of all vertices
size = 1; %cell size
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


hold on;
pts = P;
% % % Import an STL mesh, returning a PATCH-compatible face-vertex structure
[F,V,N]=stlread('AlignedScaledSymm_12.stl');
V(V<1e-6)=0;
V(:,2) = V(:,2)-1.5;
V(:,3) = V(:,3)-4.5;
patch('Vertices',V,...
    'Faces',F,...
    'FaceColor',[.7 .7 .7],...
    'LineWidth',0.25);
in = inpolyhedron(F,V,pts,'FlipNormals',false);
Part = [ ]; %initialize for part point cloud
Del = [ ]; %point cloud that is deleted to create support structure
E = [P in]; %each vertex that is the part
Pp = 1;
n = 1;
N = 1;
out   = permute(reshape(E',4,8,[]),[2,1,3]); %sort into cubes with 8 vertices
%% Get points of potential support structure locations
for I=1:length(out)
    if ismember(1,out(:,4,I))
        Del(Pp) = I;
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
Part = unique(Part,'rows');
% P(Del,:) = [ ];
x_c =0;
y_c= 0;
z_c = 0;
xmax = max(Part(:,1));
ymax = max(Part(:,2));
zmax = max(Part(:,3));
PartSing = [ ];

% % %%%%%% Final support structure list %%%%%%
%% Sort unit cells from point cloud
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



% pts4 = [S_; Sing_];
% WSing_ = [ ];
Nn = 1;
for i =1:length(Sing_)
    X = Sing_(i,1); %X coordinate of single unit cell
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
WS_ = [ ];
Nn = 1;
for i =1:length(S_)
    X = S_(i,1);
    Y = S_(i,2);
    Z = S_(i,3);
    for j =0:2
        for k = 0:2
            for I = 0:2
                WS_(Nn,1)=X+I;
                WS_(Nn,2)=Y+k;
                WS_(Nn,3)=Z+j;
                Nn=Nn+1;
            end
        end
    end
end
pts1 = WS_; % all vertices of 2mm unit cells
pts2 = WSing_; % all vertices of 1mm unit cells

plot3(pts1(:,1),pts1(:,2),pts1(:,3),'.b', 'MarkerSize',15); %add support position points to plot

plot3(pts2(:,1),pts2(:,2),pts2(:,3),'.C', 'MarkerSize',15); %add support position points to plot
