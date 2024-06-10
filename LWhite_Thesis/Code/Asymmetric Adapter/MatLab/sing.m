%% Identify full 1 mm unit cells
function Temp_ = sing(x_c, y_c, z_c,Part)
Temp_ = [ ];
Na = 1;
for z = 0:1
    for y = 0:1
        x_c = x_c;
        y_c2 = y_c+y;
        z_c2 = z_c+z;
        A = [x_c y_c2 z_c2];
        if ismember(A,Part,'rows')==1
            x_p=x_c;
            y_p=y_c+1;
            z_p=z_c2;
            B = [x_p y_p z_p];
            if ismember(B,Part,'rows')==1
                x_p=x_c+1;
                y_p=y_c2+1;
                z_p=z_c2+1;
                C = [x_p y_p z_p];
                if ismember(C,Part,'rows')==1
                    x_p=x_c+1;
                    y_p=y_c2+1;
                    z_p=z_c2;
                    D = [x_p y_p z_p];
                    if ismember(D,Part,'rows')==1
                        x_p=x_c+1;
                        y_p=y_c2;
                        z_p=z_c2;
                        E = [x_p y_p z_p];
                        if ismember(E,Part,'rows')==1
                            x_p=x_c+1;
                            y_p=y_c2;
                            z_p=z_c+1;
                            F = [x_p y_p z_p];
                            if ismember(F,Part,'rows')==1
                                x_p=x_c;
                                y_p=y_c2;
                                z_p=z_c+1;
                                G = [x_p y_p z_p];
                                if ismember(G,Part,'rows')==1
                                    x_p=x_c;
                                    y_p=y_c2+1;
                                    z_p=z_c2+1;
                                    H = [x_p y_p z_p];
                                    if ismember(H,Part,'rows')==1
                                        Temp_(Na,1) = x_c;
                                        Temp_(Na,2) = y_c2;
                                        Temp_(Na,3) = z_c2;
                                        Na = Na+1;

                                    end

                                end

                            end
 
                        end

                    end

                end

            end
        end
    end
end
end