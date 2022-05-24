% Binning based on final x and y of photon
% Initialize bins to zero
bins  = zeros(10,10);

% Choose i based on x value
if x < 0.1 && x >= 0.
    i = 1;
else 
    if x < 0.2 && x >= 0.1
        i = 2;
    else
        if x < 0.3 && x >= 0.2
            i = 3;
        else
            if x < 0.4 && x >= 0.3
                i = 4;
            else
                if x < 0.5 && x >= 0.4
                    i = 5;
                else
                    if x < 0.6 && x >= 0.5
                        i = 6;
                    else
                        if x < 0.7 && x >= 0.6
                            i = 7;
                        else
                            if x < 0.8 && x >= 0.7
                                i = 8;
                            else
                                if x < 0.9 && x >= 0.8
                                    i = 9;
                                else
                                    if x < 1.0 && x >= 0.9
                                        i = 10;
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

% Choose i based on x value
if y < -0.4 && y >= -0.5
    j = 1;
else 
    if y < -0.3 && y >= -0.4
        j = 2;
    else
        if y < -0.2 && y >= -0.3
            j = 3;
        else
            if y < -0.1 && y >= -0.2
                j = 4;
            else
                if y < 0.0 && y >= -0.1
                    j = 5;
                else
                    if y < 0.1 && y >= 0.
                        j = 6;
                    else
                        if y < 0.2 && y >= 0.1
                            j = 7;
                        else
                            if y < 0.3 && y >= 0.2
                                j = 8;
                            else
                                if y < 0.4 && y >= 0.3
                                    j = 9;
                                else
                                    if y < 0.5 && y >= 0.4
                                        j = 10;
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

bins(i,j) = bins(i,j) + 1;