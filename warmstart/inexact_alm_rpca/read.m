clear;
clc;

label = ['yaleB39_P00A+000E+00.pgm'
'yaleB39_P00A+010E-20.pgm'
'yaleB39_P00A+020E-10.pgm'
'yaleB39_P00A+025E+00.pgm'
'yaleB39_P00A+020E+10.pgm'
'yaleB39_P00A+015E+20.pgm'
'yaleB39_P00A+000E+20.pgm'
'yaleB39_P00A-015E+20.pgm'
'yaleB39_P00A-020E+10.pgm'
'yaleB39_P00A-025E+00.pgm'
'yaleB39_P00A-020E-10.pgm'
'yaleB39_P00A-010E-20.pgm'
'yaleB39_P00A+000E-20.pgm'
'yaleB39_P00A+050E-40.pgm'
'yaleB39_P00A+060E-20.pgm'
'yaleB39_P00A+070E+00.pgm'
'yaleB39_P00A+060E+20.pgm'
'yaleB39_P00A+070E+45.pgm'
'yaleB39_P00A+035E+65.pgm'
'yaleB39_P00A-035E+65.pgm'
'yaleB39_P00A-070E+45.pgm'
'yaleB39_P00A-060E+20.pgm'
'yaleB39_P00A-070E+00.pgm'
'yaleB39_P00A-060E-20.pgm'
'yaleB39_P00A-050E-40.pgm'
'yaleB39_P00A-130E+20.pgm'
'yaleB39_P00A-110E+15.pgm'
'yaleB39_P00A-120E+00.pgm'
'yaleB39_P00A-110E-20.pgm'
'yaleB39_P00A+130E+20.pgm'
'yaleB39_P00A+110E+15.pgm'
'yaleB39_P00A+120E+00.pgm'
'yaleB39_P00A+110E-20.pgm'
'yaleB39_P00A-070E-35.pgm'
'yaleB39_P00A-085E-20.pgm'
'yaleB39_P00A-095E+00.pgm'
'yaleB39_P00A-085E+20.pgm'
'yaleB39_P00A-110E+40.pgm'
'yaleB39_P00A-110E+65.pgm'
'yaleB39_P00A+000E+90.pgm'
'yaleB39_P00A+110E+65.pgm'
'yaleB39_P00A+110E+40.pgm'
'yaleB39_P00A+085E+20.pgm'
'yaleB39_P00A+095E+00.pgm'
'yaleB39_P00A+085E-20.pgm'
'yaleB39_P00A+070E-35.pgm'
'yaleB39_P00A-020E-40.pgm'
'yaleB39_P00A-035E-20.pgm'
'yaleB39_P00A-050E+00.pgm'
'yaleB39_P00A-035E+15.pgm'
'yaleB39_P00A-035E+40.pgm'
'yaleB39_P00A+000E+45.pgm'
'yaleB39_P00A+035E+40.pgm'
'yaleB39_P00A+035E+15.pgm'
'yaleB39_P00A+050E+00.pgm'
'yaleB39_P00A+035E-20.pgm'
'yaleB39_P00A+020E-40.pgm'
'yaleB39_P00A+000E-35.pgm'
'yaleB39_P00A-005E-10.pgm'
'yaleB39_P00A-010E+00.pgm'
'yaleB39_P00A-005E+10.pgm'
'yaleB39_P00A+005E+10.pgm'
'yaleB39_P00A+010E+00.pgm'
'yaleB39_P00A+005E-10.pgm'];

k = 1;
kbad = 1;
for i = 1:9
    p1 = strcat('E:\DeYuMeng''s files\Experiments\MFRealYaleB\CroppedYale\yaleB0',num2str(i));
    cd(p1);
    Tlabel = label;
    Tlabel(:,6) = '0';
    Tlabel(:,7) = num2str(i);
    for j = 1:64
        im = imread(Tlabel(j,:),'pgm');
        vim = vec(im);
        X(:,k) = vim; 
        k = k+1;
    end
end

indbad = [];
for i = [10:13 15:39]
    p1 = strcat('E:\DeYuMeng''s files\Experiments\MFRealYaleB\CroppedYale\yaleB',num2str(i));
    cd(p1);
    Tlabel = label;
    p2 = num2str(i);
    Tlabel(:,6) = p2(1);
    Tlabel(:,7) = p2(2);
    for j = 1:64
        if exist(Tlabel(j,:),'file') == 0
            p3 = strcat(Tlabel(j,:),'.bad');
            im = imread(p3,'pgm');
            vim = vec(im);
            X(:,k) = vim; 
            indbad = [indbad k];
            k = k+1;
        else
            im = imread(Tlabel(j,:),'pgm');
            vim = vec(im);
            X(:,k) = vim; 
            k = k+1;
        end
    end
end