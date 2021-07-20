clear all
clc
fp = fopen('cdf_table.txt', 'w');
d1 = 0.002;
d2 = 0.01;
i1=2:d1:3-d1;
i2 = 3:d2:5-d2;
y1 = normcdf(i1);
y2 = normcdf(i2);
for i=1:length(y1)
    fprintf(fp, '%16.15f,  ',y1(i));
    if(mod(i, 20)==0)
        fprintf(fp, '\n');
    end
end
fprintf(fp, '\n');
fprintf(fp, '\n');
for i=1:length(y2)
    fprintf(fp, '%16.15f,  ',y2(i));
    if(mod(i, 20)==0)
        fprintf(fp, '\n');
    end
end
fclose(fp);