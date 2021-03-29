%---------------static_test.m Script-----------------
% This script runs matlab demo_meanshift2.m with a given dataset and stores
% result in outputmatlab.bin. Then it runs our executable code with the
% same dataset input in binary version. Our executable code compares matlab output 
% with it own results, comparing corresponding coordinates. We assume that a
% difference <0.2 is acceptable (this is due to matlab inaccuracy in
% results). If the percentage of coordinates differing>0.2 is larger than 2%
% of total coordinates then possibly the test has failed. That's why 
% we also compare(inside matlab script) the norm of the matrix that results from the
% differnce of two outputs. If norm value is small then it's ok. 

% Author: Matsoukas Vasileios vmatsouk@auth.gr
% demo_meanshift2.m by Dimitrios Floros

%Define parameters
totalPoints=1024;            %Points of dataset 
matInput='test32';           %.mat file with vector L and dataset X
binaryInput='input32.bin';   %binary Input to testing code
h=25;                        %value of h
executable_file="serial_meanshift"; %choose between serial_meanshift, cuda_meanshift_sm, cuda_meanshift_gm
dimension=32;                %dimension of set (number of Coordinates per point) 
iter=3;                      %iterations for CUDA executable files


%run matlab meanshift script
run('demo_meanshift2.m')

%obtain matlab's output and write to a binary file
Out=y;
Out=Out';
fileID = fopen('outputmatlab.bin','w');
fwrite(fileID,Out,'double');
fclose(fileID);

%Test the executable with the binaryInput and compare 
%results with output binary files produced by matlab script
%check output in matlab command window
fprintf("Testing dataset on: %s\n",executable_file);
h=sqrt(h);
h=num2str(h);
dim=num2str(dimension);
iter=num2str(iter);
if executable_file=="serial_meanshift"
    command=['./',char(executable_file),' ',char(binaryInput),' ',char(h),' ','outputTEST.bin ',char(dim)];
else
    command=['./',char(executable_file),' ',char(binaryInput),' ',char(h),' ',char(iter),' ','outputTEST.bin ',char(dim)];
end
[status,cmdout] = system(command);

cmdout   
 
%Get output of testing code and compare with matlab results. inside matlab
fileID = fopen('outputTEST.bin','r');
B=fread(fileID,[dimension,totalPoints],'double');
B=B';
fclose(fileID);
Out=Out';
C=B-Out; %C matrix holds the difference of matlab output from testing code output
fprintf("Norm of comparison matrix C is: %f\n",norm(C));
if(norm(C)>10)
    fprintf("Considerable large value of norm(C)\n");
else
    fprintf("Accepted norm value\n"); 
end
%max(C)
fprintf("End of test\n");


