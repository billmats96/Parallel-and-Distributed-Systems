
%---------------test_meanshift.m Script-----------------
% This script generates random 2-dimensional datasets suited for clustering using
% generateData.m . Then we execute demo_meanshift.m with the random dataset
% and compare matlab output with the ouput of our testing code using the same
% dataset as input. If the results satisfy identification prerequisit (
% difference of ouputs <0.2) then test is successful.If the percentage of 
% coordinates differing>0.2 is larger than 2% of total coordinates then 
% possibly the test has failed. That's why we also compare(inside matlab script) 
%the norm of the matrix that results from the  differnce of two outputs. 
%If norm value is small then it's ok. 

% Author: Matsoukas Vasileios vmatsouk@auth.gr

%demo_meanshift.m by Dimitrios Floros
%generateData.m by Nuno Fachada

close all
clear 

%Define parameters
totalPoints=1024;           %Points of dataset 
numClusters=30;             %Clusters to be created
xClustAvgSep=4.5;           %average distance of cluster centers in x axis
yClustAvgSep=4.5;           %average distance of cluster centers in x axis
clusterFatness=0.5;         %the "radius" of cluster 
iter=15;                     %iterations for CUDA executable files
executable_file="serial_meanshift"; %choose between serial_meanshift, cuda_meanshift_sm, cuda_meanshift_gm
h=2*clusterFatness;         % choose h value.
dim="2";

[data cp idx] = generateData(0.5, 0.3, numClusters, xClustAvgSep, yClustAvgSep, 0.4, 0.1,clusterFatness,totalPoints);

labels=ones(totalPoints,1);
L=labels;
X=data;
save test_meanshift X L

%save random input as a binary file
A=data;
A=A';
fileID = fopen('randInput.bin','w');
fwrite(fileID,A,'double');
fclose(fileID);

%run matlab meanshift script

run('demo_meanshift.m')

%obtain matlab's output and write to a binary file
Out=y;
Out=Out';
fileID = fopen('outputmatlab.bin','w');
fwrite(fileID,Out,'double');
fclose(fileID);

%Test the executable with the random binary input produced above and compare 
%results with output binary files produced by matlab script
%execute and check output in matlab command window
fprintf("Testing dataset on: %s\n",executable_file);
h=sqrt(h);
h=num2str(h);
% dataset_cols=num2str(dataset_cols);
% if executable_file=="serial_meanshift"
iter=num2str(iter);
binaryInput="randInput.bin";

if executable_file=="serial_meanshift"
    command=['./',char(executable_file),' ',char(binaryInput),' ',char(h),' ','outputTEST.bin ',char(dim)];
else
    command=['./',char(executable_file),' ',char(binaryInput),' ',char(h),' ',char(iter),' ','outputTEST.bin ',char(dim)];
end
[status,cmdout] = system(command);

cmdout   
 
%Get output of testing code and compare with matlab results. inside matlab
fileID = fopen('outputTEST.bin','r');
B=fread(fileID,[2,totalPoints],'double');
B=B';
fclose(fileID);
Out=Out';
C=B-Out; %C matrix holds the difference of matlab output from testing code output
fprintf('Norm of comparison matrix C is: %f\n',norm(C));
if(norm(C)>10)
    fprintf('Considerable large value of norm(C)\n');
else
    fprintf('Accepted norm value\n'); 
end
%max(C)
fprintf('End of test\n');







