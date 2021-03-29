#!/bin/bash

clear

echo Running script..

touch results.txt


echo Denoising image1 176x117 ..

INPUT="in1_noisy.bin"
ROWS="117"
COLS="176"
OUTPUT1="out1_snlm.bin"
OUTPUT2="out1_snlm_ad.bin"
OUTPUT3="out1_nlm.bin"
OUTPUT4="out1_ad.bin"



echo in1&>>results.txt

./serialnlm $INPUT $OUTPUT1 $ROWS $COLS &>>results.txt
./serialnlm_adaptive $INPUT $OUTPUT2 $ROWS $COLS &>>results.txt
./nlm $INPUT $OUTPUT3 $ROWS $COLS &>>results.txt
./adaptive $INPUT $OUTPUT4 $ROWS $COLS &>>results.txt

echo ---------------------------------------------------------------- &>>results.txt




echo Denoising image2 360x241 ..

INPUT="in2_noisy.bin"
ROWS="241"
COLS="360"
OUTPUT1="out2_snlm.bin"
OUTPUT2="out2_snlm_ad.bin"
OUTPUT3="out2_nlm.bin"
OUTPUT4="out2_ad.bin"



echo in2&>>results.txt

./serialnlm $INPUT $OUTPUT1 $ROWS $COLS &>>results.txt
./serialnlm_adaptive $INPUT $OUTPUT2 $ROWS $COLS &>>results.txt
./nlm $INPUT $OUTPUT3 $ROWS $COLS &>>results.txt
./adaptive $INPUT $OUTPUT4 $ROWS $COLS &>>results.txt

echo --------------------------------------------------------------- &>>results.txt


echo Denoising image3 640x360 ..

INPUT="in3_noisy.bin"
ROWS="360"
COLS="640"
OUTPUT1="out3_snlm.bin"
OUTPUT2="out3_snlm_ad.bin"
OUTPUT3="out3_nlm.bin"
OUTPUT4="out3_ad.bin"



echo in3&>>results.txt

./serialnlm $INPUT $OUTPUT1 $ROWS $COLS &>>results.txt
./serialnlm_adaptive $INPUT $OUTPUT2 $ROWS $COLS &>>results.txt
./nlm $INPUT $OUTPUT3 $ROWS $COLS &>>results.txt
./adaptive $INPUT $OUTPUT4 $ROWS $COLS &>>results.txt

echo -------------------------------------------------------------- &>>results.txt


echo Denoising image4 813x460 ..

INPUT="in4_noisy.bin"
ROWS="460"
COLS="813"
OUTPUT1="out4_snlm.bin"
OUTPUT2="out4_snlm_ad.bin"
OUTPUT3="out4_nlm.bin"
OUTPUT4="out4_ad.bin"



echo in4&>>results.txt

./serialnlm $INPUT $OUTPUT1 $ROWS $COLS &>>results.txt
./serialnlm_adaptive $INPUT $OUTPUT2 $ROWS $COLS &>>results.txt
./nlm $INPUT $OUTPUT3 $ROWS $COLS &>>results.txt
./adaptive $INPUT $OUTPUT4 $ROWS $COLS &>>results.txt

echo ------------------------------------------------------------ &>>results.txt


echo Denoising image5 1000x1000 ..

INPUT="in5_noisy.bin"
ROWS="1000"
COLS="1000"
OUTPUT1="out5_snlm.bin"
OUTPUT2="out5_snlm_ad.bin"
OUTPUT3="out5_nlm.bin"
OUTPUT4="out5_ad.bin"



echo in5&>>results.txt

./serialnlm $INPUT $OUTPUT1 $ROWS $COLS &>>results.txt
./serialnlm_adaptive $INPUT $OUTPUT2 $ROWS $COLS &>>results.txt
./nlm $INPUT $OUTPUT3 $ROWS $COLS &>>results.txt
./adaptive $INPUT $OUTPUT4 $ROWS $COLS &>>results.txt

echo ------------------------------------------------------------- &>>results.txt

echo END &>>results.txt

echo End of script..
