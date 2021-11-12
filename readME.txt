SP assignment 3 - Vowel Recognition
Roll No - 214101011
Name - Atul Bunkar
-----------------------------------------------------------------

NOTE: Please ignore my other submission, and consider this one.

-----------------------------------------------------------------
How to run :

- Check that folder "test" and "reference" are empty before running the code.

- I have provided the 100 input files in "vowels" folder , so please dont use assigment-2's folder as it has a missing file which is throwing error. 

- Simply build and run the file "main.cpp" . There is only one cpp file with all the functions.
 
- No utility file is used separately. All steps are done in main.cpp itself.

----------------------------------------------------------------
Creation of Files :

- During execution of program , all files will be read from vowels and stored in 2 folders,
for this, I have created  2 empty folders, named "reference" and "test" where respective files will be created on running the program.

- reference will have 5 files, 1 for each vowel. Each file will have 60 ci values i.e 12 ci for 5 frames. 1 ci stored in 1 line. Total 60 lines.

- test will have 10 files for each vowel (Ex: For a ,test_a_1 ,test_a_2 ... so on) ,these files also has 60 lines of ci's. Read the function explantion for more detail. 

----------------------------------------------------------------

Flow of main() :
step 1. First checks ri,ai,ci for given sample test and outputs in console.
	-> order of functions called: get_samples -> calc_ri -> calc_ai -> calc_ci .
			(print for showing array)

step 2. generates all reference and then all test files. 
	->Functn called : generate_refernce() and generate_test()

step 3. Predction of vowels
	-> Functn called : vowel_recognition()
------------------------------------------------------------------

Functions used :-

1. skip_line - To skip some lines in reading a file.

2. dc_fix - Reads the dc file and calculates the dc shift, stores it in dc_bias.

3. get_samples - reads the sample test file and fills 1D samples array of 320 i.e one frame. Used only for step 1.

4. calc_ri - calculates autocorrelation for one frame then stores it in ri[] array(passed by reference)

5. calc_ai - Using Durbin , calculates ai using ri, for one frame and stores it in array ai[] .

6. calc_ci - Calculated cepstral coeffs (using ai and ri) for one frame and stores in ci[].

7. raised_sine - applies raised sine window to ci by multiplying with each ci's.

8. preprocess - Normalises the data and stores it in static array norm_data(defined at top).

9. find_start - Finds the start of stable array by getting max energy frame and moving back 2 frames.

10. generate_reference - Generates reference files , taking "1-10" from each vowel and applying above functions in given order - 
 	: preprocess (which calls dc_fix) ->  find_start ->  calc_ri -> calc_ai -> calc_ci -> raised_sine .
Then averaging ci over 10 files to get 1 reference file for 1 vowel. These are named ref_a , ref_e...so on.
(Ex: ref_a.txt is avg_ci of 214101011_a_1.txt to 214101011_a_10.txt.)
Creates 1 reference file for 10 vowel file, So 5 files created for 50 vowel file.

11. generate_test - generated test files , reads "11-20" files of each vowel , generates 12*5 ci for each file by calling same functions as above and store in test.
(Ex: test_a_1.txt is ci for file 214101011_a_11.txt and so on. )
creates 1 ci file for 1 vowel file ,so total 50 test files are created for 50 vowel files .

12. Vowel_recognition - Reads all generated test file one by one and for each file, read 5 reference files and calculate tokuharo distance and report it , and find how many accurate prediction is made based on the distance.

13. print - To print array elements in given range(used in step 1 only)

-----------------------------------------------------------------

Console Output: 

Step 1: prerequisite:
Output can be checked in console and matched with dump ri,ai,ci text files, they are same.

for step 2 : generation of files:
No output is shown in console.  Files are created which can be checked in the folders.

For step 3 : Vowel Recognition :
Tokuharo dist is shown for each test file with each reference file and how many accurate vowels are predicted out of 10 for each vowel.
----------------------------------------------------------------

