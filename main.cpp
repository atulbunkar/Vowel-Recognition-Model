//Assignment 3 : Vowel Recognition
//Roll No: 214101011
#include<stdlib.h>
#include<stdio.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#pragma warning(disable:4996)

#define p 12
#define N 320
static double norm_data[50000];

int skip_line(FILE *fp) //  function to skip lines, i.e move cursor to next line !
{
    int c;

    while (c = fgetc(fp), c != '\n' && c != EOF);

    return c;
} 

double dc_fix(FILE* fp_dc_check){ //check for dc shift and fix it
	double dc_bias = 0.0 ; 
	char line_buffer[20] ;
	int Xi=0 ;
	int total_datapoints =0  ; 
	while(!feof(fp_dc_check)){                      // Summing over all datapoints then averaging to get dc_shift if its not zero.
		fgets(line_buffer, 20, fp_dc_check);
			Xi = atoi(line_buffer);              //Converting string got from file to integer.
			dc_bias += Xi;   
			total_datapoints++;              
	}
	dc_bias = dc_bias/total_datapoints;
	return dc_bias;
}

void get_samples(FILE* fp , double *samples ){  // Read 320 samples from test.txt for step 1
	char line_buffer[20];
	double Xi; int c = 0;
	double sum =0.0; 
	int j = 0 , k =0;
	
	while(!feof(fp)){          
		fgets(line_buffer, 20, fp);
		Xi = atof(line_buffer);
		if(Xi != 0.0){
			samples[c] = Xi;
			c++;
		}
	}
}

void calc_ri(double *ri , double *samples){ // calculates autocorrelation of data and stores in ri
	int j = 0 , k =0;
	double sum = 0.0;
	for(k; k <= p ; k++){
		sum = 0.0;
		for(j=0 ; j < N-k ; j++){
			sum += samples[j]*samples[j+k]*1.0; 
		}
		ri[k] = sum;
	}

}

void calc_ai(double *ri , double *ai){ // calculates ai's using durbin algo and stores in ai
	double e[p+1] , a[p+1][p+1] ,k[p+1] ;
	int i= 0 ,j = 0 ; double sum = 0.0;
	e[0] = ri[0];
	for(i=1 ; i <= p ; i++){
		sum = 0;
		for(j = 1 ; j <= i-1 ; j++){
			sum+= a[i-1][j]*ri[i-j]*1.0;
		}
		k[i] = (ri[i]-sum)*1.0/e[i-1];
		a[i][i] = k[i];
		for(j =1 ; j <= i-1 ;j++){
			a[i][j] = a[i-1][j] - k[i]*a[i-1][i-j]*1.0;
		}
		e[i] =(1-(k[i]*k[i])) * e[i-1] *1.0;
	}
	for(i=1;  i <= p ; i++){
		ai[i] = a[p][i];
	}

}

void calc_ci(double *ri ,double *ai , double *ci){  // calculates cepstral's coefficients  ci's and stores in ci
	double sigma = ri[0]; int m ,k ; double sum = 0; 
	
	ci[0] = log(sigma*sigma) ;
	for(m = 1; m <= p ; m++){
		sum = 0;
		for(k = 1; k <= m-1 ; k++){
			sum += (k*ci[k]*ai[m-k] )*1.0/m;
		}
		ci[m] = ai[m] + sum ; 
	}
}
 
void print(double *arr , int start ,int end){		// print array elements
	int x ;
	for(x =start; x <=end; x++)
		printf("%lf \t" , arr[x]);
	printf("\n");
}

void raised_sine(double *ci){ //applying the raised sine window on ci
	int m;
	double pi = 3.141 ;

	for(m = 1 ; m <= p ; m++){
		ci[m] = ci[m] * (1.0 +  (p/2.0)* sin( pi*m*1.0/p )) ;   
	}
	
}

int preprocess(FILE* fp_raw ,double dc_bias  ){
	
	int frame_size = 100 ; double norm_mult;
	double Xi=0 ;
	double peak_Xi=0 ; 
	char line_buffer[20] ; int new_peak = 5000 ;         // new_peak will be the peak signal after normalization.
	int i =0 ; double xi = 0.0, norm_xi = 0.0;

	int line_skips = 10; 
	while(line_skips--)skip_line(fp_raw);				// skipping 10 lines.

	// Normalization                   
	while(!feof(fp_raw)){ 
		fgets(line_buffer, 20, fp_raw);						 // Here, I found the max amplitude signal from the data.
			Xi = atof(line_buffer);
			if(abs(Xi) > peak_Xi ){
				peak_Xi = abs(Xi);
			}    
	 }
	fseek(fp_raw , 0 , SEEK_SET);				// Setting the file ptr back to 0 and skipping 10 lines again.
	line_skips = 10 ;
	while(line_skips--){
		skip_line(fp_raw) ; 
	}
	norm_mult = new_peak*1.0/abs(peak_Xi) ;  // This gives us the normalization factor. 
	
	while(!feof(fp_raw)){                        // I am transforming the raw data here by subtracting from each point the dc_bias , then multiplying them with 
		fgets(line_buffer, 20, fp_raw);			 // the normalization factor.
		xi = atof(line_buffer);
		norm_xi = ((xi-dc_bias) * norm_mult*1.0 );      
		if(norm_xi != 0)                                
			norm_data[i++] = norm_xi;        // Storing the normalised data in array "norm_data".
		}
	return i;								// returns the number of data points present in normalized data !
}

int find_start(int size){ //finds stable frame 
	int start = 0;
	int st = 0;      
	double cur_Xi = 0, energy = 0;
	double max_energy = 0 ;

	while(start < size){                                // This loops the normalised data.
		cur_Xi =  norm_data[start++];

		energy += cur_Xi*cur_Xi;                     // calculating energy for current frame
		
		if(start%N == 0){							  // check end of frame.
			energy = energy*1.0/N ;					// dividing by num of data .  														      
			if(energy > max_energy){			//found max energy frame
				st = start-N ;					//start is at end of frame so st has start of frame sample point's index.
				max_energy = energy;		
			}
		}
	}
	st = st- 2*N ;		//taking 2 frame behind max energy frame as starting frame
	return st ;	

}

void genrate_reference(FILE *fp_dump, const char* vowel_start_path , double dc_bias ){
	double total_ci[11][5][p+1] = {0}; // stores ci for 10 files , 5 frames each having 12 ci values.
	double avg_ci[5][p+1] = {0} ;	// after averaging ci is stored here.
	double ri[p+1] ;   
	double ai[p+1] ;					// for storing single frame ai,ci,ri's
	double ci[p+1];
	double samples[N+1];			//stores 320 samples , i.e one frame.
	char path[40] ;				//total path of raw file
	char filenum[3] ;		   //store file num
	int i , j,frame , x =0 ;
	
	for(i =1 ; i <= 10 ; i++){				// looping over 10 files 1-10 of vowel
	
		strcpy(path , vowel_start_path );   //vowel_start_path stores start path of raw vowels
		sprintf(filenum ,"%d", i);
		strcat(path, filenum) ;
		strcat(path, ".txt") ;
		//final path looks like : "./vowels/214101011_a_i.txt" for vowel A where i is replaced by number i;

		FILE* fpp;
		fpp = fopen(path , "r");
		if (fpp == NULL){						 // checking if file has opened successfully r not.
			printf("error opening file - %d\n",i);
		}
		int size = preprocess(fpp , dc_bias);	// preprocesses the file and stores in norm_data , returns size(index) till which values are stored.
		fclose(fpp);

		int start = find_start(size);			// finds the datapoint of stable vowel from where we can analyse 5 frames.
	
		for(frame = 0; frame < 5 ; frame++){	// runs for 5 frames
			for(j =0 ; j < N; j++){				//gives one frame 
				samples[j] = norm_data[start++];	//samples store one frame samples, then for next frame stores next 320 as start value is not changed.
			}
			calc_ri(ri,samples);			//calc ri for one frame
			calc_ai(ri,ai);					//calc ai for one frame
			
			calc_ci(ri,ai,ci);			// calculating ci
			raised_sine(ci);			// applying raised sine on ci

			for(x=1; x <= p ; x++){				//filling total_ci for i'th file's current frame 
				total_ci[i-1][frame][x] = ci[x];	// so for i => 1 to 10 ,fills 0 to 9 , and fills frame 0 to 4 and fills ci's 1 to 12;
			}
		}
	}

	// looping over same frame of all 10 files and storing avg of their ci's
	for(frame = 0; frame<5; frame++){		//loop for each frame
		for(x = 1 ; x <= 12 ; x++){			//loop over ci's
			for(i =0 ; i <=9 ;i++){			//loop files
				avg_ci[frame][x] += total_ci[i][frame][x];

			}
			avg_ci[frame][x] /= 10.0;
		}	
	}
	//dumping avg_ci in reference folder.
	if (fp_dump == NULL){						 
		printf("error opening dump file \n");
		return ;
	}

	for(i = 0; i<5; i++){		
		for(x =1 ; x <=p ;x++) {
			fprintf(fp_dump,"%lf\n",avg_ci[i][x]); //stores single value in each line
		}										  //total lines 12x5 = 60.
	}
}

void genrate_test( const char* dump_test_path,  const char* vowel_start_path , double dc_bias ) {
	double ri[p+1] ;   
	double ai[p+1] ;					// for storing single frame ai,ci,ri's
	double ci[p+1];
	double final_ci[5][p+1];		//storing all calculated ci's
	double samples[N+1];			 //stores 320 samples , i.e one frame.
	char vowel_path[40] ;			//store total path of raw file
	char dump_path[40];				//store total path of dump file
	char vow_num[3] ;				//store file num
	char test_num[3] ;	
	int i , j,frame , x =0 ;
	
	for(i =1 ; i <= 10 ; i++){
		//create path to get vowel i+10 ,i.e 11 ,12 ,13 ...
		strcpy(vowel_path , vowel_start_path );   
		sprintf(vow_num ,"%d", i+10);
		strcat(vowel_path, vow_num) ;
		strcat(vowel_path, ".txt") ;	
		
		//create path to dump file numbered i ,i,e 1,2,3 ... 
		strcpy(dump_path , dump_test_path);
		sprintf(test_num ,"%d", i);
		strcat(dump_path, test_num) ;
		strcat(dump_path, ".txt") ;

		FILE* fpp;		// for reading vowel
		FILE* fp_dump;	//for dumping test file generated

		fpp = fopen(vowel_path , "r");
		if (fpp == NULL){						 // checking if file has opened successfully r not.
			printf("error opening file - %d\n",i+1);
			return ;
		}
		int size = preprocess(fpp , dc_bias);  // calls preprocess , store in norm_data and return size(index) till which values are stored.
		fclose(fpp);


		int start = find_start(size);			// finds the datapoint of stable vowel from where we can analyse 5 frames.
	
		for(frame = 0; frame < 5 ; frame++){	// runs for 5 frames
			for(j =0 ; j < N; j++){				//gives one frame 
				samples[j] = norm_data[start++];	//samples store one frame samples, then for next frame stores next 320 as start value is not changed.
			}

			calc_ri(ri,samples);			//calculate ri for one frame
			calc_ai(ri,ai);					//calculate ai for one frame
			
			calc_ci(ri,ai,ci);			 // calculating ci
			raised_sine(ci);			// applying raised sine on ci

			for(x=1; x <= p ; x++){				//filling final_ci for i'th file's current frame 
				final_ci[frame][x] = ci[x];		// so for i => 1 to 10 ,fills 0 to 9 , and fills frame 0 to 4 and fills ci's 1 to 12;
			}
		}
		//dumping in text file.
		fp_dump = fopen(dump_path , "w");
		if(fp_dump == NULL){					
			printf("error opening dump file - %d\n",i);
			return ;
		}

		for(j = 0; j<5; j++){		
			for(x =1 ; x <=p ;x++) {
				fprintf(fp_dump,"%lf\n",final_ci[j][x]);  //stores single value in each line
			}	//total lines 12x5 = 60.
		}
		//NOW WE HAVE 10 FILES for EACH Vowel WITH 12*5 CIs. 
		fclose(fp_dump);
	}	
}

void vowel_recognition(char* test_path , char vowels[5] ){ //test_vowel_ind is index of cur vowel testing . 0 for a, 1 for e , ...4 for u.
	char vow[3] ;		
	char test_num[3] ;
	char ref_path[40] ;      
	char full_test_path[40]; //stores full test path
	
	double ci_ref=0 , ci_test=0 ,sum =0,diff=0,tok_dist=0 , min_td = 99999999999; 
	int correct_pred = 0 , vow_ind=0 ;
	char line_buffer[20] ,  line_buffer2[20] ; 
	int i=0,j=0,ref=0,test =1 ;
	double tok_wt[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0}; //tokuhara weights

	FILE* fp_test ; FILE* fp_ref;
	printf("\n\n--------------------------------------------------\n\n");
	printf("PREDICTION FOR VOWEL '%c' : \n\n",test_path[12]);

	for(test = 0 ; test < 10 ; test++){				//read 10 test files in loop and check it's ci with every refernce vowel ci
		strcpy(full_test_path , test_path);
		sprintf(test_num ,"%d", test+1);
		strcat(full_test_path, test_num) ;			
		strcat(full_test_path, ".txt") ;			// Generate full test path
		//printf("%s\n",full_test_path );
		printf("\n\ntesting For file %s : \n\n", full_test_path);	

		fp_test = fopen(full_test_path , "r");
		if(!fp_test)printf("Cannot open test %d",test);
											
		min_td = 99999999999;
		for(ref =0 ; ref<5 ; ref++ ){						//loop through 5 refernce files
			sprintf(vow , "%c" ,vowels[ref]);			//get vowel char wrt ref number. a for 0 ,  e for 1 .... u for 4;
			strcpy(ref_path ,"./reference/ref_" );
			strcat(ref_path, vow);
			strcat(ref_path, ".txt") ;					//generate full ref path	

			fp_ref = fopen(ref_path , "r");
			if(!fp_ref)printf("Cannot open ref %d", ref);

			for(i = 0 ; i < 5 ; i++){						//loop 5 frames
				for(j = 0 ; j < 12 ; j++){					//loop 12 entry
					fgets(line_buffer, 20, fp_test);		//read ci_test from test file
					ci_test = atof(line_buffer);		
					
					fgets(line_buffer2, 20, fp_ref);	//read ci_ref from reference file
					ci_ref = atof(line_buffer2);
					diff =  ci_ref - ci_test ;
					sum += tok_wt[j] *(diff*diff);		// getting tokuharo distance for all 12 ci and summing it up for "all 5 frames" in sum.
				}
			}
			tok_dist = sum*1.0/5;						//averaging over 5 frames to get final tokuharo distance
			sum = 0;
			printf("tokuhara distance for test vowel %c with reference vowel %c is: %lf\n\n" ,test_path[12], vowels[ref], tok_dist  );
			if(tok_dist < min_td){ //store current vowel index if min distance found	
				min_td = tok_dist;
				vow_ind = ref;					
			}
			fclose(fp_ref);
			fseek(fp_test , 0 , SEEK_SET); //reset test file pointer to start of file.
		}

		if(test_path[12] == vowels[vow_ind])correct_pred++;
		//test_path[12] is test vowel being checked(got from parameter path),  and vowels[vow_ind] is vowel with min distance found.
		//if both equal then correct predicted for current file.
		
		printf("Actual Vowel : %c , Predicted Vowel : %c \n\n", test_path[12] , vowels[vow_ind] );
		fclose(fp_test);
	}
	correct_pred = (correct_pred/10) * 100 ;	//Getting accuracy percentage
	printf("Accuracy of predicting vowel %c is : %d percent \n\n",test_path[12] , correct_pred); //Showing accuracy of vowel prediction for 10 files
	correct_pred = 0;

}//end of VR


int main(){
	printf("\t\tAssignment 3 : Vowel Recognition ! \n \t\tSubmitted By : Atul Bunkar , Roll No: 214101011\n");

	//STEP 1 : Prerequisite task : testing from sample , if calculated ai , ri , ci match.
	printf("\n STEP 1 : Prerequisite task :- \n");

	FILE *fp ; 
	double samples[N+1] ;	//store 320 samples
	double sample_ri[p+1] ; //store 13 ri
	double sample_ai[p+1] ;	//store 12 ai
	double sample_ci[p+1];  //store 12 ci
	int x =0;

	fp = fopen("test.txt" , "r");
	if (fp == NULL) {						 // checking if file has opened successfully r not.
      fprintf(stderr, "error opening %s: %s", "test.txt", strerror(errno));
	}
	
	get_samples(fp, samples);					// gets samples from sample test file and stores in array samples[]
	fclose(fp);

	calc_ri(sample_ri , samples);				//calls function calc_ri to calculate ri from sample data, total 13 ri, stores in array sample_ri .
	printf("\nCalculated ri from sample data : \n\n");
	print(sample_ri , 0, p);					// to print array elements in console. You can check we get the correct output. 

	printf("\n\nCalculated ai from  sample data: \n\n");
	calc_ai(sample_ri , sample_ai);				   // calls function calc_ai to calculate ai of sample data, total 12 ai , stores in array sample_ai .
	print(sample_ai,1 ,p);

	printf("\n\nCalculated ci from sample data: \n\n");
	calc_ci(sample_ri , sample_ai, sample_ci );		// calls function calc_ci to calculate ci of sample data, total 12 ci , stores in array sample_ci .
	print(sample_ci ,1,p);
	printf("\n\n");
	
	//-------------------------------------------------------------------------------------------------------------
//	STEP2 : Generating reference files

	printf("\n ----------------------------------------------------------------------------------------\n ");
	printf("\nSTEP 2 - Reference and Test Generation , check the folders for this step where the test and reference files are generated !!!\n\n");

	double dc_bias = 0;
	int k =0 , i = 0 ,j =0 ; 
	//Get dc shift ,need to get this only once. Will be used in preprocess.
	FILE* fp_dc ;
	fp_dc = fopen("DCshift.txt", "r") ;
	dc_bias = dc_fix(fp_dc) ;
	fclose(fp_dc);

	FILE* fp_dump;
	char vowel_start_path[30] = "./vowels/214101011_a_";      // Start path to raw vowel txt , later will concat with filenum
	char ci_dump[30] = "./reference/ref_a.txt";                //Full path to dump ci's of vowel 
	fp_dump = fopen(ci_dump , "w");								//file to dump refernce ci.
	genrate_reference(fp_dump, vowel_start_path , dc_bias );	//args  passed are : file ptr to dump in , first half of path to raw vowel files , calculated dc_vias.
	//calls generate_reference
	fclose(fp_dump);


	//Now do same for vowel E , here path of vowel e is passed and file to dump in it is passed.
	strcpy(vowel_start_path,"./vowels/214101011_e_");	//replces path of vowel a with vowel e
	strcpy(ci_dump , "./reference/ref_e.txt");			//replaces dump file path.
	fp_dump = fopen(ci_dump , "w");								
	genrate_reference(fp_dump, vowel_start_path , dc_bias );	
	fclose(fp_dump);

	//For Vowel I
	strcpy(ci_dump , "./reference/ref_i.txt");
	strcpy(vowel_start_path,"./vowels/214101011_i_");
	fp_dump = fopen(ci_dump , "w");								
	genrate_reference(fp_dump, vowel_start_path , dc_bias );	
	fclose(fp_dump);

	//For vowel O
	strcpy(ci_dump , "./reference/ref_o.txt");
	strcpy(vowel_start_path,"./vowels/214101011_o_");
	fp_dump = fopen(ci_dump , "w");								
	genrate_reference(fp_dump, vowel_start_path , dc_bias );	
	fclose(fp_dump);

	//For vowel U
	strcpy(vowel_start_path,"./vowels/214101011_u_");
	strcpy(ci_dump , "./reference/ref_u.txt");
	fp_dump = fopen(ci_dump , "w");								
	genrate_reference(fp_dump, vowel_start_path , dc_bias );	
	fclose(fp_dump);


	//------------------------------------------------
	// Now to generate test files
	
	//generate 10 test files for vowel a_11 to a_20
	strcpy(vowel_start_path,"./vowels/214101011_a_");			
	char dump_test_path[30] = "./test/test_a_";					// start path of dump file to dump test ci's .
	genrate_test(dump_test_path,  vowel_start_path , dc_bias ); //args passed are : first half of test file path (rest will be filled in function) , 
																					// first half of raw vowel path and dc_bias.

	//Same for Vowel E
	strcpy(vowel_start_path,"./vowels/214101011_e_");			
	strcpy(dump_test_path, "./test/test_e_");  
	genrate_test(dump_test_path,  vowel_start_path , dc_bias );

	//For Vowel I
	strcpy(vowel_start_path,"./vowels/214101011_i_");
	strcpy(dump_test_path, "./test/test_i_"); 
	genrate_test(dump_test_path,  vowel_start_path , dc_bias );

	//For Vowel O
	strcpy(vowel_start_path,"./vowels/214101011_o_");
	strcpy(dump_test_path, "./test/test_o_"); 
	genrate_test(dump_test_path,  vowel_start_path , dc_bias );

	//For Vowel U
	strcpy(vowel_start_path,"./vowels/214101011_u_");
	strcpy(dump_test_path, "./test/test_u_"); 
	genrate_test(dump_test_path,  vowel_start_path , dc_bias );
	
	printf("\n\n --- Successfully Created Test and Reference Files--- \n\n");


	//-------------------------------------------------------------------------------
	//STEP 3 : Prediction/Testing 

	printf("\n ----------------------------------------------------------------------------------------\n ");
	printf("\n\n STEP 3 - Testing test files against reference Files !!!\n\n");

	char vowels[5] = {'a', 'e' , 'i' , 'o' , 'u'};	//will be used to check if correct prediction is made by program

	//checking for vowel A 
	char test_path[30] = "./test/test_a_";	//we are now testing 10 test files of a , test path has first half path , to this we append 1,2,3.. in the function to open respective file.
	//13th char in the test_path (test_path[12]) is the vowel currently tested ,will use this to check
	//if correct vowel is predicted and check accuracy of program.
	vowel_recognition(test_path , vowels);  // args passed are : first half of test file , vowels[5] array
					

	//checking for vowel E
	strcpy(test_path ,"./test/test_e_") ;	
										
	vowel_recognition(test_path , vowels); 


	//checking for vowel I 
	strcpy(test_path ,"./test/test_i_") ;											
	vowel_recognition(test_path , vowels); 


	//checking for vowel O 
	strcpy(test_path ,"./test/test_o_") ;											
	vowel_recognition(test_path , vowels); 


	//checking for vowel U 
	strcpy(test_path ,"./test/test_u_") ;											
	vowel_recognition(test_path , vowels); 
	//----------------------------------------------------------
	//Check output in console for step 3.
}

//END OF ASSIGNMENT . THANK YOU :) 
