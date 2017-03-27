
/******************************************************************************/
/* TCidd.c - source code for Turbo code simulations with QAM encoding         */
/******************************************************************************/

/******************************************************************************/
/* QAM signaling (modulation) is used.                                        */
/* The modulated data is transmitted over an AWGN channel.                    */
/******************************************************************************/

/**********************/
/* Include Statements */
/**********************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/********************/
/* Define Constants */
/********************/

/* Simulation Parameters */

#define Num_QAMbits	1	/* Number of Bits to represent QAM point */
#define QAM_Size        2	/* Number of points in QAM */
#define M               4	/* Number of memory elements in encoder */
#define Num_of_States   16	/* Number of states = 2^M in trellis */
#define Eb_Over_N0_Min	4	/* Min/Max values for Channel SNR (dB) */
#define Eb_Over_N0_Max	8.0
#define Num_Points	9	/* Number of points for the error curve */
#define Max_TC_Iter	7	/* Number of iterations for turbo decoding */
#define Max_FB_Iter	1	/* Number of iterations for feedback decoding */
#define Sim_Repeat	1	/* Number of times to repeat simulation */
#define Num_of_Blocks   10	/* Number of blocks to randomize */
#define Block_Size    2048	/* Number of symbols per block */
#define Symbol_Size     1	/* Number of bits per symbol (must be 1) */
#define Ortho_Size	2	/* 2^(Symbol_Size) */
#define Puncture_Size   16	/* Size of puncture array */
#define Repetition		2   /*Repetition parameter: 1 for no repetition, 2 for duplication */
#define LLR_summation   1  /* zero for No LLR summation, 1 for LLR summation */

const char inter_file_name[]  = "Int4096s.dat";	 /* Interleaver file name */
/* Interleaver for rate 2/3 code */
/* const char inter2_file_name[] = "Int6144s.dat"; */
/* Interleaver for rate 1/2 code */
const char inter2_file_name[] = "Int8192s.dat";

/* Interleaver for rate 1/3 code */
/* const char inter2_file_name[] = "Int12288s.dat"; */

const char stat_file_name[]   = "status-test3.dat";  /* Status file name      */
const char test_file_name[]   = "results-test3.dat"; /* Results file name     */
const char log_file_name[]    = "log-test3.dat";     /* Intermediate log file */
const char QAM_file_name[]    = "QAM16.dat";   /* QAM code file name	*/


/* Array Sizes */

#define AS0  (Block_Size*Repetition+M+1)               /* Block_Size items */
#define AS02 (AS0*3+Num_QAMbits)            /* Block_Size x max Code_Rate */
#define AS1  (AS0*Symbol_Size)              /* Block_Size x Symbol_Size */
#define AS4  (AS0*Ortho_Size)               /* Block_Size x Ortho_Size */
#define AS6  ((Num_of_Blocks+1)*Block_Size*Symbol_Size) /* random info bits */

/********************************/
/* Global Variable Declarations */
/********************************/

int gfeedback[7]={1,0,0,1,1};           /* Convolutional Encoder */
int gforward[7]={1,1,1,0,1};            /* Transfer Functions    */
                                        /* Use generators 37, 21 */


/* Both gforward and gfeedback should have M+1 elements. */
/* Note: In practical codes gforward starts with {1,...} */
/* Note: gfeedback usually starts with {0,...}           */

/* Data variables are declared as type "char" instead of the "int" */
/* because we are dealing with binary random variables             */

char data_org[AS6];       /* Input bits, random generated data  */
char data[AS6*Repetition];	/*Output of Duplicate function */
char source_bit[AS1]; /* Input bits partitioned into blocks */

double QAM_Array[QAM_Size][2];	/* Array for QAM distribution */
double eb; /* Average bit energy per transmitted bit for QAM */
double Code_Rate; /* Code rate of the transmitted code */

/* Arrays specifying puncturing pattern */
/* Puncture arrays for Code_Rate = 2/3 */
/* int puncture1[Puncture_Size] = {0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0};
int puncture2[Puncture_Size] = {0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1}; */

/* Puncture arrays for Code_Rate = 1/2 */
int puncture1[Puncture_Size] = {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};
int puncture2[Puncture_Size] = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};

/* Puncture arrays for Code_Rate = 1/3 */
/* int puncture1[Puncture_Size] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
int puncture2[Puncture_Size] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; */

int interleaver[AS0];   /* Array containing interleaver structure */
int interleaver2[AS02];	/* Array containing second interleaver structure */
char interleaved_bit[AS1]; /* Array containing interleaved information bits */
char datastream_int[AS02]; /* Array containing interleaved stream bits */

int extra_bits; /* Number of extra zeroes to add to datastream for QAM */

/* output of the encoders */
char encoder1_bit[AS1], encoder2_bit[AS1];
/* output of the decoders */
char decoder1_bit[AS1], decoder2_bit[AS1];
char detected_data_dup[AS1];

/* X and Y coordinates for QAM points */
double Xstream[AS02];
double Ystream[AS02];

/* x and y coordinates after adding AWGN */
double xstream[AS02];
double ystream[AS02];

/* Probability for exp(-x^2/N0) */
double QAM_exp[AS02][QAM_Size];

/* Received bits after adding noise */
double received_sourceprob_bit[AS0][2];
double received_sourceprob_int[AS0][2];
double received_encoder1prob_bit[AS0][2];
double received_encoder2prob_bit[AS0][2];

/* Probability distribution of decoded bits for each user */
double syst_prob[AS1][Ortho_Size];
double par1_prob[AS1][Ortho_Size];
double par2_prob[AS1][Ortho_Size];

/* Difference of LLR's in and out of constellation function for info bits */
double Syst_LLR[AS0][2];

/* error pattern of the two decoders */
char error1[AS1], error2[AS1];

/* Structure of trellis and trellis tail */
int trellis[Num_of_States][Ortho_Size][2];
int trellis_tail[Num_of_States][Ortho_Size][2];

/* Binary decomposition of trellis input/output values */
int symbol_bit[Ortho_Size][Symbol_Size+1];

/* Extrinsic information for FB iteration of decoder */
double extrinsic[AS0][Ortho_Size];

double N0;				/* Noise power coefficent */

int error_sum[Max_TC_Iter+1][3];        /* error sums */

/* File variables for input and output files */
FILE *outfile;
FILE *STATUS_File;

/*************************/
/* Function Declarations */
/*************************/

void Initialize (void);			/* Initialization */
void QAM_Init (void);			/* Builds QAM distribution */
void Trellis (void);			/* Builds the trellis diagram of code */
void Calc_Sym_Bits (void);		/* Calculates symbol_bit[i][j] values */
void Generate (void);			/* Generates random input data */
void Interleave (void);			/* Interleaves source symbols */
void Encode (char encoder_num);         /* RSC encoders */
void Interleave2 (void);		/* Interleave source and encoded bits */
void Modulate (void);			/* Orthogonal Modulation */
void Add_Noise (void);			/* Correlates noise and adds to bits */
void Decode_Init (void);		/* Initialization for decoder */
void QAM_dist (void);			/* Calculates QAM probability dist */
void Turbo_Decode (void);		/* Turbo Code Decoder */
void Duplicate (void);			/* Duplicates generated data */
void Interleave_dup(void);		/* Interleave the duplicated data */
void Deinterleave_dup(void);	/* Deinterleave the duplicated data */

/******************************************************************************/
/* main ()                                                                    */
/*      Controls all initialization                                           */
/*      Controls Main Loop of program - communication system simulation       */
/******************************************************************************/

main ()
{
  int i,j,k,l,im,iter,iteration, temp_count;
  int LOOP;
  int error[Max_FB_Iter+1][Max_TC_Iter+1][3]; 
  double error_pr[Max_FB_Iter+1][Max_TC_Iter+1][3];
  double eb_over_no, delta_snr;

  /* Global Initialization Functions */

  /* Calculate Code_Rate */
  temp_count = Puncture_Size;
  for (i=0; i< Puncture_Size; i++) /*YES*/
  {
    /*printf("first for  %d \n ",i);*/
    if (puncture1[i] == 1)
      temp_count += 1;
    if (puncture2[i] == 1)
      temp_count += 1;
  }
  Code_Rate = temp_count / (double) Puncture_Size;
/*printf("This is an updated version \n");*/
printf(" code rate : %d ", Code_Rate);
printf(" temp_count : %d ", temp_count);
printf(" Puncture_Size : %d ", Puncture_Size);

  Initialize ();	/* Initialize certain variables and files */
  QAM_Init ();		/* Build the QAM distribution */
  Trellis ();		/* Builds the trellis diagram of the code  */
  Calc_Sym_Bits ();	/* Calculates symbol_bit[i][j] values */

  /* Initialize Channel SNR Variables */

  if (Num_Points != 1)
    delta_snr = (Eb_Over_N0_Max - Eb_Over_N0_Min)/(double)(Num_Points - 1);
  else
    delta_snr = 0;

  eb_over_no = Eb_Over_N0_Min - delta_snr;
printf(" numPoints : %d \n", Num_Points);
  /* Loop for each channel SNR value */
  /* Used if want to cover a range of SNR values */
  for (LOOP=1; LOOP<=Num_Points; LOOP++)
  {
 /*printf(" Loop in Num_Points - 2nd for : %d \n", LOOP);*/
    eb_over_no += delta_snr;
    N0 = eb*(pow(10,(-(eb_over_no)/10)));
    outfile = fopen(&test_file_name[0], "a");

    fprintf(outfile,"\n10*log10(Eb/N0) = %f dB,  Eb = %f,  N0 = %f\n",
        eb_over_no, eb, N0);
    fclose(outfile);
    outfile = fopen(&log_file_name[0], "w");
    fprintf(outfile,"\n10*log10(Eb/N0) = %f dB,  Eb = %f,  N0 = %f\n",
        eb_over_no, eb, N0);
    fclose(outfile);

    for (iter=1; iter<=Max_FB_Iter; iter++)
      for (iteration=1; iteration<=Max_TC_Iter; iteration++)
      {

        /* Number of errors in decoder 1 and 2 in */
        /* a given iteration of FB and TC for each user */
        error[iter][iteration][1]=0;
        error[iter][iteration][2]=0;
        /* Overall prob. of bit error for decoder 1 and decoder 2 */
        error_pr[iter][iteration][1]=0;
        error_pr[iter][iteration][2]=0;
      }

    /********************* MAIN LOOP ***********************/

    for (im=1; im<=Sim_Repeat; im++)
    {
      //printf("%d\n",im);
      STATUS_File = fopen(&stat_file_name[0], "a");
      fprintf(STATUS_File,"%d\n",im);
      fclose(STATUS_File);
      /* Generate the random input data */
      Generate ();
      Duplicate();

      for (i=1; i<=Num_of_Blocks; i++)
      {

        /* Load current block of data */
        for (j=1; j<=Block_Size * Repetition; j++)
        {
          l = Block_Size * Repetition*(i-1)+j;
          source_bit[j] = data[l];
        }
        Interleave_dup();
        Encode (1);
        Interleave ();
        Encode (2);

        Interleave2 ();

        Modulate ();

        Add_Noise ();

      Decode_Init ();

        for (iter=1; iter<=Max_FB_Iter; iter++)
        {

          QAM_dist ();
          /* Do not add LLR to L2 of TD for first iteration */
         if (iter ==1)
            for (j=1; j<Block_Size * Repetition+M; j++)
            {
              Syst_LLR[j][0] = 0;
              Syst_LLR[j][1] = 0;
            }
          Turbo_Decode ();
          Deinterleave_dup();

          for (iteration=1; iteration<=Max_TC_Iter; iteration++)
          {
            error[iter][iteration][1] += error_sum[iteration][1];
            error[iter][iteration][2] += error_sum[iteration][2];
          }
        }
      }

      /* Output Intermediate Statistics to Log File */
      outfile = fopen(&log_file_name[0], "a+");
      for (iter=1; iter<=Max_FB_Iter; iter++)
        for (iteration=1; iteration <= Max_TC_Iter; iteration++)
        {
          error_pr[iter][iteration][1]=(double)error[iter][iteration][1]/
		(im*Num_of_Blocks*Block_Size * Repetition*Symbol_Size);
          error_pr[iter][iteration][2]=(double)error[iter][iteration][2]/
		(im*Num_of_Blocks*Block_Size * Repetition*Symbol_Size);
        }

      fprintf(outfile,"\nNumber of bits = %d\n",
			im*Num_of_Blocks*Block_Size * Repetition*Symbol_Size);
      fprintf(outfile,"Prob. of Error for each FB and TC Iteration");
      fprintf(outfile,": -log10 values");
      fprintf(outfile,"\n\nFB_Iter  TC_Iter   Dec_1           Dec_2\n");

      for (iter=1; iter<=Max_FB_Iter; iter++)
        for (iteration=1; iteration <= Max_TC_Iter; iteration++)
          fprintf(outfile,"%d\t%d\t%e\t%e\n",
			iter, iteration,
			-log10(error_pr[iter][iteration][1]),
			-log10(error_pr[iter][iteration][2]));	

      fclose(outfile);
    }

    /* Calculate Error Statistics */
    for (iter=1; iter<=Max_FB_Iter; iter++)
      for (iteration=1; iteration <= Max_TC_Iter; iteration++)
      {
        error_pr[iter][iteration][1]=(double)error[iter][iteration][1]/
		(Sim_Repeat*Num_of_Blocks*Block_Size * Repetition*Symbol_Size);
        error_pr[iter][iteration][2]=(double)error[iter][iteration][2]/
		(Sim_Repeat*Num_of_Blocks*Block_Size * Repetition*Symbol_Size);
        printf("%d \t%d \n", error[iter][iteration][1],
			error[iter][iteration][2]);
      }
    /* Output Error Statistic to File */
    outfile = fopen(&test_file_name[0], "a+");

    fprintf(outfile,"Prob. of Error for each FB and TC Iteration");
    fprintf(outfile,": -log10 values");
    fprintf(outfile,"\n\nFB_Iter  TC_Iter   Dec_1           Dec_2\n");

    for (iter=1; iter<=Max_FB_Iter; iter++)
      for (iteration=1; iteration <= Max_TC_Iter; iteration++)
        fprintf(outfile,"%d\t%d\t%e\t%e\n",
		iter, iteration,
		-log10(error_pr[iter][iteration][1]),
		-log10(error_pr[iter][iteration][2]));

    fclose(outfile);
  }

  exit(0);
}

/******************************** FUNCTIONS ***********************************/

/******************************************************************************/
/* Initialize (void)                                                          */
/*      Initializes variables                                                 */
/******************************************************************************/

void Initialize (void)
{
  int i, j, location1, location2;
  FILE *interleavefile;

  /* Initialize data_org */
 
  for (j=0; j<AS6; j++)
    data_org[j] = 0;
       
  /* Initialize source bits */

  for(j=1; j<=Block_Size * Repetition; j++)
    source_bit[j] = 0;

  /* Initialize Interleaver */

  interleavefile = fopen(&inter_file_name[0], "r");
  for (i=1; i<=Block_Size * Repetition; i++)
  {
    fscanf(interleavefile, "%d %d", &location1, &location2);
    interleaver[location1] = location2;
  }
  fclose(interleavefile);

  /* Initialize Second Interleaver */

  interleavefile = fopen(&inter2_file_name[0], "r");
  for (i=1; i<=Block_Size * Repetition*Code_Rate; i++)
  {
    fscanf(interleavefile, "%d %d", &location1, &location2);
    interleaver2[location1] = location2;
  }
  fclose(interleavefile);

  /* File Initialization */

  STATUS_File = fopen(&stat_file_name[0], "w");
  fclose(STATUS_File);

  outfile = fopen(&test_file_name[0], "w");
  fprintf(outfile,"\nBlock_Size %d, Symbol_Size %d, ",
        Block_Size, Symbol_Size);
  fprintf(outfile,"Rate = %f, ", (1/Code_Rate));
  fprintf(outfile,"M = %d, Bits %d\n", M,
        (Block_Size*Symbol_Size*Num_of_Blocks*Sim_Repeat));
  fclose(outfile);

  outfile = fopen(&log_file_name[0], "w");
  fprintf(outfile,"\nBlock_Size %d, Symbol_Size %d, ",
        Block_Size, Symbol_Size);
  fprintf(outfile,"Rate = %f, ", (1/Code_Rate));
  fprintf(outfile,"M = %d, Bits %d\n", M,
        (Block_Size*Symbol_Size*Num_of_Blocks*Sim_Repeat));
  fprintf(outfile,"\nIntermediate Error Statistics\n");
  fclose(outfile);

  printf("\nBlock_Size %d, Symbol_Size %d, ",
        Block_Size, Symbol_Size);
  printf("Rate = %f, ", (1/Code_Rate));
  printf("M = %d, Bits %d\n", M,
	(Block_Size*Symbol_Size*Num_of_Blocks*Sim_Repeat));
}

/******************************************************************************/
/* QAM_Init (void)                                                            */
/*      Initialize QAM distribution                                           */
/******************************************************************************/

void QAM_Init (void)
{
  int i,j;
  double temp, E_s;

  /* Initialize Array */
  for (i=0; i<QAM_Size; i++)
    for (j=0; j<=1; j++)
      QAM_Array[i][j] = 0;



  /* Load Points (x,y) */

    QAM_Array[0][0] = -1;
    QAM_Array[1][0] = 1;


  temp = 0;
  for (i=0; i<QAM_Size; i++)
  {
    temp += (pow(QAM_Array[i][0],2)+pow(QAM_Array[i][1],2));
  }
  E_s = temp / (double)QAM_Size;  /* Calculate average symbol energy */

  /* Calculate average energy per bit, E_s*Code_Rate = energy per info bit */
  eb = E_s * Repetition * Code_Rate / (double)Num_QAMbits;
}

/******************************************************************************/
/* Trellis (void)                                                             */
/*      Generates Trellis and Trellis tail structure                          */
/*                                                                            */
/* prev_state=(a[1] a[2] .... )                                               */
/* current_state=(a[0] a[1] .... )                                            */
/* The trellis diagram is represented by a matrix 2^M *2^M matrix,            */
/* trellis[prev_state][input][0] = current_state                              */
/* trellis[prev_state][input][1] = output value                               */
/* (no connection is represented by -1)                                       */
/******************************************************************************/

void Trellis (void)
{
  int prev_state,current_state,state;
  int a[6];
  int sum_forward, sum_feedback;
        /* sum_forward is the sum of the forward path excluding a[0] */
        /* sum_feedback is the sum of the feedback path              */
  int i,j,k,Tmp_int;
  int input, in_bit, output;
  int step_no, mask, shift;

  /* Assuming at most one transition between states */
  /* trellis[mp][i][0] = m      */
  /* trellis[mp][i][1] = output */

  /* Initialize Trellis and Trellis tail to no paths */

  for (i=0; i<Num_of_States; i++)
  {
    for (j=0;j<Ortho_Size;j++)
    {
      for (k=0; k<2; k++)
      {
         trellis[i][j][k]=-1;
         trellis_tail[i][j][k]=-1;
      }
    }
  }

  for (prev_state=0; prev_state<Num_of_States;prev_state++)
  {
    for (input=0;input<Ortho_Size;input++)
    {
      state = prev_state;
      step_no = 0; output = 0;
      while(step_no < Symbol_Size)
      {
        /* Constructs a[1],a[2],... from prev_state=(a[1],a[2],...,a[M]) */

        Tmp_int = state;
        for(j=M;j>=1;j--)
        {
          a[j] = Tmp_int & 1;
          Tmp_int = Tmp_int >> 1;
        }

        sum_forward = 0;
        sum_feedback = 0;
        for (j=1; j<=M; j++)
        {
          sum_forward = sum_forward ^ (gforward[j] & a[j]);
          sum_feedback = sum_feedback ^ (gfeedback[j] & a[j]);
        }

        /* Constructs current_state=(a[0], a[1],...,a[M-1]) */

        shift = Symbol_Size - step_no -1;
        mask = (int)pow(2,shift);
        in_bit = (input & mask) >> shift;
        a[0] = in_bit ^ sum_feedback;
        output=(2*output)+((gforward[0] & a[0]) ^ sum_forward);

        Tmp_int = 1;
        current_state=0;
        for(j=M-1;j>=0;j--)
        {
          current_state += (a[j]*Tmp_int);
          Tmp_int *= 2;
        }
        state = current_state;
        step_no++;
      }
      trellis[prev_state][input][0] = current_state;
      trellis[prev_state][input][1] = output;
    }
  }

  /* Build the trellis tail */

  for (prev_state=0; prev_state<=Num_of_States-1; prev_state++)
  {
    state=prev_state;
    step_no=0; input = 0; output = 0;
    while(step_no < Symbol_Size)
    {
      /* Constructs a[1],a[2],... from prev_state=(a[1],a[2],...,a[M]) */

      Tmp_int = state;
      for(j=M;j>=1;j--)
      {
        a[j] = Tmp_int & 1;
        Tmp_int = Tmp_int >> 1;
      }

      sum_forward = 0;
      sum_feedback = 0;
      for (j=1; j<=M; j++)
      {
        sum_forward = sum_forward ^ (gforward[j] & a[j]);
        sum_feedback = sum_feedback ^ (gfeedback[j] & a[j]);
      }

      /* Calculate current state and output */

      a[0]=0;
      input = (2*input)+sum_feedback;
      /* In the termination of trellis, sum_feedback acts as the new input */
      output = (2*output)+((gforward[0] & a[0]) ^ sum_forward);

      current_state=0;
      Tmp_int=1;
      for (j=M-1; j>=0; j--)
      {
        current_state += (a[j]*Tmp_int);
        Tmp_int *= 2;
      }
      state = current_state;
      step_no++;
    }
    trellis_tail[prev_state][input][0] = current_state;
    trellis_tail[prev_state][input][1] = output;
  }
}

/******************************************************************************/
/* Calc_Sym_Bits (void)                                                       */
/*      Decomposes i into j bits and stores in symbol_bit[i][j]               */
/******************************************************************************/

void Calc_Sym_Bits(void)
{
  int i,j;

  for (i=0; i<Ortho_Size; i++)
    for (j=1; j<=Symbol_Size; j++)
      symbol_bit[i][j] = (i & ((int)pow(2,Symbol_Size-j)))>>(Symbol_Size-j);
}

/******************************************************************************/
/* Generate (void)                                                            */
/*      Generates source data randomly                                        */
/*      A 31 bit random number generated and each bit placed in Data[m]	      */
/*      Total number of bits = Num_of_Blocks*Block_Size*Symbol_Size           */
/******************************************************************************/

void Generate (void)
{
  long ran_number, bitmask;
  int j,k,l,m,now;
double max;
time_t t;
  max = 2147483647;

  /*srandom(time(&now));*/
  srand((unsigned) time(&t));
  for (k=1; k<=(((Num_of_Blocks*Block_Size)/31)+1); k++)
  {

    bitmask = 0x40000000;
    ran_number = rand() % 2147483647;
    l = 31*(k-1);
    for (j=1; j<=31; j++)
    {
      m=l+j;
      if ((ran_number & bitmask) > 0) data_org[m] = 1;
      else data_org[m] = 0;
      // data[m]=data_org[m];
      bitmask = bitmask >> 1;
    }
  }
}

/******************************************************************************/
/* Duplicate (void)                                                           */
/*      Duplicates source data bits                                           */
/*      Dup_Distance put every duplicated bit in distance dup_distance	      */
/*      Total number of bits = 2*Num_of_Blocks*Block_Size*Symbol_Size         */
/******************************************************************************/

void Duplicate ()
{
	int i;
	if (Repetition != 1)
	{
		for (i=1; i<=AS6; i++) {
			data[2*i-1] = data_org[i];
			data[2*i] = data_org[i];
		}
	}
	else
		for (i=1; i<=AS6; i++)
			data[i] = data_org[i];
}

/******************************************************************************/
/* Interleave_dup (void)                                                      */
/*      Interleave duplicated source data bits                                */
/*                                                                  	      */
/******************************************************************************/

void Interleave_dup (void)
{
  int k;
  char source_bit_p[Block_Size*Repetition];

  for (k=1; k<=(Block_Size*Repetition); k++)
	  source_bit_p[k] = source_bit[k];

  for (k=1; k<=(Block_Size*Repetition); k++)
    source_bit[interleaver[k]] = source_bit_p[k];
}


/******************************************************************************/
/* Deinterleave_dup (void)                                                      */
/*      De-interleave duplicated source data bits                                */
/*                                                                  	      */
/******************************************************************************/

void Deinterleave_dup (void)
{
  int k;

  for (k=1; k<=(Block_Size*Repetition); k++)
    detected_data_dup[k] = decoder2_bit[interleaver[k]];
}




/******************************************************************************/
/* Encoder (char encoder_num)                                                 */
/*      Encoders source bits using 'encoder_num' encoder                      */
/******************************************************************************/

void Encode (char encoder_num)
{
  int i,k,m;
  int N;        /* Size of the block */
  char a[AS1];
  char sum_forward, sum_feedback;
        /* Used in the forward and backward paths of the encoder filter */
        /* sum_forward is the sum of the forward path excluding a[0]    */
        /* sum_feedback is the sum of the feedback path                 */

  /* Set block size 'N' for encoder_num */

   N = Block_Size * Repetition;

  /* N+M = total number of bits in block of size N with M terminating bits */

  for (m=0; m<=N+M; m++)
    a[m] = 0;

  /* k is the index of the bit within the block */
  /* a[k] is the input bit after the first adder of the encoder */
  /* ^ is the binary addition */

  for (k=1; k<=N; k++)
  {
    sum_forward = 0;
    sum_feedback = 0;

    for (i=1; i<=M; i++)
    {
      if (k>i)  /* This loop takes care of the effect of previous bits  */
                /* to the extent that they are available                */
                /* Takes care of the empty bit places before having the */
                /* first M bits                                         */
      {
        sum_forward = sum_forward ^ (gforward[i] & a[k-i]);
        sum_feedback = sum_feedback ^ (gfeedback[i] & a[k-i]);
      }
    }

    switch (encoder_num)
    {
      case 1:
        /* Add new incomming bit with sum_feedback */
        a[k] = source_bit[k] ^ sum_feedback;
        encoder1_bit[k] = (gforward[0] & a[k]) ^ sum_forward;
        break;
      case 2:
        /* Add new incomming bit with sum_feedback */
        a[k] = interleaved_bit[k] ^ sum_feedback;
        encoder2_bit[k] = (gforward[0] & a[k]) ^ sum_forward;
        break;
    }
  }

  /* This part terminates the trellis */

  for (k=(N+1); k<=(N+M+Symbol_Size-1); k++)
  {
    sum_forward = 0;
    sum_feedback = 0;

    for (i=1; i<=M; i++)
    {
      sum_forward = sum_forward ^ (gforward[i] & a[k-i]);
      sum_feedback = sum_feedback ^ (gfeedback[i] & a[k-i]);
    }

    /* Considers sum_feedback as the new input and assigns a[k] = 0 */

    switch (encoder_num)
    {
      case 1:
        a[k] = 0;
        encoder1_bit[k] = (gforward[0] & a[k]) ^ sum_forward;
        break;
      case 2:
        source_bit[k] =  sum_feedback;
        a[k] = source_bit[k] ^ sum_feedback;
        /* Is equivalent to setting a[k]=0 */
        encoder2_bit[k] = (gforward[0] & a[k]) ^ sum_forward;
        break;
    }
  }
}

/******************************************************************************/
/* Interleave (void)                                                          */
/*      Interleaves source bits                                               */
/******************************************************************************/

void Interleave (void)
{
  int k;

  for (k=1; k<=Block_Size * Repetition; k++)
    interleaved_bit[interleaver[k]] = source_bit[k];
}

/******************************************************************************/
/* Interleave2 (void)                                                         */
/*      Interleaves source and encoder bits after Turbo encoding              */
/******************************************************************************/

void Interleave2 (void)
{
  int i,j,k,l,Punc_Blocks,punc_position;
  char datastream[AS02];
  double temp_a;
  int temp_b, leftover_bits;

  /* Initialize datastream */
  for (i=1; i<=AS02; i++)
    datastream[i] = 0;

  temp_a = (Block_Size * Repetition+M) / (double) Puncture_Size;
  temp_b = (Block_Size * Repetition+M) / Puncture_Size;
  if ((temp_a - temp_b) != 0)
    Punc_Blocks = temp_b + 1;
  else
    Punc_Blocks = temp_b;

  /* Combine the source and encoder bits into one stream */
  punc_position = 1;
  for (i=1; i<=Punc_Blocks; i++)
    for (j=0; j<Puncture_Size; j++)
    {
      l = Puncture_Size*(i-1) + j + 1;
      if (l <= (Block_Size * Repetition+M))
      {
        datastream[punc_position] = source_bit[l];
        punc_position += 1;
        if (puncture1[j] == 1)
        {
          datastream[punc_position] = encoder1_bit[l];
          punc_position += 1;
        }
        if (puncture2[j] == 1)
        {
          datastream[punc_position] = encoder2_bit[l];
          punc_position += 1;
        }
      }
    }

  /* Add extra bits, 0's, for proper QAM calculation */
  temp_a = (((Block_Size * Repetition+M)*Code_Rate)/(double)Num_QAMbits);
  temp_b = (((Block_Size * Repetition+M)*Code_Rate)/Num_QAMbits);
  leftover_bits = (temp_a - temp_b) * Num_QAMbits;
//  if (leftover_bits != 0)
 //   extra_bits = Num_QAMbits - leftover_bits;
//  else
    extra_bits = 0;

  if (extra_bits != 0)
   for (i=(Block_Size * Repetition+M)*Code_Rate+1;i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits;i++)
     datastream[i] = 0;

  /* Interleave information bit stream */
  for (k=1; k<=Block_Size * Repetition*Code_Rate; k++)
    datastream_int[interleaver2[k]] = datastream[k];

  for (k=Block_Size * Repetition*Code_Rate+1; k<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; k++)
    datastream_int[k] = datastream[k];
}

/******************************************************************************/
/* Modulate (void)                                                            */
/*      Modulates source bits and encoder bits                                */
/******************************************************************************/

void Modulate (void)
{
  int i,j,counter,Point,length;
  int u[Num_QAMbits+1];

  /* Calculate length of QAM stream */
  length = ((((Block_Size * Repetition+M)*Code_Rate)+extra_bits)/Num_QAMbits);

  /* Initialize X and Y point streams */
  for (i=1; i<=length; i++)
  {
    Xstream[i] = 0;
    Ystream[i] = 0;
  }


  for (i=1; i<=((Block_Size * Repetition+M)*Code_Rate)+extra_bits; i++)
  {
      /* Set X and Y coordinate */
      Xstream[i] = (2*datastream_int[i])-1;
      Ystream[i] = 0;

  }
}

/******************************************************************************/
/* Add_Noise (void)                                                           */
/*      Adds White Gaussian noise to encoded bits                             */
/******************************************************************************/

void Add_Noise (void)
{
  int i,j,k,now,length;
  float uniform[2], normal[2][AS1], s;
  double max, var;
time_t t;
  /* Calculate length of QAM stream */
  length = ((((Block_Size * Repetition+M)*Code_Rate)+extra_bits)/Num_QAMbits);

  /* Initialization */
  max = 2147483647;
     /* maximum value that random() can generate */

  var = N0/2;           /* Set variance of AWGN */

  /* Generate normal random variables */
  for (i=0; i<1; i++)
  {

    /*srandom((i+1)*time(&now));   Seed random number generator */
	
srand((unsigned) time(&t));

    for (j=0; j<=length; j++)

    {

      do
      {
        uniform[0] = (rand() % 2147483647)/max;

        uniform[0] = 2*uniform[0]-1;

        uniform[1] = (rand() % 2147483647)/max;

        uniform[1] = 2*uniform[1]-1;

        s = (uniform[0]*uniform[0])+(uniform[1]*uniform[1]);


      } while(s>=1);

      normal[i][j] = uniform[0]*sqrt(-2*log(s)/s);

    }


  }

  /* Add uncorrelated Gaussian noise to X and Y streams */
  for (i=1; i<=length; i++)
  {
    xstream[i] = Xstream[i]+sqrt(var)*normal[0][i];
    ystream[i] = Ystream[i];//+sqrt(var)*normal[1][i];
  }
}

/******************************************************************************/
/* Decode_Init ()                                                             */
/*      Initialization for decoder                                            */
/******************************************************************************/

void Decode_Init (void)
{
  int i,j,k,l,length,temp,temp_b,punc_position,Punc_Blocks;
  double QAM_pointprob[AS02][QAM_Size];
  double bitstream_prob_int[AS02][2],bitstream_prob[AS02][2];
  double temp_sum[AS02], temp_a;

  /* Calculate length of QAM stream */
  length = ((((Block_Size * Repetition+M)*Code_Rate)+extra_bits)/Num_QAMbits);

  temp_a = (Block_Size * Repetition+M) / (double) Puncture_Size;
  temp_b = (Block_Size * Repetition+M) / Puncture_Size;
  if ((temp_a - temp_b) != 0)
    Punc_Blocks = temp_b + 1;
  else
    Punc_Blocks = temp_b;

  /* Initialize the bit probabilities */
  for (k=1; k<=Block_Size * Repetition+M; k++)
    for (j=0; j<Ortho_Size; j++)
    {
      syst_prob[k][j]=(1/(double)Ortho_Size);
      par1_prob[k][j]=(1/(double)Ortho_Size);
      par2_prob[k][j]=(1/(double)Ortho_Size);
    }

  /* Initialize the extrinsic information */
  for (k=1; k<=Block_Size * Repetition+M; k++)
    for (j=0; j<Ortho_Size; j++)
      extrinsic[k][j] = 0;

  /* Initialize QAM_exp */
  for (k=1; k<=length; k++)
    for (j=0; j<QAM_Size; j++)
      QAM_exp[k][j] = 0;

  /* Calculate exp(-(x^2+y^2)/N0) for all QAM points */
  for (k=1; k<=length; k++)
    for (j=0; j<QAM_Size; j++)
      QAM_exp[k][j] = exp(-(pow((xstream[k]-QAM_Array[j][0]),2)
			+pow((ystream[k]-QAM_Array[j][1]),2))/N0);

  for (i=1; i<=length; i++)
      temp_sum[i] = 0;

  /* Calculate updated probability of each QAM point */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
    {
      QAM_pointprob[i][j] = QAM_exp[i][j] / (double) QAM_Size;
      temp_sum[i] += QAM_pointprob[i][j];
    }

  /* Normalize QAM point probabilities */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
      QAM_pointprob[i][j] /= temp_sum[i];      

  /* Initialize individual bit probabilities */
  for (i=1; i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; i++)
    for (j=0; j<Ortho_Size; j++)
      bitstream_prob_int[i][j] = 0;

  /* Convert QAM point probability to bit probability */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
    {
      temp = 1;
      for (k=0; k<Num_QAMbits; k++)
      {
        l = (j & temp) >> k;
        bitstream_prob_int[i*Num_QAMbits-k][l] += QAM_pointprob[i][j];
        temp = temp << 1;
      }
    }  

  /* Deinterleave bit stream */
  for (i=1; i<=Block_Size * Repetition*Code_Rate; i++)
    for (j=0; j<Ortho_Size; j++)
      bitstream_prob[i][j] = bitstream_prob_int[interleaver2[i]][j];

  for (i=Block_Size * Repetition*Code_Rate+1; i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; i++)
    for (j=0; j<Ortho_Size; j++)
      bitstream_prob[i][j] = bitstream_prob_int[i][j];

  /* Destream bit probabilities into systematic and parity bit streams */
  punc_position = 1;
  for (i=1; i<=Punc_Blocks; i++)
    for (j=0; j<Puncture_Size; j++)
    {
      l = Puncture_Size*(i-1) + j + 1;
      if (l <= (Block_Size * Repetition+M))
      {
        received_sourceprob_bit[l][0] = bitstream_prob[punc_position][0];
        received_sourceprob_bit[l][1] = bitstream_prob[punc_position][1];
        punc_position += 1;
        if (puncture1[j] == 1)
        {
          received_encoder1prob_bit[l][0] = bitstream_prob[punc_position][0];
          received_encoder1prob_bit[l][1] = bitstream_prob[punc_position][1];
          punc_position += 1;
        }
        else
        {
          received_encoder1prob_bit[l][0] = par1_prob[l][0];
          received_encoder1prob_bit[l][1] = par1_prob[l][1];
        }
        if (puncture2[j] == 1)
        {
          received_encoder2prob_bit[l][0] = bitstream_prob[punc_position][0];
          received_encoder2prob_bit[l][1] = bitstream_prob[punc_position][1];
          punc_position += 1;
        }
        else
        {
          received_encoder2prob_bit[l][0] = par2_prob[l][0];
          received_encoder2prob_bit[l][1] = par2_prob[l][1];
        }
      }
    }

  /* Interleave systematic probability for second RSC decoder of TD */
  for (i=1; i<=Block_Size * Repetition; i++)
    for (j=0; j<Ortho_Size; j++)
      received_sourceprob_int[interleaver[i]][j]=received_sourceprob_bit[i][j];

  for (i=Block_Size * Repetition+1; i<=Block_Size * Repetition+M; i++)
    for (j=0; j<Ortho_Size; j++)
      received_sourceprob_int[i][j]=received_sourceprob_bit[i][j];
}

/******************************************************************************/
/* QAM_dist ()                                                                */
/*      Calculates pdf of received bits to be used in turbo decoder           */
/******************************************************************************/

void QAM_dist (void)
{
  int i,j,k,l,temp,counter,length,par_counter;
  int Punc_Blocks,temp_b,punc_position;
  double a_priori[2];
  double prob_stream[AS02][2], prob_stream_int[AS02][2];
  double p[Num_QAMbits+1][2], syst_bit_prob[AS0][2];
  double QAM_APPprob[AS02][QAM_Size];
  double QAM_pointprob[AS02][QAM_Size];
  double bitstream_prob_int[AS02][2], bitstream_prob[AS02][2];
  double temp_sum[AS02], temp_a;

  /* Calculate length of QAM stream */
  length = ((((Block_Size * Repetition+M)*Code_Rate)+extra_bits)/Num_QAMbits);

  temp_a = (Block_Size * Repetition+M) / (double) Puncture_Size;
  temp_b = (Block_Size * Repetition+M) / Puncture_Size;
  if ((temp_a - temp_b) != 0)
    Punc_Blocks = temp_b + 1;
  else
    Punc_Blocks = temp_b;

  /* Combine systematic and parity bit probabilities from Turbo Decoder*/
  punc_position = 1;
  for (i=1; i<=Punc_Blocks; i++)
    for (j=0; j<Puncture_Size; j++)
    {
      l = Puncture_Size*(i-1) + j + 1;
      if (l <= (Block_Size * Repetition+M))
      {
        prob_stream[punc_position][0] = syst_prob[l][0];
        prob_stream[punc_position][1] = syst_prob[l][1];
        punc_position += 1;
        if (puncture1[j] == 1)
        {
          prob_stream[punc_position][0] = par1_prob[l][0];
          prob_stream[punc_position][1] = par1_prob[l][1];
          punc_position += 1;
        }
        if (puncture2[j] == 1)
        {
          prob_stream[punc_position][0] = par2_prob[l][0];
          prob_stream[punc_position][1] = par2_prob[l][1];
          punc_position += 1;
        }
      }
    }

  /* Add extra bits, known to be zeroes, if necessary */
  if (extra_bits != 0)
   for (i=(Block_Size * Repetition+M)*Code_Rate+1;i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits;i++)
   {
     prob_stream[i][0] = 1;
     prob_stream[i][1] = 0;
   }

  /* Interleave APP probabilities */
  for (i=1; i<=Block_Size * Repetition*Code_Rate; i++)
    for (j=0; j<Ortho_Size; j++)
      prob_stream_int[interleaver2[i]][j] = prob_stream[i][j];

  for (i=Block_Size * Repetition*Code_Rate+1; i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; i++)
    for (j=0; j<Ortho_Size; j++)
      prob_stream_int[i][j] = prob_stream[i][j];
      
  /* Initialize QAM point APP probabilities */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
      QAM_APPprob[i][j] = 1;

  counter = 0;
  for (i=1; i<=Num_QAMbits; i++)
    for (j=0; j<Ortho_Size; j++)
      p[i][j] = 0;

  /* Calculate QAM Point APP Probabilities */
  for (i=1; i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; i++)
  {
    counter += 1;
    for (j=0; j<Ortho_Size; j++)
      p[counter][j] = prob_stream_int[i][j];
    if (counter == Num_QAMbits)
    {
      for (j=0; j<QAM_Size; j++)
      {
        temp = 1;
        for (k=0; k<Num_QAMbits; k++)
        {
          l = (j & temp) >> k;
          QAM_APPprob[i/Num_QAMbits][j] *= p[Num_QAMbits-k][l];
          temp = temp << 1;
        }
      }
      counter = 0;
      for (j=1; j<=Num_QAMbits; j++)
        for (k=0; k<Ortho_Size; k++)
          p[j][k] = 0;
    }
  }

  for (i=1; i<=length; i++)
      temp_sum[i] = 0;

  /* Calculate updated probability of each QAM point */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
    {     
     QAM_pointprob[i][j] = QAM_exp[i][j];
     temp_sum[i] += QAM_pointprob[i][j];
    }

  /* Normalize QAM point probabilities */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
      QAM_pointprob[i][j] /= temp_sum[i];      

  /* Initialize individual bit probabilities */
  for (i=1; i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; i++)
    for (j=0; j<Ortho_Size; j++)
      bitstream_prob_int[i][j] = 0;

  /* Convert QAM point probability to bit probability */
  for (i=1; i<=length; i++)
    for (j=0; j<QAM_Size; j++)
    {
      temp = 1;
      for (k=0; k<Num_QAMbits; k++)
      {
        l = (j & temp) >> k;
        a_priori[l] = QAM_APPprob[i][j]/prob_stream_int[i*Num_QAMbits-k][l];
        bitstream_prob_int[i*Num_QAMbits-k][l] += (QAM_pointprob[i][j]*a_priori[l]);
        temp = temp << 1;
      }
    }  

  /* Deinterleave bit stream */
  for (i=1; i<=Block_Size * Repetition*Code_Rate; i++)
    for (j=0; j<Ortho_Size; j++)
      bitstream_prob[i][j] = bitstream_prob_int[interleaver2[i]][j];

  for (i=Block_Size * Repetition*Code_Rate+1; i<=(Block_Size * Repetition+M)*Code_Rate+extra_bits; i++)
    for (j=0; j<Ortho_Size; j++)
      bitstream_prob[i][j] = bitstream_prob_int[i][j];

/* Destream bit probabilities into systematic bit stream */
  punc_position = 1;
  for (i=1; i<=Punc_Blocks; i++)
    for (j=0; j<Puncture_Size; j++)
    {
      l = Puncture_Size*(i-1) + j + 1;
      if (l <= (Block_Size * Repetition+M))
      {        
        received_sourceprob_bit[l][0] = bitstream_prob[punc_position][0];
        received_sourceprob_bit[l][1] = bitstream_prob[punc_position][1];
        punc_position += 1;
        if (puncture1[j] == 1)
          punc_position += 1;
        if (puncture2[j] == 1)
          punc_position += 1;
      }
    }
/* Interleave systematic probability for second RSC decoder of TD */
 for (i=1; i<=Block_Size * Repetition; i++)
    for (j=0; j<Ortho_Size; j++)
      received_sourceprob_int[interleaver[i]][j]=received_sourceprob_bit[i][j];

  for (i=Block_Size * Repetition+1; i<=Block_Size * Repetition+M; i++)
    for (j=0; j<Ortho_Size; j++)
      received_sourceprob_int[i][j]=received_sourceprob_bit[i][j];
 

 /* Calculate the difference of the LLR's coming in and out */
  for (i=1; i<=Block_Size * Repetition; i++)
  {
    Syst_LLR[i][0] = 0;
    Syst_LLR[i][1] = log(received_sourceprob_bit[i][1]*syst_prob[i][0])-
                     log(received_sourceprob_bit[i][0]*syst_prob[i][1]);
  }

  for (i=Block_Size * Repetition+1; i<=Block_Size * Repetition+M; i++)
  {
    Syst_LLR[i][0] = 0;
    Syst_LLR[i][1] = 0;
  }
}

/******************************************************************************/
/* Turbo_Decode ()                                                            */
/*      Decodes received source bits using encoder bits                       */
/******************************************************************************/

void Turbo_Decode (void)
{
  int i,j,k,m,mp;
  int iteration, max_ortho;
  int num_of_dec1_errors, num_of_dec2_errors;
  double llr[AS0][Ortho_Size], app[AS0][Ortho_Size], app_sum;
  double gamma1[AS0][Num_of_States][Ortho_Size];
  double gamma2[AS0][Num_of_States][Ortho_Size];
  double alpha_temp[Num_of_States], alpha_denominator;
  double alpha1[AS0][Num_of_States], alpha2[AS0][Num_of_States];
  double beta_num[Num_of_States], beta_den[AS0], beta_denominator;
  double beta1[AS0][Num_of_States], beta2[AS0][Num_of_States];
  double delta_temp[Ortho_Size], Delta_int[AS0][Ortho_Size];
  double delta_par1[Ortho_Size], delta_par2[Ortho_Size];
  double Delta1[AS0][Ortho_Size], Delta2[AS0][Ortho_Size];
  double L1[AS0][Ortho_Size], L2[AS0][Ortho_Size]; /* extrinsic info. */
  double L2_enhanced[AS0][Ortho_Size];
  double L2_deinterleaved[AS0][Ortho_Size];
  double max_delta, temp_tot;
  double syst_sum, syst_prob_int[AS1][Ortho_Size];
  double syst_prob_temp[AS1][Ortho_Size];
  double par1_sum, par1_prob_temp[AS1][Ortho_Size];
  double par2_sum, par2_prob_temp[AS1][Ortho_Size];
  int max_llr_log, min_llr_log;


  /*   Initialization   */

  max_llr_log = 150;
  min_llr_log = -150;

  alpha1[0][0] = 1;
  alpha2[0][0] = 1;
  beta1[Block_Size * Repetition+M][0] = 1;
  beta2[Block_Size * Repetition+M][0] = 1;

  for (m=1; m<Num_of_States; m++)
  {
    alpha1[0][m] = 0;
    alpha2[0][m] = 0;
    beta1[Block_Size * Repetition+M][m] = 0;
    beta2[Block_Size * Repetition+M][m] = 0;
  }

  /* Initializing the extrinsic info for dec1 */
  for (k=1; k<=Block_Size * Repetition+M; k++)
    for (i=0; i<Ortho_Size; i++)
      L2[k][i] = extrinsic[k][i];
      
  /* Initialize error_sums */
  for (iteration=1; iteration<=Max_TC_Iter; iteration++)
  {
    error_sum[iteration][1]=0;
    error_sum[iteration][2]=0;
  }

  /* Calculating contant portion of gamma values */
  for (k=1; k<=Block_Size * Repetition; k++)
    for (mp=0; mp<Num_of_States; mp++)
      for (i=0; i<Ortho_Size; i++)
      {
        gamma1[k][mp][i] = 1.0;
        gamma2[k][mp][i] = 1.0;
      }

  /* Calculate gamma1 */
  for (k=1; k<=Block_Size * Repetition; k++)
    for (mp=0; mp<Num_of_States; mp++)
      for (i=0; i<Ortho_Size; i++)
      {
        gamma1[k][mp][i]*=received_sourceprob_bit[k][symbol_bit[i][1]];
        gamma1[k][mp][i]*=received_encoder1prob_bit[k]
                                        [symbol_bit[trellis[mp][i][1]][1]];
      }

  /* The tail bits for gamma1 */
  for(k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
    for (mp=0; mp<Num_of_States; mp++)
      for (i=0; i<Ortho_Size; i++)
        if (trellis_tail[mp][i][0] != -1)
        {
          gamma1[k][mp][i] = 1.00;
          /* Second encoder terminated */
          /* gamma1[k][mp][i]*=received_sourceprob_bit[k]
                                                       [symbol_bit[i][1]]; */
          gamma1[k][mp][i]*=received_encoder1prob_bit[k]
                                       [symbol_bit[trellis[mp][i][1]][1]];
        }
        else
          gamma1[k][mp][i] = 0;

  /* Calculate gamma2 */
  for (k=1; k<=Block_Size * Repetition; k++)
    for (mp=0; mp<Num_of_States; mp++)
      for (i=0; i<Ortho_Size; i++)
      {
        gamma2[k][mp][i]*=received_sourceprob_int[k][symbol_bit[i][1]];
        gamma2[k][mp][i]*=received_encoder2prob_bit[k]
                                        [symbol_bit[trellis[mp][i][1]][1]];
      }

  /* The tail bits for gamma2 */
  for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
    for (mp=0; mp<Num_of_States; mp++)
      for (i=0; i<Ortho_Size; i++)
        if (trellis_tail[mp][i][0] != -1)
        {
          gamma2[k][mp][i] = 1.00;
          gamma2[k][mp][i]*=received_sourceprob_int[k][symbol_bit[i][1]];
          gamma2[k][mp][i]*=received_encoder2prob_bit[k]
                                          [symbol_bit[trellis[mp][i][1]][1]];
        }
        else
          gamma2[k][mp][i] = 0;

  /* Begin Decoding Loop for each iteration */
  for (iteration=1; iteration<=Max_TC_Iter; iteration++)
  {
    /* Initialize Error counts */
    num_of_dec1_errors = 0;
    num_of_dec2_errors = 0;

    /* Receive LLR values from previous iteration */
    for (k=1; k<=Block_Size * Repetition+M; k++)
      for(i=0; i<Ortho_Size; i++)
        llr[k][i] = L2[k][i];

    /* Calculate app values from llr values */
    for (k=1; k<=Block_Size * Repetition+M; k++)
    {
      app_sum = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        app[k][i] = exp(llr[k][i]);
        app_sum += app[k][i];
      }
      for (i=0; i<Ortho_Size; i++)
        app[k][i] /= app_sum;
    }

    /* Calculate alpha values */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      for (m=0; m<Num_of_States; m++)
        alpha_temp[m] = 0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          alpha_temp[trellis[mp][i][0]]+=gamma1[k][mp][i]*app[k][i]
						*alpha1[k-1][mp];

      alpha_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        alpha_denominator += alpha_temp[m];

      for (m=0; m<Num_of_States; m++)
        alpha1[k][m] = alpha_temp[m]/alpha_denominator;
    }

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
    {
      for (m=0; m<Num_of_States; m++)
        alpha_temp[m] = 0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[mp][i][0] != -1)
            alpha_temp[trellis_tail[mp][i][0]]+=gamma1[k][mp][i]
					*app[k][i]*alpha1[k-1][mp];

      alpha_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        alpha_denominator += alpha_temp[m];

      for (m=0; m<Num_of_States; m++)
        alpha1[k][m] = alpha_temp[m]/alpha_denominator;
    }

    /* Calculate beta values */
    for (k=(Block_Size * Repetition+M-1); k>=Block_Size * Repetition; k--)
    {
      for (m=0; m<Num_of_States; m++)
        beta_num[m] = 0;

      for (m=0; m<Num_of_States; m++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[m][i][0] != -1)
            beta_num[m]+=gamma1[k+1][m][i]*app[k+1][i]
				*beta1[k+1][trellis_tail[m][i][0]];

      for (m=0; m<Num_of_States; m++)
        beta_den[m] =0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[m][i][0] != -1)
            beta_den[trellis_tail[mp][i][0]]+=gamma1[k+1][mp][i]
                                *app[k+1][i]*alpha1[k][mp];

      beta_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        beta_denominator += beta_den[m];

      for (m=0; m<Num_of_States; m++)
        beta1[k][m] = beta_num[m]/beta_denominator;
    }

    for (k=(Block_Size * Repetition-1); k>=1; k--)
    {
      for (m=0; m<Num_of_States; m++)
        beta_num[m] = 0;

      for (m=0; m<Num_of_States; m++)
        for (i=0; i<Ortho_Size; i++)
          beta_num[m]+=gamma1[k+1][m][i]*app[k+1][i]
				*beta1[k+1][trellis[m][i][0]];

      for (m=0; m<Num_of_States; m++)
        beta_den[m] =0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          beta_den[trellis[mp][i][0]]+=gamma1[k+1][mp][i]
                                *app[k+1][i]*alpha1[k][mp];

      beta_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        beta_denominator += beta_den[m];

      for (m=0; m<Num_of_States; m++)
        beta1[k][m] = beta_num[m]/beta_denominator;
    }

    /* Calculate delta values */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      par1_sum = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        delta_temp[i] = 0;
        delta_par1[i] = 0;
      }

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
        {
          delta_temp[i] += gamma1[k][mp][i]*app[k][i]
    		            *alpha1[k-1][mp]*beta1[k][trellis[mp][i][0]];
          delta_par1[trellis[mp][i][1]] += gamma1[k][mp][i]*app[k][i]
    		            *alpha1[k-1][mp]*beta1[k][trellis[mp][i][0]];
        }

      for (i=0; i<Ortho_Size; i++)
      {
        delta_temp[i] = log(delta_temp[i]);
        par1_sum += delta_par1[i];
      }

      Delta_int[k][0] = 0;
      Delta_int[k][1] = delta_temp[1]-delta_temp[0];
      if (Delta_int[k][1] < min_llr_log)
        Delta_int[k][1] = min_llr_log;
      if (Delta_int[k][1] > max_llr_log)
        Delta_int[k][1] = max_llr_log;

      for (i=0; i<Ortho_Size; i++)
        par1_prob_temp[k][i] = (delta_par1[i]/par1_sum);
    }

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
    {
      par1_sum = 0;
      for (i=0; i<Ortho_Size; i++)
        delta_par1[i] = 0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[mp][i][0] != -1)
            delta_par1[trellis_tail[mp][i][1]] += gamma1[k][mp][i]
              *app[k][i]*alpha1[k-1][mp]*beta1[k][trellis_tail[mp][i][0]];

      for (i=0; i<Ortho_Size; i++)
        par1_sum += delta_par1[i];

      for (i=0; i<Ortho_Size; i++)
        par1_prob_temp[k][i] = (delta_par1[i]/par1_sum);
    }

    /* Copy Delta_int values into Delta1 */
    for (k=1; k<=Block_Size * Repetition; k++)
      for (i=0; i<Ortho_Size; i++)
        Delta1[k][i] = Delta_int[k][i];

    /* Calculate L1[k][i] */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      L1[k][0] = 0;
      L1[k][1] = Delta1[k][1] - L2[k][1];
      L1[k][1] -= log(received_sourceprob_bit[k][1])
                 -log(received_sourceprob_bit[k][0]);
                /* subtracting the systematic info.*/
      if (L1[k][1] > max_llr_log)
        L1[k][1] = max_llr_log;
      if (L1[k][1] < min_llr_log)
        L1[k][1] = min_llr_log;
    }

    for (k=Block_Size * Repetition+1; k<=Block_Size+M; k++)
      for (i=0; i<Ortho_Size; i++)
        L1[k][i] = 0;

    /* Decode using Delta1[k][i] and check for errors */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      /* Find largest LLR value */
      max_delta = 0;
      max_ortho = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        if (Delta1[k][i] > max_delta)
        {
          max_delta = Delta1[k][i];
          max_ortho = i;
        }
      }

      /* Convert largest LLR into a data symbol */
      decoder1_bit[k] = max_ortho;

      /* Error checking and tabulation */
      if (decoder1_bit[k] != source_bit[k])
      {
        error1[k] = 1;
        num_of_dec1_errors++;
      }
      else
        error1[k] = 0;
    }


    /* DEC 2 */

    /* interleave L1[k] */
    for (k=1; k<=Block_Size * Repetition; k++)
      for (i=0; i<Ortho_Size; i++)
        llr[interleaver[k]][i] = L1[k][i];

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
      for (i=0; i<Ortho_Size; i++)
        llr[k][i] = L1[k][i];

    /* Calculate app values from llr values */
    for (k=1; k<=Block_Size * Repetition+M; k++)
    {
      app_sum = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        app[k][i] = exp(llr[k][i]);
        app_sum += app[k][i];
      }
      for (i=0; i<Ortho_Size; i++)
        app[k][i] /= app_sum;
    }

    /* Calculate alpha values */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      for (m=0; m<Num_of_States; m++)
        alpha_temp[m] = 0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          alpha_temp[trellis[mp][i][0]]+=gamma2[k][mp][i]*app[k][i]
						*alpha2[k-1][mp];

      alpha_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        alpha_denominator += alpha_temp[m];

      for (m=0; m<Num_of_States; m++)
        alpha2[k][m] = alpha_temp[m]/alpha_denominator;
    }

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
    {
      for (m=0; m<Num_of_States; m++)
        alpha_temp[m] = 0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[mp][i][0] != -1)
            alpha_temp[trellis_tail[mp][i][0]]+=gamma2[k][mp][i]
					*app[k][i]*alpha2[k-1][mp];

      alpha_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        alpha_denominator += alpha_temp[m];

      for (m=0; m<Num_of_States; m++)
        alpha2[k][m] = alpha_temp[m]/alpha_denominator;
    }

    /* Calculate beta values */
    for (k=(Block_Size * Repetition+M-1); k>=Block_Size * Repetition; k--)
    {
      for (m=0; m<Num_of_States; m++)
        beta_num[m] = 0;

      for (m=0; m<Num_of_States; m++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[m][i][0] != -1)
            beta_num[m]+=gamma2[k+1][m][i]*app[k+1][i]
			*beta2[k+1][trellis_tail[m][i][0]];

      for (m=0; m<Num_of_States; m++)
        beta_den[m] =0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[m][i][0] != -1)
            beta_den[trellis_tail[mp][i][0]]+=gamma2[k+1][mp][i]
                                *app[k+1][i]*alpha2[k][mp];

      beta_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        beta_denominator += beta_den[m];

      for (m=0; m<Num_of_States; m++)
        beta2[k][m] = beta_num[m]/beta_denominator;
    }

    for (k=(Block_Size * Repetition-1); k>=1; k--)
    {
      for (m=0; m<Num_of_States; m++)
        beta_num[m] = 0;

      for (m=0; m<Num_of_States; m++)
        for (i=0; i<Ortho_Size; i++)
          beta_num[m]+=gamma2[k+1][m][i]*app[k+1][i]
			*beta2[k+1][trellis[m][i][0]];

      for (m=0; m<Num_of_States; m++)
        beta_den[m] =0;

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          beta_den[trellis[mp][i][0]]+=gamma2[k+1][mp][i]
                                *app[k+1][i]*alpha2[k][mp];

      beta_denominator = 0;
      for (m=0; m<Num_of_States; m++)
        beta_denominator += beta_den[m];

      for (m=0; m<Num_of_States; m++)
        beta2[k][m] = beta_num[m]/beta_denominator;
    }


    /* Calculate delta values */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      syst_sum = 0;
      par2_sum = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        delta_temp[i] = 0;
        delta_par2[i] = 0;
      }

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
        {
          delta_temp[i] += gamma2[k][mp][i]*app[k][i]
				*alpha2[k-1][mp]*beta2[k][trellis[mp][i][0]];
          delta_par2[trellis[mp][i][1]] += gamma2[k][mp][i]*app[k][i]
				*alpha2[k-1][mp]*beta2[k][trellis[mp][i][0]];
        }

      for (i=0; i<Ortho_Size; i++)
      {
        syst_prob_int[k][i] = delta_temp[i];
        syst_sum += syst_prob_int[k][i];
        delta_temp[i] = log(delta_temp[i]);
        par2_sum += delta_par2[i];
      }

      Delta_int[k][0] = 0;
      Delta_int[k][1] = delta_temp[1]-delta_temp[0];
      if (Delta_int[k][1] < min_llr_log)
        Delta_int[k][1] = min_llr_log;
      if (Delta_int[k][1] > max_llr_log)
        Delta_int[k][1] = max_llr_log;
      for (i=0; i<Ortho_Size; i++)
      {
        syst_prob_int[k][i] /= syst_sum;
        par2_prob_temp[k][i] = (delta_par2[i]/par2_sum);
      }
    }

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
    {
      syst_sum = 0;
      par2_sum = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        delta_temp[i] = 0;
        delta_par2[i] = 0;
      }

      for (mp=0; mp<Num_of_States; mp++)
        for (i=0; i<Ortho_Size; i++)
          if (trellis_tail[mp][i][0] != -1)
          {
            delta_temp[i] += gamma2[k][mp][i]*app[k][i]
                             *alpha2[k-1][mp]*beta2[k][trellis_tail[mp][i][0]];
            delta_par2[trellis_tail[mp][i][1]] += gamma2[k][mp][i]
              *app[k][i]*alpha2[k-1][mp]*beta2[k][trellis_tail[mp][i][0]];
          }

      for (i=0; i<Ortho_Size; i++)
      {
        syst_prob_int[k][i] = delta_temp[i];
        syst_sum += syst_prob_int[k][i];
        par2_sum += delta_par2[i];
      }

      for (i=0; i<Ortho_Size; i++)
      {
        syst_prob_int[k][i] /= syst_sum;
        par2_prob_temp[k][i] = (delta_par2[i]/par2_sum);
      }
    }

    /* De-interleaving the output of dec2 */
    for (k=1; k<=Block_Size * Repetition; k++)
      for (i=0; i<Ortho_Size; i++)
      {
        Delta2[k][i] = Delta_int[interleaver[k]][i];
        syst_prob_temp[k][i] = syst_prob_int[interleaver[k]][i];
      }

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
      for (i=0; i<Ortho_Size; i++)
        syst_prob_temp[k][i] = syst_prob_int[k][i];


    /* Calculate L2[k][i] */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      L2[k][0] = 0;
      L2[k][1] = Delta2[k][1] - L1[k][1];
      L2[k][1] -= log(received_sourceprob_bit[k][1])
                 -log(received_sourceprob_bit[k][0]);
                  /* subtracting the systematic info.*/
    }




    /********************************************************************
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     **********************************************************************/




if (LLR_summation == 1)
{
    /* De-interleaving the output of dec2 */
    for (k=1; k<=Block_Size*Repetition; k++)
      for (i=0; i<Ortho_Size; i++)
      {
        L2_deinterleaved[k][i] = L2[interleaver[k]][i];
      }


    for (k=Block_Size*Repetition+1; k<=Block_Size*Repetition+M; k++)
      for (i=0; i<Ortho_Size; i++)
        L2_deinterleaved[k][i] = L2[k][i];

    /* Add LLRs*/
    for (k=1; k<=Block_Size; k++)
      for (i=0; i<Ortho_Size; i++)
      {
    	  L2_enhanced[2*k-1][i] = (L2_deinterleaved[2*k-1][i]+L2_deinterleaved[2*k][i]);
    	  L2_enhanced[2*k][i] = (L2_deinterleaved[2*k-1][i]+L2_deinterleaved[2*k][i]);
      }

    for (k=Block_Size*Repetition+1; k<=Block_Size*Repetition+M; k++)
      for (i=0; i<Ortho_Size; i++)
    	  L2_enhanced[k][i] = L2_deinterleaved[k][i];

    /* interleave L2_enhanced[k] */
    for (k=1; k<=Block_Size*Repetition; k++)
      for (i=0; i<Ortho_Size; i++)
        L2[interleaver[k]][i] = L2_enhanced[k][i];

    for (k=Block_Size*Repetition+1; k<=Block_Size*Repetition+M; k++)
      for (i=0; i<Ortho_Size; i++)
        L2[k][i] = L2_enhanced[k][i];


  }


    /*************************************************************************
     *
     ***********************************************************************88*/






    for (k=1; k<=Block_Size * Repetition; k++)
    {
      if (L2[k][1] > max_llr_log)
        L2[k][1] = max_llr_log;
      if (L2[k][1] < min_llr_log)
        L2[k][1] = min_llr_log;
    }

    for (k=Block_Size * Repetition+1; k<=Block_Size * Repetition+M; k++)
      for (i=0; i<Ortho_Size; i++)
        L2[k][i] = 0;

    /* Decode using Delta2 and check for errors */
    for (k=1; k<=Block_Size * Repetition; k++)
    {
      /* Find largest LLR value */
      max_delta = 0;
      max_ortho = 0;
      for (i=0; i<Ortho_Size; i++)
      {
        if (Delta2[k][i] > max_delta)
        {
          max_delta = Delta2[k][i];
          max_ortho = i;
        }
      }

      /* Convert largest LLR into data symbol */
      for (j=1; j<=Symbol_Size; j++)
        decoder2_bit[k] = max_ortho;

      /* Error checking and tabulation */
      if (decoder2_bit[k] != source_bit[k])
      {
        error2[k] = 1;
        num_of_dec2_errors++;
      }
      else
        error2[k] = 0;
    }

    error_sum[iteration][1]=error_sum[iteration][1]+ num_of_dec1_errors;
    error_sum[iteration][2]=error_sum[iteration][2]+ num_of_dec2_errors;
  }

  /* Store L2[k][i] for next FB iteration */
  for (k=1; k<=Block_Size * Repetition+M; k++)
    for (i=0; i<Ortho_Size; i++)
      extrinsic[k][i] = L2[k][i];

  /* Subtract LLR of input to TD from the output LLR of TD */
  for (k=1; k<=Block_Size * Repetition+M; k++)
  {
    temp_tot = 0;
    syst_prob[k][0] = received_sourceprob_bit[k][1]*syst_prob_temp[k][0];
    syst_prob[k][1] = received_sourceprob_bit[k][0]*syst_prob_temp[k][1];
    temp_tot = syst_prob[k][0] + syst_prob[k][1];
    syst_prob[k][0] /= temp_tot;
    syst_prob[k][1] /= temp_tot;
    temp_tot = 0;
    par1_prob[k][0] = received_encoder1prob_bit[k][1]*par1_prob_temp[k][0];
    par1_prob[k][1] = received_encoder1prob_bit[k][0]*par1_prob_temp[k][1];
    temp_tot = par1_prob[k][0] + par1_prob[k][1];
    par1_prob[k][0] /= temp_tot;
    par1_prob[k][1] /= temp_tot;
    temp_tot = 0;
    par2_prob[k][0] = received_encoder2prob_bit[k][1]*par2_prob_temp[k][0];
    par2_prob[k][1] = received_encoder2prob_bit[k][0]*par2_prob_temp[k][1];
    temp_tot = par2_prob[k][0] + par2_prob[k][1];
    par2_prob[k][0] /= temp_tot;
    par2_prob[k][1] /= temp_tot;
  }

}

