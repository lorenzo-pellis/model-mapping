//
//  typedecls.h
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

/////////////////////////////* Type declarations *///////////////////////////////////

#pragma once
/////////////////////////////* Precompiler switches *///////////////////////////
#define CMDL_INPUT_SOCIAL_STRUCTURE_NETWORK_TYPE 1
#define CMDL_INPUT_Rg 1
#define CMDL_INPUT_Rh 0  // Remember Rh is true transmission in schools
#define CMDL_INPUT_Rw 1
#define CMDL_INPUT_SIGMA 1 // Relative susceptibility of children VS adults
#define CMDL_INPUT_RHO 1 // Relativs infectivity of children VS adults
#define CMDL_INPUT_PHI 1 // Global ssortativity
#define CMDL_INPUT_INIT_INF_A 1 // Global ssortativity
#define CMDL_INPUT_INIT_INF_C 1 // Global ssortativity
#define CMDL_INPUT_GAMMA 1 // Relativs infectivity of children VS adults

#define SAME_SIGMA 1 // Same relative susceptibility of children VS adults in the community and in the households?
#define SAME_RHO 1 // Same relative infectivity of children VS adults in the community and in the households?
#define SAME_PHI 0 // Same assortativity in the community and in the households?
#define RANDOM_MIXING_IN_HOUSEHOLDS 1 // Should I use random mixing in households? Definitely yes, otherwise I run into troubles


#define SWITCH_SCHOOL_CLOSURE 0
#define SWITCH_SCHOOL_REOPEN 0
#define SWITCH_CLOSE_EACH_SCHOOL 0
#define SWITCH_CLOSE_ALL_SCHOOLS 0
#define SWITCH_PREVALENCE_TRIGGER 1 // If 0, I use cumulative incidence to trigger closure
#define CMDL_INPUT_DURATION_CLOSURE 0
#define CMDL_INPUT_WITHIN_SCHOOL_PREV_TRIGGER 0
#define CMDL_INPUT_PERC_SCHOOLS_INFECTED_TRIGGER 0
#define CMDL_INPUT_DELAY 0
#define SPLIT_DELAY 0 // One big delay or break it in components

/* Define what to write in the name of the output file */
#define OUT_NETWORK_TYPE 1
#define OUT_Rg 1
#define OUT_Rh 0
#define OUT_Rw 1
#define OUT_SIGMA 1 // Relative susceptibility of children VS adults
#define OUT_RHO 1 // Relative infectivity of children VS adults
#define OUT_PHI 1 // Assortativity
#define OUT_GAMMA 1 // Ratio of contact rates of adults VS children (in each environment
#define OUT_WITHIN_SCHOOL_PREV_TRIGGER 0
#define OUT_PERC_SCHOOLS_INFECTED_TRIGGER 0
#define OUT_DELAY 0
#define OUT_DURATION 0
#define OUT_INF_MEAN 0
#define OUT_INF_ALPHA 0
#define OUT_GEN_TIME 0
#define OUT_GEN_ALPHA 0

/////////////////////////////* Typedefs *///////////////////////////////////

typedef struct parameters {
    /* Model settings */
    int MODEL_TYPE;		/*	MODEL_TYPE = 0 is the sSIR model (possibly with latent period)
                         MODEL_TYPE = 1 is the time-since-infection model with truncated gamma
                         MODEL_TYPE != 0 and 1 is the model with recovery delay... */
    int LATENCE;		/*	LATENT = 0 is the SIR, LATENT = 1 is the SEIR	*/
    double DUR_PRODROMAL;	/* Time lag from the infection (before latent period)
                             to when children feel sick and stay at home */
    double PROB_ASYMPT;		/* Probability of children being asymptomatic and
                             rematining in school even when infectious */
    double FRACTION_OF_Rh_IN_SCHOOL_IF_SYMPT; // If children are symptomatic and stay at home when sick, Rh is reduced to Rh*FRACTION
    double DUR_SICK;		/* Time from the infection at which individual is not
                             infectious anymore. This is the case when MODEL_TYPE = 1	*/
    double REC_DELAY;		/*	Time at which individuals recover after the last contact
                             in the time-since-infection model ( MODEL_TYPE != 0,1 ) */
    int SOCIAL_STRUCTURE_NETWORK_TYPE_LEN;
    char* SOCIAL_STRUCTURE_NETWORK_TYPE;
    
    double DELAY_INF_ONSET;
    double DELAY_ONSET_DET;
    double DELAY_DET_CONFIRM;
    int DELAY_VEC_LEN;
    double *DELAY_VEC;
    double DELAY; // ( DELAY_INF_ONSET + DELAY_ONSET_DET + DELAY_DET_CONFIRM )
    
    int WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN;
    int *WITHIN_SCHOOL_PREV_TRIGGER_VEC;
    int WITHIN_SCHOOL_PREV_TRIGGER;	// prevalence level in the school to trigger school closure
    int PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN;
    int *PERC_SCHOOLS_INFECTED_TRIGGER_VEC;
    int PERC_SCHOOLS_INFECTED_TRIGGER;
    int n_schools_flagged;
    int DURATION_CLOSURE_VEC_LEN;
    int *DURATION_CLOSURE_VEC;
    int DURATION_CLOSURE;
    
    /* Number of runs */
    int HOW_MANY;	// Number of epidemics to run
    
    /* Population */
    int TOT_A; /* Total population size */
    int TOT_C; /* Total population size */
    int TOT;
    double F_A;
    double F_C;
    int N_INIT_INF_A; /* Initial number of infectives, starting the inf period at time 0 */
    int N_INIT_INF_C; /* Initial number of infectives, starting the inf period at time 0 */
    int N_INIT_INF;
    int N_INIT_IMM_A; /* Initial number of immunes */
    int N_INIT_IMM_C; /* Initial number of immunes */
    int N_INIT_IMM;
    int N_NON_SUSC_A; // (N_INIT_INF_A+N_INIT_IMM_A)
    int N_NON_SUSC_C; // (N_INIT_INF_C+N_INIT_IMM_C)
    int N_NON_SUSC;
    
    /* Infection process for the sSIR model */
    int Rg_LEN;
    int Rh_LEN;
    int Rw_LEN;
    double *Rg; /* Average number of contacts from an infective */
    double *Rh; /* Average number of contacts from an infective */
    double *Rw; /* Average number of contacts from an infective */
    /* Rs distribution, with average Rh */
    int Rs_DISTR;	/* The distribution for the length of the infectious period */
    /* 0 --> Dirac delta-function, concentrated on MEAN */
    /* 1 --> exp, with parameter lambda=1/MEAN */
    /* 2 --> gamma, with alpha=ALPHA and lambda=ALPHA/MEAN */
    int Rs_ALPHA; /* Only for the gamma distribution, and only integer */
    double FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A;
    double FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C;
    double FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A;
    double FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C;
    //double CANCEL_RATIO_HOUSEHOLD_TO_A;
    //double CANCEL_RATIO_HOUSEHOLD_TO_C;
    //double CANCEL_RATIO_COMM_TO_A;
    //double CANCEL_RATIO_COMM_TO_C;
    
    int INF_LENGTH_DISTR;	/* The distribution for the length of the infectious period */
    /* 0 --> Dirac delta-function, concentrated on MEAN */
    /* 1 --> exp, with parameter lambda=1/MEAN */
    /* 2 --> gamma, with alpha=ALPHA and lambda=ALPHA/MEAN */
    int INF_MEAN_LEN;
    double *INF_MEAN_VEC;
    double INF_MEAN;
    int INF_ALPHA_LEN;
    int *INF_ALPHA_VEC; /* Only for the gamma distribution, and only integer */
    int INF_ALPHA; /* Only for the gamma distribution, and only integer */
    double BETAg; // (Rg/INF_MEAN) /* Average number of contacts per unit of time */
    double BETAh; // (Rh/INF_MEAN) /* Average number of contacts per unit of time */
    double BETAw; // (Rw/INF_MEAN) /* Average number of contacts per unit of time */
    
    double REL_SUSC_C_VS_A_IN_HOUSEHOLDS; // relative susceptibility of children versus adults
    double REL_INF_C_VS_A_IN_HOUSEHOLDS; // relative infectivity of children versus adults
    double ASS_IN_HOUSEHOLDS; // assortativity in households
    
    int REL_SUSC_C_VS_A_IN_COMM_VEC_LEN;
    double *REL_SUSC_C_VS_A_IN_COMM_VEC;
    double REL_SUSC_C_VS_A_IN_COMM;	// prevalence level in the school to trigger school closure
    int REL_INF_C_VS_A_IN_COMM_VEC_LEN;
    double *REL_INF_C_VS_A_IN_COMM_VEC;
    double REL_INF_C_VS_A_IN_COMM;	// prevalence level in the school to trigger school closure
    int ASS_IN_COMM_VEC_LEN;
    double *ASS_IN_COMM_VEC;
    double ASS_IN_COMM;	// prevalence level in the school to trigger school closure
    int RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC_LEN;
    double *RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC;
    double RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS;
    int RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC_LEN;
    double *RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC;
    double RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM;
    //double REL_SUSC_C_VS_A_IN_COMM; // same thing but in the community. They are considered equal at the beginning
    //double REL_INF_C_VS_A_IN_COMM; // same thing but in the community. They are considered equal at the beginning
    //double ASS_IN_COMM; // assortativity in the community
    //double FRAC_CONT_WITH_C_IN_HOUSEHOLDS; // calculated directly from the rel_susc_C
    //double FRAC_CONT_WITH_C_IN_COMM;
    double FACTOR_REDUCING_GLOBAL_C_TO_C_DUE_TO_SCHOOLS;
    int STRING_DEN_FREQUENCY_DEPENDENT_LEN;
    char *STRING_DEN_FREQUENCY_DEPENDENT; // String telling if frequency dependent transmission in households has n or n-1 as denominator
    double EXPONENT_FREQUENCY_DEPENDENT; // Exponent for the scaling with the household size
    
    /* Latent period */
    int LAT_LENGTH_DISTR;	/* The distribution for the length of the infectious period */
    /* 0 --> Dirac delta-function, concentrated on MEAN */
    /* 1 --> exp, with parameter lambda=1/MEAN */
    /* 2 --> gamma, with alpha=ALPHA and lambda=ALPHA/MEAN */
    double LAT_MEAN;
    int LAT_ALPHA; /* Only for the gamma distribution, and only integer */
    
    /* Time-since-infection model */
    int GEN_TIME_DISTR;	/* The shape of infectivity over time, normalised (R0 specified above) */
    /* 0 --> Dirac delta-function, concentrated on MEAN (don't know if doable)*/
    /* 1 --> exp, with parameter lambda=1/MEAN */
    /* 2 --> gamma, with alpha=ALPHA and lambda=ALPHA/MEAN */
    int GEN_TIME_LEN;
    double *GEN_TIME_VEC;
    double GEN_TIME;
    int GEN_ALPHA_LEN;
    int *GEN_ALPHA_VEC;
    int GEN_ALPHA; /* Only for the gamma distribution, and only integer */
    
    /* Graphics as function of t */
    double dt; /* time step */
    double TIME_STOP; /* top limit for the time (not sharp) */
    int TOT_STEPS; // ( int ) ( TIME_STOP / dt )
    
    /* Generation analysis */
    int MAX_GEN; /* maximum number of generations considered */
    
    double r_MEAN;
    int r_GAP; // ( int ) ( r_MEAN / dt )
    
    /* Random access to the event list */
    double CONST; // A constant of proportionality between dh and 1/TOT
    double dh; // ( CONST / 1000 ) // The time step between two cells of the helping structure
    int TOT_h_STEPS; // ( int ) ( TIME_STOP / dh ) // The length of the helping structure
    
    /* Initial guess for the length of the event list */
    double MAX_EV_FACTOR;
    int MAX_EV;
    
    /* Indicator of the separation between a large and a small epidemic */
    double CUTOFF;
    
    ////////////// Other global variables
    long *idum;
    long seed;
    int epid_index; // This is a global counter for the epidemic running. It is initialised in the Main loop
    int ev_count; /* Counts the number of events used */
    int n_large;
    int n_r, r_length;
    int n_H, n_W;
    int MAX_H_SIZE, MAX_W_SIZE, MAX_WA_SIZE, MAX_WC_SIZE;
    
    char *extinction_flag;//[ HOW_MANY ];
    int *final_size;//[ HOW_MANY ];
    int *final_size_H;//[ HOW_MANY ];
    int *final_size_W;//[ HOW_MANY ];
    double *inc_peak;
    int *t_index_peak;//[ HOW_MANY ];
    double *t_last_event;//[ HOW_MANY ];
    int *n_sync_large;//[ TOT_STEPS+1 ];
    char base_name[64]; // That initial part of the name containing the parameters values
    int cum_I_counter;
    double *t_school_closure;//[ HOW_MANY ];
    double *t_school_reopen;
    double TIME_CLOSURE;
    double TIME_REOPEN;
    int *n_school_closure;//[ TOT_STEPS+1 ];
    double *vector_mean_freq_h_sizes;
    double *vector_sd_freq_h_sizes;
    int *table_h_fs;
    int **table_h_fs_ptr_intermediate;
    int ***table_h_fs_ptr;
    double *table_mean_fs_strat;
    double **table_mean_fs_ptr;
    double *table_sd_fs_strat;
    double **table_sd_fs_ptr;
    int *table_freq_h_sizes;
    int **table_freq_h_sizes_ptr;
    
    
    
    int *gen_cum_large;//[ MAX_GEN+1 ];
    
    int *cum_r;//[ TOT_STEPS+1 ];
    int *r_100;//[ HOW_MANY ];
    int *r_10ghosts;//[ HOW_MANY ];
    
} PAR;
extern PAR P;

typedef struct individual { /* individual-based model */
    int index;		/* from 0 to TOT-1, corresponding to the actual index in the vector */
    char age;		/* A or C */
    int H_index;
    int W_index;
    char is;		/* infectious state: S, I or R */
    double ti;		/* time of infection */
    double tl;		/* time of moving from latent to infectious (if no latency, it is equal to ti */
    double ta;		/* time at which children only are absent from school */
    double tb;		/* time at which children only go back to school */
    double tr;		/* time of recovery */
    int gen;		/* generation of infection */
    int G_conts;	/* number of contacts the individuals makes */
    int G_infs;		/* number of infections the individual generates */
    int H_infs;
    int W_infs;
    int n_conts;	/* G_conts + H_infs + W_infs, only to see when there are global collisions */
    int n_infs;		/* G_infs + H_infs + W_infs */
} IND;
extern IND *ind;	/* from 0 to TOT-1 */

typedef struct event {	/* I create a list of events, adding new events for each infection */
    char heap;			/* Keep track of what needs to be freed */
    double t;			/* time of the event */
    char type;			/* infection (I), from latent to infectious (E) or recovery (R) */
    int i_sub;	/* index of the individual causing the event (subject) */
    int i_obj;	/* index of the individual hit by the event (object) */
    /* for a recovery sub=obj */
    char inf_type;		/* When it is an infection, this stores the type of infection */
    struct event* next;	/* self reference useful to organize the list */
} EVENT;
extern EVENT *ev;

typedef struct helper {
    double t;
    EVENT* p; // pointer to the first event in this h-step
} HELPER;
extern HELPER *h;
extern HELPER *ptr_h;

// When I construct the infectiosu life of a new case, I look for who they infect. This structure helps doing that.
typedef struct obj_type_and_time {
    int index_obj;
    char type; // H, W or G
    double time;
} INFECTEDS;

// Recall that "household" refers to schools, and workplace to households
typedef struct household {
    int index;			/* from 0 to n_H, using calloc */
    double Rs;			/* This is the R0 in schools, which is irrelevant now */
    int size;			/* it's incremented while people are added: the last value is true */
    int *people;
    char is;			/* S, I or R */
    double ti;
    double tc;			/* time when schools close */
    double to;			/* time when schools reopen again */
    double tm;          /* first time of multiple reintroduction */
    double tr;				/* First time at which the household is infected again from outside */
    int preval;			/* Prevalence: counts how many infs there are in the household */
    int abs;			/* Absenteeism level (i.e. prevalence of symptomatic at home for sickness) */
    int ufs;			/* Ultimate final size */
    //	char inf_type;		/* Locally or Globally infected */
    //int counter;		/* Counts up the number of W infected locally and then counts down how many of them recovers */
    //int final_W;		/* Stores only the final value of the number of W infected locally */
    //	int H_local_infs;	/* Counts the number of H infected locally in two steps */
    //	int H_global_infs;	/* Counts the number of H infected globally in one step */
    //int W_infs_min;
    //int W_infs_max;
    //int W_conts_max;
    //int G_infs_min;
    //int G_infs_max;
    //int G_conts_max;
    //int L_infs_min_min;
    //int L_infs_max_max;
    //int L_conts_max_max;
    double NGM_H_AA;
    double NGM_H_AC;
    double NGM_H_CA;
    double NGM_H_CC;
} HOUSEHOLD;
extern HOUSEHOLD* H_list;
extern int *H_people;

typedef struct workplace {
    int index;		/* from 0 to n_W, using calloc */
    int sizeA;			/* it's incremented while people are added: the last value is true */
    int sizeC;			/* it's incremented while people are added: the last value is true */
    int *peopleA;
    int *peopleC;
    char is;			/* S, I or R */
    double ti;
    double tm;          /* first time of multiple reintroduction */
    double tr;				/* First time at which the household is infected again from outside */
    int preval;			/* Prevalence: counts how many infs there are in the household */
    //int fs_max;			/* Final size computed when the local epidemic finishes */
    //int fs_min;			/* Final size computed when the first collision from outside occur */
    int ufs;			/* Ultimate final size */
    //	char inf_type;		/* Locally or Globally infected */
    //int H_infector;		// If of type L, stores the index of the infector household
    // If of type G, stores the index of the households of the first case
    //	int H_local_infs;	/* Counts the number of H infected locally */
    //int H_infs_min;
    //int H_infs_max;
    //int H_conts_max;
    double NGM_W_AA;
    double NGM_W_AC;
    double NGM_W_CA;
    double NGM_W_CC;
} WORKPLACE;
extern WORKPLACE* W_list;
extern int *W_people;

typedef struct generations {
    int number;
    int gen_cum_conts_large;
    int gen_cum_infs_large;
    int C;		/* Number of children	*/
    int A;		/* Number of adults		*/
    int I;		/* Total number of cases in this generation = A+C */
    double TransmToC;
    double TransmToA;
    double propC;
    double propA;
    double R0;
    double Rt;
    int colls;
} GEN;
extern GEN *gen_large;

typedef struct changes {
    int S;
    int I;
    int R;
    int C;			// used to count the cumulative number of children
    int A;
//    int ABSsick;		// Number of children absent from school because sick
//    int ABSclosed;		// Number of children absent from school because schools closed
    int n_conts;	// Only "proj" type
    int n_infs;		// Only "proj" type
    int G_infs;		// Only "proj" type
    int H_infs;		// Only "proj" type
    int W_infs;		// Only "proj" type
    int colls;		// Only "proj" type
    int H_S;
    int H_I;
    int H_C;
    int H_R;
    //int H_fs_min;
    //int H_fs_max;
    int W_S;
    //int W_fs_min;
    //int W_fs_max;
    //int G_infs_min;
    //int G_infs_max;
    //int G_conts_max;
    //int L_infs_min_min;
    //int L_infs_max_max;
    //int L_conts_max_max;
    //int n_H_infs_min;
    //int n_H_infs_max;
    //int n_H_conts_max;
} CHANGES;
extern CHANGES *ch;

typedef struct system_state {	/* to draw S, I and R over time */
    int S;			/* S(t) */
    int I;			/* I(t) */
    int R;			/* R(t) */
    int inc;		/* incidence */
    int inc_C;		/* incidence children */
    int inc_A;		/* incidence adults */
    int cum_inc;	/* cumulative incidence, useful (but not necessary), to use the "weighted" type */
    int cum_inc_C;
    int cum_inc_A;
//    int ABSsick;		// Number of children absent from school because sick
//    int ABSclosed;		// Number of children absent from school because schools closed
    int n_conts; // Note: they are only of "projected" type, as suggested by Christophe
    int n_infs;  // Note: they are only of "projected" type, as suggested by Christophe
    int G_infs;  // Note: they are only of "projected" type, as suggested by Christophe
    int H_infs;  // Note: they are only of "projected" type, as suggested by Christophe
    int W_infs;  // Note: they are only of "projected" type, as suggested by Christophe
    int colls;	 // Note: they are only of "projected" type, as suggested by Christophe
    int H_S;
    int H_I;
    int H_C;
    int H_R;
    int H_inc;
    int H_cum_inc;
    //int H_fs_min;
    //int H_fs_max;
    int W_cum_inc;
    //int W_fs_min;
    //int W_fs_max;
    //int G_infs_min;
    //int G_infs_max;
    //int G_conts_max;
    //int L_infs_min_min;
    //int L_infs_max_max;
    //int L_conts_max_max;
    //int n_H_infs_min;
    //int n_H_infs_max;
    //int n_H_conts_max;
} SYS;
extern SYS *sys_matrix;	/* A row for each run */
extern SYS **sys_line;
extern SYS **sys;


typedef struct summary {
    double t;
    double S;
    double I;
    double R;
    double inc;
    double inc_C;
    double inc_A;
    double cum_inc;
    double cum_inc_C;
    double cum_inc_A;
//    double ABSsick;		// Number of children absent from school because sick
//    double ABSclosed;	// Number of children absent from school because schools closed
    double R0;	// Only of "weighted proj" type
    double Rt;	// Only of "weighted proj" type
    double R_in_G;	// Only of "weighted proj" type
    double R_in_H;	// Only of "weighted proj" type
    double R_in_W;	// Only of "weighted proj" type
    double r;
    double colls;	// Only of "proj" type
    double H_S;
    double H_I;
    double H_C;
    double H_R;
    double H_inc;
    double H_cum_inc;
    //double muH_min;
    //double muH_max;
    double H_r;
    double W_cum_inc;
    //double muW_min;
    //double muW_max;
    //double RG_min;
    //double RG_obs;
    //double RG_max;
    //double RL_min;
    //double RL_obs;
    //double RL_max;
    //double RH_min;
    //double RH_obs;	// This is what can be really observed
    //double RH_max;
    double cum_inc_gap;
    double R0gap;
    double Rtgap;
} SUM;
extern SUM *sync_large; // only sync and large //, global[ TOT_STEPS+1 ], large[ TOT_STEPS+1 ], sync_global[ TOT_STEPS+1 ];

extern double Rg;
extern double Rh;
extern double Rw;
extern double NGM_G_AA;
extern double NGM_G_AC;
extern double NGM_G_CA;
extern double NGM_G_CC;

///////////////////////// Functions declarations /////////////////////////

void epidemic( double a, double b, double c );
void print_output( void ); //char *base_name );

void prepare_model( void );
void load_network( void );
char* create_file_name( char* name, char* beg, char* type );
void read_size_cum( void );
void read_size_distr( void );
void create_cum_distr( void );
int draw_H_size( void );
int draw_W_size( void );
void assign_households( void );
void assign_workplaces( void );
void print_network( void );
void initialyse_basic_quantities( void );


void reset_epidemic_quantities( void );
void initial_infectives( void );
void random_initial_infectives( void );
void infectious_life( int i, double t );
void create_and_place_event( double t, char type, int i_sub, int i_obj, char inf_type );
void place_node( EVENT* p_ev );
void infection_process( void );
void analyse_epidemic( void );
void print_epidemic( void );
void clean_and_free( void );

void store_epidemic( void );
void get_final_size( void );
void create_cum_gen_R0( void );
void create_cum_r( void );

void prepare_summary( void );
void print_final_size( void ); //const char name[] );
void print_generations( void );
void print_real_time( void ); //char* pre_name );
void print_epidemic_curves( void );
void print_r( void );

double ran2( long* idum );
double rnd_time( int distr, double mean, int alpha );
double rnd_unif( double interval );
int rnd_Poiss( double mean );
double expdev(long *idum);
double gamdev(int ia, long *idum);
double poidev( double xm, long *idum );
double bnldev(double pp, int n, long *idum);
double gammln(double k);
int my_round( double x );

////////////////////////////
void effect_of_infection_on_H( EVENT* p );
void effect_of_infection_on_W( EVENT* p );
void effect_of_recovery_on_H( EVENT* p );
void effect_of_recovery_on_W( EVENT* p );

int read_param( int, char** );
int	GetInputParameter( FILE*, char*, char*, void*, int );
char* my_int2str( char* , int, int );
char* my_float2str( char* , double, int, int );

void print_epidemiological_outputs( void );

