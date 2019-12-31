//
//  read.cpp
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

#include "headers.h"

// This function reads the external parameter file "param.txt". It uses low-level functions (maybe deprecated) but it's fast.
// Different types of inputs are read differently by the function "GetInputParameter" below.
// If I read a vector, "GetInputParameter" places what read in a vector and returns it length
// (in case I want to run multiple epidemics with different parameter values).
// Some parameter (the ones that change in every iteration of the Matlab code) are inputted via the command line
// (see CMDL_INPUT_name in "typedecls.h" to decide whether "name" is read from command line or "param.txt" file).
int read_param( int n_cmdl_arg, char **cmdl_arg_list )
{
    int cmdl_arg_index;
    FILE *par;
    
    cmdl_arg_index = 1;
    if (!(par = fopen(cmdl_arg_list[cmdl_arg_index],"r")))
    {
        printf("\nUnable to open the parameter file: %s", cmdl_arg_list[cmdl_arg_index]);
        exit(1);
    }
    else
    {
        cmdl_arg_index = 2;
        if ( cmdl_arg_index > n_cmdl_arg )
            fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        /* Model */
        GetInputParameter( par, "MODEL_TYPE", "%d", ( void* ) &P.MODEL_TYPE, 0 );
        if ( P.MODEL_TYPE == 0 ) // If model is SIR or SEIR, but not a time-since-infection model
        {
            // Quantities useful only for P.MODEL_TYPE == 1 or 2 (i.e. not now)
            P.GEN_TIME_LEN = 1; // This is just to run through a cycle in the main file
            P.GEN_TIME_VEC = ( double* ) malloc( sizeof( double ) );
            *P.GEN_TIME_VEC = 0;
            P.GEN_ALPHA_LEN = 1; // This is just to run through a cycle in the main file
            P.GEN_ALPHA_VEC = ( int* ) malloc( sizeof( int ) );
            *P.GEN_ALPHA_VEC = 0;
            // Quantities useful now. Distributions are always: 0 for delta-peaked, 1 for exponential, 2 for gamma
            GetInputParameter( par, "INF_LENGTH_DISTR", "%d", ( void* ) &P.INF_LENGTH_DISTR, 0 );
            P.INF_MEAN_LEN = GetInputParameter( par, "INF_MEAN", "%f", ( void* ) &P.INF_MEAN_VEC, 1 );
            P.INF_ALPHA_LEN = GetInputParameter( par, "INF_ALPHA", "%d", ( void* ) &P.INF_ALPHA_VEC, 1 );
            GetInputParameter( par, "LATENCE", "%d", ( void* ) &P.LATENCE, 0 );
            if ( P.LATENCE )
            {
                GetInputParameter( par, "LAT_LENGTH_DISTR", "%d", ( void* ) &P.LAT_LENGTH_DISTR, 0 );
                GetInputParameter( par, "LAT_MEAN", "%f", ( void* ) &P.LAT_MEAN, 0 );
                GetInputParameter( par, "LAT_ALPHA", "%d", ( void* ) &P.LAT_ALPHA, 0 );
            }
        }
        else
        {
            // Quantities useful only for P.MODEL_TYPE == 1 or 2 (i.e. not now)
            P.INF_MEAN_LEN = 1; // This is just to run through a cycle in the main file
            P.INF_MEAN_VEC = ( double* ) malloc( sizeof( double ) );
            *P.INF_MEAN_VEC = 0;
            P.INF_ALPHA_LEN = 1; // This is just to run through a cycle in the main file
            P.INF_ALPHA_VEC = ( int* ) malloc( sizeof( int ) );
            *P.INF_ALPHA_VEC = 0;
            // Quantities useful now. Distributions are always: 0 for delta-peaked, 1 for exponential, 2 for gamma
            GetInputParameter( par, "GEN_TIME_DISTR", "%d", ( void* ) &P.GEN_TIME_DISTR, 0 );
            P.GEN_TIME_LEN = GetInputParameter( par, "GEN_TIME", "%f", ( void* ) &P.GEN_TIME_VEC, 1 );
            P.GEN_ALPHA_LEN = GetInputParameter( par, "GEN_ALPHA", "%d", ( void* ) &P.GEN_ALPHA_VEC, 1 );
            GetInputParameter( par, "REC_DELAY", "%f", ( void* ) &P.REC_DELAY, 0 );
            GetInputParameter( par, "DUR_SICK", "%f", ( void* ) &P.DUR_SICK, 0 );
            P.DUR_PRODROMAL = P.DUR_SICK;
            GetInputParameter( par, "DUR_PRODROMAL", "%f", ( void* ) &P.DUR_PRODROMAL, 0 );
            if ( P.DUR_PRODROMAL > P.DUR_SICK ) {
                P.DUR_PRODROMAL = P.DUR_SICK;
                printf( "Prodromal period can't be longer than duration of sickness (from infection to recovery): value modified.\n" );
                printf( "DUR_PRODROMAL: %f\n", P.DUR_PRODROMAL );
            }
            P.PROB_ASYMPT = 0;
            GetInputParameter( par, "PROB_ASYMPT", "%f", ( void* ) &P.PROB_ASYMPT, 0 );
//            P.FRACTION_OF_Rh_IN_SCHOOL_IF_SYMPT = 1;
//            P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A = 1;
//            P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C = 1;
//            P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A = 1;
//            P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C = 1;
//            if ( SWITCH_SCHOOL_CLOSURE )
//            {
//                GetInputParameter( par, "FRACTION_OF_Rh_IN_SCHOOL_IF_SYMPT", "%f", ( void* ) &P.FRACTION_OF_Rh_IN_SCHOOL_IF_SYMPT, 0 );
//                GetInputParameter( par, "FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A", "%f", ( void* ) &P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A, 0 );
//                GetInputParameter( par, "FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C", "%f", ( void* ) &P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C, 0 );
//                GetInputParameter( par, "FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A", "%f", ( void* ) &P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A, 0 );
//                GetInputParameter( par, "FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C", "%f", ( void* ) &P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C, 0 );
//            }
            //P.CANCEL_RATIO_HOUSEHOLD_TO_A = P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A / ( 1 + P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A );
            //P.CANCEL_RATIO_HOUSEHOLD_TO_C = P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C / ( 1 + P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C );
            //P.CANCEL_RATIO_COMM_TO_A = P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A / ( 1 + P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A );
            //P.CANCEL_RATIO_COMM_TO_C = P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C / ( 1 + P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C );
//            P.Rs_DISTR = 0;
//            GetInputParameter( par, "Rs_DISTR", "%d", ( void* ) &P.Rs_DISTR, 0 );
//            P.Rs_ALPHA = 1;
//            GetInputParameter( par, "Rs_ALPHA", "%d", ( void* ) &P.Rs_ALPHA, 0 );
//        }
        
        /* Population */
        //GetInputParameter( par, "TOT_A", "%d", ( void* ) &P.TOT_A, 0 );
        //GetInputParameter( par, "TOT_C", "%d", ( void* ) &P.TOT_C, 0 );
        //P.TOT = P.TOT_A + P.TOT_C;
        GetInputParameter( par, "TOT", "%d", ( void* ) &P.TOT, 0 );
        printf( "TOT: %d\n", P.TOT );
            
        /* Read in the name of the population structure file, to read later on */
        if ( CMDL_INPUT_SOCIAL_STRUCTURE_NETWORK_TYPE )
        {
            P.SOCIAL_STRUCTURE_NETWORK_TYPE = cmdl_arg_list[cmdl_arg_index];
            P.SOCIAL_STRUCTURE_NETWORK_TYPE_LEN = ( int ) strlen(P.SOCIAL_STRUCTURE_NETWORK_TYPE);
            printf( "SOCIAL_STRUCTURE_NETWORK_TYPE: [%d] %s\n", P.SOCIAL_STRUCTURE_NETWORK_TYPE_LEN, P.SOCIAL_STRUCTURE_NETWORK_TYPE );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.SOCIAL_STRUCTURE_NETWORK_TYPE_LEN = GetInputParameter( par, "SOCIAL_STRUCTURE_NETWORK_TYPE", "%s", ( void* ) &P.SOCIAL_STRUCTURE_NETWORK_TYPE, 0 );
        
        /* Epidemic */
        GetInputParameter( par, "HOW_MANY", "%d", ( void* ) &P.HOW_MANY, 0 ); // Total number of simulations
        if ( CMDL_INPUT_Rg ) // Read Rg from command line?
        {
            P.Rg = ( double* ) malloc( sizeof( double ) );
            *P.Rg = atof(cmdl_arg_list[cmdl_arg_index]);
            P.Rg_LEN = 1;
            printf( "Rg: [%d] %f\n", P.Rg_LEN, *P.Rg );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.Rg_LEN = GetInputParameter( par, "Rg", "%f", ( void* ) &P.Rg, 1 );
        if ( CMDL_INPUT_Rh ) // Read the household transmission from command line? (Remember this would be my school)
        {
            P.Rh = ( double* ) malloc( sizeof( double ) );
            *P.Rh = atof(cmdl_arg_list[cmdl_arg_index]);
            P.Rh_LEN = 1;
            printf( "Rh: [%d] %f\n", P.Rh_LEN, *P.Rh );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.Rh_LEN = GetInputParameter( par, "Rh", "%f", ( void* ) &P.Rh, 1 );
        if ( CMDL_INPUT_Rw ) // Read the workplace transmission from command line? (Remember this is my household)
        {
            P.Rw = ( double* ) malloc( sizeof( double ) );
            *P.Rw = atof(cmdl_arg_list[cmdl_arg_index]);
            P.Rw_LEN = 1;
            printf( "Rw: [%d] %f\n", P.Rw_LEN, *P.Rw );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.Rw_LEN = GetInputParameter( par, "Rw", "%f", ( void* ) &P.Rw, 1 );
        
        /* Mixing matrices */
        if ( CMDL_INPUT_SIGMA ) // Read the relative susceptibility from command line?
        {
            P.REL_SUSC_C_VS_A_IN_COMM_VEC = ( double* ) malloc( sizeof( double ) );
            *P.REL_SUSC_C_VS_A_IN_COMM_VEC = atof(cmdl_arg_list[cmdl_arg_index]);
            P.REL_SUSC_C_VS_A_IN_COMM_VEC_LEN = 1;
            printf( "REL_SUSC_C_VS_A_IN_COMM: [%d] %f\n", P.REL_SUSC_C_VS_A_IN_COMM_VEC_LEN, *P.REL_SUSC_C_VS_A_IN_COMM_VEC );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.REL_SUSC_C_VS_A_IN_COMM_VEC_LEN = GetInputParameter( par, "REL_SUSC_C_VS_A_IN_COMM", "%f", ( void* ) &P.REL_SUSC_C_VS_A_IN_COMM_VEC, 1 );
        
        if ( SAME_SIGMA ) // If we assume the relative susceptility is the same in households and the community
            printf( "REL_SUSC_C_VS_A_IN_HOUSEHOLDS: same as for the community\n" );
        else
            GetInputParameter( par, "REL_SUSC_C_VS_A_IN_HOUSEHOLDS", "%f", ( void* ) &P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS, 0 );
        
        if ( CMDL_INPUT_RHO ) // Read the relative infectivity from command line?
        {
            P.REL_INF_C_VS_A_IN_COMM_VEC = ( double* ) malloc( sizeof( double ) );
            *P.REL_INF_C_VS_A_IN_COMM_VEC = atof(cmdl_arg_list[cmdl_arg_index]);
            P.REL_INF_C_VS_A_IN_COMM_VEC_LEN = 1;
            printf( "REL_INF_C_VS_A_IN_COMM: [%d] %f\n", P.REL_INF_C_VS_A_IN_COMM_VEC_LEN, *P.REL_INF_C_VS_A_IN_COMM_VEC );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.REL_INF_C_VS_A_IN_COMM_VEC_LEN = GetInputParameter( par, "REL_INF_C_VS_A_IN_COMM", "%f", ( void* ) &P.REL_INF_C_VS_A_IN_COMM_VEC, 1 );
        
        if ( SAME_RHO ) // If we assume the relative infectivity is the same in households and the community
            printf( "REL_INF_C_VS_A_IN_HOUSEHOLDS: same as for the community\n" );
        else
            GetInputParameter( par, "REL_INF_C_VS_A_IN_HOUSEHOLDS", "%f", ( void* ) &P.REL_INF_C_VS_A_IN_HOUSEHOLDS, 0 );
        
        if ( CMDL_INPUT_PHI ) // Read the global assortativity from command line?
        {
            P.ASS_IN_COMM_VEC = ( double* ) malloc( sizeof( double ) );
            *P.ASS_IN_COMM_VEC = atof(cmdl_arg_list[cmdl_arg_index]);
            P.ASS_IN_COMM_VEC_LEN = 1;
            printf( "ASS_IN_COMM: [%d] %f\n", P.ASS_IN_COMM_VEC_LEN, *P.ASS_IN_COMM_VEC );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            P.ASS_IN_COMM_VEC_LEN = GetInputParameter( par, "ASS_IN_COMM", "%f", ( void* ) &P.ASS_IN_COMM_VEC, 1 );
        
        if ( RANDOM_MIXING_IN_HOUSEHOLDS ) // Assume random mixing in household? Yes, in general...
        {
            printf( "Not reading ASS_IN_HOUSEHOLDS: use random mixing\n" );
            P.ASS_IN_HOUSEHOLDS = -1;
        }
        else // I really want to input a value for the assortativity in households, despite being difficult to define independently of the household size
        {
            if ( SAME_PHI )
                printf( "ASS_IN_HOUSEHOLDS: same as for the community\n" );
            else
                GetInputParameter( par, "ASS_IN_HOUSEHOLDS", "%f", ( void* ) &P.ASS_IN_HOUSEHOLDS, 0 );
        }
        
        if ( CMDL_INPUT_GAMMA ) // Read the ratios of contact between adults and children from command line?
        {
            P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC = ( double* ) malloc( sizeof( double ) );
            *P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC = atof(cmdl_arg_list[cmdl_arg_index]);
            P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC_LEN = 1;
            printf( "RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM: [%d] %f\n", P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC_LEN, *P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
            P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC = ( double* ) malloc( sizeof( double ) );
            *P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC = atof(cmdl_arg_list[cmdl_arg_index]);
            P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC_LEN = 1;
            printf( "RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS: [%d] %f\n", P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC_LEN, *P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
        {
            GetInputParameter( par, "RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS", "%f", ( void* ) &P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS, 0 );
            GetInputParameter( par, "RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM", "%f", ( void* ) &P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM, 0 );
        }
        
        if ( CMDL_INPUT_INIT_INF_A ) // Read the initial number of infected adults from command line?
        {
            //P.N_INIT_INF_A = ( int ) malloc( sizeof( int ) );
            P.N_INIT_INF_A = atoi( cmdl_arg_list[cmdl_arg_index] );
            printf( "N_INIT_INF_A: %d\n", P.N_INIT_INF_A );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            GetInputParameter( par, "N_INIT_INF_A", "%d", ( void* ) &P.N_INIT_INF_A, 0 );
        if ( CMDL_INPUT_INIT_INF_C ) // Read the initial number of infected children from command line?
        {
            //P.N_INIT_INF_C = ( int ) malloc( sizeof( int ) );
            P.N_INIT_INF_C = atoi( cmdl_arg_list[cmdl_arg_index] );
            printf( "N_INIT_INF_C: %d\n", P.N_INIT_INF_C );
            cmdl_arg_index++;
            if ( cmdl_arg_index > n_cmdl_arg )
                fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
        }
        else
            GetInputParameter( par, "N_INIT_INF_C", "%d", ( void* ) &P.N_INIT_INF_C, 0 );
        P.N_INIT_INF = P.N_INIT_INF_A + P.N_INIT_INF_C;
        printf( "N_INIT_INF: %d\n", P.N_INIT_INF );
        GetInputParameter( par, "N_INIT_IMM_A", "%d", ( void* ) &P.N_INIT_IMM_A, 0 );
        GetInputParameter( par, "N_INIT_IMM_C", "%d", ( void* ) &P.N_INIT_IMM_C, 0 );
        P.N_INIT_IMM = P.N_INIT_IMM_A + P.N_INIT_IMM_C;
        printf( "N_INIT_IMM: %d\n", P.N_INIT_IMM );
        P.N_NON_SUSC_A = P.N_INIT_IMM_A + P.N_INIT_INF_A;
        printf( "N_NON_SUSC_A: %d\n", P.N_NON_SUSC_A );
        P.N_NON_SUSC_C = P.N_INIT_IMM_C + P.N_INIT_INF_C;
        printf( "N_NON_SUSC_C: %d\n", P.N_NON_SUSC_C );
        P.N_NON_SUSC = P.N_INIT_IMM + P.N_INIT_INF;
        printf( "N_NON_SUSC: %d\n", P.N_NON_SUSC );
        //P.FRAC_CONT_WITH_C_IN_HOUSEHOLDS = P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS / ( P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS + 1 );
        //P.FRAC_CONT_WITH_C_IN_COMM = P.REL_SUSC_C_VS_A_IN_COMM / ( P.REL_SUSC_C_VS_A_IN_COMM + 1 );
//        GetInputParameter( par, "FACTOR_REDUCING_GLOBAL_C_TO_C_DUE_TO_SCHOOLS", "%f", ( void* ) &P.FACTOR_REDUCING_GLOBAL_C_TO_C_DUE_TO_SCHOOLS, 0 );
        P.STRING_DEN_FREQUENCY_DEPENDENT = NULL;
        P.STRING_DEN_FREQUENCY_DEPENDENT_LEN = GetInputParameter( par, "STRING_DEN_FREQUENCY_DEPENDENT", "%s", ( void* ) &P.STRING_DEN_FREQUENCY_DEPENDENT, 0 );
        GetInputParameter( par, "EXPONENT_FREQUENCY_DEPENDENT", "%f", ( void* ) &P.EXPONENT_FREQUENCY_DEPENDENT, 0 );
        
//        /* Control */
//        if ( SWITCH_SCHOOL_CLOSURE )
//        {
//            if ( CMDL_INPUT_WITHIN_SCHOOL_PREV_TRIGGER )
//            {
//                P.WITHIN_SCHOOL_PREV_TRIGGER_VEC = ( int* ) malloc( sizeof( int ) );
//                *P.WITHIN_SCHOOL_PREV_TRIGGER_VEC = atoi(cmdl_arg_list[cmdl_arg_index]);
//                P.WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN = 1;
//                printf( "WITHIN_SCHOOL_PREV_TRIGGER: [%d] %d\n", P.WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN, *P.WITHIN_SCHOOL_PREV_TRIGGER_VEC );
//                cmdl_arg_index++;
//                if ( cmdl_arg_index > n_cmdl_arg )
//                    fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//            }
//            else
//                P.WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN = GetInputParameter( par, "WITHIN_SCHOOL_PREV_TRIGGER", "%d", ( void* ) &P.WITHIN_SCHOOL_PREV_TRIGGER_VEC, 1 );
//            if ( CMDL_INPUT_PERC_SCHOOLS_INFECTED_TRIGGER )
//            {
//                P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC = ( int* ) malloc( sizeof( int ) );
//                *P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC = atoi(cmdl_arg_list[cmdl_arg_index]);
//                P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN = 1;
//                printf( "PERC_SCHOOLS_INFECTED_TRIGGER: [%d] %d\n", P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN, *P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC );
//                cmdl_arg_index++;
//                if ( cmdl_arg_index > n_cmdl_arg )
//                    fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//            }
//            else
//                P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN = GetInputParameter( par, "PERC_SCHOOLS_INFECTED_TRIGGER", "%d", ( void* ) &P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC, 1 );
//            if ( CMDL_INPUT_DELAY )
//            {
//                if ( SPLIT_DELAY )
//                {
//                    P.DELAY_INF_ONSET = atoi(cmdl_arg_list[cmdl_arg_index]);
//                    printf( "DELAY_INF_ONSET: %f\n", P.DELAY_INF_ONSET );
//                    cmdl_arg_index++;
//                    if ( cmdl_arg_index > n_cmdl_arg )
//                        fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//                    P.DELAY_ONSET_DET = atoi(cmdl_arg_list[cmdl_arg_index]);
//                    printf( "DELAY_ONSET_DET: %f\n", P.DELAY_ONSET_DET );
//                    cmdl_arg_index++;
//                    if ( cmdl_arg_index > n_cmdl_arg )
//                        fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//                    P.DELAY_DET_CONFIRM = atoi(cmdl_arg_list[cmdl_arg_index]);
//                    printf( "DELAY_DET_CONFIRM: %f\n", P.DELAY_DET_CONFIRM );
//                    cmdl_arg_index++;
//                    if ( cmdl_arg_index > n_cmdl_arg )
//                        fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//                    P.DELAY = P.DELAY_INF_ONSET + P.DELAY_ONSET_DET + P.DELAY_DET_CONFIRM;
//                    printf( "DELAY: %f\n", P.DELAY );
//                }
//                else
//                {
//                    P.DELAY_VEC = ( double* ) malloc( sizeof( double ) );
//                    *P.DELAY_VEC = atof(cmdl_arg_list[cmdl_arg_index]);
//                    P.DELAY_VEC_LEN = 1;
//                    printf( "DELAY: [%d] %f\n", P.DELAY_VEC_LEN, *P.DELAY_VEC );
//                    cmdl_arg_index++;
//                    if ( cmdl_arg_index > n_cmdl_arg )
//                        fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//                }
//            }
//            else
//            {
//                if ( SPLIT_DELAY )
//                {
//                    GetInputParameter( par, "DELAY_INF_ONSET", "%f", ( void* ) &P.DELAY_INF_ONSET, 0 );
//                    GetInputParameter( par, "DELAY_ONSET_DET", "%f", ( void* ) &P.DELAY_ONSET_DET, 0 );
//                    GetInputParameter( par, "DELAY_DET_CONFIRM", "%f", ( void* ) &P.DELAY_DET_CONFIRM, 0 );
//                    P.DELAY = P.DELAY_INF_ONSET + P.DELAY_ONSET_DET + P.DELAY_DET_CONFIRM;
//                    printf( "DELAY: %f\n", P.DELAY );
//                }
//                else
//                    P.DELAY_VEC_LEN = GetInputParameter( par, "DELAY", "%f", ( void* ) &P.DELAY_VEC, 1 );
//            }
//            if ( SWITCH_SCHOOL_REOPEN )
//            {
//                if ( CMDL_INPUT_DURATION_CLOSURE )
//                {
//                    P.DURATION_CLOSURE_VEC = ( int* ) malloc( sizeof( int ) );
//                    *P.DURATION_CLOSURE_VEC = atoi(cmdl_arg_list[cmdl_arg_index]);
//                    P.DURATION_CLOSURE_VEC_LEN = 1;
//                    printf( "DURATION_CLOSURE: [%d] %d\n", P.DURATION_CLOSURE_VEC_LEN, *P.DURATION_CLOSURE_VEC );
//                    cmdl_arg_index++;
//                    if ( cmdl_arg_index > n_cmdl_arg )
//                        fprintf( stderr, "Critical error: attepting to read too many inputs from the command line!\n" );
//                }
//                else
//                    P.DURATION_CLOSURE_VEC_LEN = GetInputParameter( par, "DURATION_CLOSURE", "%d", ( void* ) &P.DURATION_CLOSURE_VEC, 1 );
//            }
//            else // schools do not reopen
//            {
//                printf( "No school reopen...\n\t...skipping parameter duration (default = TIME_STOP)..." );
//                P.DURATION_CLOSURE_VEC_LEN = 1;
//                P.DURATION_CLOSURE_VEC = ( int* ) malloc( sizeof( int ) );
//                P.DURATION_CLOSURE_VEC[ 0 ] = P.TIME_STOP;
//                printf( "DURATION_CLOSURE: [%d] %d\n", P.DURATION_CLOSURE_VEC_LEN, *P.DURATION_CLOSURE_VEC );
//            }
//        }
//        else // No school closure
//        {
//            printf( "No school closure...\n\t...skipping all other parameters (default = TOT or TIME_STOP)...\n\n" );
//            P.WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN = 1;
//            P.WITHIN_SCHOOL_PREV_TRIGGER_VEC = ( int* ) malloc( sizeof( int ) );
//            P.WITHIN_SCHOOL_PREV_TRIGGER_VEC[ 0 ] = P.TOT;
//            printf( "WITHIN_SCHOOL_PREV_TRIGGER: [%d] %d\n", P.WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN, *P.WITHIN_SCHOOL_PREV_TRIGGER_VEC );
//            P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN = 1;
//            P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC = ( int* ) malloc( sizeof( int ) );
//            P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC[ 0 ] = P.TOT;
//            printf( "PERC_SCHOOLS_INFECTED_TRIGGER: [%d] %d\n", P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN, *P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC );
//            P.DELAY_VEC_LEN = 1;
//            P.DELAY_VEC = ( double* ) malloc( sizeof( double ) );
//            P.DELAY_VEC[ 0 ] = P.TIME_STOP;
//            printf( "DELAY: [%d] %f\n", P.DELAY_VEC_LEN, *P.DELAY_VEC );
//            P.DURATION_CLOSURE_VEC_LEN = 1;
//            P.DURATION_CLOSURE_VEC = ( int* ) malloc( sizeof( int ) );
//            P.DURATION_CLOSURE_VEC[ 0 ] = P.TIME_STOP;
//            printf( "DURATION_CLOSURE: [%d] %d\n", P.DURATION_CLOSURE_VEC_LEN, *P.DURATION_CLOSURE_VEC );
        }
        
        
        /* Code structure */
        // The epidemic is exact in time, but it's then post-processed for output/plot in discrete time steps dt.
        // TIME_STOP is the largest time in post processing: if the epidemic runs beyond it, it remains exact, but becomes more and more inefficient.
        GetInputParameter( par, "dt", "%f", ( void* ) &P.dt, 0 ); // Time step of epidemic post-processing
        GetInputParameter( par, "TIME_STOP", "%f", ( void* ) &P.TIME_STOP, 0 ); // Final time of epidemic post-processing
        P.TOT_STEPS = ( int ) ( P.TIME_STOP / P.dt );
        printf( "TOT_STEPS: %d\n", P.TOT_STEPS );
        
        // I use a hash table to store the list of events, i.e. I split the list in many small sublists, one for each fixed time step of length dh.
        // CONST is a constant of proportionality between the total population size and the step dh, with decent value chosen but trial and errors.
        GetInputParameter( par, "CONST", "%f", ( void* ) &P.CONST, 0 );
        P.dh = P.CONST / P.TOT;
        printf( "dh: %f\n", P.dh );
        P.TOT_h_STEPS = ( int ) ( P.TIME_STOP / P.dh );
        printf( "TOT_h_STEPS: %d\n", P.TOT_h_STEPS );
            
        // How large should the memory for all events be? I don't know, so I choose MAX_EV_FACTOR * population size. If I go beyond, I need to allocate more memory dynamically.
        GetInputParameter( par, "MAX_EV_FACTOR", "%f", ( void* ) &P.MAX_EV_FACTOR, 0 );
        P.MAX_EV = ( int ) ( P.MAX_EV_FACTOR * P.TOT );
        printf( "MAX_EV: %d\n", P.MAX_EV );
        GetInputParameter( par, "MAX_GEN", "%d", ( void* ) &P.MAX_GEN, 0 );
        
        /* Real-time growth rate */
        if ( P.MODEL_TYPE == 0 )
        {
        	if ( P.LATENCE == 0 )
        		P.r_MEAN = P.INF_MEAN;
        	else
        		P.r_MEAN = ( P.LAT_MEAN + P.INF_MEAN );
        }
        else
        	P.r_MEAN = P.GEN_TIME;
        printf( "r_MEAN: %f\n", P.r_MEAN );
        if ( P.r_MEAN < P.dt ) {
        	printf( "Careful! Some results may be altered because dt is too large\n" );
        	P.r_GAP = P.dt;
        }
        else
        	P.r_GAP = ( int ) ( P.r_MEAN / P.dt );
        printf( "r_GAP: %d\n", P.r_GAP );
        
        // What is considered the cutoff between a small and a large epidemic?
        GetInputParameter( par, "CUTOFF", "%f", ( void* ) &P.CUTOFF, 0 );
    }
    return 0;
}

int	GetInputParameter( FILE *par, char *StringItemName, char *ItemType, void *ItemPtr, int dim )
{
    char NewWord[1024]="", line[1024]="", *ReadItem; //ReadItemName[1024]="", *PtrReadItemName = ReadItemName, ItemName[1024], ;
    int NumItemsRead = 0;
    int i, len, count;
    char *beg, *end;
    double *fReadItem;
    int *iReadItem;
    char *cReadItem;
    char *sReadItem;
    
    fseek( par, 0, 0 );
    while ( ( !NumItemsRead ) && ( !feof( par ) ) )
    {
        fgets( line, 1024, par );
        i = 0;
        while ( ( line[i] == ' ' ) || ( line[i] == '\t' ) )
            i++;
        if ( line[i] == '[' )
        {
            beg = end = line;
            i++;
            while ( ( line[i] == ' ' ) || ( line[i] == '\t' ) )
                i++;
            beg = line + i;
            while ( ( line[i] != '\n' ) && ( line[i] != '%' ) && ( line[i] != ']' ) )
                i++;
            if ( line[i] == ']' )
            {
                i--;
                while ( ( line[i] == ' ' ) || ( line[i] == '\t' ) )
                    i--;
                i++;
                end = line + i;
                len = end - beg;
                if ( !strncmp( StringItemName, beg, len ) )
                {
                    // Starts reading
                    if ( dim == 0 )
                    {
                        if ( !strcmp( ItemType, "%d" ) )
                        {
                            NumItemsRead = fscanf( par, "%d", ( int* ) ItemPtr );
                            printf( "%s: %d\n", StringItemName, *( ( int* ) ItemPtr ) );
                        }
                        else if ( !strcmp( ItemType, "%c" ) )
                        {
                            NumItemsRead = fscanf( par, "%c", ( char* ) ItemPtr );
                            printf( "%s: %c\n", StringItemName, *( ( char* ) ItemPtr ) );
                        }
                        else if ( !strcmp( ItemType, "%f" ) )
                        {
                            NumItemsRead = fscanf( par, "%lf", ( double* ) ItemPtr );
                            printf( "%s: %f\n", StringItemName, *( ( double* ) ItemPtr ) );
                        }
                        else if ( !strcmp( ItemType, "%s" ) )
                        {
                            fgets( line, 1024, par );
                            ReadItem = strtok( line, "% \t\n" );
                            NumItemsRead = strlen( ReadItem );
                            sReadItem = ( char* ) calloc( NumItemsRead, sizeof( char ) );
                            sReadItem = strcpy( sReadItem, ReadItem );
                            *( ( char** ) ItemPtr ) = sReadItem;
                            //							ItemPtr = calloc( NumItemsRead, sizeof( char ) );
                            printf( "%s: [%d] %s\n", StringItemName, NumItemsRead, *( ( char** ) ItemPtr ) );
                        }
                        else
                            fprintf( stderr, "Type of %s does not match\n", StringItemName );
                    }
                    else if ( dim == 1 )
                    {
                        fgets( line, 1024, par );
                        count = 0;
                        i = 0;
                        do
                        {
                            if ( ( line[i] == ' ' ) || ( line[i] == '\t' ) )
                                i++;
                            else if ( ( line[i] != '%' ) && ( line[i] != '\n' ) )
                            {
                                count++;
                                while ( ( line[i] != ' ' ) && ( line[i] != '\t' ) && ( line[i] != '%' ) && ( line[i] != '\0' ) && ( line[i] != '\n' ) )
                                    i++;
                            }
                        } while ( ( line[i] != '%' ) && ( line[i] != '\n' ) && ( line[i] != '\0' ) );
                        if ( !strcmp( ItemType, "%d" ) )
                        {
                            iReadItem = ( int* ) calloc( count, sizeof( int ) );
                            printf( "%s: [%d]", StringItemName, count );
                            ReadItem = strtok( line, " \t%" );
                            for ( i = 0; i < count; i++ )
                            {
                                NumItemsRead += sscanf( ReadItem, "%d", ( iReadItem ) + i ); 
                                printf( " %d", *( iReadItem + i ) );
                                //( ( double* ) ItemPtr )++;
                                ReadItem = strtok( NULL, " \t%" );
                            }
                            *( ( int** ) ItemPtr ) = iReadItem; 
                            printf( "\n" );
                        }
                        else if ( !strcmp( ItemType, "%c" ) )
                        {
                            cReadItem = ( char* ) calloc( count, sizeof( char ) );
                            printf( "%s: [%d]", StringItemName, count );
                            ReadItem = strtok( line, " \t%" );
                            for ( i = 0; i < count; i++ )
                            {
                                NumItemsRead += sscanf( ReadItem, "%c", ( cReadItem ) + i ); 
                                printf( " %c", *( cReadItem + i ) );
                                //( ( double* ) ItemPtr )++;
                                ReadItem = strtok( NULL, " \t%" );
                            }
                            *( ( char** ) ItemPtr ) = cReadItem; 
                            printf( "\n" );
                        }
                        else if ( !strcmp( ItemType, "%f" ) )
                        {
                            fReadItem = ( double* ) calloc( count, sizeof( double ) );
                            printf( "%s: [%d]", StringItemName, count );
                            ReadItem = strtok( line, " \t%" );
                            for ( i = 0; i < count; i++ )
                            {
                                NumItemsRead += sscanf( ReadItem, "%lf", ( fReadItem ) + i ); 
                                printf( " %f", *( fReadItem + i ) );
                                //( ( double* ) ItemPtr )++;
                                ReadItem = strtok( NULL, " \t%" );
                            }
                            *( ( double** ) ItemPtr ) = fReadItem; 
                            printf( "\n" );
                        }
                        else
                            fprintf( stderr, "Type of %s does not match\n", StringItemName );
                    }
                }
            }
        }
    }
    return NumItemsRead;
}


