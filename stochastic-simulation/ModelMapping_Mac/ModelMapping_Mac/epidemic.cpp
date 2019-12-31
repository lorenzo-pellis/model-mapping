//
//  epidemic.cpp
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

#include "headers.h"


void epidemic( double a, double b, double c )
{
    Rg = a;
    Rh = b;
    Rw = c;
    NGM_G_AA = Rg * ( P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM - ( 1 - P.ASS_IN_COMM ) * ( double ) P.TOT_C / P.TOT_A );
    if ( NGM_G_AA < 0 )
    {
        if ( NGM_G_AA > -0.0001 )
        {
            printf("Assortativity in imprecise but still viable!");
            NGM_G_AA = 0;
        } else {
            printf( "Problem, specified assortativity in the community is incompatible with population structure and contact patterns.\n" );
            printf( "Revert back to random mixing to avoid disasters!\n" );
            P.ASS_IN_COMM = ( double ) 1 / ( 1 + P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM * P.TOT_A / P.TOT_C );
            NGM_G_AA = Rg * ( P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM - ( 1 - P.ASS_IN_COMM ) * ( double ) P.TOT_C / P.TOT_A );
        }
    }
    NGM_G_CA = Rg * P.REL_SUSC_C_VS_A_IN_COMM * ( 1 - P.ASS_IN_COMM ) * ( double ) P.TOT_C / P.TOT_A;
    NGM_G_AC = Rg * ( 1 - P.ASS_IN_COMM ) * P.REL_INF_C_VS_A_IN_COMM;
    NGM_G_CC = Rg * P.REL_SUSC_C_VS_A_IN_COMM * P.REL_INF_C_VS_A_IN_COMM * P.ASS_IN_COMM;
    //    NGM_G_AC = P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A * Rg * ( 1 - P.ASS_IN_COMM ) * P.REL_INF_C_VS_A_IN_COMM;
    //    NGM_G_CC = P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C * Rg * P.REL_SUSC_C_VS_A_IN_COMM * P.REL_INF_C_VS_A_IN_COMM * P.ASS_IN_COMM;

    reset_epidemic_quantities();
    //initial_infectives(); // Non-random: choose either this or the next one (next one is better)
    random_initial_infectives();
    infection_process();
    //print_epidemic(); // Print the list of events, for debugging purposes
    analyse_epidemic(); 
}

void reset_epidemic_quantities( void ) /* Initialization of the individual and population states */
{
    int i, s, e, k, u, w;
    int size = 0;

    for ( s = 0; s < P.n_H; s++ ) {
        H_list[ s ].Rs = rnd_time( P.Rs_DISTR, Rh, P.Rs_ALPHA );
    }
    printf( "\n" );
    
    P.ev_count = 0;
    P.cum_I_counter = P.N_INIT_INF;
    P.TIME_CLOSURE = P.TIME_STOP;
    P.TIME_REOPEN = P.TIME_STOP;
    P.n_schools_flagged = 0;
    
    
    for ( i = 0; i < P.TOT; i++ ) {
        ind[i].index = i;
        ind[i].is = 'S';
        ind[i].ti = -1;
        ind[i].tl = -1;
//        ind[i].ta = -1;
//        ind[i].tb = -1;
        ind[i].tr = -1;
        ind[i].gen = -1;
        ind[i].G_conts = 0;
        ind[i].G_infs = 0;
        ind[i].H_infs = 0;
        ind[i].W_infs = 0;
    }
    
    /* Reset vector changes */
    for ( k = 0; k < P.TOT_STEPS; k++ ) {
        ch[ k ].S = 0;
        ch[ k ].I = 0;
        ch[ k ].R = 0;
        ch[ k ].C = 0;
        ch[ k ].A = 0;
//        ch[ k ].ABSsick = 0;
//        ch[ k ].ABSclosed = 0;
        ch[ k ].n_conts = 0;
        ch[ k ].n_infs = 0;
        ch[ k ].G_infs = 0;
        ch[ k ].H_infs = 0;
        ch[ k ].W_infs = 0;
        ch[ k ].colls = 0;
        ch[ k ].H_S = 0;
        ch[ k ].H_I = 0;
        ch[ k ].H_C = 0;
        ch[ k ].H_R = 0;
        ch[ k ].W_S = 0;
    }
    
    /* Reset event list */
    for ( e = 0; e < P.MAX_EV; e++ ) {
        // the heap 'N' is set in main
        ev[ e ].t = -1;
        ev[ e ].type = 0;
        ev[ e ].i_sub = -1;
        ev[ e ].i_obj = -1;
        ev[ e ].inf_type = 0;
        ev[ e ].next = NULL;
    }
    
    /* Household list */
    for ( u = 0; u < P.n_H; u++ ) {
        H_list[ u ].is = 'S';
        H_list[ u ].ti = -1;
//        H_list[ u ].tc = -1;
//        H_list[ u ].to = -1;
        H_list[ u ].tm = -1;
        H_list[ u ].tr = -1;
        H_list[ u ].preval = 0;
//        H_list[ u ].abs = 0;
        H_list[ u ].ufs = 0;
    }
    
    /* Workplace list */
    for ( w = 0; w < P.n_W; w++ ) {
        W_list[ w ].is = 'S';
        W_list[ w ].ti = -1;
        W_list[ w ].tm = -1;
        W_list[ w ].tr = -1;
        W_list[ w ].preval = 0;
        W_list[ w ].ufs = 0;
        size = W_list[ w ].sizeA + W_list[ w ].sizeC;
        if ( RANDOM_MIXING_IN_HOUSEHOLDS )
        {
            P.ASS_IN_HOUSEHOLDS = ( double ) ( W_list[ w ].sizeC - 1 ) / ( ( W_list[ w ].sizeC - 1 ) + P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS * W_list[ w ].sizeA );
        }
        
        if ( size <= 1 ) { // There is only 1 adult or 1 children
            W_list[ w ].NGM_W_AA = 0;
            W_list[ w ].NGM_W_AC = 0;
            W_list[ w ].NGM_W_CA = 0;
            W_list[ w ].NGM_W_CC = 0;
        } else { // if there are at least 2 people...
            if ( W_list[ w ].sizeA == 0 ) {
                W_list[ w ].NGM_W_AA = 0;
                W_list[ w ].NGM_W_AC = 0;
                W_list[ w ].NGM_W_CA = 0;
                // W_list[ w ].NGM_W_CC = P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C * Rw * P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS * P.REL_INF_C_VS_A_IN_HOUSEHOLDS * P.ASS_IN_HOUSEHOLDS;
                W_list[ w ].NGM_W_CC = Rw * P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS * P.REL_INF_C_VS_A_IN_HOUSEHOLDS * P.ASS_IN_HOUSEHOLDS;
            } else if ( W_list[ w ].sizeC == 0 ) {
                W_list[ w ].NGM_W_AA = Rw * ( P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS - ( 1 - P.ASS_IN_HOUSEHOLDS ) * ( double ) W_list[ w ].sizeC / W_list[ w ].sizeA );
                W_list[ w ].NGM_W_AC = 0;
                W_list[ w ].NGM_W_CA = 0;
                W_list[ w ].NGM_W_CC = 0;
            } else { // if there are at least one adult and one children
                // W_list[ w ].NGM_W_AC = P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A * Rw * ( 1 - P.ASS_IN_HOUSEHOLDS ) * P.REL_INF_C_VS_A_IN_HOUSEHOLDS;
                W_list[ w ].NGM_W_AC = Rw * ( 1 - P.ASS_IN_HOUSEHOLDS ) * P.REL_INF_C_VS_A_IN_HOUSEHOLDS;
                W_list[ w ].NGM_W_CA = Rw * P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS * ( 1 - P.ASS_IN_HOUSEHOLDS ) * ( double ) W_list[ w ].sizeC / W_list[ w ].sizeA;
                if ( W_list[ w ].sizeA == 1 )
                    W_list[ w ].NGM_W_AA = 0;
                else
                    W_list[ w ].NGM_W_AA = Rw * ( P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS - ( 1 - P.ASS_IN_HOUSEHOLDS ) * ( double ) W_list[ w ].sizeC / W_list[ w ].sizeA );
                if ( W_list[ w ].sizeC == 1 )
                    W_list[ w ].NGM_W_CC = 0;
                else {
                    // W_list[ w ].NGM_W_CC = P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C * Rw * P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS * P.REL_INF_C_VS_A_IN_HOUSEHOLDS * P.ASS_IN_HOUSEHOLDS;
                    W_list[ w ].NGM_W_CC = Rw * P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS * P.REL_INF_C_VS_A_IN_HOUSEHOLDS * P.ASS_IN_HOUSEHOLDS;
                }
            }
            if ( P.EXPONENT_FREQUENCY_DEPENDENT != 1 )
            {
                if ( strcmp( P.STRING_DEN_FREQUENCY_DEPENDENT, "n-1" ) ) // The string is NOT n-1, so the denominator is put to n by default
                {
                    W_list[ w ].NGM_W_AA *= W_list[ w ].sizeA / pow( W_list[ w ].sizeA, P.EXPONENT_FREQUENCY_DEPENDENT );
                    W_list[ w ].NGM_W_AC *= W_list[ w ].sizeA / pow( W_list[ w ].sizeA, P.EXPONENT_FREQUENCY_DEPENDENT );
                    W_list[ w ].NGM_W_CA *= W_list[ w ].sizeC / pow( W_list[ w ].sizeC, P.EXPONENT_FREQUENCY_DEPENDENT );
                    W_list[ w ].NGM_W_CC *= W_list[ w ].sizeC / pow( W_list[ w ].sizeC, P.EXPONENT_FREQUENCY_DEPENDENT );
                } else { // denominator is n-1. The external factor remains sizeA or sizeC, because then I select the person to contact among all, including the infector
                    W_list[ w ].NGM_W_AA *= W_list[ w ].sizeA / pow( W_list[ w ].sizeA - 1, P.EXPONENT_FREQUENCY_DEPENDENT );
                    W_list[ w ].NGM_W_AC *= W_list[ w ].sizeA / pow( W_list[ w ].sizeA, P.EXPONENT_FREQUENCY_DEPENDENT );
                    W_list[ w ].NGM_W_CA *= W_list[ w ].sizeC / pow( W_list[ w ].sizeC, P.EXPONENT_FREQUENCY_DEPENDENT );
                    W_list[ w ].NGM_W_CC *= W_list[ w ].sizeC / pow( W_list[ w ].sizeC - 1, P.EXPONENT_FREQUENCY_DEPENDENT );
                }
            }
        }
    }
}


void place_node( EVENT* ptr_ev ) {
    int h_index;
    EVENT *curr, *prev, *temp;
    
    h_index = ( int ) ( ptr_ev->t / P.dh ) ;
    if ( h_index > P.TOT_h_STEPS )
        h_index = P.TOT_h_STEPS;
    if ( h[ h_index ].p == NULL )
        h[ h_index ].p = ptr_ev;
    else {
        if ( h[ h_index ].p->t > ptr_ev->t ) {
            temp = h[ h_index ].p;
            h[ h_index ].p = ptr_ev;
            ptr_ev->next = temp;
        } else {
            curr = h[ h_index ].p->next;
            prev = h[ h_index ].p;
            while ( curr != NULL ) {
                if ( curr->t < ptr_ev->t ) {
                    prev = curr;
                    curr = curr->next;
                } else {
                    prev->next = ptr_ev;
                    ptr_ev->next = curr;
                    return;
                }
            }
            assert( curr == NULL );
            prev->next = ptr_ev;
        }
    }
}

void random_initial_infectives( void ) // Efficient only for almost virgin populations
{
    int index, c = 0;
    
    while ( c < P.N_INIT_INF_C ) {
        index = ( int )( P.TOT_C * ran2( P.idum ) ); /* Initial cases, randomly picked */
        if ( ind[ index ].is == 'S' ) {
            c++;
            if ( P.MODEL_TYPE == 0 && P.LATENCE == 1 )
                ind[ index ].is = 'E';
            else
                ind[ index ].is = 'I';
            ind[ index ].gen = 0;
            ind[ index ].ti = 0;
            ind[ index ].tl = 0; /* Changed later if there is a latent period in the model */
            infectious_life( index, 0 );
            // Household stuff
            H_list[ ind[ index ].H_index ].is = 'I';
            H_list[ ind[ index ].H_index ].ti = 0;
            H_list[ ind[ index ].H_index ].preval++;
            H_list[ ind[ index ].H_index ].ufs++;
            // Workplaces stuff
            W_list[ ind[ index ].W_index ].is = 'I';
            W_list[ ind[ index ].W_index ].ti = 0;
            W_list[ ind[ index ].W_index ].preval++;
            W_list[ ind[ index ].W_index ].ufs++;
        }
    }
    c = 0;
    while ( c < P.N_INIT_IMM_C ) {
        index = ( int )( P.TOT_C * ran2( P.idum ) ); /* Initial immunes, randomly picked */
        if ( ind[ index ].is == 'S' ) {
            c++;
            ind[ index ].is = 'R';
            ind[ index ].gen = -1;
            ind[ index ].ti = -2; // Just to distinguish them from people who escape infection
            ind[ index ].tl = -2; // Just to distinguish them from people who escape infection
            ind[ index ].tr = -2; // Just to distinguish them from people who escape infection
        }
    }
    /* Adults */
    while ( c < P.N_INIT_INF_A ) {
        index = P.TOT_C + ( int )( P.TOT_A * ran2( P.idum ) ); /* Initial cases, randomly picked */
        assert( ind[ index ].age == 'A' );
        if ( ind[ index ].is == 'S' ) {
            c++;
            if ( P.MODEL_TYPE == 0 && P.LATENCE == 1 )
                ind[ index ].is = 'E';
            else
                ind[ index ].is = 'I';
            ind[ index ].gen = 0;
            ind[ index ].ti = 0;
            ind[ index ].tl = 0; /* Changed later if there is a latent period in the model */
            infectious_life( index, 0 );
            // Adults don't go to school (remember households are schools!)
            //H_list[ ind[ index ].H_index ].is = 'I';
            //H_list[ ind[ index ].H_index ].ti = 0;
            //H_list[ ind[ index ].H_index ].preval++;
            //H_list[ ind[ index ].H_index ].ufs++;
            // Workplaces stuff
            W_list[ ind[ index ].W_index ].is = 'I';
            W_list[ ind[ index ].W_index ].ti = 0;
            W_list[ ind[ index ].W_index ].preval++;
            W_list[ ind[ index ].W_index ].ufs++;
        }
    }
    c = 0;
    while ( c < P.N_INIT_IMM_A ) {
        index = P.TOT_C + ( int )( P.TOT_A * ran2( P.idum ) ); /* Initial immunes, randomly picked */
        if ( ind[ index ].is == 'S' ) {
            c++;
            ind[ index ].is = 'R';
            ind[ index ].gen = -1;
            ind[ index ].ti = -2; // Just to distinguish them from people who escape infection
            ind[ index ].tl = -2; // Just to distinguish them from people who escape infection
            ind[ index ].tr = -2; // Just to distinguish them from people who escape infection
            // Think about immune inds in households and workplaces
        }
    }
}

void initial_infectives( void ) // Non-random initialization
{
    int c;
    
    for ( c = 0; c < P.N_INIT_INF_C; c++ ) {
        assert( ind[ c ].is == 'S' );
        if ( P.MODEL_TYPE == 0 && P.LATENCE == 1 )
            ind[ c ].is = 'E';
        else
            ind[ c ].is = 'I';
        ind[ c ].gen = 0;
        ind[ c ].ti = 0;
        ind[ c ].tl = 0; /* Changed later if there is a latent period in the model */
        infectious_life( c, 0 );
        // Household stuff
        H_list[ ind[ c ].H_index ].is = 'I';
        H_list[ ind[ c ].H_index ].ti = 0;
        H_list[ ind[ c ].H_index ].preval++;
        H_list[ ind[ c ].H_index ].ufs++;
        // Workplaces stuff
        W_list[ ind[ c ].W_index ].is = 'I';
        W_list[ ind[ c ].W_index ].ti = 0;
        W_list[ ind[ c ].W_index ].preval++;
        W_list[ ind[ c ].W_index ].ufs++;
        //W_list[ ind[ c ].W_index ].H_infector = ind[ c ].H_index;
    }
    for ( ; c < P.N_NON_SUSC_C; c++ ) {
        assert( ind[ c ].is == 'S' );
        ind[ c ].is = 'R';
        ind[ c ].gen = -2;
        ind[ c ].ti = -2; // Just to distinguish them from people who escape infection
        ind[ c ].tl = -2; // Just to distinguish them from people who escape infection
        ind[ c ].tr = -2; // Just to distinguish them from people who escape infection
    }
    for ( c = P.TOT_C; c < P.TOT_C + P.N_INIT_INF_A; c++ ) {
        assert( ind[ c ].is == 'S' );
        assert( ind[ c ].age == 'A' );
        if ( P.MODEL_TYPE == 0 && P.LATENCE == 1 )
            ind[ c ].is = 'E';
        else
            ind[ c ].is = 'I';
        ind[ c ].gen = 0;
        ind[ c ].ti = 0;
        ind[ c ].tl = 0; /* Changed later if there is a latent period in the model */
        infectious_life( c, 0 );
        // Adults don't go to school (remember households are schools!)
        //H_list[ ind[ c ].H_index ].is = 'I';
        //H_list[ ind[ c ].H_index ].ti = 0;
        //H_list[ ind[ c ].H_index ].preval++;
        //H_list[ ind[ c ].H_index ].ufs++;
        // Workplaces stuff
        W_list[ ind[ c ].W_index ].is = 'I';
        W_list[ ind[ c ].W_index ].ti = 0;
        W_list[ ind[ c ].W_index ].preval++;
        W_list[ ind[ c ].W_index ].ufs++;
        //W_list[ ind[ c ].W_index ].H_infector = ind[ c ].H_index;
    }
}

// This is one of the key functions: as soon as an individual is infected, I project their future infectious life, create all infectious events (or recovery, etc.) and store them in the list of events.
void infectious_life( int index_sub, double basic_time )
{
    double inf_length, lat_length = 0;
    double hitting_time;
    int G_num_conts, index_ind_hit, x, y, k = 0;
    int H_num_conts, H_ind, H_index_ind_hit, W_num_conts, W_ind, W_index_ind_hit;
    int W_num_conts_to_C, W_num_conts_to_A, G_num_conts_to_C, G_num_conts_to_A;
    INFECTEDS* infecteds;
    double last_inf = 0; // Just to decide which is the last infectious contact, for rec in tsi model
    double lat_time = 0, rec_time;
    
    H_ind = ind[ index_sub ].H_index;
    W_ind = ind[ index_sub ].W_index;
    
    if ( P.MODEL_TYPE == 0 ) { // If we use the SIR or SEIR model (MODEL_TYPE == 0)
        if ( P.LATENCE == 1 ) {
            lat_length = rnd_time( P.LAT_LENGTH_DISTR, P.LAT_MEAN, P.LAT_ALPHA ); /* function drawing the length of the latent period from an already known distribution */
            lat_time = basic_time+lat_length;
            create_and_place_event( lat_time, 'E', index_sub, -1 , '0' );
        }
        inf_length = rnd_time( P.INF_LENGTH_DISTR, P.INF_MEAN, P.INF_ALPHA ); /* function drawing the length of the infectious period from an already known distribution */
        rec_time = basic_time+lat_length+inf_length;
        
        // Place the recovery in the event list
        create_and_place_event( rec_time, 'R', index_sub, -2, '0' );
        
        // Decide how many infectious contacts the individual makes
        if ( ind[ index_sub ].age == 'A' ) { // The subject is an adult
            // Workplace i.e. true households
            if ( W_list[ W_ind ].NGM_W_AA == 0 )
                W_num_conts_to_A = 0;
            else
                W_num_conts_to_A = rnd_Poiss( inf_length * W_list[ W_ind ].NGM_W_AA / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            if ( W_list[ W_ind ].NGM_W_CA == 0 )
                W_num_conts_to_C = 0;
            else
                W_num_conts_to_C = rnd_Poiss( inf_length * W_list[ W_ind ].NGM_W_CA / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            // Community
            if ( NGM_G_AA == 0 )
                G_num_conts_to_A = 0;
            else
                G_num_conts_to_A = rnd_Poiss( inf_length * NGM_G_AA / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            if ( NGM_G_CA == 0 )
                G_num_conts_to_C = 0;
            else
                G_num_conts_to_C = rnd_Poiss( inf_length * NGM_G_CA / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            
            // Household, i.e. true school (adults don't have a school)
            assert( H_ind == -1 );
            H_num_conts = 0;
        
        } else { // The subject is a child
            assert( ind[ index_sub ].age == 'C' );
            // Workplace i.e. true households
            if ( W_list[ W_ind ].NGM_W_AC == 0 )
                W_num_conts_to_A = 0;
            else
                W_num_conts_to_A = rnd_Poiss( inf_length * W_list[ W_ind ].NGM_W_AC / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            if ( W_list[ W_ind ].NGM_W_CC == 0 )
                W_num_conts_to_C = 0;
            else
                W_num_conts_to_C = rnd_Poiss( inf_length * W_list[ W_ind ].NGM_W_CC / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            // Community
            if ( NGM_G_AC == 0 )
                G_num_conts_to_A = 0;
            else
                G_num_conts_to_A = rnd_Poiss( inf_length * NGM_G_AC / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            if ( NGM_G_CC == 0 )
                G_num_conts_to_C = 0;
            else
                G_num_conts_to_C = rnd_Poiss( inf_length * NGM_G_CC / P.INF_MEAN ); /* Number of contacts	is Poisson distributed with random mean */
            
            // Household, i.e. true school (children do have a school)
            if ( Rh == 0 )
                H_num_conts = 0;
            else {
                if ( ( H_list[ H_ind ].tc > 0 ) && ( basic_time > H_list[ H_ind ].tc ) && ( basic_time < H_list[ H_ind ].to ) )
                    H_num_conts = 0;
                else
                    H_num_conts = rnd_Poiss( inf_length * Rh / P.INF_MEAN );
            }
        }
        
    } else { // Time-since-infection model
        if ( ind[ index_sub ].age == 'A' ) { // The subject is an adult
            // Workplace i.e. true households
            if ( W_list[ W_ind ].NGM_W_AA == 0 )
                W_num_conts_to_A = 0;
            else
                W_num_conts_to_A = rnd_Poiss( W_list[ W_ind ].NGM_W_AA ); /* Number of contacts	is Poisson distributed with random mean */
            if ( W_list[ W_ind ].NGM_W_CA == 0 )
                W_num_conts_to_C = 0;
            else
                W_num_conts_to_C = rnd_Poiss( W_list[ W_ind ].NGM_W_CA ); /* Number of contacts	is Poisson distributed with random mean */
            // Community
            if ( NGM_G_AA == 0 )
                G_num_conts_to_A = 0;
            else
                G_num_conts_to_A = rnd_Poiss( NGM_G_AA ); /* Number of contacts	is Poisson distributed with random mean */
            if ( NGM_G_CA == 0 )
                G_num_conts_to_C = 0;
            else
                G_num_conts_to_C = rnd_Poiss( NGM_G_CA ); /* Number of contacts	is Poisson distributed with random mean */
            
            // Household, i.e. true school (adults don't have a school)
            assert( H_ind == -1 );
            H_num_conts = 0;
            
        } else { // The subject is a child
            assert( ind[ index_sub ].age == 'C' );
            // Workplace i.e. true households
            if ( W_list[ W_ind ].NGM_W_AC == 0 )
                W_num_conts_to_A = 0;
            else
                W_num_conts_to_A = rnd_Poiss( W_list[ W_ind ].NGM_W_AC ); /* Number of contacts	is Poisson distributed with random mean */
            if ( W_list[ W_ind ].NGM_W_CC == 0 )
                W_num_conts_to_C = 0;
            else
                W_num_conts_to_C = rnd_Poiss( W_list[ W_ind ].NGM_W_CC ); /* Number of contacts	is Poisson distributed with random mean */
            // Community
            if ( NGM_G_AC == 0 )
                G_num_conts_to_A = 0;
            else
                G_num_conts_to_A = rnd_Poiss( NGM_G_AC ); /* Number of contacts	is Poisson distributed with random mean */
            if ( NGM_G_CC == 0 )
                G_num_conts_to_C = 0;
            else
                G_num_conts_to_C = rnd_Poiss( NGM_G_CC ); /* Number of contacts	is Poisson distributed with random mean */

            // Household, i.e. true school (children do have a school)
            if ( Rh == 0 )
                H_num_conts = 0;
 ////// This is possibly wrong: the next line should be uncommented. However, that changes the random numbers so I'm not sure I want to correct it. It makes no difference - just draws a random number when not needed.
            //else {
            //	if ( ( H_list[ H_ind ].tc > 0 ) && ( basic_time > H_list[ H_ind ].tc ) && ( basic_time < H_list[ H_ind ].to ) )
            //		H_num_conts = 0;
            //	else
            //	{
            if ( ran2(P.idum) > P.PROB_ASYMPT ) // The child is symptomatic
            {
                ind[ index_sub ].ta = basic_time + P.DUR_PRODROMAL;
                create_and_place_event( ind[ index_sub ].ta, 'A', index_sub, index_sub, '0' );
                ind[ index_sub ].tb = basic_time + P.DUR_SICK;
                H_num_conts = rnd_Poiss( H_list[ H_ind ].Rs * P.FRACTION_OF_Rh_IN_SCHOOL_IF_SYMPT );
            } else { // Asymptomatic: they don't even stay at home
                H_num_conts = rnd_Poiss( H_list[ H_ind ].Rs );
                //	}
                //}
            }
        }
    }
    
    G_num_conts = G_num_conts_to_C + G_num_conts_to_A;
    W_num_conts = W_num_conts_to_C + W_num_conts_to_A;
    
    /* We now study the individual's infectivity in H, W and G (global). There is a unique list
     of infecteds since he could infect the same ind both in H and W or G... */
    infecteds = ( INFECTEDS* ) calloc( H_num_conts + W_num_conts + G_num_conts , sizeof( INFECTEDS ) ); // k is the index used to move here
    
    /* Add the H infections to the infected list */
    for ( x = 1; x <= H_num_conts; x++ ) {
        do { // infectivity is towards other individuals, but the households has a list with all members. If the individual infect him/herself, draw another individual - this might be inefficient for very small sizes
            H_index_ind_hit = ( int ) ( H_list[ H_ind ].size * ran2(P.idum) );
            index_ind_hit = H_list[ H_ind ].people[ H_index_ind_hit ];
        } while ( index_ind_hit == index_sub );
        if ( ind[ index_ind_hit ].is  == 'S' ) {
            ind[ index_ind_hit ].is = 'X'; // Mark individual "destined to be infected"
            infecteds[k].index_obj = index_ind_hit;
            infecteds[k].type = 'H'; // It's a household infectious contact
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                infecteds[k].time = rnd_unif( inf_length );
            else // Time-since-infection models
                if ( ind[ index_sub ].ta > 0 ) // Symptomatic case
                    do {
                        infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                    } while ( infecteds[k].time > P.DUR_PRODROMAL );
                else // Asymptomatic
                    if ( P.MODEL_TYPE == 1 ) // Truncated gen time distr
                        do {
                            infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                        } while ( infecteds[k].time > P.DUR_SICK );
                    else // Time-since-infection model with full gen time distr
                        infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            k++;
        }
        else if ( ind[ index_ind_hit ].is == 'X' ) { // If the individual was marked already, store new time only if sooner
            /* If the ind hit is X, it means that he was already infected before TIME_CLOSURE
             and therefore:	1) if they try to infect him later, I don't store it;
             2) if they want to infect him before, it is automatically before TIME_CLOSURE
             So, no need to do anything */
            y = 0;
            while ( infecteds[y].index_obj != index_ind_hit )
                y++;
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                hitting_time = rnd_unif( inf_length );
            else // Time-since-infection models
                //if ( P.MODEL_TYPE == 1 ) // Truncated gen time distr
                //	do {
                //		hitting_time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                //	} while ( hitting_time > P.DUR_SICK );
                //else // Time-since-infection model with full gen time distr
                // // This part above is useless, as I store the new time of infection
                // // only when it is smaller than another one (already smaller than DUR_SICK)
                hitting_time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            if ( hitting_time < infecteds[y].time ) {
                infecteds[y].time = hitting_time;
                infecteds[y].type = 'H';
            }
        }
        else
            assert( ind[ index_ind_hit ].is == 'E' || ind[ index_ind_hit ].is == 'I' || ind[ index_ind_hit ].is == 'R' );
    }
    /* Add the W infections to children to the infected list */
    for ( x = 1; x <= W_num_conts_to_C; x++ ) {
        do { // infectivity is towards other children, but the workplace has a list with all members. If the individual infect him/herself, draw another child - this might be inefficient for very small sizes
            W_index_ind_hit = ( int ) ( W_list[ W_ind ].sizeC * ran2(P.idum) );
            index_ind_hit = W_list[ W_ind ].peopleC[ W_index_ind_hit ];
            assert( ind[ index_ind_hit ].age == 'C' );
        } while ( index_ind_hit == index_sub );
        if ( ind[ index_ind_hit ].is  == 'S' ) {
            ind[ index_ind_hit ].is = 'X'; // Mark individual "destined to be infected"
            infecteds[k].index_obj = index_ind_hit;
            infecteds[k].type = 'W';
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                infecteds[k].time = rnd_unif( inf_length );
            else // Time-since-infection models
                if ( P.MODEL_TYPE == 1 ) // Truncated gen time distr
                    do {
                        infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                    } while ( infecteds[k].time > P.DUR_SICK );
                else // Time-since-infection model with full gen time distr
                    infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            k++;
        }
        else if ( ind[ index_ind_hit ].is == 'X' ) { // If the individual was marked already, store new time only if sooner
            y = 0;
            while ( infecteds[y].index_obj != index_ind_hit )
                y++;
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                hitting_time = rnd_unif( inf_length );
            else // Time-since-infection models (both, see above)
                hitting_time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            if ( hitting_time < infecteds[y].time ) {
                infecteds[y].time = hitting_time;
                infecteds[y].type = 'W';
            }
        }
        else
            assert( ind[ index_ind_hit ].is == 'E' || ind[ index_ind_hit ].is == 'I' || ind[ index_ind_hit ].is == 'R' );
    }
    /* Add the W infections to adults to the infected list */
    for ( x = 1; x <= W_num_conts_to_A; x++ ) {
        do { // infectivity is towards other adults, but the workplace has a list with all members. If the individual infect him/herself, draw another child - this might be inefficient for very small sizes
            W_index_ind_hit = ( int ) ( W_list[ W_ind ].sizeA * ran2(P.idum) );
            index_ind_hit = W_list[ W_ind ].peopleA[ W_index_ind_hit ];
            assert( ind[ index_ind_hit ].age == 'A' );
        } while ( index_ind_hit == index_sub );
        if ( ind[ index_ind_hit ].is  == 'S' ) {
            ind[ index_ind_hit ].is = 'X';
            infecteds[k].index_obj = index_ind_hit;
            infecteds[k].type = 'W';
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                infecteds[k].time = rnd_unif( inf_length );
            else // Time-since-infection models
                if ( P.MODEL_TYPE == 1 ) // Truncated gen time distr
                    do {
                        infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                    } while ( infecteds[k].time > P.DUR_SICK );
                else // Time-since-infection model with full gen time distr
                    infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            k++;
        }
        else if ( ind[ index_ind_hit ].is == 'X' ) { // If the individual was marked already, store new time only if sooner
            y = 0;
            while ( infecteds[y].index_obj != index_ind_hit )
                y++;
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                hitting_time = rnd_unif( inf_length );
            else // Time-since-infection models (both, see above)
                hitting_time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            if ( hitting_time < infecteds[y].time ) {
                infecteds[y].time = hitting_time;
                infecteds[y].type = 'W';
            }
        }
        else
            assert( ind[ index_ind_hit ].is == 'E' || ind[ index_ind_hit ].is == 'I' || ind[ index_ind_hit ].is == 'R' );
    }
    
    /* Add G infections to children to the infected list */
    for ( x = 1; x <= G_num_conts_to_C; x++ ) {
        do { // infectivity is towards other children. If the individual infect him/herself, draw another child.
            index_ind_hit = ( int ) ( P.TOT_C * ran2(P.idum) );
            assert( ind[ index_ind_hit ].age == 'C' );
        } while ( index_ind_hit == index_sub );
        if ( ind[ index_ind_hit ].is  == 'S' ) {
            ind[ index_ind_hit ].is = 'X';
            infecteds[k].index_obj = index_ind_hit;
            infecteds[k].type = 'G';
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                infecteds[k].time = rnd_unif( inf_length );
            else // Time-since-infection models
                if ( P.MODEL_TYPE == 1 ) // Truncated gen time distr
                    do {
                        infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                    } while ( infecteds[k].time > P.DUR_SICK );
                else // Time-since-infection model with full gen time distr
                    infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            k++;
        }
        else if ( ind[ index_ind_hit ].is == 'X' ) { // If the individual was marked already, store new time only if sooner
            y = 0;
            while ( infecteds[y].index_obj != index_ind_hit )
                y++;
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                hitting_time = rnd_unif( inf_length );
            else // Time-since-infection models (both, see above)
                hitting_time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            if ( hitting_time < infecteds[y].time ) {
                infecteds[y].time = hitting_time;
                infecteds[y].type = 'G';
            }
        }
        else
            assert( ind[ index_ind_hit ].is == 'E' || ind[ index_ind_hit ].is == 'I' || ind[ index_ind_hit ].is == 'R' );
    }
    /* Add G infections to adults to the infected list */
    for ( x = 1; x <= G_num_conts_to_A; x++ ) {
        do { // infectivity is towards other adult. If the individual infect him/herself, draw another child.
            index_ind_hit = P.TOT_C + ( int ) ( P.TOT_A * ran2(P.idum) );
            assert( ind[ index_ind_hit ].age == 'A' );
        } while ( index_ind_hit == index_sub );
        if ( ind[ index_ind_hit ].is  == 'S' ) {
            ind[ index_ind_hit ].is = 'X';
            infecteds[k].index_obj = index_ind_hit;
            infecteds[k].type = 'G';
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                infecteds[k].time = rnd_unif( inf_length );
            else // Time-since-infection models
                if ( P.MODEL_TYPE == 1 ) // Truncated gen time distr
                    do {
                        infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
                    } while ( infecteds[k].time > P.DUR_SICK );
                else // Time-since-infection model with full gen time distr
                    infecteds[k].time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            k++;
        }
        else if ( ind[ index_ind_hit ].is == 'X' ) { // If the individual was marked already, store new time only if sooner
            y = 0;
            while ( infecteds[y].index_obj != index_ind_hit )
                y++;
            if ( P.MODEL_TYPE == 0 ) // SIR and SEIR models
                hitting_time = rnd_unif( inf_length );
            else // Time-since-infection models (both, see above)
                hitting_time = rnd_time( P.GEN_TIME_DISTR, P.GEN_TIME, P.GEN_ALPHA );
            if ( hitting_time < infecteds[y].time ) {
                infecteds[y].time = hitting_time;
                infecteds[y].type = 'G';
            }
        }
        else
            assert( ind[ index_ind_hit ].is == 'E' || ind[ index_ind_hit ].is == 'I' || ind[ index_ind_hit ].is == 'R' );
    }
    
    /* Now clean all the X that we added */
    for ( x = 0; x < k; x++ ) {
        ind[ infecteds[x].index_obj ].is = 'S';
        create_and_place_event( basic_time+lat_length+infecteds[x].time, 'I', index_sub, infecteds[x].index_obj, infecteds[x].type );
        if ( P.MODEL_TYPE != 0 ) { // Time-since-infection model: keep track of last infection time
            if ( infecteds[x].time > last_inf )
                last_inf = infecteds[x].time;
        }
    }
    if ( P.MODEL_TYPE != 0 ) // Time-since-infection model: place recovery event
    {
        if ( P.MODEL_TYPE == 1 )
            create_and_place_event( basic_time+P.DUR_SICK, 'R', index_sub, index_sub, '0' );
        else // Non-truncated time-since-infection model: here recovery is placed a fixed time REC_DELAY after last infection
            create_and_place_event( basic_time+last_inf+P.REC_DELAY, 'R', index_sub, index_sub, '0' );
    }
    free( infecteds );
}

// Only for SIR or SEIR models, i.e. MODEL_TYPE == 0
double rnd_time( int distr, double mean, int alpha ) /* Return the duration of infectious period */
{
    switch (distr) {
        case 0: /* Dirac distribution */
            return mean;
            break;
        case 1:	/* Exp distr */
            return mean*expdev(P.idum);
            break;
        case 2:	/* Gamma distr */
            return mean*gamdev(alpha,P.idum)/alpha;
            break;
        default:
            printf("Error in rnd_time!");
            return mean;
            break;
    }
}

int rnd_Poiss( double mean )
{
    return ( int ) poidev( mean, P.idum );
}

double rnd_unif( double interval )
{
    return interval*ran2( P.idum );
}

void create_and_place_event( double t, char type, int i1, int i2, char inf_type )
{
    EVENT* ptr;
    
    // Prepare memory
    if ( P.ev_count >= P.MAX_EV ) { // If we have used the whole vector of events that was preallocated in the heap
        // We need to add more using malloc
        ptr =  ( EVENT* ) malloc( sizeof( EVENT ) );
        ptr->heap = 'Y'; // This signals that the event must be freed later on
    } else
        ptr = &ev[ P.ev_count ];
    
    // Fill the event
    ptr->t = t;
    ptr->type = type;
    ptr->i_sub = i1; // Same subject and object if it is a recovery
    ptr->i_obj = i2;
    ptr->inf_type = inf_type;
    ptr->next = NULL;
    
    // Place the event
    P.ev_count++; // keep track, just to have an idea about how many events we need
    place_node( ptr );
}

/* The core of the programme. It follows the event list and adds new infections */
void infection_process( void )
{
    int n, s;
    int i_s, i_o, H_s, W_s, H_o, W_o;
    double tc_s, to_s;
    EVENT* bookmark;
    double t_ev = 0;
    
    for ( n = 0; n <= P.TOT_h_STEPS; n++ ) {
        /* Just a form of control: I don't want TIME_STOP to occur before the epidemic is finished */
        if ( ( n == P.TOT_h_STEPS ) && ( h[ n ].p != NULL ) )
            printf( "\nCareful!!! Inefficiency may arise due to early TIME_STOP!!!\n" ); // The sublist of events for the last element of the helper may be long
        
        if ( h[ n ].p != NULL ) {
            bookmark = h[ n ].p;
            while ( bookmark != NULL ) {
                t_ev = bookmark->t;
                i_s = bookmark->i_sub;
                i_o = bookmark->i_obj;
                H_s = ind[ i_s ].H_index;
                W_s = ind[ i_s ].W_index;
                H_o = ind[ i_o ].H_index;
                W_o = ind[ i_o ].W_index;
//                tc_s = H_list[ H_s ].tc; // Time of closure of the school of the subject
//                to_s = H_list[ H_s ].to; // Time of reopening of the school of the subject
                //printf( "%f\n", t_ev ); // For debugging purposes
                if ( bookmark->type == 'I' )
                { /* If this event is an infectious contact... */
                    /* Subject AND WHEN TO IGNORE THE EVENT */
//                    if ( SWITCH_SCHOOL_CLOSURE )
//                    {
//                        if ( bookmark->inf_type == 'H' )
//                        {
//                            assert( ind[ i_s ].age == 'C' );
//                            if ( ( t_ev > tc_s ) && ( t_ev < to_s ) ) // if schools are closed, this infection does not exists
//                            {
//                                bookmark->inf_type = 'F';
//                                goto jump;
//                            }
//                            else // if it exists
//                            {
//                                if ( ind[ i_o ].is == 'S' )
//                                    ind[ i_s ].H_infs++;
//                                else
//                                {
//                                    assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                    bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                }
//                            }
//                        }
//                        else if ( bookmark->inf_type == 'W' )
//                        {
//                            if ( ind[ i_s ].age == 'C' )
//                            {
//                                if ( !( ( t_ev > tc_s ) && ( t_ev < to_s ) ) )
//                                { // schools are open. I need to reduce the extra transmission due to compensatory behaviour to the correct one
//                                    if ( ind[ i_o ].age == 'A' )
//                                    {
//                                        if ( ran2( P.idum ) > ( 1 / P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_A ) ) // This infection does not exists
//                                        {
//                                            bookmark->inf_type = 'N';
//                                            goto jump;
//                                        }
//                                        else // If it exists
//                                        {
//                                            if ( ind[ i_o ].is == 'S' )
//                                                ind[ i_s ].W_infs++;
//                                            else
//                                            {
//                                                assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                                bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                            }
//                                        }
//                                    }
//                                    else
//                                    {
//                                        assert( ind[ i_o ].age == 'C' );
//                                        if ( ran2( P.idum ) > ( 1 / P.FACTOR_COMP_BEHAVIOUR_IN_HOUSEHOLDS_TO_C ) ) // If this infection does not exists
//                                        {
//                                            bookmark->inf_type = 'N';
//                                            goto jump;
//                                        }
//                                        else // If it exists
//                                        {
//                                            if ( ind[ i_o ].is == 'S' )
//                                                ind[ i_s ].W_infs++;
//                                            else
//                                            {
//                                                assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                                bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                            }
//                                        }
//                                    }
//                                }
//                                else // If schools are closed, all infectious attempts occur
//                                {
//                                    if ( ind[ i_o ].is == 'S' )
//                                        ind[ i_s ].W_infs++;
//                                    else
//                                    {
//                                        assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                        bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                    }
//                                }
//                                
//                            }
//                            else // The subject is an adult, the infectious attempts exists always
//                            {
//                                assert( ind[ i_s ].age == 'A' );
//                                if ( ind[ i_o ].is == 'S' )
//                                    ind[ i_s ].W_infs++;
//                                else
//                                {
//                                    assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                    bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                }
//                            }
//                        }
//                        else // It's a global infection
//                        {
//                            assert( bookmark->inf_type == 'G' );
//                            if ( ind[ i_s ].age == 'C' )
//                            {
//                                if ( !( ( t_ev > tc_s ) && ( t_ev < to_s ) ) )
//                                { // schools are open. I need to reduce the extra transmission due to compensatory behaviour to the correct one
//                                    if ( ind[ i_o ].age == 'A' )
//                                    {
//                                        if ( ran2( P.idum ) > ( 1 / P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_A ) ) // This infection does not exists
//                                        {
//                                            bookmark->inf_type = 'M';
//                                            goto jump;
//                                        }
//                                        else // If it does, instead
//                                        {
//                                            ind[ i_s ].G_conts++;
//                                            if ( ind[ i_o ].is == 'S' )
//                                                ind[ i_s ].G_infs++;
//                                            else
//                                            {
//                                                assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                                bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                            }
//                                        }
//                                    }
//                                    else
//                                    {
//                                        assert( ind[ i_o ].age == 'C' );
//                                        if ( ran2( P.idum ) > ( 1 / P.FACTOR_COMP_BEHAVIOUR_IN_COMM_TO_C ) ) // This infection does not exists
//                                        {
//                                            bookmark->inf_type = 'M';
//                                            goto jump;
//                                        }
//                                        else // if it does exists
//                                        {
//                                            ind[ i_s ].G_conts++;
//                                            if ( ind[ i_o ].is == 'S' )
//                                                ind[ i_s ].G_infs++;
//                                            else
//                                            {
//                                                assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                                bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                            }
//                                        }
//                                    }
//                                }
//                                else // If schools are closed, all infectious attempts occur
//                                {
//                                    ind[ i_s ].G_conts++;
//                                    if ( ind[ i_o ].is == 'S' )
//                                        ind[ i_s ].G_infs++;
//                                    else
//                                    {
//                                        assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                        bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                assert( ind[ i_s ].age == 'A' );
//                                ind[ i_s ].G_conts++;
//                                if ( ind[ i_o ].is == 'S' )
//                                    ind[ i_s ].G_infs++;
//                                else
//                                {
//                                    assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
//                                    bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
//                                }
//                            }
//                        }
//                    }
//                    else // There is no school closure, i.e. there's no infectivity to correct or events to ignore
//                    {	 // Therefore, compensatory behaviour is automatically set to 1 (to avoid generating extra events)
                        if ( bookmark->inf_type == 'H' )
                            ind[ i_s ].H_infs++;
                        else if ( bookmark->inf_type == 'W' )
                            ind[ i_s ].W_infs++;
                        else // It's a global infection
                        {
                            ind[ i_s ].G_conts++;
                            if ( ind[ i_o ].is == 'S' )
                                ind[ i_s ].G_infs++;
                            else
                            {
                                assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
                                bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
                            }
                        }
//                    }
                    
                    // All the rest is happening only if the infection actually occurs
                    /* Object */
                    if ( ind[ i_o ].is == 'S' )
                    {
                        if ( P.MODEL_TYPE == 0 && P.LATENCE == 1 )
                            ind[ i_o ].is = 'E';
                        else
                            ind[ i_o ].is = 'I';
                        ind[ i_o ].ti = t_ev;
                        if ( ind[ i_s ].gen < P.MAX_GEN )
                            ind[ i_o ].gen = ind[ i_s ].gen + 1;
                        else
                        {
                            assert( ind[ i_s ].gen == P.MAX_GEN );
                            ind[ i_o ].gen = P.MAX_GEN;
                        }
                        infectious_life( i_o , t_ev );
                        if ( H_o != -1 ) // If this is a child, i.e. it has a household (a true school)
                        {
                            assert( ind[ i_o ].age == 'C' );
                            /* Household (of the object) */
                            effect_of_infection_on_H( bookmark );
                        }
                        /* Workplace (of the object) */
                        effect_of_infection_on_W( bookmark );
                    }
                    else
                    {
                        assert( ind[ i_o ].is == 'E' || ind[ i_o ].is == 'I' || ind[ i_o ].is == 'R' );
                        bookmark->inf_type = 'X'; // To signal it's an attempt of hitting an already infected person
                    }
                }
                else if ( bookmark->type == 'E' )
                { /* If the latent period finishes (if there is any) */
                    assert ( P.MODEL_TYPE == 0 && P.LATENCE == 1 );
                    assert ( ind[ i_s ].is == 'E' );
                    ind[ i_s ].is = 'I';
                    ind[ i_s ].tl = t_ev;
                }
                else if ( bookmark->type == 'R' )
                { /* If the event is a recovery... */
                    assert( ind[ i_s ].is == 'I' );
                    ind[ i_s ].is = 'R';
                    ind[ i_s ].tr = t_ev;
                    if ( H_s != -1 ) // If this is a child, i.e. it has a household (a true school)
                    {
                        assert( ind[ i_s ].age == 'C' );
                        effect_of_recovery_on_H( bookmark );
//                        if ( (  ind[ i_s ].ta > 0 ) && ( SWITCH_SCHOOL_CLOSURE ) )// only if symptomatic
//                        {
//                            assert( H_list[ H_s ].abs > 0 );
//                            H_list[ H_s ].abs--;
//                        }
                    }
                    effect_of_recovery_on_W( bookmark );
                }
                else
                {
                    assert( bookmark->type == 'A' );
//                    //////////////////////////////////
//                    if ( SWITCH_SCHOOL_CLOSURE )
//                    {
//                        /* Check if schools are closed or not */
//                        if ( H_o != -1 ) //NB recall that i_obj = i_sub in this case
//                        {
//                            // If SWITCH_PREVALENCE_TRIGGER = 1, I flag schools when prevalence crosses a trigger, otherwise I use cum inc: note that e1?e2:e3 has value e2 if e1 is true, otherwise e3
//                            //if ( H_list[ ind[ bookmark->i_obj ].H_index ].ufs >= P.WITHIN_SCHOOL_PREV_TRIGGER - 1 )
//                            H_list[ H_o ].abs++;
//                            if ( ( ( SWITCH_PREVALENCE_TRIGGER == 1 ) ? H_list[ H_o ].abs : H_list[ H_o ].ufs ) >= P.WITHIN_SCHOOL_PREV_TRIGGER )
//                            {	// If the school needs to be flagged
//                                if ( H_list[ H_o ].tc == -1 ) { // It has not been set yet
//                                    P.n_schools_flagged++;
//                                    if ( SWITCH_CLOSE_EACH_SCHOOL ) {
//                                        H_list[ H_o ].tc = t_ev + P.DELAY;
//                                        if ( SWITCH_SCHOOL_REOPEN )
//                                            H_list[ H_o ].to = H_list[ H_o ].tc + P.DURATION_CLOSURE;
//                                    }
//                                    else // If I wait to close all schools in one go
//                                        H_list[ H_o ].tc = -2; // To signal the school is flagged but not closed
//                                    if ( ( SWITCH_CLOSE_ALL_SCHOOLS ) && ( ( P.n_schools_flagged * 100 / P.n_H ) >= P.PERC_SCHOOLS_INFECTED_TRIGGER ) ) // Careful, integer division, and percentage is also an integer
//                                    {
//                                        P.TIME_CLOSURE = H_list[ H_o ].tc = t_ev + P.DELAY;
//                                        P.TIME_REOPEN = H_list[ H_o ].tc + P.DURATION_CLOSURE;
//                                        P.t_school_closure[ P.epid_index ] = P.TIME_CLOSURE;
//                                        P.t_school_reopen[ P.epid_index ] = P.TIME_REOPEN;
//                                        // close all schools
//                                        for ( s = 0; s < P.n_H; s++ ) {
//                                            H_list[ s ].tc = P.TIME_CLOSURE;
//                                            H_list[ s ].to = P.TIME_REOPEN;
//                                        }
//                                    }
//                                }
//                                //if ( bookmark->inf_type == 'H' )
//                                //{
//                                //	if ( ( H_list[ ind[ bookmark->i_sub ].H_index ].tc > 0 ) && ( t_ev > H_list[ ind[ bookmark->i_sub ].H_index ].tc ) && ( t_ev < H_list[ ind[ bookmark->i_sub ].H_index ].to ) )
//                                //	{	// not only ts != -1 because it could also be -2 (flagged but not closed)
//                                //		bookmark->inf_type = 'F';	// Just ignore this event
//                                //		goto jump; // Skip everything (as if this event wouldn't exists
//                                //	}
//                                //}
//                            }
//                        }
//                    }
                }
            jump:
                bookmark = bookmark->next;
            }
        }
    }
    printf("");
}

void effect_of_infection_on_H( EVENT* ptr )
{
    int H_obj, H_sub, W_sub;
    
    H_obj = ind[ ptr->i_obj ].H_index;
    H_sub = ind[ ptr->i_sub ].H_index;
    W_sub = ind[ ptr->i_sub ].W_index;
    
    /* First infection of the household */
    if ( H_list[ H_obj ].is == 'S' ) {
        assert( H_list[ H_obj ].preval == 0 );
        assert( H_list[ H_obj ].ufs == 0 );
        assert( H_list[ H_obj ].ti == -1 );
        assert( H_list[ H_obj ].tm == -1 );
        assert( H_list[ H_obj ].tr == -1 );
        H_list[ H_obj ].is = 'I';
        H_list[ H_obj ].ti = ptr->t;
        H_list[ H_obj ].preval++;
        H_list[ H_obj ].ufs++;
        //// [[[
        if ( ptr->inf_type == 'H' )			//// This case can never happen in theory
            assert( 0 );					////
        else if ( ptr->inf_type == 'W' ) {	//// A new household has been infected locally
            //////		assert( W_list[ W_sub ].is == 'I' );
            assert( W_list[ W_sub ].preval > 0 );
            assert( W_list[ W_sub ].ti >= 0 );
        } else {							//// A new household has been infected globally
            assert( ptr->inf_type == 'G' );
        }
        //// ]]]
    } else if ( H_list[ H_obj ].is == 'I' ) {
        /* The household has already been infected */
        assert( H_list[ H_obj ].ti >= 0 );
        assert( H_list[ H_obj ].tr < 0 ); // hasn't recovered yet
        assert( H_list[ H_obj ].preval > 0 );
        assert( H_list[ H_obj ].ufs > 0 );
        if ( ptr->inf_type == 'H' ) {
            // Everything is regular, because we are in a household epidemic
            assert( H_obj == ind[ ptr->i_sub ].H_index ); // check
            H_list[ H_obj ].preval++;
            H_list[ H_obj ].ufs++; // This is the only parameter that never gets frozen. I use it as a counter
        } else { // Reintroduction from outside (W or G)
            assert( ptr->inf_type == 'W' || ptr->inf_type == 'G' );
            if ( H_list[ H_obj ].tm == -1 ) { // If it is the first reintroduction
                H_list[ H_obj ].tm = ptr->t;
                H_list[ H_obj ].preval++;
                H_list[ H_obj ].ufs++;
            } else { // It is not the first reintroduction
                assert( H_list[ H_obj ].tm >= 0 ); // be careful: two initial cases in the same household may be a problem
                H_list[ H_obj ].preval++;
                H_list[ H_obj ].ufs++;
            }
        }
    } else { // The household has recovered, but gets reinfected
        assert( H_list[ H_obj ].ufs > 0 );
        assert( H_list[ H_obj ].is == 'R' );
        assert( H_list[ H_obj ].tr > 0 );
        if ( H_list[ H_obj ].tm == -1 ) { // If it is the first reintroduction
            assert( H_list[ H_obj ].preval == 0 );
            H_list[ H_obj ].tm = ptr->t;
            H_list[ H_obj ].preval++;
            H_list[ H_obj ].ufs++;
        } else { // It is not the first reintroduction
            assert( H_list[ H_obj ].tm >= 0 ); // be careful: two initial cases in the same household may be a problem
            H_list[ H_obj ].preval++;
            H_list[ H_obj ].ufs++;
        }
    }
}

void effect_of_infection_on_W( EVENT* ptr )
{
    int W_obj, H_sub, W_sub;
    
    W_obj = ind[ ptr->i_obj ].W_index;
    H_sub = ind[ ptr->i_sub ].H_index;
    W_sub = ind[ ptr->i_sub ].W_index;
    
    /* First infection of the workplace */
    if ( W_list[ W_obj ].is == 'S' ) {
        assert( W_list[ W_obj ].preval == 0 );
        assert( W_list[ W_obj ].ufs == 0 );
        assert( W_list[ W_obj ].ti == -1 );
        assert( W_list[ W_obj ].tm == -1 );
        assert( W_list[ W_obj ].tr == -1 );
        W_list[ W_obj ].is = 'I';
        W_list[ W_obj ].ti = ptr->t;
        W_list[ W_obj ].preval++;
        W_list[ W_obj ].ufs++;
        //// [[[
        if ( ptr->inf_type == 'H' )	{		//// A new workplace has been infected locally
            assert( H_list[ H_sub ].preval > 0 );
            assert( H_list[ H_sub ].ti >= 0 );
            assert( ind[ ptr->i_obj ].H_index == H_sub );
        } else if ( ptr->inf_type == 'W' ) {//// This case can never happen in theory
            assert( 0 );					////
        } else {							//// A new household has been infected globally
            assert( ptr->inf_type == 'G' );
        }
        //// ]]]
    } else if ( W_list[ W_obj ].is == 'I' ) {
        /* The workplace has already been infected */
        assert( W_list[ W_obj ].ti >= 0 );
        assert( W_list[ W_obj ].tr < 0 ); // hasn't recovered yet
        assert( W_list[ W_obj ].preval > 0 );
        assert( W_list[ W_obj ].ufs > 0 );
        if ( ptr->inf_type == 'W' ) {
            // Everything is regular, because we are in a workplace epidemic
            assert( W_obj == ind[ ptr->i_sub ].W_index ); // check
            W_list[ W_obj ].preval++;
            W_list[ W_obj ].ufs++; // This is the only parameter that never gets frozen. I use it as a counter
        } else { // Reintroduction from outside (H or G)
            assert( ptr->inf_type == 'H' || ptr->inf_type == 'G' );
            if ( W_list[ W_obj ].tm == -1 ) { // If it is the first reintroduction
                W_list[ W_obj ].tm = ptr->t;
                W_list[ W_obj ].preval++;
                W_list[ W_obj ].ufs++;
            } else { // It is not the first reintroduction
                assert( W_list[ W_obj ].tm >= 0 ); // be careful: two initial cases in the same workplace may be a problem
                W_list[ W_obj ].preval++;
                W_list[ W_obj ].ufs++;
            }
        }
    } else { // The workplace has recovered, but gets reinfected
        assert( W_list[ W_obj ].ufs > 0 );
        assert( W_list[ W_obj ].is == 'R' );
        assert( W_list[ W_obj ].tr > 0 );
        if ( W_list[ W_obj ].tm == -1 ) { // If it is the first reintroduction
            assert( W_list[ W_obj ].preval == 0 );
            W_list[ W_obj ].tm = ptr->t;
            W_list[ W_obj ].preval++;
            W_list[ W_obj ].ufs++;
        } else { // It is not the first reintroduction
            assert( W_list[ W_obj ].tm >= 0 ); // be careful: two initial cases in the same workplace may be a problem
            W_list[ W_obj ].preval++;
            W_list[ W_obj ].ufs++;
        }
    }
}

void effect_of_recovery_on_H( EVENT* ptr )
{
    int H_sub;
    
    H_sub = ind[ ptr->i_sub ].H_index;
    
    assert( H_list[ H_sub ].is != 'S' );
    assert( H_list[ H_sub ].ti >= 0 );
    assert( H_list[ H_sub ].preval > 0 );
    assert( H_list[ H_sub ].ufs > 0 );
    H_list[ H_sub ].preval--;
    if ( H_list[ H_sub ].preval == 0 ) { // there are no more infectives in the household
        if ( H_list[ H_sub ].is == 'I' ) { // The household fully recovers for the first time
            assert( H_list[ H_sub ].tr < 0 );
            H_list[ H_sub ].is = 'R';
            H_list[ H_sub ].tr = ptr->t;
        }
    }
}

void effect_of_recovery_on_W( EVENT* ptr )
{
    int W_sub;
    
    W_sub = ind[ ptr->i_sub ].W_index;
    
    assert( W_list[ W_sub ].is != 'S' );
    assert( W_list[ W_sub ].ti >= 0 );
    assert( W_list[ W_sub ].preval > 0 );
    assert( W_list[ W_sub ].ufs > 0 );
    W_list[ W_sub ].preval--;
    if ( W_list[ W_sub ].preval == 0 ) { // there are no more infectives in the workplace
        if ( W_list[ W_sub ].is == 'I' ) { // The workplace fully recovers for the first time
            assert( W_list[ W_sub ].tr < 0 );
            W_list[ W_sub ].is = 'R';
            W_list[ W_sub ].tr = ptr->t;
        }
    }
}



void print_epidemic( void )
{
    FILE* ep;
    int l;
    EVENT* bmk;
    
    ep = fopen( "epidemic.dat", "w" );
    
    fprintf( ep, "Full epidemic:\n\n" );
    for ( l = 0; l <= P.TOT_h_STEPS; l++ ) {
        bmk = h[ l ].p;
        while ( bmk != NULL ) {
            if ( bmk->type == 'I' ) 
                fprintf( ep, "%.3f\t%c\t%d\t%d\t%c\n", bmk->t, bmk->type, bmk->i_sub, bmk->i_obj, bmk->inf_type );
            else
                fprintf( ep, "%.3f\t%c\t%d\t%d\t%c\n", bmk->t, bmk->type, bmk->i_sub, bmk->i_obj, bmk->inf_type );
            bmk = bmk->next;
        }
    }
    fclose( ep );
}

void clean_and_free( void )
{
    int l;
    EVENT *bmk, *temp;
    
    for ( l = 0; l <= P.TOT_h_STEPS; l++ ) {
        bmk = h[ l ].p;
        while ( bmk != NULL ) {
            temp = bmk->next;
            if ( bmk->heap == 'Y' ) { // Out of the vector: need to free
                free( bmk );
            } else {
                assert( bmk->heap == 'N' ); // In the vector: just need to clean 
                bmk->t = -1;
                bmk->type = 0;
                bmk->i_sub = -1;
                bmk->i_obj = -1;
                bmk->next = NULL;
            }
            bmk = temp;
        }
        h[ l ].p = NULL;
    }
}


