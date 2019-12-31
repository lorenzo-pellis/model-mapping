//
//  analyse_epidemic.cpp
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

#include "headers.h"

// Here I analyse the epidemic by converting a list of events into discrete-time epidemic curves.
// After that, I still have to reprocess them to synchronise them at the peak
void analyse_epidemic( void )
{
    
    store_epidemic();	// I can't create cumulative SIR straight away, because I want to synchronize
    get_final_size();	// store_epidemic calculates also the final size: use only one or the other (otherwise n_large is twice what it should be)

    // These functions are less relevant, but are used to study how R0 and r change over time
    create_cum_gen_R0();
    create_cum_r();

    // Finally, I clean all the extra memory allocated and I re-initialise values used (this function is in "epidemic.cpp")
    clean_and_free();
}

void get_final_size( void )
{
    int i, u, w;
    int fs = 0, fs_H = 0, fs_W = 0;
    
    for ( i = 0; i < P.TOT; i++ ) {
        if ( ind[ i ].ti >= 0 ) {	// Only those that have a true infectious period
            assert( ind[ i ].is == 'R' );	// It shouldn't be anything else
            assert( ind[ i ].gen > -1 );	// Check that he took part in the epidemic
            fs++; // Keep track of the final size, so to check if it is a large epidemic or not
        }
    }
    P.final_size[ P.epid_index ] = fs;
    /* Is it a large epidemic? */
    if ( fs > P.CUTOFF*P.TOT ) {
        P.n_large++;
        P.extinction_flag[ P.epid_index ] = 'L';
    } else
        P.extinction_flag[ P.epid_index ] = 'S';
    
    for ( u = 0; u < P.n_H; u++ ) {
        if ( H_list[ u ].ti >= 0 ) {	// Only those that have a true infectious period
            assert( H_list[ u ].is == 'R' );	// It shouldn't be anything else
            assert( H_list[ u ].tr > 0 );	// Check that he took part in the epidemic
            fs_H++; // Keep track of the final size, so to check if it is large or not
        }
    }
    P.final_size_H[ P.epid_index ] = fs_H;
    
    for ( w = 0; w < P.n_W; w++ ) {
        if ( W_list[ w ].ti >= 0 ) {	// Only those that have a true infectious period
            assert( W_list[ w ].is == 'R' );	// It shouldn't be anything else
            assert( W_list[ w ].tr > 0 );	// Check that he took part in the epidemic
            fs_W++; // Keep track of the final size, so to check if it is large or not
        }
    }
    P.final_size_W[ P.epid_index ] = fs_W;
    
    for ( w = 0; w < P.n_W; w++ )
    {
        P.table_freq_h_sizes_ptr[ P.epid_index ][ ( W_list[ w ].sizeA + W_list[ w ].sizeC ) - 1 ]++;
        P.table_h_fs_ptr[ P.epid_index ][ ( W_list[ w ].sizeA + W_list[ w ].sizeC ) - 1 ][ W_list[ w ].ufs ]++;
    }
    
}

void store_epidemic( void )
{
    int i, k, t_index;
    int I, cum_inc, R, n_conts, n_infs, G_infs, H_infs, W_infs;
    int cum_inc_C, cum_inc_A;
    int ABSsick, ABSclosed;
    int peak, t_peak;
    double last_t = 0;
    int u, w;
    int H_I, H_cum_inc, H_R, H_C;
    int W_cum_inc;
    
    n_conts = 0;
    n_infs = 0;
    G_infs = 0;
    H_infs = 0;
    W_infs = 0;
    H_I = 0;
    H_C = 0;
    H_R = 0;
    W_cum_inc = 0; // I try to do this in order not to add too many variables
    for ( i = 0; i < P.TOT; i++ ) {
        // First deal with the initial infectives, then all the initial susceptibles (ignore immunes)
        if ( ind[ i ].ti == 0 ) { // Only the initial infectives
            assert( ind[ i ].is == 'R' );
            n_conts += ind[ i ].G_conts + ind[ i ].H_infs + ind[ i ].W_infs;
            n_infs += ind[ i ].G_infs + ind[ i ].H_infs + ind[ i ].W_infs;
            G_infs += ind[ i ].G_infs;
            H_infs += ind[ i ].H_infs;
            W_infs += ind[ i ].W_infs;
            sys[ P.epid_index ][ 0 ].colls += ind[ i ].G_conts - ind[ i ].G_infs;
            
            if ( ind[ i ].tr > last_t )
                last_t = ind[ i ].tr;
            t_index = ( int ) ( ind[ i ].tr / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].I--;
            ch[ t_index ].R++;
        } else if ( ind[ i ].ti > 0 ) { // In this way we exclude the initial infectives and immunes
            assert( ind[ i ].is != 'S' );
            
            t_index = ( int ) ( ind[ i ].ti / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].S--;
            ch[ t_index ].I++;
            //////////////////////////
            if ( ind[ i ].age == 'C' )
                ch[ t_index ].C++;
            else {
                assert( ind[ i ].age == 'A' );
                ch[ t_index ].A++;
            }
            ///////////////////////////
            ch[ t_index ].n_conts += ind[ i ].G_conts + ind[ i ].H_infs + ind[ i ].W_infs;
            ch[ t_index ].n_infs += ind[ i ].G_infs + ind[ i ].H_infs + ind[ i ].W_infs;
            ch[ t_index ].G_infs += ind[ i ].G_infs;
            ch[ t_index ].H_infs += ind[ i ].H_infs;
            ch[ t_index ].W_infs += ind[ i ].W_infs;
            ch[ t_index ].colls += ind[ i ].G_conts - ind[ i ].G_infs;
            
            if ( ind[ i ].tr > last_t )
                last_t = ind[ i ].tr;
            t_index = ( int ) ( ind[ i ].tr / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].I--;
            ch[ t_index ].R++;
        } else {
            assert( ind[ i ].ti < 0 && ind[ i ].tr < 0 );
            if ( P.N_INIT_IMM == 0 ) {
                assert( ind[ i ].is == 'S' );
                assert( ind[ i ].ti == -1 && ind[ i ].tr == -1 );
            }
        }
//        /* Sickness and school closure */
//        if ( ind[ i ].age == 'C' )
//        {
//            if ( ( ind[ i ].ta < 0 ) && ( H_list[ ind[ i ].H_index ].tc < 0 ) );
//            else if ( ( ind[ i ].ta < 0 ) && ( H_list[ ind[ i ].H_index ].tc > 0 ) )
//            {
//                t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                if ( t_index > P.TOT_STEPS-1 )
//                    t_index = P.TOT_STEPS-1;
//                ch[ t_index ].ABSclosed++;
//                if ( H_list[ ind[ i ].H_index ].to > 0 )
//                {
//                    t_index = ( int ) ( H_list[ ind[ i ].H_index ].to / P.dt );
//                    if ( t_index > P.TOT_STEPS-1 )
//                        t_index = P.TOT_STEPS-1;
//                    ch[ t_index ].ABSclosed--;
//                }
//            }
//            else if ( ( ind[ i ].ta > 0 ) && ( H_list[ ind[ i ].H_index ].tc < 0 ) )
//            {
//                t_index = ( int ) ( ind[ i ].ta / P.dt );
//                if ( t_index > P.TOT_STEPS-1 )
//                    t_index = P.TOT_STEPS-1;
//                ch[ t_index ].ABSsick++;
//                if ( ind[ i ].tb > 0 )
//                {
//                    t_index = ( int ) ( ind[ i ].tb / P.dt );
//                    if ( t_index > P.TOT_STEPS-1 )
//                        t_index = P.TOT_STEPS-1;
//                    ch[ t_index ].ABSsick--;
//                }
//            }
//            else // both events occurred: s/he gets sick and his/her school is closed
//            {
//                if ( ( ind[ i ].tb < 0 ) && ( H_list[ ind[ i ].H_index ].to < 0 ) )
//                {
//                    if ( ind[ i ].ta < H_list[ ind[ i ].H_index ].tc )
//                    {
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                    }
//                    else
//                    {
//                        assert( ind[ i ].ta > H_list[ ind[ i ].H_index ].tc );
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed++;
//                        
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                        ch[ t_index ].ABSclosed--;
//                    }
//                }
//                else if ( ( ind[ i ].tb < 0 ) && ( H_list[ ind[ i ].H_index ].to > 0 ) )
//                {
//                    if ( ind[ i ].ta < H_list[ ind[ i ].H_index ].tc )
//                    {
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                    }
//                    else if ( ( ind[ i ].ta > H_list[ ind[ i ].H_index ].tc ) && ( ind[ i ].ta < H_list[ ind[ i ].H_index ].to ) )
//                    {
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed++;
//                        
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                        ch[ t_index ].ABSclosed--;
//                    }
//                    else
//                    {
//                        assert( ind[ i ].ta > H_list[ ind[ i ].H_index ].to );
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed++;
//                        
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].to / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed--;
//                        
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                    }
//                }
//                else if ( ( ind[ i ].tb > 0 ) && ( H_list[ ind[ i ].H_index ].to < 0 ) )
//                {
//                    if ( H_list[ ind[ i ].H_index ].tc < ind[ i ].ta )
//                    {
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed++;
//                        
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                        ch[ t_index ].ABSclosed--;
//                        
//                        t_index = ( int ) ( ind[ i ].tb / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick--;
//                        ch[ t_index ].ABSclosed++;
//                    }
//                    else if ( ( H_list[ ind[ i ].H_index ].tc > ind[ i ].ta ) && ( H_list[ ind[ i ].H_index ].tc < ind[ i ].tb ) )
//                    {
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                        
//                        t_index = ( int ) ( ind[ i ].tb / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick--;
//                        ch[ t_index ].ABSclosed++;
//                    }
//                    else
//                    {
//                        assert( H_list[ ind[ i ].H_index ].tc > ind[ i ].tb );
//                        
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                        
//                        t_index = ( int ) ( ind[ i ].tb / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick--;
//                        
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed++;
//                    }
//                }
//                else
//                {
//                    assert( ( ind[ i ].tb > 0 ) && ( H_list[ ind[ i ].H_index ].to > 0 ) );
//                    if ( ind[ i ].ta < H_list[ ind[ i ].H_index ].tc )
//                    {
//                        t_index = ( int ) ( ind[ i ].ta / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSsick++;
//                        if ( ind[ i ].tb < H_list[ ind[ i ].H_index ].tc )
//                        {
//                            t_index = ( int ) ( ind[ i ].tb / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSsick--;
//                            
//                            t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSclosed++;
//                            
//                            t_index = ( int ) ( H_list[ ind[ i ].H_index ].to / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSclosed--;
//                        }
//                        else if ( ( ind[ i ].tb > H_list[ ind[ i ].H_index ].tc ) && ( ind[ i ].tb < H_list[ ind[ i ].H_index ].to ) )
//                        {
//                            t_index = ( int ) ( ind[ i ].tb / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSsick--;
//                            ch[ t_index ].ABSclosed++;
//                            
//                            t_index = ( int ) ( H_list[ ind[ i ].H_index ].to / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSclosed--;
//                        }
//                        else
//                        {
//                            assert( ind[ i ].tb > H_list[ ind[ i ].H_index ].to );
//                            t_index = ( int ) ( ind[ i ].tb / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSsick--;
//                        }
//                    }
//                    else
//                    {
//                        assert( ind[ i ].ta >= H_list[ ind[ i ].H_index ].tc );
//                        t_index = ( int ) ( H_list[ ind[ i ].H_index ].tc / P.dt );
//                        if ( t_index > P.TOT_STEPS-1 )
//                            t_index = P.TOT_STEPS-1;
//                        ch[ t_index ].ABSclosed++;
//                        if ( ind[ i ].ta < H_list[ ind[ i ].H_index ].to )
//                        {
//                            t_index = ( int ) ( ind[ i ].ta / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSsick++;
//                            ch[ t_index ].ABSclosed--;
//                            if ( ind[ i ].tb < H_list[ ind[ i ].H_index ].to )
//                            {
//                                t_index = ( int ) ( ind[ i ].tb / P.dt );
//                                if ( t_index > P.TOT_STEPS-1 )
//                                    t_index = P.TOT_STEPS-1;
//                                ch[ t_index ].ABSsick--;
//                                ch[ t_index ].ABSclosed++;
//                                
//                                t_index = ( int ) ( H_list[ ind[ i ].H_index ].to / P.dt );
//                                if ( t_index > P.TOT_STEPS-1 )
//                                    t_index = P.TOT_STEPS-1;
//                                ch[ t_index ].ABSclosed--;
//                            }
//                            else
//                            {
//                                assert( ind[ i ].tb > H_list[ ind[ i ].H_index ].to );
//                                t_index = ( int ) ( ind[ i ].tb / P.dt );
//                                if ( t_index > P.TOT_STEPS-1 )
//                                    t_index = P.TOT_STEPS-1;
//                                ch[ t_index ].ABSsick--;
//                            }
//                        }
//                        else
//                        {
//                            assert( ind[ i ].ta >= H_list[ ind[ i ].H_index ].to );
//                            t_index = ( int ) ( H_list[ ind[ i ].H_index ].to / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSclosed--;
//                            
//                            t_index = ( int ) ( ind[ i ].ta / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSsick++;
//                            
//                            t_index = ( int ) ( ind[ i ].tb / P.dt );
//                            if ( t_index > P.TOT_STEPS-1 )
//                                t_index = P.TOT_STEPS-1;
//                            ch[ t_index ].ABSsick--;
//                        }
//                    }
//                }
//            }
//        }
        
    }

    /* Households */
    for ( u = 0; u < P.n_H; u++ ) {
        if ( H_list[ u ].ti == 0 ) { // Initially infected households (assume for the moment that there is only 1 initial infective)
            assert( H_list[ u ].is == 'R' );
            //assert( H_list[ u ].tm > 0 ); // Just for the moment, no double initial cases
            assert( H_list[ u ].tr > 0 );
            
            H_I++;
            
            t_index = ( int ) ( H_list[ u ].tr / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].H_I--;
            ch[ t_index ].H_R++;
        }
        if ( H_list[ u ].ti > 0 ) { // Only non-initially infected households
            assert( H_list[ u ].is == 'R' );
            //assert( H_list[ u ].tm > 0 ); // Just for the moment, no double initial cases
            assert( H_list[ u ].tr > 0 );
            
            t_index = ( int ) ( H_list[ u ].ti / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].H_S--;
            ch[ t_index ].H_I++;
            
            t_index = ( int ) ( H_list[ u ].tr / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].H_I--;
            ch[ t_index ].H_R++;
        }
//        if ( SWITCH_SCHOOL_CLOSURE )
//        {
//            if ( H_list[ u ].tc > 0 )
//            {
//                t_index = ( int ) ( H_list[ u ].tc / P.dt );
//                if ( t_index > P.TOT_STEPS-1 )
//                    t_index = P.TOT_STEPS-1;
//                ch[ t_index ].H_C++;
//                //ch[ t_index ].ABS += H_list[ u ].size;
//            }
//            if ( SWITCH_SCHOOL_REOPEN )
//            {
//                if ( H_list[ u ].to > 0 )
//                {
//                    t_index = ( int ) ( H_list[ u ].to / P.dt );
//                    if ( t_index > P.TOT_STEPS-1 )
//                        t_index = P.TOT_STEPS-1;
//                    ch[ t_index ].H_C--;
//                    //ch[ t_index ].ABS -= H_list[ u ].size;
//                }
//            }
//        }
        
    }
    
    /* Workplaces */
    for ( w = 0; w < P.n_W; w++ ) {
        if ( W_list[ w ].ti == 0 ) { // Initially infected households (assume for the moment that there is only 1 initial infective)
            assert( W_list[ w ].is == 'R' );
            //assert( W_list[ w ].tm > 0 ); // Just for the moment, no double initial cases
            assert( W_list[ w ].tr > 0 );
            
            W_cum_inc++;
        }
        if ( W_list[ w ].ti > 0 ) { // Only non-initially infected households
            assert( W_list[ w ].is == 'R' );
            //assert( W_list[ w ].tm > 0 ); // Just for the moment, no double initial cases
            assert( W_list[ w ].tr > 0 );
            
            t_index = ( int ) ( W_list[ w ].ti / P.dt );
            if ( t_index > P.TOT_STEPS-1 )
                t_index = P.TOT_STEPS-1;
            ch[ t_index ].W_S--;
        }
    }
    
    /* store the epidemic */
    I = P.N_INIT_INF;
    cum_inc_C = P.N_INIT_INF_C;
    cum_inc_A = P.N_INIT_INF_A;
    R = P.N_INIT_IMM;
    cum_inc = I;
    peak = P.N_INIT_INF;
    t_peak = 0;
    ABSsick = 0;
    ABSclosed = 0;
    sys[ P.epid_index ][ 0 ].S = P.TOT - P.N_NON_SUSC;
    sys[ P.epid_index ][ 0 ].I = I;
    sys[ P.epid_index ][ 0 ].R = R;
    sys[ P.epid_index ][ 0 ].inc = 0;
    sys[ P.epid_index ][ 0 ].inc_C = 0;
    sys[ P.epid_index ][ 0 ].inc_A = 0;
    sys[ P.epid_index ][ 0 ].cum_inc = cum_inc;
    sys[ P.epid_index ][ 0 ].cum_inc_C = cum_inc_C;
    sys[ P.epid_index ][ 0 ].cum_inc_A = cum_inc_A;
//    sys[ P.epid_index ][ 0 ].ABSsick = ABSsick;
//    sys[ P.epid_index ][ 0 ].ABSclosed = ABSclosed;
    sys[ P.epid_index ][ 0 ].n_conts = n_conts;
    sys[ P.epid_index ][ 0 ].n_infs = n_infs;
    sys[ P.epid_index ][ 0 ].G_infs = G_infs;
    sys[ P.epid_index ][ 0 ].H_infs = H_infs;
    sys[ P.epid_index ][ 0 ].W_infs = W_infs;
    sys[ P.epid_index ][ 0 ].colls = n_conts - n_infs;
    
    H_cum_inc = H_I;
    sys[ P.epid_index ][ 0 ].H_S = P.n_H - H_I;
    sys[ P.epid_index ][ 0 ].H_I = H_I;
    sys[ P.epid_index ][ 0 ].H_R = H_R;
    sys[ P.epid_index ][ 0 ].H_inc = 0;
    sys[ P.epid_index ][ 0 ].H_cum_inc = H_cum_inc;
    sys[ P.epid_index ][ 0 ].H_C = H_C;
    
    sys[ P.epid_index ][ 0 ].W_cum_inc = W_cum_inc;
    for ( k = 1; k <= P.TOT_STEPS; k++ ) {
        sys[ P.epid_index ][ k ].S = sys[ P.epid_index ][ k-1 ].S + ch[ k-1 ].S;
        I += ch[ k-1 ].I;
        sys[ P.epid_index ][ k ].I = I;
        R += ch[ k-1 ].R;
        sys[ P.epid_index ][ k ].R = R;
        sys[ P.epid_index ][ k ].inc = -ch[ k-1 ].S;
        if ( sys[ P.epid_index ][ k ].inc > peak ) {
            peak = sys[ P.epid_index ][ k ].inc;
            t_peak = k;
        }
        sys[ P.epid_index ][ k ].inc_C = ch[ k-1 ].C;
        sys[ P.epid_index ][ k ].inc_A = ch[ k-1 ].A;
        cum_inc -= ch[ k-1 ].S;
        cum_inc_C += ch[ k-1 ].C;
        cum_inc_A += ch[ k-1 ].A;
        sys[ P.epid_index ][ k ].cum_inc = cum_inc;
        sys[ P.epid_index ][ k ].cum_inc_C = cum_inc_C;
        sys[ P.epid_index ][ k ].cum_inc_A = cum_inc_A;
//        ABSsick += ch[ k-1 ].ABSsick;
//        sys[ P.epid_index ][ k ].ABSsick = ABSsick;
//        ABSclosed += ch[ k-1 ].ABSclosed;
//        sys[ P.epid_index ][ k ].ABSclosed = ABSclosed;
        n_conts += ch[ k-1 ].n_conts;
        sys[ P.epid_index ][ k ].n_conts = n_conts;
        n_infs += ch[ k-1 ].n_infs;
        G_infs += ch[ k-1 ].G_infs;
        H_infs += ch[ k-1 ].H_infs;
        W_infs += ch[ k-1 ].W_infs;
        sys[ P.epid_index ][ k ].n_infs = n_infs;
        sys[ P.epid_index ][ k ].G_infs = G_infs;
        sys[ P.epid_index ][ k ].H_infs = H_infs;
        sys[ P.epid_index ][ k ].W_infs = W_infs;
        sys[ P.epid_index ][ k ].colls = sys[ P.epid_index ][ k-1 ].colls + ch[ k-1 ].colls;
        
        sys[ P.epid_index ][ k ].H_S = sys[ P.epid_index ][ k-1 ].H_S + ch[ k-1 ].H_S;
        H_I += ch[ k-1 ].H_I;
        sys[ P.epid_index ][ k ].H_I = H_I;
        H_C += ch[ k-1 ].H_C;
        sys[ P.epid_index ][ k ].H_C = H_C;
        H_R += ch[ k-1 ].H_R;
        sys[ P.epid_index ][ k ].H_R = H_R;
        sys[ P.epid_index ][ k ].H_inc = -ch[ k-1 ].H_S;
        H_cum_inc -= ch[ k-1 ].H_S;
        sys[ P.epid_index ][ k ].H_cum_inc = H_cum_inc;
        
        assert( ( sys[ P.epid_index ][ k ].S >= 0 ) && ( sys[ P.epid_index ][ k ].S <= P.TOT - P.N_NON_SUSC ) );
        assert( ( I >= 0 ) && ( I <= P.TOT ) );
        assert( ( R >= 0 ) && ( R <= P.TOT ) );
        assert( sys[ P.epid_index ][ k ].colls >= 0 );
    }
    P.inc_peak[ P.epid_index ] = ( double ) peak / P.dt;
    P.t_index_peak[ P.epid_index ] = t_peak;
    P.t_last_event[ P.epid_index ] = last_t;
}

void create_cum_gen_R0( void )
{
    int i;
    
    for ( i = 0; i < P.TOT; i++ ) {
        if ( ind[ i ].gen > -1 ) {	// Only individuals that took part in the epidemic (no those escaping, no initial immunes)
            assert( ind[ i ].is != 'S' ); 
            if ( P.extinction_flag[ P.epid_index ] == 'L' ) {
                gen_large[ ind[ i ].gen ].I++;
                if ( ind[ i ].age == 'C' )
                    gen_large[ ind[ i ].gen ].C++;
                else {
                    assert( ind[ i ].age == 'A' );
                    gen_large[ ind[ i ].gen ].A++;
                }
                gen_large[ ind[ i ].gen ].colls += ind[ i ].G_conts - ind[ i ].G_infs;
                gen_large[ ind[ i ].gen ].gen_cum_conts_large += ind[ i ].G_conts + ind[ i ].H_infs + ind[ i ].W_infs;
                gen_large[ ind[ i ].gen ].gen_cum_infs_large += ind[ i ].G_infs + ind[ i ].H_infs + ind[ i ].W_infs;
                P.gen_cum_large[ ind[ i ].gen ]++;
            }
        }
    }
}

void create_cum_r( void ) 
{
    int k, k_in, k_fin, flag = 0;
    int I, colls;
    
    if ( P.TOT*P.CUTOFF < 100 )
        printf( "\nEpidemic %4d: Careful!!! Not enough people for a nice exponential phase!\n", P.epid_index );
    else {
        if ( P.extinction_flag[ P.epid_index ] == 'L' ) {
            I = P.N_INIT_INF;
            colls = 0;
            for ( k = 0; k <= P.TOT_STEPS; k++ ) {
                if ( I > 100 ) {
                    flag = 1;
                    break;
                } else {
                    I += ch[ k ].I;
                    colls += ch[ k ].colls;
                }
            }
            if ( flag == 0 ) {
                printf( "Epidemic %4d: not even reached 100 individuals!\n", P.epid_index+1 );
                return;
            } else {
                if ( colls > 10 ) {
                    printf( "Epidemic %4d: 100 individuals infected, but 10 colls already occurred!\n", P.epid_index+1 );
                    return;
                } else {
                    k_in = k;
                    P.r_100[ P.epid_index ] = k_in;
                    P.n_r++;
                    for ( ; k <= P.TOT_STEPS; k++ ) {
                        if ( colls > 10 ) {
                            flag = 1;
                            break;
                        }
                        P.cum_r[ k - k_in ] += I;
                        I += ch[ k ].I;
                        colls += ch[ k ].colls;
                    }
                    assert( flag == 1 );
                    k_fin = k;
                    P.r_10ghosts[ P.epid_index ] = k_fin;
                    if (  k_fin - k_in < P.r_length )
                        P.r_length = k_fin - k_in;
                }
            }
        }
    }
}

