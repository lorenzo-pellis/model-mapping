//
//  print_output.cpp
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

#include "headers.h"

// In order to print the final epidemic curve, I need to prepare the summary (i.e. to post-process the epidemic - I don't print the list of events itself).
// The output can be a file with final size and other summary information (print_final_size), or a file with the numbers over time (print_real_time).
void print_output( void ) //const char *base_name )
{
    prepare_summary();
    print_final_size(); // Use this function to print the output needed in the model mapping procedure
    //print_real_time(); // Use this function to print the mean epidemic curve (Fig 2 of main paper)

    /* Other functions: */
    //print_generations();
    //print_epidemiological_outputs();
    //print_epidemic_curves();
    //print_correlations();
    //print_r();
}

// The main part of the preparation of the output is to synchronise epidemics at the peak
void prepare_summary( void )
{
    int max_t_peak;
    int z, k, k_rel;

    max_t_peak = 0;
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.t_index_peak[ z ] > max_t_peak )
            max_t_peak = P.t_index_peak[ z ]; // Finding the time of the peak for each epidemic
    }
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        for ( k = max_t_peak - P.t_index_peak[ z ]; k <= P.TOT_STEPS; k++ ) {
            k_rel = k - max_t_peak + P.t_index_peak[ z ]; // This is the point where epidemics are synchronised at the peak
            if ( P.extinction_flag[ z ] == 'L' ) { // If this is a large epidemic...
                P.n_sync_large[ k ]++;
//                if ( k_rel > P.t_school_closure[ z ] / P.dt )
//                    P.n_school_closure[ k ]++;
                sync_large[ k ].S += sys[ z ][ k_rel ].S;
                sync_large[ k ].I += sys[ z ][ k_rel ].I;
                sync_large[ k ].R += sys[ z ][ k_rel ].R;
                sync_large[ k ].inc += sys[ z ][ k_rel ].inc;
                sync_large[ k ].inc_C += sys[ z ][ k_rel ].inc_C;
                sync_large[ k ].inc_A += sys[ z ][ k_rel ].inc_A;
                sync_large[ k ].cum_inc += sys[ z ][ k_rel ].cum_inc;
                sync_large[ k ].cum_inc_C += sys[ z ][ k_rel ].cum_inc_C;
                sync_large[ k ].cum_inc_A += sys[ z ][ k_rel ].cum_inc_A;
//                sync_large[ k ].ABSsick += sys[ z ][ k_rel ].ABSsick;
//                sync_large[ k ].ABSclosed += sys[ z ][ k_rel ].ABSclosed;
                sync_large[ k ].R0 += sys[ z ][ k_rel ].n_conts;
                sync_large[ k ].Rt += sys[ z ][ k_rel ].n_infs;
                sync_large[ k ].R_in_G += sys[ z ][ k_rel ].G_infs;
                sync_large[ k ].R_in_H += sys[ z ][ k_rel ].H_infs;
                sync_large[ k ].R_in_W += sys[ z ][ k_rel ].W_infs;
                sync_large[ k ].colls += sys[ z ][ k_rel ].colls;
                sync_large[ k ].H_S += sys[ z ][ k_rel ].H_S;
                sync_large[ k ].H_I += sys[ z ][ k_rel ].H_I;
                sync_large[ k ].H_C += sys[ z ][ k_rel ].H_C;
                sync_large[ k ].H_R += sys[ z ][ k_rel ].H_R;
                sync_large[ k ].H_inc += sys[ z ][ k_rel ].H_inc;
                sync_large[ k ].H_cum_inc += sys[ z ][ k_rel ].H_cum_inc;
                sync_large[ k ].W_cum_inc += sys[ z ][ k_rel ].W_cum_inc;
            }
        }
    }
    for ( k = 0; k < P.r_GAP; k++ ) {
        sync_large[ k ].cum_inc_gap += sync_large[ k ].cum_inc;
        sync_large[ k ].R0gap += sync_large[ k ].R0;
        sync_large[ k ].Rtgap += sync_large[ k ].Rt;
    }
    for ( k = P.r_GAP; k <= P.TOT_STEPS; k++ ) {
        sync_large[ k ].r = ( log( sync_large[ k ].I > 0 ? sync_large[ k ].I : 1 ) - log( sync_large[ k - P.r_GAP ].I > 0 ? sync_large[ k - P.r_GAP ].I : 1 ) ) / ( P.r_GAP * P.dt );
        sync_large[ k ].H_r = ( log( sync_large[ k ].I > 0 ? sync_large[ k ].I : 1 ) - log( sync_large[ k - P.r_GAP ].I > 0 ? sync_large[ k - P.r_GAP ].I : 1 ) ) / ( P.r_GAP * P.dt );
        sync_large[ k ].cum_inc_gap += sync_large[ k ].cum_inc - sync_large[ k - P.r_GAP ].cum_inc;
        sync_large[ k ].R0gap += sync_large[ k ].R0 - sync_large[ k - P.r_GAP ].R0;
        sync_large[ k ].Rtgap += sync_large[ k ].Rt - sync_large[ k - P.r_GAP ].Rt;
        if ( sync_large[ k ].cum_inc_gap == 0 )
        {
            sync_large[ k ].cum_inc_gap = 1;
        }
    }
    if ( P.n_large == 0 )
        P.n_large = 1;
    for ( k = 0; k <= P.TOT_STEPS; k++ ) {
        sync_large[ k ].S /= P.n_sync_large[ k ];
        sync_large[ k ].I /= P.n_sync_large[ k ];
        sync_large[ k ].R /= P.n_sync_large[ k ];
        sync_large[ k ].inc /= P.n_sync_large[ k ];
        sync_large[ k ].inc_C /= P.n_sync_large[ k ];
        sync_large[ k ].inc_A /= P.n_sync_large[ k ];
//        sync_large[ k ].ABSsick /= P.n_sync_large[ k ];
//        sync_large[ k ].ABSclosed /= P.n_sync_large[ k ];
        sync_large[ k ].R0 /= sync_large[ k ].cum_inc;
        sync_large[ k ].Rt /= sync_large[ k ].cum_inc;
        sync_large[ k ].R_in_G /= sync_large[ k ].cum_inc;
        sync_large[ k ].R_in_H /= sync_large[ k ].cum_inc;
        sync_large[ k ].R_in_W /= sync_large[ k ].cum_inc;
        sync_large[ k ].R0gap /= sync_large[ k ].cum_inc_gap;
        sync_large[ k ].Rtgap /= sync_large[ k ].cum_inc_gap;
        sync_large[ k ].colls /= P.n_sync_large[ k ];
        sync_large[ k ].H_S /= P.n_sync_large[ k ];
        sync_large[ k ].H_I /= P.n_sync_large[ k ];
        sync_large[ k ].H_C /= P.n_sync_large[ k ];
        sync_large[ k ].H_R /= P.n_sync_large[ k ];
        sync_large[ k ].H_inc /= P.n_sync_large[ k ];
        sync_large[ k ].cum_inc /= P.n_sync_large[ k ];
        sync_large[ k ].cum_inc_C /= P.n_sync_large[ k ];
        sync_large[ k ].cum_inc_A /= P.n_sync_large[ k ];
    }
}

void print_final_size( void )
{
    int z, test;
    double temp_large, mean_large, var_large, mean_inc, var_inc;
    double mean_t_peak, var_t_peak, mean_t_last_event, var_t_last_event;
//    double mean_t_school_closure, var_t_school_closure;
    double temp_large_H, mean_large_H, var_large_H;
    double temp_large_W, mean_large_W, var_large_W;
    FILE *p_av;
    char name[ 1024 ] = "";
    char *p_name = name;
    
    /* Final size */
    temp_large = 0;
    temp_large_H = 0;
    temp_large_W = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += P.final_size[ z ];
            temp_large_H += P.final_size_H[ z ];
            temp_large_W += P.final_size_W[ z ];
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    mean_large = ( double ) temp_large / P.n_large;
    mean_large_H = ( double ) temp_large_H / P.n_large;
    mean_large_W = ( double ) temp_large_W / P.n_large;
    
    temp_large = 0;
    temp_large_H = 0;
    temp_large_W = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += ( P.final_size[ z ] - mean_large ) * ( P.final_size[ z ] - mean_large );
            temp_large_H += ( P.final_size_H[ z ] - mean_large_H ) * ( P.final_size_H[ z ] - mean_large_H );
            temp_large_W += ( P.final_size_W[ z ] - mean_large_W ) * ( P.final_size_W[ z ] - mean_large_W );
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    var_large = ( double ) temp_large / P.n_large;
    var_large_H = ( double ) temp_large_H / P.n_large;
    var_large_W = ( double ) temp_large_W / P.n_large;
    
    /* Incidence */
    temp_large = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += P.inc_peak[ z ];
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    mean_inc = ( double ) temp_large / P.n_large;
    
    temp_large = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += ( P.inc_peak[ z ] - mean_inc ) * ( P.inc_peak[ z ] - mean_inc );
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    var_inc = ( double ) temp_large / P.n_large;
    
    /* Time of Incidence */
    temp_large = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += P.t_index_peak[ z ]*P.dt;
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    mean_t_peak = ( double ) temp_large / P.n_large;
    
    temp_large = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += ( P.t_index_peak[ z ]*P.dt - mean_t_peak ) * ( P.t_index_peak[ z ]*P.dt - mean_t_peak );
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    var_t_peak = ( double ) temp_large / P.n_large;
    
//    /* Time of School closure */
//    temp_large = 0;
//    test = 0;
//    
//    for ( z = 0; z < P.HOW_MANY; z++ ) {
//        if ( P.extinction_flag[ z ] == 'L' ) {
//            temp_large += P.t_school_closure[ z ];
//            test++;
//        }
//    }
//    if ( test == 0 )
//        test = 1;
//    assert( test == P.n_large );
//    mean_t_school_closure = ( double ) temp_large / P.n_large;
//    
//    temp_large = 0;
//    test = 0;
//    
//    for ( z = 0; z < P.HOW_MANY; z++ ) {
//        if ( P.extinction_flag[ z ] == 'L' ) {
//            temp_large += ( P.t_school_closure[ z ] - mean_t_school_closure ) * ( P.t_school_closure[ z ] - mean_t_school_closure );
//            test++;
//        }
//    }
//    if ( test == 0 )
//        test = 1;
//    assert( test == P.n_large );
//    var_t_school_closure = ( double ) temp_large / P.n_large;
    
    /* Time of last event */
    temp_large = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += P.t_last_event[ z ];
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    mean_t_last_event = ( double ) temp_large / P.n_large;
    
    temp_large = 0;
    test = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        if ( P.extinction_flag[ z ] == 'L' ) {
            temp_large += ( P.t_last_event[ z ] - mean_t_last_event ) * ( P.t_last_event[ z ] - mean_t_last_event );
            test++;
        }
    }
    if ( test == 0 )
        test = 1;
    assert( test == P.n_large );
    var_t_last_event = ( double ) temp_large / P.n_large;
    
    printf( "Average final size (cond on a large epidemic):\t%7.3f\n", mean_large );
    printf( "Variance (cond on a large epidemic):\t%7.3f\n\n", var_large );
    printf( "Average final size for HOUSEHOLDS (cond on a large epidemic):\t%7.3f\n", mean_large_H );
    printf( "Variance (cond on a large epidemic):\t%7.3f\n\n", var_large_H );
    printf( "Average final size for WORKPLACES (cond on a large epidemic):\t%7.3f\n", mean_large_W );
    printf( "Variance (cond on a large epidemic):\t%7.3f\n\n", var_large_W );
    
    sprintf(p_name,"%saverages.dat",P.base_name);
    p_av = fopen(p_name, "w" );
    fprintf( p_av, "HOW_MANY\tn_large\tProb no ext\tAverage final size(%%)\tStandard error\tAverage inc peak(%%)\tStandard error inc peak\tAv time peak (in gens)\tSd time peak (in gens)\tAv time last event\tSd time last event\n" );
    fprintf( p_av, "%7d\t%7d\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n", P.HOW_MANY, P.n_large, ( ( double ) 100 ) * P.n_large/P.HOW_MANY, 100*mean_large/P.TOT, 100*sqrt( var_large )/P.TOT, 100*mean_inc/P.TOT, 100*sqrt( var_inc )/P.TOT, mean_t_peak / P.GEN_TIME, sqrt(var_t_peak) / P.GEN_TIME, mean_t_last_event, sqrt(var_t_last_event) );
    fclose( p_av );
//    fprintf( p_av, "HOW_MANY\tn_large\tProb no ext\tAverage final size(%%)\tStandard error\tAverage inc peak(%%)\tStandard error inc peak\tAv time peak (in gens)\tSd time peak (in gens)\tAv time school closure\tSd time school closure\tAv time last event\tSd time last event\tAv time school reopen\tGap closure-peak\tGap reopen-peak\n" );
//    fprintf( p_av, "%7d\t%7d\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n", P.HOW_MANY, P.n_large, ( ( double ) 100 ) * P.n_large/P.HOW_MANY, 100*mean_large/P.TOT, 100*sqrt( var_large )/P.TOT, 100*mean_inc/P.TOT, 100*sqrt( var_inc )/P.TOT, mean_t_peak / P.GEN_TIME, sqrt(var_t_peak) / P.GEN_TIME, mean_t_school_closure, sqrt(var_t_school_closure), mean_t_last_event, sqrt(var_t_last_event), mean_t_school_closure+P.DURATION_CLOSURE, mean_t_school_closure-mean_t_peak, mean_t_school_closure+P.DURATION_CLOSURE-mean_t_peak );
//    fclose( p_av );
}

void print_real_time( void ) //char* pre_name )
{
    int k;
    FILE *p_large;
    
    char name[ 1024 ] = "";
    char *p_name = name;
    
    p_name = strcpy( p_name, P.base_name );
    sprintf(p_name,"%ssync_large.dat",P.base_name);
    p_large = fopen( p_name, "w" );
    // This is what is printed:
    // t    S   I   R   inc   cum_inc   r    R0    Rt   collisions  #epidemic_started_by_now   incA   incC   cum_incA   cum_incC   prop_incA   prop_incC
    
    for ( k = 0; k <= P.TOT_STEPS; k++ )
        fprintf( p_large, "%3.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%-2.3f\t%2.3f\t%2.3f\t%8.2f\t%3d\t\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%0.3f\t%0.3f\n",
                sync_large[ k ].t, sync_large[ k ].S, sync_large[ k ].I, sync_large[ k ].R, sync_large[ k ].inc, sync_large[ k ].cum_inc, sync_large[ k ].r, sync_large[ k ].R0, sync_large[ k ].Rt, sync_large[ k ].colls, P.n_sync_large[ k ], sync_large[ k ].inc_A, sync_large[ k ].inc_C, sync_large[ k ].cum_inc_A, sync_large[ k ].cum_inc_C, sync_large[ k ].inc_A / ( ( sync_large[ k ].inc == 0 ) ? 1 : sync_large[ k ].inc ), sync_large[ k ].inc_C / ( ( sync_large[ k ].inc == 0 ) ? 1 : sync_large[ k ].inc ) );
//    for ( k = 0; k <= P.TOT_STEPS; k++ )
//        fprintf( p_large, "%3.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%-2.3f\t%2.3f\t%2.3f\t%8.2f\t%3d\t%3d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%0.3f\t%0.3f\n",
//                sync_large[ k ].t, sync_large[ k ].S, sync_large[ k ].I, sync_large[ k ].R, sync_large[ k ].inc, sync_large[ k ].cum_inc, sync_large[ k ].r, sync_large[ k ].R0, sync_large[ k ].Rt, sync_large[ k ].colls, P.n_sync_large[ k ], P.n_school_closure[ k ], sync_large[ k ].inc_A, sync_large[ k ].inc_C, sync_large[ k ].cum_inc_A, sync_large[ k ].cum_inc_C, sync_large[ k ].inc_A / ( ( sync_large[ k ].inc == 0 ) ? 1 : sync_large[ k ].inc ), sync_large[ k ].inc_C / ( ( sync_large[ k ].inc == 0 ) ? 1 : sync_large[ k ].inc ) );
    fclose( p_large );
    }

//void print_generations( void )
//{
//    FILE *p_large; //*p_global,
//    int g;
//    char name[ 1024 ] = "";
//    char *p_name = name;
//    
//    for ( g = 0; g < P.MAX_GEN; g++ ) {
//        //gen_global[ g ].Rt = ( double ) gen_cum_global[ g+1 ] / gen_cum_global[ g ];
//        gen_large[ g ].R0 = ( double ) gen_large[ g ].gen_cum_conts_large / P.gen_cum_large[ g ];
//        gen_large[ g ].Rt = ( double ) gen_large[ g ].gen_cum_infs_large / P.gen_cum_large[ g ];
//    }
//    for ( g = 1; g < P.MAX_GEN; g++ ) {
//        gen_large[ g ].TransmToC = ( double ) gen_large[ g ].C / ( ( gen_large[ g-1 ].I == 0 ) ? 1 : gen_large[ g-1 ].I );
//        gen_large[ g ].TransmToA = ( double ) gen_large[ g ].A / ( ( gen_large[ g-1 ].I == 0 ) ? 1 : gen_large[ g-1 ].I );
//        gen_large[ g ].propC = ( double ) gen_large[ g ].C / ( ( gen_large[ g ].I == 0 ) ? 1 : gen_large[ g ].I );
//        gen_large[ g ].propA = ( double ) gen_large[ g ].A / ( ( gen_large[ g ].I == 0 ) ? 1 : gen_large[ g ].I );
//        //printf( "%2.3f\t%2.3f\n", gen_large[ g ].TransmToC, gen_large[ g ].TransmToA );
//    }
//    sprintf(p_name,"%sgenerations.dat",P.base_name);
//    p_large = fopen( p_name, "w" );
//    for ( g = 0; g < P.MAX_GEN; g++ )
//        fprintf( p_large, "%3d\t\t%2.2f\t%2.2f\t%8d\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n", gen_large[ g ].number, gen_large[ g ].R0, gen_large[ g ].Rt, gen_large[ g ].colls, gen_large[ g ].TransmToC, gen_large[ g ].TransmToA, gen_large[ g ].propA, gen_large[ g ].propC );
//    fclose( p_large );
//}

//void print_epidemic_curves( void )
//{
//    FILE *p_curves, *p_curves_synchro;
//    FILE *p_H_curves, *p_H_curves_synchro;
//    int max_t_peak = 0;
//    int z, k, k_rel;
//    char name[ 100 ] = "";
//    char *p_name = name;
//    
//    
//    for ( z = 0; z < P.HOW_MANY; z++ )
//        if ( P.t_index_peak[ z ] > max_t_peak )
//            max_t_peak = P.t_index_peak[ z ];
//    
//    p_name = strcpy( p_name, P.base_name );
//    p_curves = fopen( strcat( p_name, "curves.dat" ), "w" );
//    for ( k = 0; k <= P.TOT_STEPS; k++ ) {
//        fprintf( p_curves, "%2.2f", ( double ) k*P.dt );
//        for ( z = 0; z < P.HOW_MANY; z++ )
//            fprintf( p_curves, "\t%8d", sys[ z ][ k ].I );
//        fprintf( p_curves, "\n" );
//    }
//    fclose( p_curves );
//    /* Households */
//    p_name = strcpy( p_name, P.base_name );
//    p_H_curves = fopen( strcat( p_name, "H_curves.dat" ), "w" );
//    for ( k = 0; k <= P.TOT_STEPS; k++ ) {
//        fprintf( p_H_curves, "%2.2f", ( double ) k*P.dt );
//        for ( z = 0; z < P.HOW_MANY; z++ )
//            fprintf( p_H_curves, "\t%8d", sys[ z ][ k ].H_I );
//        fprintf( p_H_curves, "\n" );
//    }
//    fclose( p_H_curves );
//    
//    p_name = strcpy( p_name, P.base_name );
//    p_curves_synchro = fopen( strcat( p_name, "curves_synchro.dat" ), "w" );
//    for ( k = 0; k <= P.TOT_STEPS; k++ ) {
//        fprintf( p_curves_synchro, "%2.2f", ( double ) k*P.dt );
//        for ( z = 0; z < P.HOW_MANY; z++ ) {
//            k_rel = k - max_t_peak + P.t_index_peak[ z ];
//            if ( k_rel < 0 )
//                fprintf( p_curves_synchro, "\t%8d", 0 );
//            else
//                fprintf( p_curves_synchro, "\t%8d", sys[ z ][ k_rel ].I );
//        }
//        fprintf( p_curves_synchro, "\n" );
//    }
//    fclose( p_curves_synchro );
//    /* Households */
//    p_name = strcpy( p_name, P.base_name );
//    p_H_curves_synchro = fopen( strcat( p_name, "H_curves_synchro.dat" ), "w" );
//    for ( k = 0; k <= P.TOT_STEPS; k++ ) {
//        fprintf( p_H_curves_synchro, "%2.2f", ( double ) k*P.dt );
//        for ( z = 0; z < P.HOW_MANY; z++ ) {
//            k_rel = k - max_t_peak + P.t_index_peak[ z ];
//            if ( k_rel < 0 )
//                fprintf( p_H_curves_synchro, "\t%8d", 0 );
//            else
//                fprintf( p_H_curves_synchro, "\t%8d", sys[ z ][ k_rel ].H_I );
//        }
//        fprintf( p_H_curves_synchro, "\n" );
//    }
//    fclose( p_H_curves_synchro );
//}

//void print_r( void )
//{
//    printf( "\nGood exp phase seen in %4d epidemics among %4d large ones (total of %4d)\n", P.n_r, P.n_large, P.HOW_MANY );
//    printf( "The minimum time interval for exponential phase is: %3.2f\n", ( double ) P.r_length * P.dt );
//    printf( "Best guess: r = %1.4f\n", log( ( double ) P.cum_r[ P.r_length - 1 ] / P.cum_r[ 0 ] ) / ( P.r_length * P.dt ) );
//}

//void print_epidemiological_outputs( void )
//{
//    int z, k;
//    int s, fs; // size and final size
//    FILE *p_strat, *p_env;
//    char name[ 100 ] = "";
//    char *p_name = name;
//    
//    //table = calloc( P.MAX_W_SIZE * ( P.MAX_W_SIZE + 1 ), sizeof( int ) );
//    //table_ptr = calloc( P.MAX_W_SIZE, sizeof( int* ) );
//    //for ( s = 0; s < P.MAX_W_SIZE; s++ )
//    //	table_ptr[ s ] = &table[ s * ( P.MAX_W_SIZE + 1 ) ];
//    //sizes_vector = calloc( P.MAX_W_SIZE, sizeof( int ) );
//    //for ( w = 0; w < P.n_W; w++ )
//    //{
//    //	sizes_vector[ ( W_list[ w ].sizeA + W_list[ w ].sizeC ) - 1 ]++;
//    //	table_ptr[ ( W_list[ w ].sizeA + W_list[ w ].sizeC ) - 1 ][ W_list[ w ].ufs ]++;
//    //}
//    
//    for ( s = 0; s < P.MAX_W_SIZE; s++ )
//    {
//        for ( fs = 0; fs <= s+1; fs++ )
//        {
//            for ( z = 0; z < P.HOW_MANY; z++ )
//                if ( P.extinction_flag[ z ] == 'L' )
//                    P.table_mean_fs_ptr[ s ][ fs ] += P.table_h_fs_ptr[ z ][ s ][ fs ];
//            P.table_mean_fs_ptr[ s ][ fs ] /= P.n_large;
//            for ( z = 0; z < P.HOW_MANY; z++ )
//                if ( P.extinction_flag[ z ] == 'L' )
//                    P.table_sd_fs_ptr[ s ][ fs ] += ( P.table_h_fs_ptr[ z ][ s ][ fs ] - P.table_mean_fs_ptr[ s ][ fs ] ) * ( P.table_h_fs_ptr[ z ][ s ][ fs ] - P.table_mean_fs_ptr[ s ][ fs ] );
//            P.table_sd_fs_ptr[ s ][ fs ] = sqrt( P.table_sd_fs_ptr[ s ][ fs ] / P.n_large );
//        }
//        for ( z = 0; z < P.HOW_MANY; z++ )
//            if ( P.extinction_flag[ z ] == 'L' )
//                P.vector_mean_freq_h_sizes[ s ] += P.table_freq_h_sizes_ptr[ z ][ s ];
//        P.vector_mean_freq_h_sizes[ s ] /= P.n_large;
//        for ( z = 0; z < P.HOW_MANY; z++ )
//            if ( P.extinction_flag[ z ] == 'L' )
//                P.vector_sd_freq_h_sizes[ s ] += ( P.table_freq_h_sizes_ptr[ z ][ s ] - P.vector_mean_freq_h_sizes[ s ] ) * ( P.table_freq_h_sizes_ptr[ z ][ s ] - P.vector_mean_freq_h_sizes[ s ] );
//        P.vector_sd_freq_h_sizes[ s ] = sqrt( P.vector_sd_freq_h_sizes[ s ] / P.n_large );
//        
//    }
//    
//    sprintf(p_name,"%sfs_strat.dat",P.base_name);
//    p_strat = fopen(p_name, "w" );
//    for ( s = 0; s < P.MAX_W_SIZE; s++ )
//    {
//        fprintf( p_strat, "%7.2f\t\t", P.vector_mean_freq_h_sizes[ s ] );
//        for ( fs = 0; fs < P.MAX_W_SIZE; fs++ )
//            fprintf( p_strat, "%7.2f\t", P.table_mean_fs_ptr[ s ][ fs ] );
//        fprintf( p_strat, "%7.2f\n", P.table_mean_fs_ptr[ s ][ P.MAX_W_SIZE ] );
//    }
//    fprintf( p_strat, "\n\n" );
//    for ( s = 0; s < P.MAX_W_SIZE; s++ )
//    {
//        fprintf( p_strat, "%7.2f\t\t", P.vector_sd_freq_h_sizes[ s ] );
//        for ( fs = 0; fs < P.MAX_W_SIZE; fs++ )
//            fprintf( p_strat, "%7.2f\t", P.table_sd_fs_ptr[ s ][ fs ] );
//        fprintf( p_strat, "%7.2f\n", P.table_sd_fs_ptr[ s ][ P.MAX_W_SIZE ] );
//    }
//    fclose( p_strat );
//    
//    // Printing the proportion of infections occurring in different environments
//    sprintf(p_name,"%sprop_env.dat",P.base_name);
//    p_env = fopen(p_name, "w" );
//    for ( k = 0; k <= P.TOT_STEPS; k++ )
//        fprintf( p_env, "%3.2f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n", sync_large[ k ].t, sync_large[ k ].Rt, sync_large[ k ].R_in_G / sync_large[ k ].Rt, sync_large[ k ].R_in_H / sync_large[ k ].Rt, sync_large[ k ].R_in_W / sync_large[ k ].Rt );
//    fclose( p_env );
//}




