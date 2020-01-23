//
//  main.cpp
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

#include "headers.h"
#include "typeconstr.h"
using namespace std;


int main(int argc, char** argv)
{
    int a, b, c, m1, m2, m3, m4, m5, aa, bb, cc, dd;
//    int d, e, f, g;
    string name;
    char help[ 64 ];
    
    //blah = "hi!\n"; cout << blah;
    clock_t time = clock();
    
    read_param( argc, argv );
    time = clock() - time;
    printf( "Time needed to read the parameters:\t%4.5f seconds\n\n", ( double ) time / CLOCKS_PER_SEC );
    
    prepare_model();
    time = clock() - time;
    printf( "Time needed to load the network:\t%4.5f seconds\n\n", ( double ) time / CLOCKS_PER_SEC );
    
    for ( a = 0; a < P.Rg_LEN; a++ )
    {
        for ( b = 0; b < P.Rh_LEN; b++ )
        {
            for ( c = 0; c < P.Rw_LEN; c++ )
            {
                for ( m1 = 0; m1 < P.REL_SUSC_C_VS_A_IN_COMM_VEC_LEN; m1++ )
                {
                    P.REL_SUSC_C_VS_A_IN_COMM = P.REL_SUSC_C_VS_A_IN_COMM_VEC[ m1 ];
                    if ( SAME_SIGMA )
                        P.REL_SUSC_C_VS_A_IN_HOUSEHOLDS = P.REL_SUSC_C_VS_A_IN_COMM;
                    for ( m2 = 0; m2 < P.REL_INF_C_VS_A_IN_COMM_VEC_LEN; m2++ )
                    {
                        P.REL_INF_C_VS_A_IN_COMM = P.REL_INF_C_VS_A_IN_COMM_VEC[ m2 ];
                        if ( SAME_RHO )
                            P.REL_INF_C_VS_A_IN_HOUSEHOLDS = P.REL_INF_C_VS_A_IN_COMM;
                        for ( m3 = 0; m3 < P.ASS_IN_COMM_VEC_LEN; m3++ )
                        {
                            P.ASS_IN_COMM = P.ASS_IN_COMM_VEC[ m3 ];
                            if ( SAME_PHI )
                                P.ASS_IN_HOUSEHOLDS = P.ASS_IN_COMM;
                            for ( m4 = 0; m4 < P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC_LEN; m4++ )
                            {
                                P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM = P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM_VEC[ m3 ];
                                for ( m5 = 0; m5 < P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC_LEN; m5++ )
                                {
                                    P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS = P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS_VEC[ m3 ];
//                                    for ( d = 0; d < P.WITHIN_SCHOOL_PREV_TRIGGER_VEC_LEN; d++ )
//                                    {
//                                        P.WITHIN_SCHOOL_PREV_TRIGGER = P.WITHIN_SCHOOL_PREV_TRIGGER_VEC[ d ];
//                                        for ( e = 0; e < P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC_LEN; e++ )
//                                        {
//                                            P.PERC_SCHOOLS_INFECTED_TRIGGER =P.PERC_SCHOOLS_INFECTED_TRIGGER_VEC[ e ];
//                                            for ( f = 0; f < P.DELAY_VEC_LEN; f++ )
//                                            {
//                                                P.DELAY = P.DELAY_VEC[ f ];
//                                                for ( g = 0; g < P.DURATION_CLOSURE_VEC_LEN; g++ )
//                                                {
//                                                    P.DURATION_CLOSURE = P.DURATION_CLOSURE_VEC[ g ];
                                                    for ( aa = 0; aa < P.INF_MEAN_LEN; aa++ )
                                                    {
                                                        P.INF_MEAN = P.INF_MEAN_VEC[ aa ];
                                                        for ( bb = 0; bb < P.INF_ALPHA_LEN; bb++ )
                                                        {
                                                            P.INF_ALPHA = P.INF_ALPHA_VEC[ bb ];
                                                            for ( cc = 0; cc < P.GEN_TIME_LEN; cc++ )
                                                            {
                                                                P.GEN_TIME = P.GEN_TIME_VEC[ cc ];
                                                                for ( dd = 0; dd < P.GEN_ALPHA_LEN; dd++ )
                                                                {
                                                                    P.GEN_ALPHA = P.GEN_ALPHA_VEC[ dd ];
                                                                    if ( P.MODEL_TYPE == 0 )
                                                                        P.r_MEAN = P.INF_MEAN;
                                                                    else if ( P.MODEL_TYPE == 1 )
                                                                        P.r_MEAN = P.LAT_MEAN + P.INF_MEAN;
                                                                    else
                                                                        P.r_MEAN = P.GEN_TIME;
                                                                    P.r_GAP = ( int ) P.r_MEAN / P.dt;
                                                                    
                                                                    initialyse_basic_quantities();
                                                                    time = clock() - time;
                                                                    printf( "Time needed to initialyse:\t%4.5f seconds\n\n", ( double ) time / CLOCKS_PER_SEC );
                                                                    
                                                                    for ( P.epid_index = 0; P.epid_index < P.HOW_MANY; P.epid_index++ )
                                                                    {
                                                                        epidemic( P.Rg[ a ], P.Rh[ b ], P.Rw[ c ] );
                                                                        // printf( "Epidemic %3d of %3d finished. Final size: %7d. Peak at time: %2.2f.\nCount: %7d. Last recovery: %2.2f. School closed at time %2.2f and reopened at time %2.2f\n", P.epid_index+1, P.HOW_MANY, P.final_size[ P.epid_index ], ( double ) P.t_index_peak[ P.epid_index ] * P.dt, P.ev_count, ( double ) P.t_last_event[ P.epid_index ], P.TIME_CLOSURE, P.TIME_REOPEN );
                                                                        printf( "Epidemic %3d of %3d finished. Final size: %7d. Peak at time: %2.2f.\nCount: %7d. Last recovery: %2.2f\n", P.epid_index+1, P.HOW_MANY, P.final_size[ P.epid_index ], ( double ) P.t_index_peak[ P.epid_index ] * P.dt, P.ev_count, ( double ) P.t_last_event[ P.epid_index ] );
                                                                     }
                                                                    time = clock() - time;
                                                                    printf( "Time needed to run the epidemics:\t%4.5f seconds\n\n", ( double ) time / CLOCKS_PER_SEC );
                                                                    
                                                                    name = "";
                                                                    if ( OUT_NETWORK_TYPE )
                                                                        name = name + P.SOCIAL_STRUCTURE_NETWORK_TYPE + "_";
                                                                    if ( OUT_Rg ) {
                                                                        sprintf( help, "Rg%.3f_", P.Rg[ a ] );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_Rh ) {
                                                                        sprintf( help, "Rh%.3f_", P.Rh[ b ] );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_Rw ) {
                                                                        sprintf( help, "Rw%.3f_", P.Rw[ b ] );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_SIGMA ) {
                                                                        sprintf( help, "sigma%.4f_", P.REL_SUSC_C_VS_A_IN_COMM );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_RHO ) {
                                                                        sprintf( help, "rho%.1f_", P.REL_INF_C_VS_A_IN_COMM );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_PHI ) {
                                                                        sprintf( help, "ass%.3f_", P.ASS_IN_COMM );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_GAMMA ) {
                                                                        sprintf( help, "gammaG%.2f_", P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_COMM );
                                                                        name = name + help;
                                                                        sprintf( help, "H%.2f_", P.RATIO_CONTACT_RATE_ADULTS_VS_CHILDREN_IN_HOUSEHOLDS );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_WITHIN_SCHOOL_PREV_TRIGGER ) {
                                                                        sprintf( help, "trsch%d_", P.WITHIN_SCHOOL_PREV_TRIGGER );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_PERC_SCHOOLS_INFECTED_TRIGGER ) {
                                                                        sprintf( help, "trall%d_", P.PERC_SCHOOLS_INFECTED_TRIGGER );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_DELAY ) {
                                                                        sprintf( help, "delay%.1f_", P.DELAY );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_DURATION ) {
                                                                        sprintf( help, "dur%d_", P.DURATION_CLOSURE );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_INF_MEAN ) {
                                                                        sprintf( help, "infmean%.2f_", P.INF_MEAN );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_INF_ALPHA ) {
                                                                        sprintf( help, "infalpha%d_", P.INF_ALPHA );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_GEN_TIME ) {
                                                                        sprintf( help, "Tg%.2f_", P.GEN_TIME );
                                                                        name = name + help;
                                                                    }
                                                                    if ( OUT_GEN_ALPHA ) {
                                                                        sprintf( help, "genalpha%d_", P.GEN_ALPHA );
                                                                        name = name + help;
                                                                    }
                                                                    sprintf( P.base_name, "%s_", name.c_str() );

                                                                    print_output();
                                                                    time = clock() - time;
                                                                    printf( "Time needed to print:\t%4.5f seconds\n\n", ( double ) time / CLOCKS_PER_SEC );
                                                                }
                                                            }
                                                        }
//                                                    }
//                                                }
//                                            }
//                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}


