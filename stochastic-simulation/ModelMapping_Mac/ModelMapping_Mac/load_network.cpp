//
//  load_network.cpp
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//  All rights reserved.
//

#include "headers.h"

// This function sets up the simulation by loading the network of the social structure, allocating all memory that can be set up to start with and setting the seed for the random number generator.
void prepare_model( void )
{
    int z, e;
    
    // Set the random number generator
    // P.seed = time( NULL );
    P.seed = -7;
    P.idum = &P.seed;
    
    ind = ( IND* ) calloc( P.TOT, sizeof( IND ) );
    ev = ( EVENT* ) calloc( P.MAX_EV, sizeof( EVENT ) );
    h = ( HELPER* ) calloc( P.TOT_h_STEPS+1, sizeof( HELPER ) );
    gen_large = ( GEN* ) calloc( P.MAX_GEN, sizeof( GEN ) );
    ch = ( CHANGES* ) calloc( P.TOT_STEPS, sizeof( CHANGES ) );
    sys_matrix = ( SYS* ) calloc( P.HOW_MANY * ( P.TOT_STEPS + 1 ), sizeof( SYS ) );
    sys_line = ( SYS** ) calloc( P.HOW_MANY, sizeof( SYS* ) );
    for ( z = 0; z < P.HOW_MANY; z++ )
        sys_line[ z ] = &sys_matrix[ z*( P.TOT_STEPS+1 ) ];
    sys = sys_line;
    sync_large = ( SUM* ) calloc( P.TOT_STEPS + 1, sizeof( SUM ) );
    
    load_network(); // I need the value of MAX_W_SIZE
    P.F_A = ( double ) P.TOT_A / P.TOT;
    P.F_C = ( double ) P.TOT_C / P.TOT;
    printf( "F_A: %f\n", P.F_A );
    printf( "F_C: %f\n", P.F_C );
    
    
    P.extinction_flag = ( char* ) calloc( P.HOW_MANY + 1, sizeof( char ) );
    P.final_size = ( int* ) calloc( P.HOW_MANY, sizeof( int ) );
    P.final_size_H = ( int* ) calloc( P.HOW_MANY, sizeof( int ) );
    P.final_size_W = ( int* ) calloc( P.HOW_MANY, sizeof( int ) );
    P.inc_peak = ( double* ) calloc( P.HOW_MANY, sizeof( double ) );
    P.t_index_peak = ( int* ) calloc( P.HOW_MANY, sizeof( int ) );
    P.t_last_event = ( double* ) calloc( P.HOW_MANY, sizeof( double ) );
//    P.t_school_closure = ( double* ) calloc( P.HOW_MANY, sizeof( double ) );
//    P.t_school_reopen = ( double* ) calloc( P.HOW_MANY, sizeof( double ) );
    P.table_h_fs = ( int* ) calloc( P.HOW_MANY * P.MAX_W_SIZE * ( P.MAX_W_SIZE + 1 ), sizeof( int ) );
    P.table_h_fs_ptr_intermediate = ( int** ) calloc( P.HOW_MANY * P.MAX_W_SIZE, sizeof( int* ) );
    P.table_h_fs_ptr = ( int*** ) calloc( P.HOW_MANY, sizeof( int** ) );
    for ( e = 0; e < P.HOW_MANY; e++ )
    {
        P.table_h_fs_ptr[ e ] = &P.table_h_fs_ptr_intermediate[ e * P.MAX_W_SIZE ];
        for ( z = 0; z < P.MAX_W_SIZE; z++ )
            P.table_h_fs_ptr[ e ][ z ] = &P.table_h_fs[ e * P.MAX_W_SIZE * ( P.MAX_W_SIZE + 1 ) + z * ( P.MAX_W_SIZE + 1 ) ];
    }
    P.table_mean_fs_strat = ( double* ) calloc( P.MAX_W_SIZE * ( P.MAX_W_SIZE + 1 ), sizeof( double ) );
    P.table_mean_fs_ptr = ( double** ) calloc( P.MAX_W_SIZE, sizeof( double* ) );
    for ( z = 0; z < P.MAX_W_SIZE; z++ )
        P.table_mean_fs_ptr[ z ] = &P.table_mean_fs_strat[ z * ( P.MAX_W_SIZE + 1 ) ];
    P.table_sd_fs_strat = ( double* ) calloc( P.MAX_W_SIZE * ( P.MAX_W_SIZE + 1 ), sizeof( double ) );
    P.table_sd_fs_ptr = ( double** ) calloc( P.MAX_W_SIZE, sizeof( double* ) );
    for ( z = 0; z < P.MAX_W_SIZE; z++ )
        P.table_sd_fs_ptr[ z ] = &P.table_sd_fs_strat[ z * ( P.MAX_W_SIZE + 1 ) ];
    P.table_freq_h_sizes = ( int* ) calloc( P.HOW_MANY * P.MAX_W_SIZE, sizeof( int ) );
    P.table_freq_h_sizes_ptr = ( int** ) calloc( P.HOW_MANY, sizeof( int* ) );
    for ( e = 0; e < P.HOW_MANY; e++ )
        P.table_freq_h_sizes_ptr[ e ] = &P.table_freq_h_sizes[ e * P.MAX_W_SIZE ];
    P.vector_mean_freq_h_sizes = ( double* ) calloc( P.MAX_W_SIZE, sizeof( double ) );
    P.vector_sd_freq_h_sizes = ( double* ) calloc( P.MAX_W_SIZE, sizeof( double ) );
    
    P.gen_cum_large = ( int* ) calloc( P.MAX_GEN + 1, sizeof( int ) );

    P.n_sync_large = ( int* ) calloc( P.TOT_STEPS + 1, sizeof( int ) );
//    P.n_school_closure = ( int* ) calloc( P.TOT_STEPS + 1, sizeof( int ) );
    
    P.cum_r = ( int* ) calloc( P.TOT_STEPS + 1, sizeof( int ) );
    
    P.r_100 = ( int* ) calloc( P.HOW_MANY, sizeof( int ) );
    P.r_10ghosts = ( int* ) calloc( P.HOW_MANY, sizeof( int ) );
    
    return;
}

// Load the network from 2 files called H_GB_5.dat and W_GB_5.dat for network named GB and a population size of 10^5
void load_network( void ) /* This is the main function of this file */
{
    FILE *H_ptr, *W_ptr;
    char name_array[ 1024 ] = "", *name = name_array;
    int checkC = 0;
    int u, v, size, sizeA, sizeC, z, i, j = 0;
    
    /* Individuals */
    for ( i = 0; i < P.TOT; i++ ) {
        ind[ i ].H_index = -1;
        ind[ i ].W_index = -1;
    }
    
    name = create_file_name( name, "H_", P.SOCIAL_STRUCTURE_NETWORK_TYPE );
    if ( !( H_ptr = fopen( name, "r" ) ) )
    {
        printf("\nUnable to open the network file: %s", name);
    }
    fscanf( H_ptr, "%d%d%d", &P.TOT_C, &P.n_H, &P.MAX_H_SIZE );
    printf( "Total number of children: TOT_C = %d\n", P.TOT_C );
    H_list = ( HOUSEHOLD* ) calloc( P.n_H, sizeof( HOUSEHOLD ) );
    H_people = ( int* ) calloc( P.n_H*P.MAX_H_SIZE, sizeof( int ) );
    for ( u = 0; u < P.n_H; u++ )
    {
        fscanf( H_ptr, "%d%d", &v, &size );
        assert( u == v );
        H_list[ u ].index = v;
        H_list[ u ].size = size;
        //printf( "%d\t%d\t", H_list[ u ].index, H_list[ u ].size );
        H_list[ u ].people = &H_people[ u*P.MAX_H_SIZE ];
        for ( z = 0; z < size; z++ )
        {
            fscanf( H_ptr, "%d", &i);
            H_people[ u*P.MAX_H_SIZE+z ] = i;
            ind[ i ].H_index = u;
            j++; // Just to count how many times I store an individual, so I can check in the end
            //printf( "%d\t", H_list[ u ].people[ z ] );
        }
        //assert( fscanf( H_ptr, "\n" ) );
        for ( ; z < P.MAX_H_SIZE; z++ )
        {
            H_people[ u*P.MAX_H_SIZE+z ] = -1;
            //printf( "%d\t", H_list[ u ].people[ z ] );
        }
        //printf( "\n" );
    }
    assert( P.TOT_C == j );
    j = 0;
    //printf( "\n" );
    
    name = create_file_name( name, "W_", P.SOCIAL_STRUCTURE_NETWORK_TYPE );
    if ( !( W_ptr = fopen( name, "r" ) ) )
    {
        printf("\nUnable to open the network file: %s", name);
        exit(1);
    }
    W_ptr = fopen( name, "r" );
    fscanf( W_ptr, "%d%d%d%d%d", &P.TOT_A, &checkC, &P.n_W, &P.MAX_WA_SIZE, &P.MAX_WC_SIZE );
    if ( ( P.TOT_C != checkC ) || ( P.TOT != P.TOT_A + P.TOT_C ) )
    {
        printf( "Error in loading the populations size and composition!" );
        exit( 1 );
    }
    else
    {
        printf( "The total population size has been checked correctly: TOT_A = %d\n", P.TOT_A );
        printf( "The total population size has been checked correctly: TOT_C = %d\n", P.TOT_C );
        W_list = ( WORKPLACE* ) calloc( P.n_W, sizeof( WORKPLACE ) );
        P.MAX_W_SIZE = P.MAX_WA_SIZE + P.MAX_WC_SIZE;
        W_people = ( int* ) calloc( P.n_W*P.MAX_W_SIZE, sizeof( int ) );
        for ( u = 0; u < P.n_W; u++ )
        {
            fscanf( W_ptr, "%d%d%d", &v, &sizeA, &sizeC );
            assert( u == v );
            W_list[ u ].index = v;
            W_list[ u ].sizeA = sizeA;
            W_list[ u ].sizeC = sizeC;
            //printf( "%d\t%d\t", W_list[ u ].index, W_list[ u ].size );
            W_list[ u ].peopleA = &W_people[ u*P.MAX_W_SIZE ];
            W_list[ u ].peopleC = &W_people[ u*P.MAX_W_SIZE + P.MAX_WA_SIZE ];
            for ( z = 0; z < P.MAX_WA_SIZE; z++ )
            {
                fscanf( W_ptr, "%d", &i);
                W_people[ u*P.MAX_W_SIZE+z ] = i;
                if ( i != -1 ) {
                    //ind[ i ].H_index = -1; // Just to place each adult in a different workplace of size 1
                    ind[ i ].W_index = u;
                    j++; // Just to count how many times I store an individual, so I can check in the end
                    //printf( "%d\t", W_list[ u ].peopleA[ z ] );
                }
            }
            for ( z = 0; z < P.MAX_WC_SIZE; z++ )
            {
                fscanf( W_ptr, "%d", &i);
                W_people[ u*P.MAX_W_SIZE+P.MAX_WA_SIZE+z ] = i;
                if ( i != -1 ) {
                    ind[ i ].W_index = u;
                    j++; // Just to count how many times I store an individual, so I can check in the end
                    //printf( "%d\t", W_list[ u ].peopleC[ z ] );
                }
            }
        }
        assert( P.TOT == j );
    }
    return;
}

char* create_file_name( char *name, char *beg, char* type )
{
    char suf_array[ 1024 ] = "", s_esp_array[ 1024 ] = "";
    char *suf = suf_array, *s_esp = s_esp_array;
    int n_esp, r;
    
    n_esp = my_round( log( ( double ) P.TOT ) / log( 10.0 ) );
    r = sprintf( s_esp, "%d", n_esp );
    suf = strcpy( suf, type );
    printf("\n%s\n", suf );
    suf = strcat( suf, "_" );
    //	printf("\n%s\n", suf );
    suf = strcat( suf, s_esp );
    //	printf("\n%s\n", suf );
    suf = strcat( suf, ".dat" );
    //	printf("\n%s\n", suf );
    
    name = strcpy( name, beg );
    //	printf("\n%s\n", name );
    name = strcat( name, suf );
    //	printf("\n%s\n", name );
    return name;
}

int my_round( double x )
{
    int c;
    
    c = ( int ) x;
    
    return ( x - c ) < 0.5 ? c : c + 1;
}

void initialyse_basic_quantities( void )
{
    int i, l, e, z, g, s, fs;
    
    /* Individuals */
    for ( i = 0; i < P.TOT_C; i++ )
        ind[ i ].age = 'C';
    for ( ; i < P.TOT; i++ )
        ind[ i ].age = 'A';
    
    /* Helper stuff */
    for ( l = 0; l <= P.TOT_h_STEPS; l++ ) {
        h[ l ].p = NULL;
        h[ l ].t = l*P.dh;
    }
    
    /* Event list (only the flag to signal not to free it) */
    for ( e = 0; e < P.MAX_EV; e++ )
        ev[ e ].heap = 'N';
    
    /* Household stuff */
    // In HW_random_network
    
    /* Epidemic stuff */
    P.n_large = 0;
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        P.extinction_flag[ z ] = 'U'; // Unknown...
        P.final_size[ z ] = 0;
        P.final_size_H[ z ] = 0;
        P.final_size_W[ z ] = 0;
        P.inc_peak[ z ] = 0;
        P.t_index_peak[ z ] = 0;
        P.t_last_event[ z ] = 0;
//        P.t_school_closure[ z ] = 0;
//        P.t_school_reopen[ z ] = 0;
    }
    P.extinction_flag[ P.HOW_MANY ] = '\0';
    
//    /* School closure stuff */
//    P.n_schools_flagged = 0;
    
    /* Generation stuff */
    for ( g = 0; g < P.MAX_GEN; g++ ) {
        P.gen_cum_large[ g ] = 0;
        gen_large[ g ].number = g+1;
        gen_large[ g ].R0 = 0;
        gen_large[ g ].Rt = 0;
        gen_large[ g ].colls = 0;
        gen_large[ g ].A = 0;
        gen_large[ g ].C = 0;
        gen_large[ g ].I = 0;
        gen_large[ g ].TransmToC = 0;
        gen_large[ g ].TransmToA = 0;
        gen_large[ g ].propC = 0;
        gen_large[ g ].propA = 0;
        gen_large[ g ].gen_cum_conts_large = 0;
        gen_large[ g ].gen_cum_infs_large = 0;
    }
    P.gen_cum_large[ P.MAX_GEN ] = 0;
    
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        for ( s = 0; s <= P.TOT_STEPS; s++ ) {
            sys[ z ][ s ].S = 0;
            sys[ z ][ s ].I = 0;
            sys[ z ][ s ].R = 0;
            sys[ z ][ s ].inc = 0;
            sys[ z ][ s ].inc_C = 0;
            sys[ z ][ s ].inc_A = 0;
            sys[ z ][ s ].cum_inc = 0;
            sys[ z ][ s ].cum_inc_C = 0;
            sys[ z ][ s ].cum_inc_A = 0;
            sys[ z ][ s ].n_conts = 0;
            sys[ z ][ s ].n_infs = 0;
            sys[ z ][ s ].G_infs = 0;
            sys[ z ][ s ].H_infs = 0;
            sys[ z ][ s ].W_infs = 0;
            sys[ z ][ s ].colls = 0;
            sys[ z ][ s ].H_S = 0;
            sys[ z ][ s ].H_I = 0;
            sys[ z ][ s ].H_R = 0;
            sys[ z ][ s ].H_inc = 0;
            sys[ z ][ s ].H_cum_inc = 0;
            sys[ z ][ s ].W_cum_inc = 0;
        }
    }
    for ( s = 0; s <= P.TOT_STEPS; s++ ) {
        P.n_sync_large[ s ] = 0;
//        P.n_school_closure[ s ] = 0;
        sync_large[ s ].t = s * P.dt;
        sync_large[ s ].S = 0;
        sync_large[ s ].I = 0;
        sync_large[ s ].R = 0;
        sync_large[ s ].inc = 0;
        sync_large[ s ].inc_C = 0;
        sync_large[ s ].inc_A = 0;
        sync_large[ s ].cum_inc = 0;
        sync_large[ s ].cum_inc_C = 0;
        sync_large[ s ].cum_inc_A = 0;
        sync_large[ s ].r = 0;
        sync_large[ s ].R0 = 0;
        sync_large[ s ].Rt = 0;
        sync_large[ s ].R_in_G = 0;
        sync_large[ s ].R_in_H = 0;
        sync_large[ s ].R_in_W = 0;
        sync_large[ s ].colls = 0;
        sync_large[ s ].H_S = 0;
        sync_large[ s ].H_I = 0;
        sync_large[ s ].H_R = 0;
        sync_large[ s ].H_inc = 0;
        sync_large[ s ].H_cum_inc = 0;
        sync_large[ s ].H_r = 0;
        sync_large[ s ].W_cum_inc = 0;
        sync_large[ s ].cum_inc_gap = 0;
        sync_large[ s ].R0gap = 0;
        sync_large[ s ].Rtgap = 0;
    }
    
    /* r stuff */
    P.n_r = 0;
    P.r_length = P.TOT_STEPS;
    for ( s = 0; s <= P.TOT_STEPS; s++ )
        P.cum_r[ s ] = 0; // Can be put in the previous cycle...
    for ( z = 0; z < P.HOW_MANY; z++ ) {
        P.r_100[ z ] = 0; // Can be put in the previous cycle...
        P.r_10ghosts[ z ] = 0; // Can be put in the previous cycle...
    }
    
    for ( s = 0; s < P.MAX_W_SIZE; s++ )
    {
        P.vector_mean_freq_h_sizes[ s ] = 0;
        P.vector_sd_freq_h_sizes[ s ] = 0;
        for ( fs = 0; fs <= s; fs++ )
        {
            P.table_mean_fs_ptr[ s ][ fs ] = 0;
            P.table_sd_fs_ptr[ s ][ fs ] = 0;
            for ( z = 0; z < P.HOW_MANY; z++ )
                P.table_h_fs_ptr[ z ][ s ][ fs ] = 0;
        }
        for ( z = 0; z < P.HOW_MANY; z++ )
            P.table_freq_h_sizes_ptr[ z ][ s ] = 0;
    }
    
}


