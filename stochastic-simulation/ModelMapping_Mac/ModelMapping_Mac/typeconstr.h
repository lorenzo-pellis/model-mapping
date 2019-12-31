//
//  typeconstr.h
//  ModelMapping_Mac
//
//  Created by Lorenzo Pellis on 26/05/2019.
//  Copyright © 2019 University of Manchester
//              2017 University of Warwick
//              2012 Imperial College London
//

/////////////////////////////* Type construction *///////////////////////////////////

PAR P;

IND *ind;	/* from 0 to TOT-1 */

EVENT *ev;

HELPER *h;
HELPER *ptr_h;

HOUSEHOLD* H_list;
int *H_people;

WORKPLACE* W_list;
int *W_people;

GEN *gen_large;

CHANGES *ch;

SYS *sys_matrix;	/* A row for each run */
SYS **sys_line;
SYS **sys;

SUM *sync_large; // only sync and large //, global[ TOT_STEPS+1 ], large[ TOT_STEPS+1 ], sync_global[ TOT_STEPS+1 ];

/////////////////////////// Variables declarations /////////////////////////

double Rg;
double Rh;
double Rw;
double NGM_G_AA;
double NGM_G_AC;
double NGM_G_CA;
double NGM_G_CC;

