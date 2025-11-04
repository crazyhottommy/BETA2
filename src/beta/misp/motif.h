#ifndef MIS_MOTIF_H
#define MIS_MOTIF_H
#define PWM_MAX_COL 50
#define PVAL_DP_MULTIPLIER 1000.0
/* #define PWM_BASE_A 0 */
/* #define PWM_BASE_C 1 */
/* #define PWM_BASE_G 2 */
/* #define PWM_BASE_T 3 */

#define PWM_BASE_T 0
#define PWM_BASE_C 1
#define PWM_BASE_G 2
#define PWM_BASE_A 3

#define PWM_BASE_N 4
#define PSUEDO_E 0.1
#include "kseq.h"
#include <zlib.h>
#include <stdio.h>



struct pwm_matrix 
{
	int **weight;
	/* weight[kinds][position] */
	int len;
	short kinds;
	/* char *comment; */
};

struct pssm_matrix 
{
	double **score;
	/* weight[kinds][position] */
	int len;
	short kinds;
	char *name;
	/* char *comment; */
	struct pssm_matrix *next;
};

struct match_doublet 
{
	int position;
	float score;
	struct match_doublet *next;
};

struct order_s {
  int pos;
  double good;
};

  
int bg_counter(unsigned long *acnt, unsigned long *ccnt, unsigned long *other, char *seq);

int pwm_reader(FILE *fp, struct pwm_matrix *pm);

void pssm_reader(FILE *fp, struct pssm_matrix *pm);
	
void pssm2logodd(struct pssm_matrix *pm, double c_p);
	
void display_pwm(struct pwm_matrix *pm);

void display_pssm(struct pssm_matrix *pm);

/* Calculates a threshold for a scoring matrix from a given p value */
double threshold_fromP (struct pssm_matrix *pm, double c_p, double p);

int init_array(double *array, size_t size, double value);

int counts2Logfodds(struct pwm_matrix *pm_in, struct pssm_matrix *pm_out,  double c_p, double ps);



struct pssm_matrix *use_pssm(struct pssm_matrix *pm, char* name);

void base2code(char *seq, short *code);
#endif
