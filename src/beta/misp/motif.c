#define _GNU_SOURCE
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "motif.h"

#define PSUDO_SMALL 0.0001

int bg_counter(unsigned long *acnt, unsigned long *ccnt, unsigned long *other, char *seq)
{
	int i;
	for (i=0; seq[i] != 0; ++i) {
		switch (seq[i]) {
		case 'A': case 'T': case 'a': case 't':
			(*acnt)++;
			break;
		case 'C': case 'G': case 'c': case 'g':
			(*ccnt)++;
			break;
		default:
			(*other)++;
		}
	}
	return 1;
}



void pssm_reader(FILE *fp, struct pssm_matrix *pm)
/* allocate the memory for pm before use this function */
{
	char * read;
	int i;
	int j;
	int k;
	int first;

	int len = 1000;
	char str[1000];
	char *pch = NULL;
	if (fp == NULL){
		perror("Error opening the pwm file");
		exit(-1);
	}

	first = 1;
	pm->kinds = 4;
	pm->len = 0;
	pm->score = (double **)malloc(sizeof(double *) * pm->kinds);
	for (i = 0; i < 4; i++)
		pm->score[i] = (double *)malloc(sizeof(double) * PWM_MAX_COL);


	read = fgets (str, len, fp);
	i = -1;
	j = 0;
	while (read != NULL && i <= pm->kinds) {
		assert(j < PWM_MAX_COL);
		if (j == 0) {
			pch = strtok(str, " ");
		} else {
			pch = strtok(NULL, " ");
		}

		if (pch == NULL) {
			if (i == -1) {} /* before reading motif */
			else if (i ==0 ) {
				pm->len = j;
				j = 0;
				i++;
			}
			else if (i == pm->kinds -1) {
				i = -1;
				j = 0;
			} else {
				assert(pm->len == j);
				j = 0;
				i++;
			}

			read = fgets (str, len, fp);
			continue;
		}
		if (pch != NULL) {
			if (i == -1) {
				if (strncmp(pch,"\n", 1) != 0) { /* If it's a '\n', skip malloc new space until new line read */
					if (!first) {
						pm->next = malloc(sizeof(struct pssm_matrix)); /* init next motif */
						pm = pm->next;
						pm->kinds = 4;
						pm->len = 0;
						pm->score = (double **)malloc(sizeof(double *) * pm->kinds);
					}

					for (k = 0; k < 4; k++)
						pm->score[k] = (double *)malloc(sizeof(double) * PWM_MAX_COL);

					pm->name = (char *)malloc(strlen(pch) + 1);
					strcpy(pm->name, pch);
					i = 0;
					first = 0;
				}
				read = fgets (str, len, fp);
			} else {
				pm->score[i][j++] = atof(pch);
			}
		}
	}
	pm->next=NULL;

}

void pssm2logodd(struct pssm_matrix *pm, double c_p)
{
	int i;
	int j;
	for(; pm->next != NULL; pm = pm->next) {
		for (i = 0; i < pm->kinds; ++i) {
			for (j = 0; j < pm->len; ++j) {
				if (i==PWM_BASE_G || i==PWM_BASE_C)
					pm->score[i][j] = log(2 * pm->score[i][j] / c_p + PSUDO_SMALL);
				else
					pm->score[i][j] = log(2* pm->score[i][j] / (1 - c_p) + PSUDO_SMALL);
			}
		}
	}
}

struct pssm_matrix *use_pssm(struct pssm_matrix *pm, char* name)
{
	for(; pm->next != NULL; pm = pm->next) {
		if (strcasecmp(pm->name, name) == 0)
			return pm;
	}
	printf("Can't find the motif id: %s\n",name);
	exit(1);
}




void display_pwm(struct pwm_matrix *pm)
{
	int i;
	int j;
	for (i = 0; i < pm->kinds; ++i) {
		for (j = 0; j < pm->len; ++j) {
			printf("%d\t", pm->weight[i][j]);
		}
		putchar('\n');
	}
}

void display_pssm(struct pssm_matrix *pm)
{
	int i;
	int j;
	for(; pm->next != NULL; pm = pm->next) {
		printf("%s\n",pm->name);
		for (i = 0; i < pm->kinds; ++i) {
			for (j = 0; j < pm->len; ++j) {
				printf("%.3f\t", pm->score[i][j]);
			}
			putchar('\n');
		}
	}
}



double threshold_fromP (struct pssm_matrix *pm, double c_p, double p)
/* c_p means ``half'' of the CG percent */
{
	int i, j;
	int max, min;
	int maxT, minV;
	assert(pm->kinds > 0 && pm->len > 0);

	int **mat;
	mat = (int **)malloc(sizeof(int *) * pm->kinds);
	for (i = 0; i < 4; i++)
		mat[i] = (int *)malloc(sizeof(int) * pm->len);


	int tmp;
	/* printf("pm->len(thred) is %d\n", pm->len);	 */

	/* init */
	for (i = 0; i < pm->len; ++i) {
		for (j = 0; j < pm->kinds; ++j) {
			if (pm->score[j][i] > 0.0)
				mat[j][i] = (int)(PVAL_DP_MULTIPLIER * pm->score[j][i] + 0.5);
			else
				mat[j][i] = (int)(PVAL_DP_MULTIPLIER * pm->score[j][i] - 0.5);

			/* printf("mat[%d,%d] is %d\n", j, i, mat[j][i]); */
		}
	}


	maxT = 0;
	minV = INT_MAX;

	for (i = 0; i < pm->len; ++i) {
		max = mat[0][i];
		min = max;
		for (j = 1; j < pm->kinds; ++j) {
			tmp = mat[j][i];
			if (max < tmp)
				max = tmp;
			else if (min > tmp)
				min = tmp;
		}
		maxT += max;
		if (minV > min)
			minV = min;
	}

	tmp = maxT - pm->len * minV;
	/* printf("tmp is %d\n" , tmp); */
	double table0[tmp+1];
	double table1[tmp+1];
	int r,s;
	double sum;
	init_array(table0, tmp+1, 0.0);
	init_array(table1, tmp+1, 0.0);

	for (j = 0; j < pm->kinds; ++j) {
		if (j == PWM_BASE_C || j == PWM_BASE_G) {
			table0[mat[j][0] - minV] += c_p;
		}
		else if (j == PWM_BASE_A || j == PWM_BASE_T) {
			table0[mat[j][0] - minV] += 0.5 - c_p;
		} else {
			perror("Error");
			exit(-1);
		}
	}


	for (i = 1; i < pm->len; ++i) {
		for (j = 0; j< pm->kinds; ++j) {
			s = mat[j][i] - minV;
			for (r = s; r <= tmp; ++r) {
				if (j == PWM_BASE_C || j == PWM_BASE_G) {
					table1[r] += c_p * table0[r-s];
				}
				else if (j == PWM_BASE_A || j == PWM_BASE_T) {
					table1[r] += (0.5 - c_p) * table0[r-s];
				} else {
					perror("Error");
					exit(-1);
				}
			}

		}

		for (r = 0; r <= tmp; ++r) {
			table0[r] = table1[r];
			table1[r] = 0.0;
		}
	}
	sum = 0.0;

	for (r = tmp; r >=0; --r) {
		sum += table0[r];
		/* printf("sum is %f r is %d\n", sum,r);	 */
		if (sum > p){
			return (double) ((r + (pm->len) * minV + 1) / PVAL_DP_MULTIPLIER);
		}
	}

	/* table0 and table1 are VLAs (stack-allocated), no need to free */
	return (double) (((pm->len) * minV) / PVAL_DP_MULTIPLIER);
}


int init_array(double *array, size_t size, double value)
{
	size_t i;
	for (i = 0; i < size; ++i)
		array[i] = value;
	return 0;
}

void base2code(char *seq, short *code) {
	srand(time(NULL));
	int i;
	for (i=0; seq[i] != 0; ++i) {
		switch (seq[i]) {
		case 'A': case 'a':
			code[i] = PWM_BASE_A;
			break;
		case 'C': case 'c':
			code[i] = PWM_BASE_C;
			break;
		case 'T': case 't':
			code[i] = PWM_BASE_T;
			break;
		case 'G': case 'g':
			code[i] = PWM_BASE_G;
			break;
		default:
			code[i] = PWM_BASE_N;
			break;
		}
	}
	code[i] = -1;
	/* the end mark */


}

int pwm_reader(FILE *fp, struct pwm_matrix *pm)
/* allocate the memory for pm before use this function */
{
	char * read;
	int i;
	int j;
	int len = 1000;
	char str[1000];
	char *pch = NULL;
	if (fp == NULL){
		perror("Error opening the pwm file");
		exit(-1);
	}
	pm->kinds = 4;
	pm->len = 0;
	pm->weight = (int **)malloc(sizeof(int *) * pm->kinds);
	for (i = 0; i < 4; i++)
		pm->weight[i] = (int *)malloc(sizeof(int) * PWM_MAX_COL);

	read = fgets (str, len, fp);
	/* read == -1: error or EOF*/
	i = 0;
	j = 0;

	while(read != NULL && i < 4) {
		assert(j < PWM_MAX_COL);

		if (j == 0)
			pch = strtok(str, " ");
		else
			pch = strtok(NULL, " ");
		if (j == 0 && pch == NULL) {
			perror("PSSM file error");
			exit(-1);
		}
		else if (j != 0 && pch == NULL) {
			pm->weight[i++][j] = -1;
			if (pm->len)
			{
				assert(pm->len == j);
			}
			/* check whether all line has same count of elements */
			else
				pm->len = j;
			j = 0;	/* at the end of a line */
			read = fgets (str, len, fp);
		} else {	/* pch has value */
			assert((pm->weight[i][j++] = atoi(pch)) >= 0);
		}

	}
	printf("pm->len is %d\n", pm->len);
	return 0;
}

int counts2Logfodds(struct pwm_matrix *pm_in, struct pssm_matrix *pm_out,  double c_p, double ps)
/* allocate memory for pm_in and pm_out before using */
{
	int i, j;
	int count;
	double freq;
	double a_p = 0.5 - c_p;
	assert(a_p > 0 && c_p > 0);

	pm_out->kinds = pm_in->kinds;
	pm_out->len = pm_in->len;

	pm_out->score = (double **)malloc(sizeof(double *) * pm_out->kinds);
	for (i = 0; i < 4; i++)
		pm_out->score[i] = (double *)malloc(sizeof(double) * pm_out->len);
	for (i = 0; i < pm_out->len; ++i) {
		count = 0;
		for (j = 0; j < pm_out->kinds; ++j)
			count += pm_in->weight[j][i];
		for (j = 0; j < pm_out->kinds; ++j) {
			if (j == PWM_BASE_C || j == PWM_BASE_G) {
				freq = (pm_in->weight[j][i] + c_p) / (count + ps);
				pm_out->score[j][i] = log(freq) - log(c_p);
			}
			else if (j == PWM_BASE_A || j == PWM_BASE_T) {
				freq = (pm_in->weight[j][i] + a_p) / (count + ps);
				pm_out->score[j][i] = log(freq) - log(a_p);
			} else {
				perror("Error");
				exit(-1);
			}
		}
	}

	return 0;
}
