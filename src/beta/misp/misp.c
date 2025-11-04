#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "kseq.h"
#include "motif.h"
KSEQ_INIT(gzFile, gzread, gzrewind)
#include "debug.h"

#define BITSHIFT 2
#define LOGODD_MIN -256
#define KINDS 4			/* A, T, C, G */
float cg_percent(kseq_t* kseq);
void expected_difference(struct pssm_matrix *pm, float c_p, double* ediff);
void lookahead_filter(int q, kseq_t *kseq, struct pssm_matrix *pm, float c_p, double tol, struct match_doublet *md, FILE *output, int mini);

const char *byte_to_binary(int x)
{
	static char b[11];
	b[0] = '\0';

	int z;
	for (z = 512; z > 0; z >>= 1)
	{
		strcat(b, ((x & z) == z) ? "1" : "0");
	}

	return b;
}

uint_fast32_t flip_reverse2(uint_fast32_t value)
{
	int i;
	uint32_t new_value = 0;
	for (i = 0; i < 16; i++)
	{
		new_value <<= 2;
		new_value |= (value & 0x3);
		value >>= 2;
	}
	return ((~new_value)>>22) & 0x3ff; /* only ok when q is 5 (22 == 32 - 2*5)*/
}



int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	FILE* pssm_fp;
	FILE* output;
	struct pssm_matrix *pssm;
	char of_name[2000];
	struct match_doublet *dm;
	float p_cg;
	float p_value;

	pssm = malloc(sizeof(struct pssm_matrix));
	dm = malloc(sizeof(struct match_doublet));

	if (argc < 5) {
		fprintf(stderr, "Usage: %s <in.seq> <in.db> <p-value> <motif-id> <output-prefix>\n", argv[0]);
		return 1;
	}


	pssm_fp = fopen(argv[2],"r");
	pssm_reader(pssm_fp, pssm);
	fclose(pssm_fp);

	fp = gzopen(argv[1], "r");
	p_value = atof(argv[3]);

	seq = kseq_init(fp);
	p_cg=cg_percent(seq);
	pssm2logodd(pssm, p_cg);

	if (strcasecmp("all", argv[4]) == 0) {
		sprintf(of_name,"%s_all",argv[5]);
		output = fopen(of_name,"w");
		fprintf(output, "# Parameter List:\n");
		fprintf(output, "# Input sequence: %s\n",argv[1]);
		fprintf(output, "# Output path: %s\n",argv[5]);
		fprintf(output, "# P value: %.3f\n",p_value);
		for(; pssm->next != NULL; pssm = pssm->next) {
			kseq_rewind(seq);
			lookahead_filter(5, seq, pssm, p_cg/2.0, threshold_fromP(pssm, p_cg/2.0, p_value), dm, output,1);
		}
		fclose(output);

	} else {
		for(; pssm->next != NULL; pssm = pssm->next) {
			if (strncasecmp(pssm->name, argv[4], strlen(pssm->name) - 1) == 0)
				break;
		}

		if (strncasecmp(pssm->name, argv[4], strlen(pssm->name) -1) != 0) {
			printf("Can't find motif id %s", argv[4]);
			return 0;
		}
		kseq_rewind(seq);
		pssm->name[strlen(pssm->name) - 1] = '\0';
		sprintf(of_name,"%s_%s",argv[5],pssm->name);
		output = fopen(of_name,"w");
		fprintf(output, "# Parameter List:\n");
		fprintf(output, "# Input sequence: %s\n",argv[1]);
		fprintf(output, "# Output path: %s\n",argv[5]);
		fprintf(output, "# P value: %.3f\n",p_value);
		lookahead_filter(5, seq, pssm, p_cg/2.0, threshold_fromP(pssm, p_cg/2.0, p_value), dm, output, 0);
		fclose(output);
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
float cg_percent(kseq_t* kseq)
{
	int l;
	unsigned long a, c, other;
	a = 0;
	c = 0;
	other = 0;
	while ((l = kseq_read(kseq)) >= 0)
		bg_counter(&a, &c, &other, kseq->seq.s);
	printf("Calculating Background..\n");
	printf("A or T: %ld\n", a);
	printf("C or G: %ld\n", c);
	printf("Others: %ld\n", other);
	if (a+c<50) {
		printf("Your input sequence has too few lines, the bg score will be set to 0.5\n");
		return 0.5;
	} else  return ((double)c)/(a+c);
}


void lookahead_filter(int q, kseq_t *kseq, struct pssm_matrix *pm, float c_p, double tol, struct match_doublet* md, FILE *output, int mini)
{
	int i, j;
	int window_pos;
	struct order_s **order;

	double *good;
	double *scores;
	int pos_max;
	double hit_max;
	double tmp, tmp_max;
	int bufsize;
	short *seq_code;
	short *seq_code_gc;
	uint_fast32_t *sA;
	uint_fast32_t code;
	const uint_fast32_t size = 1 << (BITSHIFT * q);;

	int compare (const void * a, const void * b);

	assert(pm->kinds == KINDS);	/* Only for DNA */
	/* Find the window */
	if (q > pm->len) {
		return lookahead_filter(pm->len, kseq, pm, c_p, tol, md, output,mini);
	}
	sA = calloc(q, sizeof(uint_fast32_t));
	good = malloc(sizeof(double) * pm->len);
	tmp = 0;
	expected_difference(pm, c_p, good);
	for (i = 0; i < q; ++i)
		tmp += good[i];
	tmp_max = tmp;
	window_pos = 0;

	for (i = 0; i < pm->len - q; ++i) {
		tmp -= good[i];
		tmp += good[i+q];
		if (tmp > tmp_max) {
			tmp_max = tmp;
			window_pos = i + 1;
		}
	}

	/* printf("%f\n", tmp); */
	/* Arrange matrix indeces remained */
	/* NOTICE: start from 1 !!*/
	order = (struct order_s **) calloc(pm->len - q, sizeof(struct order_s*));

	for (i = 0; i < window_pos; ++i) {
		order[i] = (struct order_s *) malloc(sizeof(struct order_s));
		order[i]->pos = i;
		order[i]->good = good[i];
	}
	for (i = window_pos + q; i < pm->len; ++i)
	{
		order[i-q] = (struct order_s *) malloc(sizeof(struct order_s));
		order[i-q]->pos = i;
		order[i-q]->good = good[i];
	}
	qsort(order, pm->len - q, sizeof(struct order_s *), compare);
	/* Lookahead array for indeces outside the window */
	good = malloc(sizeof(double) * (pm->len - q));
	for (j = pm->len - q; j > 0; --j) {
		tmp_max = LOGODD_MIN;
		for (i = 0; i < pm->kinds; ++i) {
			if (tmp_max < pm->score[i][order[j-1]->pos])
				tmp_max = pm->score[i][order[j-1]->pos];
		}
		good[j-1] = good[j] + tmp_max;
	}



	/* Scores in the window for all possible strings */

	scores = malloc(sizeof(double) * (size));
	tmp = 0;
	for (j = 0; j < q; ++j) {
		tmp = tmp + pm->score[0][window_pos+j];

	}
	scores[0] = tmp;
	/* printf("window_pos is %d\n", window_pos); */
	/* printf("q is %d\n", q); */
	/* printf("score zero is %lf\n", scores[0]); */
	/* printf("%f\n",scores[0]); */
	j = 0;
	code = 0;
	while(1) {
		++code;
		if (code >= size)
			break;

		scores[code] = scores[code - 1];
		j = q - 1;	/* j is pos now */
		while (j >= 0) {
			if (sA[j] <  KINDS - 1) {
				scores[code] +=				\
					pm->score[sA[j] + 1][j + window_pos] \
					-				\
					pm->score[sA[j]][j + window_pos];
				++sA[j];
				break;
			} else {
				scores[code] +=				\
					pm->score[0][j + window_pos]	\
					-				\
					pm->score[KINDS - 1][j + window_pos];
				sA[j] = 0;
				--j;
			}

		}
		/* printf("%d's(%s) score is %f\n", code, byte_to_binary(code),scores[code]); */
	}

	int jb;
	float jb_zuiad = -1000;
	for(jb=0; jb<size; jb++) {
		if (scores[jb] > jb_zuiad) {
			jb_zuiad = scores[jb];
		}
		/* printf("score is %f for code %s\n\n",scores[jb], byte_to_binary(jb)); */
	}
	/* Debug: longest score is jb_zuiad */
	/* Actual scanning */
	tmp = 0;
	tmp_max = tol - good[0];


	if (mini == 0) {
		fprintf(output, "# CG percent: %.2f\n", c_p*2);
		fprintf(output, "# tolerance: %.2f\n", tol);
		fprintf(output, "# factor ID: %s\n", pm->name);
		fprintf(output, "sequence_name\tsequence_length\thits_score\thits_position\tsequence\n");
	} else {
		fprintf(output, "\n\n@ factor:%s# tolerance %.2f\n", pm->name, tol);
	}

	while ((bufsize = kseq_read(kseq)) >= 0) {
		if (mini == 0)
			fprintf(output,"%s\t%d\t", kseq->name.s, bufsize);

		if (bufsize < pm->len) {
			if (mini == 0)
				fprintf(output, "*\t*\t*\n");
			continue;
		}

		seq_code = malloc(sizeof(short) * (bufsize+1));
		seq_code_gc = seq_code;

		base2code(kseq->seq.s, seq_code);

		int index_N=0;

		hit_max = -10000;
		pos_max = -1;

		int i_coming_in = 0;
		int order_coming_in = 0;
		/* for(j=0; j< pm->len; j++) */
		/* 	for(jb = 0; jb< 4; jb++) */
				/* printf("[%d,%d] <= %f\n", j, jb,pm->score[jb][j]); */

		for(i = 0; i + pm->len <= bufsize; ++i) {
			/* Detection in Positive Strand */

			if(i == 0) { /* Init the code */
				code = 0;
				for (j = i+ window_pos; j < i + window_pos + q - 1; ++j)  {
					if (seq_code[j] == 4) index_N = j;
					code = (code << BITSHIFT) + seq_code[j];
				}
				if(index_N) continue;
			}

			i_coming_in = i + window_pos + q - 1;
			code = ((code << BITSHIFT) + seq_code[i_coming_in]) & (size - 1);


			if (index_N) {
				index_N --;
				continue;
			}

			if (seq_code[i_coming_in] == 4) {
				index_N = i_coming_in - i;
				continue;
			}
			/* code stores the seq info of [i+window_pos, i + window_pos + q -1]  */

			tmp = scores[code];

			if (tmp < tmp_max) /* tmp_max is the threadhood */
				continue;
			for (j = 0; j < pm->len - q; ++j) {
				order_coming_in = i + order[j]->pos;

				if (tmp + good[j] < tol) break;
				if (seq_code[order_coming_in] == 4) {
					tmp = tol -1;
					break;
				}
				/* printf("order coming %d <- %d <=  %f\n", order_coming_in, seq_code[order_coming_in], pm->score [seq_code[order_coming_in]][order[j]->pos] /\* position in the pssm *\/); */
				tmp += pm->score [seq_code[order_coming_in]] /* nucleotide */
					[order[j]->pos]; /* position in the pssm */
			}

			if (tmp > hit_max) {
				hit_max = tmp;
				pos_max = i;
			}
		}
		for(i = bufsize - 1; i + 1 - pm->len >= 0; --i){
			/* Detection in Reverse Strand */

			if(i == bufsize - 1) {
				index_N = 0;
				code = 0;
				for(j = i - window_pos; j > i - window_pos -q + 1; --j) {
					if (seq_code[j] == 4) index_N = -j;
					code = (code << BITSHIFT) + (3 - seq_code[j]);
				}
				if(index_N) continue;
			}

			i_coming_in = i - window_pos - q + 1;
			code = ((code << BITSHIFT) + (3 - seq_code[i_coming_in])) & (size - 1);
			/* printf("code is reverse <= %s\n", byte_to_binary(code)); */
			if (index_N) {
				index_N ++;
				continue;
			}

			if (seq_code[i_coming_in] == 4) {
				index_N = i_coming_in - i;
				continue;
			}

			tmp = scores[code];
			if (tmp < tmp_max) continue;
			/* printf("WOcao Neg %s!!!!\n", kseq->seq.s); */
			for (j = 0; j < pm->len - q; ++j) {
				order_coming_in = i - order[j]->pos;
				if (tmp + good[j] < tol) break;
				if (seq_code[order_coming_in] == 4) {
					tmp = tol -1;
					break;
				}
				tmp += pm->score [3 - seq_code[order_coming_in]] /* nucleotide */
					[order[j]->pos]; /* position in the pssm */
			}
			if (tmp > hit_max) {
				hit_max = tmp;
				pos_max = -i;
			}

		}
		/* printf("%f %f", tmp ,hit_max); */
		free(seq_code_gc);
		if (mini == 0) {
			if (hit_max > tol){
				char temp_char;
				int e_idx, s_idx; /* index of end and start position */
				if (pos_max>=0)
					e_idx = pos_max + pm->len;
				else
					e_idx = - pos_max + 1;
				s_idx = e_idx - pm->len;
				temp_char = *(kseq->seq.s + e_idx);
				*(kseq->seq.s + e_idx) = '\0';
				if (pos_max>=0)
					fprintf(output, "\t%.2f\t%d\t%s\n",hit_max - tol, pos_max, kseq->seq.s + s_idx);
				else
					fprintf(output, "\t%.2f\t%d(-)\t%s\n",hit_max - tol, -pos_max , kseq->seq.s + s_idx);
				*(kseq->seq.s + e_idx) = temp_char;
			}
			else
				fprintf(output, "*\t0\t*\n");
		} else {
			if (hit_max > tol)
				fprintf(output, "%.2f(%d),", hit_max - tol, pos_max);
			else
				fprintf(output, "0,");
		}

	}
}



int compare (const void * a, const void * b)
{
	if (((*(const struct order_s **)a)->good - (*(const struct order_s **)b)->good) > 0)
		return 1;
	else
		return -1;
}

void expected_difference(struct pssm_matrix *pm, float c_p, double* ediff)
/* Please allocate enough memory for ediff before call me! */
{
	int i,j;
	float max;

	for (i = 0; i < pm->len; ++i) {
		max = LOGODD_MIN;
		for (j = 0; j < pm->kinds; ++j) {
			if (max < pm->score[j][i])
				max = pm->score[j][i];
		}

		ediff[i] = max;

		for (j = 0; j < pm->kinds; ++j) {
			if (j == PWM_BASE_C || j == PWM_BASE_G)
				ediff[i] -= c_p * pm->score[j][i];
			else if (j == PWM_BASE_A || j == PWM_BASE_T)
				ediff[i] -= (0.5 - c_p) * pm->score[j][i];
			else {
				perror("Error");
				exit(1);
			}
		}
	}
}
