#include <stdio.h>
#include <math.h>
#include "hydicedat.h"

#define SIZE 2768896
#define NUM_ROWS  64
#define NUM_COLS  64
#define BODUAN  169
#define WIN_OUT  11 // outer window
#define WIN_IN  3 // inner window
#define W1 WIN_OUT * WIN_OUT
#define T  (WIN_OUT/2)  ///////////
#define T1  (WIN_IN/2)    ///////////

#pragma DATA_SECTION(dat_float,".MY_MEM");
float dat_float[SIZE / 4];  /*float����*/
#pragma DATA_SECTION(X,".MY_MEM");
float X[BODUAN][NUM_ROWS * NUM_COLS];/*��ά����*/
//#pragma DATA_SECTION(rows_sum,".MY_MEM");
//float rows_sum[BODUAN]; /*�����ξ�ֵ*/
#pragma DATA_SECTION(sigma,".MY_MEM");
float sigma[BODUAN][BODUAN]; /*X��Э�������*/
#pragma DATA_SECTION(d2,".MY_MEM");
float d2[NUM_ROWS * NUM_COLS]; /*rx�㷨�ó�������*/
#pragma DATA_SECTION(inv_sigma,".MY_MEM");
float inv_sigma[BODUAN][BODUAN];/*sigma�������*/
//#pragma DATA_SECTION(datatest,".MY_MEM");
#pragma DATA_SECTION(data_1wei,".MY_PAD");
float data_1wei[NUM_ROWS * NUM_COLS * BODUAN];/*һάԭʼ����*/
#pragma DATA_SECTION(datatest,".MY_PAD");
float datatest[3 * NUM_ROWS][3 * NUM_COLS][BODUAN]; /*��άֱ�۽ṹ����*/
#pragma DATA_SECTION(buff_2wei,".MY_PAD");
float buff_2wei[BODUAN][WIN_OUT*WIN_OUT-WIN_IN*WIN_IN];/*��άȥ��NaN����*/


#pragma DATA_SECTION(buff,".MY_MEM");
float buff[WIN_OUT][WIN_OUT][BODUAN];/*buff����������ά*/
#pragma DATA_SECTION(buff_1wei,".MY_MEM");
float buff_1wei[WIN_OUT*WIN_OUT*BODUAN];/*buff������άתһά*/
#pragma DATA_SECTION(Y,".MY_MEM");
float Y[BODUAN];/*ԭʼ���ݵ�ĳһ����*/

#pragma DATA_SECTION(d1,".MY_MEM");
float d1[BODUAN];/*��rx����ʱ���м�����*/







int flag=0;
//fix();

void Read_IMG(const unsigned short *p);
void opera(void);
void padding(void);
void local_rx(void);
void corr(void);
void inv(void);
void swap(float*, float*);
void suanzi_rx(void);


int main(void)
{
	int i;
	Read_IMG(dat_hex_h);
	opera();
	padding();
    local_rx();

//	inv();
//	rx();
	for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
	printf("%f,", d2[i]);

	printf("end\n");
	return 0;
}

void Read_IMG(const unsigned short *p)
{
	int i, x;

	printf("���ڶ�ȡͼ��...\n");

	// ��ȡͼ������
	for ( i = 0; i < SIZE / 4; i++ )/////////////////
	{
		x = ((*p) << 24) | ((*(p + 1)) << 16) | ((*(p + 2)) << 8) | (*(p + 3));
		p += 4;
		dat_float[i] = *(float*)&x;
	}

	printf("ͼ���ȡ���!\n");
}

void opera()
{
	int i, j, flag1 = 0;

	/*�任�ɶ�ά����*/
	printf("����ת����ά����...\n");
	for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
		for (j = 0; j < BODUAN; j++)
		{
			X[j][i] = dat_float[flag1++];
		}
	printf("��ά����ת���ɹ�!\n");
}


void padding()
{
	int i, j, k;
	//float data[NUM_ROWS][NUM_COLS][BODUAN];

	//float * px;
	//float * pdatatest;
	//float * pdata;
	float * pdata1wei;


	printf("���ڴ���padding����...\n");
	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)
			for (j = 0; j < (NUM_ROWS * NUM_COLS); j++)
				*(pdata1wei++) = X[i][j];
	pdata1wei = data_1wei;
//
//	for (i = 0; i < BODUAN; i++)
//		for (j = 0; j < NUM_ROWS; j++)
//			for (k = 0; k < NUM_COLS; k++)
//				data[k][j][i] = *(px++);  //data_1weiһά����ά

	for (i = 0; i < BODUAN; i++)							//padding��������
			for (j = 0; j < (3 * NUM_COLS); j++)
				for (k = 0; k < (3 * NUM_ROWS); k++)
					datatest[j][k][i] = 0;

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//��
		for (j = NUM_COLS; j < (2 * NUM_COLS); j++)
			for (k = NUM_ROWS; k < (2 * NUM_ROWS); k++)
				datatest[k][j][i] = *(pdata1wei++);
	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//��
		for (k = NUM_ROWS; k < (2 * NUM_ROWS); k++)
			for (j = (NUM_COLS - 1); j > -1; j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//��
		for (k = NUM_ROWS; k < (2 * NUM_ROWS); k++)
			for (j = (3 * NUM_COLS - 1); j > (2 * NUM_COLS - 1); j--)
				datatest[j][k][i] = *(pdata1wei++);


	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//��
		for (k = NUM_ROWS - 1; k > -1; k--)
			for (j = NUM_COLS; j < (2 * NUM_COLS); j++)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//��
		for (k = (3 * NUM_ROWS - 1); k > (2 * NUM_ROWS - 1); k--)
			for (j = NUM_COLS; j < (2 * NUM_COLS); j++)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//����
		for (k = (NUM_ROWS - 1); k > -1; k--)
			for (j = (NUM_COLS - 1); j > -1; j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//����
		for (k = (NUM_ROWS - 1); k > -1; k--)
			for (j = (3 * NUM_ROWS - 1); j > (2 * NUM_ROWS - 1); j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//����
		for (k = (3 * NUM_ROWS - 1); k > (2 * NUM_ROWS - 1); k--)
			for (j = (NUM_COLS - 1); j > -1; j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//����
		for (k = (3 * NUM_ROWS - 1); k > (2 * NUM_ROWS - 1); k--)
			for (j = (3 * NUM_COLS - 1); j > (2 * NUM_COLS - 1); j--)
				datatest[j][k][i] = *(pdata1wei++);
	printf("padding���鴴�����!\n");//////////////
}


void local_rx(void)
{
	int i,j,k,m,n,m1,n1;//m1,n1����buff��ʼ��
	float * pbuff1wei;
	printf("���ڼ���localrx����...\n");
	pbuff1wei = buff_1wei;
	//for(i=0;i<BODUAN:i++)
		for(k=NUM_ROWS;k<(2*NUM_ROWS-1);k++)  //��������ֱ�����������ƶ�
			for(j=NUM_COLS;j<(2*NUM_COLS-1);j++)
			{
				for(i=0;i<BODUAN;i++)    //�ⴰ���鸳ֵ����������
				for(m=(j-T),m1=0; m<(j+T+1); m++,m1++)
					for(n=(k-T),n1=0; n<(k+T+1); n++,n1++)
						buff[m1][n1][i] = datatest[m][n][i];

				for(i=0;i<BODUAN;i++)	//��rx�����õ������в��ε�ĳһ���ص�
				Y[i]=datatest[j][k][i];

				for(i=0;i<BODUAN;i++)   //NaNֵ���
					for(m=(T-T1);m<(T+T1+1);m++)
						for(n=(T-T1);n<(T+T1+1);n++)
							buff[m][n][i] = 0;

				for(i=0;i<BODUAN;i++)  //��������һά��
					for(m=0;m<WIN_OUT;m++)
						for(n=0;n<WIN_OUT;n++)
							*(pbuff1wei++) = buff[n][m][i];
				pbuff1wei = buff_1wei;

				for(i=0;i<BODUAN;i++)  //һά����תȥ���ڴ���Ķ�ά����169*112
					for(m=0;m<(WIN_OUT*WIN_OUT-WIN_IN*WIN_IN);m++)
						{
						if(*pbuff1wei == 0)
						{
							pbuff1wei++;
							m--;
						}
						else
							buff_2wei[i][m] = *(pbuff1wei++);
						}
				pbuff1wei = buff_1wei;
				corr();
				inv();
				suanzi_rx();
			}
		printf("local_rx���Ӽ�����ɣ�\n");
}


void corr(void)
{
	int i,j,k;
	for (i = 0; i < BODUAN; i++)
			for (j = 0; j < BODUAN; j++)
				sigma[i][j] = 0;
		//printf("������Э����...\n");
		for (i = 0; i < BODUAN; i++)
			for (j = 0; j < BODUAN; j++)
			{
				for (k = 0; k < (WIN_OUT*WIN_OUT-WIN_IN*WIN_IN); k++)
					sigma[i][j] += ((buff_2wei[i][k]) * (buff_2wei[j][k]));
				//sigma[i][j] = sigma[i][j] / (NUM_ROWS * NUM_COLS);

			}
		//printf("��Э������ɣ�\n");
}



void inv(void)
{

	int i, j, k;
	float temp, fmax;
	int is[BODUAN] = {0};
	int js[BODUAN] = {0};
	//printf("��������...\n");
	for (i = 0; i < BODUAN; i++)
		for (j = 0; j < BODUAN; j++)
			inv_sigma[i][j] = sigma[i][j];

	for (k = 0; k < BODUAN; k++)
	{
		fmax = 0.0;
		for (i = k; i < BODUAN; i++)
		{
			for (j = k; j < BODUAN; j++)
			{
				temp = fabs(inv_sigma[i][j]); /*�����ֵ*/
				if (temp > fmax)
				{
					fmax = temp;
					is[k] = i; js[k] = j;
				}
			}
		}
		/*if((fmax+1.0)==1.0)
		{
		    return 0;
		}*///////
		if ((i = is[k]) != k)
		{
			for (j = 0; j < BODUAN; j++)
				swap(&inv_sigma[k][j], &inv_sigma[i][j]); /*��k�к͵�i�н���*/
		}
		if ((j = js[k]) != k)
		{
			for (i = 0; i < BODUAN; i++)
				swap(&inv_sigma[i][k], &inv_sigma[i][j]); /*��k�к͵�j�н���*/
		}
		inv_sigma[k][k] = 1.0 / inv_sigma[k][k]; /*��*/
		for (j = 0; j < BODUAN; j++)
		{
			if (j != k)
				inv_sigma[k][j] *= inv_sigma[k][k]; /*��k�е�ÿһ��Ԫ�ض��������Խ����ϵ�Ԫ��*/
		}
		for (i = 0; i < BODUAN; i++)
		{
			if (i != k)
			{
				for (j = 0; j < BODUAN; j++)
					if (j != k)
						inv_sigma[i][j] = inv_sigma[i][j] - inv_sigma[i][k] * inv_sigma[k][j];
			}

		}
		for (i = 0; i < BODUAN; i++)
		{
			if (i != k)
				inv_sigma[i][k] *= -inv_sigma[k][k];
		}
	}
	for (k = BODUAN - 1; k >= 0; k--)
	{
		if ((j = js[k]) != k)
		{
			for (i = 0; i < BODUAN; i++)
			{
				swap(&inv_sigma[j][i], &inv_sigma[k][i]); /*��j�����k�н���*/
			}
		}
		if ((i = is[k]) != k)
		{
			for (j = 0; j < BODUAN; j++)
			{
				swap(&inv_sigma[j][i], &inv_sigma[j][k]); /*��i�����k�н���*/
			}
		}
	}
	//printf("�������!\n");
}

void swap(float *a1, float *b1)
{
	float c1;
	c1 = *a1;
	*a1 = *b1;
	*b1 = c1;
}

void suanzi_rx()
{
	int j, k;

	//float * py;
	//py = Y;


		//d1[flag] = 0;

		d2[flag] = 0;
	//printf("������rx����...\n");


		for (j = 0; j < BODUAN; j++)
		{
			d1[j] = 0;
			//py = Y;
			for (k = 0; k < BODUAN; k++)
				d1[j] += (Y[k]* inv_sigma[k][j]);
		}

		for (j = 0; j < BODUAN; j++)
			d2[flag] += (d1[j] * Y[j]);

		flag++;
	//printf("rx�������ɹ���\n");
}



