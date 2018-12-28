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
float dat_float[SIZE / 4];  /*float数据*/
#pragma DATA_SECTION(X,".MY_MEM");
float X[BODUAN][NUM_ROWS * NUM_COLS];/*二维数组*/
//#pragma DATA_SECTION(rows_sum,".MY_MEM");
//float rows_sum[BODUAN]; /*各波段均值*/
#pragma DATA_SECTION(sigma,".MY_MEM");
float sigma[BODUAN][BODUAN]; /*X的协方差矩阵*/
#pragma DATA_SECTION(d2,".MY_MEM");
float d2[NUM_ROWS * NUM_COLS]; /*rx算法得出的算子*/
#pragma DATA_SECTION(inv_sigma,".MY_MEM");
float inv_sigma[BODUAN][BODUAN];/*sigma的逆矩阵*/
//#pragma DATA_SECTION(datatest,".MY_MEM");
#pragma DATA_SECTION(data_1wei,".MY_PAD");
float data_1wei[NUM_ROWS * NUM_COLS * BODUAN];/*一维原始数据*/
#pragma DATA_SECTION(datatest,".MY_PAD");
float datatest[3 * NUM_ROWS][3 * NUM_COLS][BODUAN]; /*三维直观结构数据*/
#pragma DATA_SECTION(buff_2wei,".MY_PAD");
float buff_2wei[BODUAN][WIN_OUT*WIN_OUT-WIN_IN*WIN_IN];/*二维去掉NaN矩阵*/


#pragma DATA_SECTION(buff,".MY_MEM");
float buff[WIN_OUT][WIN_OUT][BODUAN];/*buff遍历数组三维*/
#pragma DATA_SECTION(buff_1wei,".MY_MEM");
float buff_1wei[WIN_OUT*WIN_OUT*BODUAN];/*buff数组三维转一维*/
#pragma DATA_SECTION(Y,".MY_MEM");
float Y[BODUAN];/*原始数据的某一波段*/

#pragma DATA_SECTION(d1,".MY_MEM");
float d1[BODUAN];/*求rx算子时的中间数组*/







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

	printf("正在读取图像...\n");

	// 读取图像数组
	for ( i = 0; i < SIZE / 4; i++ )/////////////////
	{
		x = ((*p) << 24) | ((*(p + 1)) << 16) | ((*(p + 2)) << 8) | (*(p + 3));
		p += 4;
		dat_float[i] = *(float*)&x;
	}

	printf("图像读取完成!\n");
}

void opera()
{
	int i, j, flag1 = 0;

	/*变换成二维数组*/
	printf("正在转化二维数组...\n");
	for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
		for (j = 0; j < BODUAN; j++)
		{
			X[j][i] = dat_float[flag1++];
		}
	printf("二维数组转化成功!\n");
}


void padding()
{
	int i, j, k;
	//float data[NUM_ROWS][NUM_COLS][BODUAN];

	//float * px;
	//float * pdatatest;
	//float * pdata;
	float * pdata1wei;


	printf("正在创建padding数组...\n");
	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)
			for (j = 0; j < (NUM_ROWS * NUM_COLS); j++)
				*(pdata1wei++) = X[i][j];
	pdata1wei = data_1wei;
//
//	for (i = 0; i < BODUAN; i++)
//		for (j = 0; j < NUM_ROWS; j++)
//			for (k = 0; k < NUM_COLS; k++)
//				data[k][j][i] = *(px++);  //data_1wei一维变三维

	for (i = 0; i < BODUAN; i++)							//padding矩阵置零
			for (j = 0; j < (3 * NUM_COLS); j++)
				for (k = 0; k < (3 * NUM_ROWS); k++)
					datatest[j][k][i] = 0;

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//中
		for (j = NUM_COLS; j < (2 * NUM_COLS); j++)
			for (k = NUM_ROWS; k < (2 * NUM_ROWS); k++)
				datatest[k][j][i] = *(pdata1wei++);
	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//西
		for (k = NUM_ROWS; k < (2 * NUM_ROWS); k++)
			for (j = (NUM_COLS - 1); j > -1; j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//东
		for (k = NUM_ROWS; k < (2 * NUM_ROWS); k++)
			for (j = (3 * NUM_COLS - 1); j > (2 * NUM_COLS - 1); j--)
				datatest[j][k][i] = *(pdata1wei++);


	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//北
		for (k = NUM_ROWS - 1; k > -1; k--)
			for (j = NUM_COLS; j < (2 * NUM_COLS); j++)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//南
		for (k = (3 * NUM_ROWS - 1); k > (2 * NUM_ROWS - 1); k--)
			for (j = NUM_COLS; j < (2 * NUM_COLS); j++)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//西北
		for (k = (NUM_ROWS - 1); k > -1; k--)
			for (j = (NUM_COLS - 1); j > -1; j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//东北
		for (k = (NUM_ROWS - 1); k > -1; k--)
			for (j = (3 * NUM_ROWS - 1); j > (2 * NUM_ROWS - 1); j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//西南
		for (k = (3 * NUM_ROWS - 1); k > (2 * NUM_ROWS - 1); k--)
			for (j = (NUM_COLS - 1); j > -1; j--)
				datatest[j][k][i] = *(pdata1wei++);

	pdata1wei = data_1wei;
	for (i = 0; i < BODUAN; i++)							//东南
		for (k = (3 * NUM_ROWS - 1); k > (2 * NUM_ROWS - 1); k--)
			for (j = (3 * NUM_COLS - 1); j > (2 * NUM_COLS - 1); j--)
				datatest[j][k][i] = *(pdata1wei++);
	printf("padding数组创建完成!\n");//////////////
}


void local_rx(void)
{
	int i,j,k,m,n,m1,n1;//m1,n1控制buff初始化
	float * pbuff1wei;
	printf("正在计算localrx算子...\n");
	pbuff1wei = buff_1wei;
	//for(i=0;i<BODUAN:i++)
		for(k=NUM_ROWS;k<(2*NUM_ROWS-1);k++)  //控制整体直观数组层面的移动
			for(j=NUM_COLS;j<(2*NUM_COLS-1);j++)
			{
				for(i=0;i<BODUAN;i++)    //外窗数组赋值给缓冲数组
				for(m=(j-T),m1=0; m<(j+T+1); m++,m1++)
					for(n=(k-T),n1=0; n<(k+T+1); n++,n1++)
						buff[m1][n1][i] = datatest[m][n][i];

				for(i=0;i<BODUAN;i++)	//求rx算子用到的所有波段的某一像素点
				Y[i]=datatest[j][k][i];

				for(i=0;i<BODUAN;i++)   //NaN值替代
					for(m=(T-T1);m<(T+T1+1);m++)
						for(n=(T-T1);n<(T+T1+1);n++)
							buff[m][n][i] = 0;

				for(i=0;i<BODUAN;i++)  //缓冲数组一维化
					for(m=0;m<WIN_OUT;m++)
						for(n=0;n<WIN_OUT;n++)
							*(pbuff1wei++) = buff[n][m][i];
				pbuff1wei = buff_1wei;

				for(i=0;i<BODUAN;i++)  //一维数组转去掉内窗后的二维数组169*112
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
		printf("local_rx算子计算完成！\n");
}


void corr(void)
{
	int i,j,k;
	for (i = 0; i < BODUAN; i++)
			for (j = 0; j < BODUAN; j++)
				sigma[i][j] = 0;
		//printf("正在求协方差...\n");
		for (i = 0; i < BODUAN; i++)
			for (j = 0; j < BODUAN; j++)
			{
				for (k = 0; k < (WIN_OUT*WIN_OUT-WIN_IN*WIN_IN); k++)
					sigma[i][j] += ((buff_2wei[i][k]) * (buff_2wei[j][k]));
				//sigma[i][j] = sigma[i][j] / (NUM_ROWS * NUM_COLS);

			}
		//printf("求协方差完成！\n");
}



void inv(void)
{

	int i, j, k;
	float temp, fmax;
	int is[BODUAN] = {0};
	int js[BODUAN] = {0};
	//printf("正在求逆...\n");
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
				temp = fabs(inv_sigma[i][j]); /*找最大值*/
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
				swap(&inv_sigma[k][j], &inv_sigma[i][j]); /*第k行和第i列交换*/
		}
		if ((j = js[k]) != k)
		{
			for (i = 0; i < BODUAN; i++)
				swap(&inv_sigma[i][k], &inv_sigma[i][j]); /*第k列和第j列交换*/
		}
		inv_sigma[k][k] = 1.0 / inv_sigma[k][k]; /*求倒*/
		for (j = 0; j < BODUAN; j++)
		{
			if (j != k)
				inv_sigma[k][j] *= inv_sigma[k][k]; /*第k行的每一个元素都除以主对角线上的元素*/
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
				swap(&inv_sigma[j][i], &inv_sigma[k][i]); /*第j行与第k行交换*/
			}
		}
		if ((i = is[k]) != k)
		{
			for (j = 0; j < BODUAN; j++)
			{
				swap(&inv_sigma[j][i], &inv_sigma[j][k]); /*第i列与第k列交换*/
			}
		}
	}
	//printf("求逆完成!\n");
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
	//printf("正在求rx算子...\n");


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
	//printf("rx算子求解成功！\n");
}



