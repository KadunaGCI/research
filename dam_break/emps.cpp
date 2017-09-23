#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

#define IN_FILE "C:\\Users\\wang\\Desktop\\solid-fluid single\\dampar\\dam\\dambreak.prof"

#define PCL_DST 0.005					//平均粒子間距離
#define MIN_X  (0.0 - PCL_DST*3)	//解析領域のx方向の最小値
#define MIN_Y  (0.0 - PCL_DST*3)	//解析領域のy方向の最小値
#define MIN_Z  (0.0 - PCL_DST*3)	//解析領域のz方向の最小値
#define MAX_X  (3.5 + PCL_DST*3)	//解析領域のx方向の最大値
#define MAX_Y  (0.5 + PCL_DST*3)	//解析領域のy方向の最大値
#define MAX_Z  (0.75 + PCL_DST*10)	//解析領域のz方向の最大値

#define GST -1			//計算対象外粒子の種類番号
#define FLD 0				//流体粒子の種類番号
#define WLL  1			//壁粒子の種類番号
#define GATE  2
#define RIGID0 3
#define RIGID1 4
#define RIGID2 5

#define NumberofRigidBodies 3
#define NUM_TYP  6		//粒子の種類数

#define DNS_FLD 1000		//流体粒子の密度
#define DNS_WLL 1000		//壁粒子の密度
#define DNS_RIGID0 464	
#define DNS_RIGID1 111	
#define DNS_RIGID2 111	

#define DT 0.0002			//時間刻み幅
#define dt_inv   double (1/DT)	
#define FIN_TIM 12		//時間の上限
#define SND 22.0			//音速
#define OPT_FQC 200		//出力間隔を決める反復数

#define KNM_VSC 0.000001	//動粘性係数
#define DIM 3				//次元数
#define CRT_NUM 0.1		//クーラン条件数
#define COL_RAT 0.2		//接近した粒子の反発率
#define DST_LMT_RAT 0.9	//これ以上の粒子間の接近を許さない距離の係数
#define G_X 0.0			//重力加速度のx成分
#define G_Y 0.0			//重力加速度のy成分
#define G_Z -9.8			//重力加速度のz成分
#define WEI(dist, re) ((re/dist) - 1.0)	//重み関数
#define kt 2.0
#define pt 0.5

FILE* fp;
char filename[256];
int iLP,iF;
double TIM;
int nP,nr0,nr1,nr2;
double *Acc,*Pos,*Vel,*Prs,*pav;
double *PrePos;
double qs0, qx0, qy0, qz0;
double qs1, qx1, qy1, qz1;
double qs2, qx2, qy2, qz2;
int *Typ;
double r,r2;
double DB,DB2,DBinv;
int nBx,nBy,nBz,nBxy,nBxyz;
int *bfst,*blst,*nxt;
double n0,lmd,A1,A2,A3,rlim,rlim2,COL;
double Dns[NUM_TYP],invDns[NUM_TYP];
double InertiaTensorInv0[3][3];
double InertiaTensorInv1[3][3];
double InertiaTensorInv2[3][3];

void ChkPcl(int i){
	if(	Pos[i*3  ]>MAX_X || Pos[i*3  ]<MIN_X ||
		Pos[i*3+1]>MAX_Y || Pos[i*3+1]<MIN_Y ||
		Pos[i*3+2]>MAX_Z || Pos[i*3+2]<MIN_Z)
	{
		Typ[i] = GST;
		Prs[i]=Vel[i*3]=Vel[i*3+1]=Vel[i*3+2]=0.0;
	}
}

void RdDat(void) {
	fopen_s(&fp, IN_FILE, "r");
	fscanf_s(fp, "%d %d %d %d", &nP, &nr0, &nr1, &nr2);
	printf("nP: %d\n",nP);
	Acc = (double*)malloc(sizeof(double)*nP*3);	//粒子の加速度
	Pos = (double*)malloc(sizeof(double)*nP*3);	//粒子の座標
	PrePos = (double*)malloc(sizeof(double)*nP*3);//粒子の座標
	Vel = (double*)malloc(sizeof(double)*nP*3);	//粒子の速度
	Prs = (double*)malloc(sizeof(double)*nP);	//粒子の圧力
	pav = (double*)malloc(sizeof(double)*nP);	//時間平均された粒子の圧力
	Typ = (int*)malloc(sizeof(int)*nP);			//粒子の種類
	for(int i=0;i<nP;i++) {
		int a[2];
		double b[8];
		fscanf_s(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&a[0],&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		Typ[i]=a[1];
		Pos[i*3]=b[0];	Pos[i*3+1]=b[1];	Pos[i*3+2]=b[2];
		Vel[i*3]=b[3];	Vel[i*3+1]=b[4];	Vel[i*3+2]=b[5];
		Prs[i]=b[6];		pav[i]=b[7];
	}
	fclose(fp);
	for(int i=0;i<nP;i++) {ChkPcl(i);}
	for(int i=0;i<nP*3;i++) {Acc[i]=0.0;}
	for (int i = 0; i<nP * 3; i++) { PrePos[i] = Pos[i]; }
}

void WrtDat(void) {
	char outout_filename[256];
	sprintf_s(outout_filename, "output%05d.prof",iF);
	fopen_s(&fp, outout_filename, "w");
	fprintf(fp,"%d\n",nP);
	for(int i=0;i<nP;i++) {
		int a[2];
		double b[8];
		a[0]=i;	a[1]=Typ[i];
		b[0]=Pos[i*3];	b[1]=Pos[i*3+1];	b[2]=Pos[i*3+2];
		b[3]=Vel[i*3];	b[4]=Vel[i*3+1];	b[5]=Vel[i*3+2];
		b[6]=Prs[i];		b[7]=pav[i]/OPT_FQC;
		if (Typ[i] == GST)continue;
		fprintf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",a[0],a[1],b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7]);
		pav[i]=0.0;
	}
	fclose(fp);
	iF++;
}

void AlcBkt(void) {
	r = PCL_DST*3.1;		//影響半径
	r2 = r*r;
	DB = r*(1.0+CRT_NUM);	//バケット1辺の長さ
	DB2 = DB*DB;
	DBinv = 1.0/DB;
	nBx = (int)((MAX_X - MIN_X)*DBinv) + 3;//解析領域内のx方向のバケット数
	nBy = (int)((MAX_Y - MIN_Y)*DBinv) + 3;//解析領域内のy方向のバケット数
	nBz = (int)((MAX_Z - MIN_Z)*DBinv) + 3;//解析領域内のz方向のバケット数
	nBxy = nBx*nBy;
	nBxyz = nBx*nBy*nBz;		//解析領域内のバケット数
	printf("nBx:%d  nBy:%d  nBz:%d  nBxy:%d  nBxyz:%d\n",nBx,nBy,nBz,nBxy,nBxyz);
	bfst = (int*)malloc(sizeof(int) * nBxyz);	//バケットに格納された先頭の粒子番号
	blst = (int*)malloc(sizeof(int) * nBxyz);	//バケットに格納された最後尾の粒子番号
	nxt  = (int*)malloc(sizeof(int) * nP);		//同じバケット内の次の粒子番号
}

void SetPara(void){
	double tn0 = 0.0;
	double tlmd = 0.0;
	for (int ix = -4; ix < 5; ix++){
		for (int iy = -4; iy < 5; iy++){
			for (int iz = -4; iz < 5; iz++){
				double x = PCL_DST* (double)ix;
				double y = PCL_DST* (double)iy;
				double z = PCL_DST* (double)iz;
				double dist2 = x*x + y*y + z*z;
				if (dist2 <= r2){
					if (dist2 == 0.0)continue;
					double dist = sqrt(dist2);
					tn0 += WEI(dist, r);
					tlmd += dist2 * WEI(dist, r);
				}
			}
		}
	}
	n0 = tn0;			//初期粒子数密度
	lmd = tlmd / tn0;	//ラプラシアンモデルの係数λ
	A1 = 2.0*KNM_VSC*DIM / n0 / lmd;//粘性項の計算に用いる係数
	A2 = SND*SND / n0;				//圧力の計算に用いる係数
	A3 = -DIM / n0;					//圧力勾配項の計算に用いる係数
	Dns[FLD] = DNS_FLD;
	Dns[WLL] = DNS_WLL;
	Dns[RIGID0] = DNS_RIGID0;
	Dns[RIGID1] = DNS_RIGID1;
	Dns[RIGID2] = DNS_RIGID2;
	invDns[FLD] = 1.0 / DNS_FLD;
	invDns[WLL] = 1.0 / DNS_WLL;
	invDns[RIGID0] = 1.0/ DNS_RIGID0;
	invDns[RIGID1] = 1.0 / DNS_RIGID1;
	invDns[RIGID2] = 1.0 / DNS_RIGID2;
	rlim = PCL_DST * DST_LMT_RAT;//これ以上の粒子間の接近を許さない距離
	rlim2 = rlim*rlim;
	COL = 1.0 + COL_RAT;
	iLP = 0;			//反復数
	iF = 0;			//ファイル番号
	TIM = 0.0;		//時刻
}


void inMatrix(double a[DIM][DIM], int n){
	int i, j, k;
	double sum;
	double b[DIM][DIM];
	double d[DIM];

	for (i = 0; i<n; i++)
	{
		for (j = i; j<n; j++)
		{
			sum = 0.0;
			for (k = 0; k<i; k++) sum += a[i][k] * a[k][j];
			a[i][j] -= sum;
			if (j>i) a[j][i] = a[i][j] / a[i][i];
		}
	}
	/*Inverse of Lower triangular matrix solution*/
	for (j = 0; j<n; j++)
	{
		d[j] = a[j][j];
		b[j][j] = 1.0;
		for (i = j + 1; i<n; i++)
		{
			b[i][j] = -a[i][j];
			for (k = j + 1; k<i; k++) b[i][j] -= a[i][k] * b[k][j];
		}
	}
	/*Multiply (l)^T x d x l to obtain inverse of A*/
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			a[i][j] = 0.0;
			for (k = i>j ? i : j; k<n; k++) a[i][j] += b[k][i] * b[k][j] / d[k];
		}
	}
}



void init_rigid0(void){

	double CenterofGravity[3] = { 0.0, 0.0, 0.0 };
	
	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID0){
			CenterofGravity[0] += Pos[i * 3];
			CenterofGravity[1] += Pos[i * 3 + 1];
			CenterofGravity[2] += Pos[i * 3 + 2];			
		}
	}
	CenterofGravity[0] = CenterofGravity[0] / nr0;
	CenterofGravity[1] = CenterofGravity[1] / nr0;
	CenterofGravity[2] = CenterofGravity[2] / nr0;
	qs0 = 1.0; qx0 = 0.0; qy0 = 0.0; qz0 = 0.0;
	double InertiaTensor[3][3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID0){
			double  v3_tmp0 = Pos[3 * i] - CenterofGravity[0];
			double  v3_tmp1 = Pos[3 * i + 1] - CenterofGravity[1];
			double  v3_tmp2 = Pos[3 * i + 2] - CenterofGravity[2];			
			double v3_tmp_squaredLength = v3_tmp0*v3_tmp0 + v3_tmp1*v3_tmp1 + v3_tmp2*v3_tmp2;
			double dd0 = Dns[Typ[i]] * PCL_DST*PCL_DST*PCL_DST;
			double dd1 = v3_tmp_squaredLength + PCL_DST*PCL_DST*0.1666666666666666666666666;			

			InertiaTensor[0][0] -= (v3_tmp0*v3_tmp0 - dd1) * dd0;
			InertiaTensor[0][1] -= (v3_tmp0*v3_tmp1) * dd0;
			InertiaTensor[0][2] -= (v3_tmp0*v3_tmp2) * dd0;
			InertiaTensor[1][0] -= (v3_tmp1*v3_tmp0) * dd0;
			InertiaTensor[1][1] -= (v3_tmp1*v3_tmp1 - dd1) * dd0;
			InertiaTensor[1][2] -= (v3_tmp1*v3_tmp2) * dd0;
			InertiaTensor[2][0] -= (v3_tmp2*v3_tmp0) * dd0;
			InertiaTensor[2][1] -= (v3_tmp2*v3_tmp1) * dd0;
			InertiaTensor[2][2] -= (v3_tmp2*v3_tmp2 - dd1) * dd0;
		}
	}

	double buf;
	double A00, A01, A02,A10, A11, A12,A20, A21, A22;
	A00 = 1.0;		A01 = 0.0;		A02 = 0.0;
	A10 = 0.0;		A11 = 1.0;		A12 = 0.0;
	A20 = 0.0;		A21 = 0.0;		A22 = 1.0;
	double B00, B01, B02, B10, B11, B12, B20, B21, B22;
	B00 = InertiaTensor[0][0];	B01 = InertiaTensor[0][1];	B02 = InertiaTensor[0][2];
	B10 = InertiaTensor[1][0];	B11 = InertiaTensor[1][1];	B12 = InertiaTensor[1][2];
	B20 = InertiaTensor[2][0];	B21 = InertiaTensor[2][1];	B22 = InertiaTensor[2][2];

	buf = 1.0 / B00;
	B01 *= buf;	B02 *= buf;
	A00 *= buf;	A01 *= buf;	A02 *= buf;
	B11 -= B01*B10;	B12 -= B02*B10;
	B21 -= B01*B20;	B22 -= B02*B20;
	A10 -= A00*B10;	A11 -= A01*B10;	A12 -= A02*B10;
	A20 -= A00*B20;	A21 -= A01*B20;	A22 -= A02*B20;

	buf = 1.0 / B11;
	B12 *= buf;
	A10 *= buf;	A11 *= buf;	A12 *= buf;
	B02 -= B12*B01;	B22 -= B12*B21;
	A00 -= A10*B01;	A01 -= A11*B01;	A02 -= A12*B01;
	A20 -= A10*B21;	A21 -= A11*B21;	A22 -= A12*B21;

	buf = 1.0 / B22;
	A20 *= buf;	A21 *= buf;	A22 *= buf;
	A00 -= A20*B02;	A01 -= A21*B02;	A02 -= A22*B02;
	A10 -= A20*B12;	A11 -= A21*B12;	A12 -= A22*B12;
	InertiaTensorInv0[0][0] = A00; InertiaTensorInv0[0][1] = A01; InertiaTensorInv0[0][2] = A02;
	InertiaTensorInv0[1][0] = A10; InertiaTensorInv0[1][1] = A11; InertiaTensorInv0[1][2] = A12;
	InertiaTensorInv0[2][0] = A20; InertiaTensorInv0[2][1] = A21; InertiaTensorInv0[2][2] = A22;
}
/*
void init_rigid1(void){
	double CenterofGravity[3] = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID1){
			CenterofGravity[0] += Pos[i * 3];
			CenterofGravity[1] += Pos[i * 3 + 1];
			CenterofGravity[2] += Pos[i * 3 + 2];
		}
	}
	CenterofGravity[0] = CenterofGravity[0] / nr1;
	CenterofGravity[1] = CenterofGravity[1] / nr1;
	CenterofGravity[2] = CenterofGravity[2] / nr1;
	qs1 = 1.0; qx1 = 0.0; qy1 = 0.0; qz1 = 0.0;
	double InertiaTensor[3][3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID1){
			double  v3_tmp0 = Pos[3 * i] - CenterofGravity[0];
			double  v3_tmp1 = Pos[3 * i + 1] - CenterofGravity[1];
			double  v3_tmp2 = Pos[3 * i + 2] - CenterofGravity[2];
			double v3_tmp_squaredLength = v3_tmp0*v3_tmp0 + v3_tmp1*v3_tmp1 + v3_tmp2*v3_tmp2;
			double dd0 = Dns[Typ[i]] * PCL_DST*PCL_DST*PCL_DST;
			double dd1 = v3_tmp_squaredLength + PCL_DST*PCL_DST*0.1666666666666666666666666;

			InertiaTensor[0][0] -= (v3_tmp0*v3_tmp0 - dd1) * dd0;
			InertiaTensor[0][1] -= (v3_tmp0*v3_tmp1) * dd0;
			InertiaTensor[0][2] -= (v3_tmp0*v3_tmp2) * dd0;
			InertiaTensor[1][0] -= (v3_tmp1*v3_tmp0) * dd0;
			InertiaTensor[1][1] -= (v3_tmp1*v3_tmp1 - dd1) * dd0;
			InertiaTensor[1][2] -= (v3_tmp1*v3_tmp2) * dd0;
			InertiaTensor[2][0] -= (v3_tmp2*v3_tmp0) * dd0;
			InertiaTensor[2][1] -= (v3_tmp2*v3_tmp1) * dd0;
			InertiaTensor[2][2] -= (v3_tmp2*v3_tmp2 - dd1) * dd0;
		}
	}

	double buf;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	A00 = 1.0;		A01 = 0.0;		A02 = 0.0;
	A10 = 0.0;		A11 = 1.0;		A12 = 0.0;
	A20 = 0.0;		A21 = 0.0;		A22 = 1.0;
	double B00, B01, B02, B10, B11, B12, B20, B21, B22;
	B00 = InertiaTensor[0][0];	B01 = InertiaTensor[0][1];	B02 = InertiaTensor[0][2];
	B10 = InertiaTensor[1][0];	B11 = InertiaTensor[1][1];	B12 = InertiaTensor[1][2];
	B20 = InertiaTensor[2][0];	B21 = InertiaTensor[2][1];	B22 = InertiaTensor[2][2];

	buf = 1.0 / B00;
	B01 *= buf;	B02 *= buf;
	A00 *= buf;	A01 *= buf;	A02 *= buf;
	B11 -= B01*B10;	B12 -= B02*B10;
	B21 -= B01*B20;	B22 -= B02*B20;
	A10 -= A00*B10;	A11 -= A01*B10;	A12 -= A02*B10;
	A20 -= A00*B20;	A21 -= A01*B20;	A22 -= A02*B20;

	buf = 1.0 / B11;
	B12 *= buf;
	A10 *= buf;	A11 *= buf;	A12 *= buf;
	B02 -= B12*B01;	B22 -= B12*B21;
	A00 -= A10*B01;	A01 -= A11*B01;	A02 -= A12*B01;
	A20 -= A10*B21;	A21 -= A11*B21;	A22 -= A12*B21;

	buf = 1.0 / B22;
	A20 *= buf;	A21 *= buf;	A22 *= buf;
	A00 -= A20*B02;	A01 -= A21*B02;	A02 -= A22*B02;
	A10 -= A20*B12;	A11 -= A21*B12;	A12 -= A22*B12;
	InertiaTensorInv1[0][0] = A00; InertiaTensorInv1[0][1] = A01; InertiaTensorInv1[0][2] = A02;
	InertiaTensorInv1[1][0] = A10; InertiaTensorInv1[1][1] = A11; InertiaTensorInv1[1][2] = A12;
	InertiaTensorInv1[2][0] = A20; InertiaTensorInv1[2][1] = A21; InertiaTensorInv1[2][2] = A22;
}


void init_rigid2(void){
	double CenterofGravity[3] = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID2){
			CenterofGravity[0] += Pos[i * 3];
			CenterofGravity[1] += Pos[i * 3 + 1];
			CenterofGravity[2] += Pos[i * 3 + 2];
		}
	}
	CenterofGravity[0] = CenterofGravity[0] / nr2;
	CenterofGravity[1] = CenterofGravity[1] / nr2;
	CenterofGravity[2] = CenterofGravity[2] / nr2;
	qs2 = 1.0; qx2 = 0.0; qy2 = 0.0; qz2 = 0.0;
	double InertiaTensor[3][3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID2){
			double  v3_tmp0 = Pos[3 * i] - CenterofGravity[0];
			double  v3_tmp1 = Pos[3 * i + 1] - CenterofGravity[1];
			double  v3_tmp2 = Pos[3 * i + 2] - CenterofGravity[2];
			double v3_tmp_squaredLength = v3_tmp0*v3_tmp0 + v3_tmp1*v3_tmp1 + v3_tmp2*v3_tmp2;
			double dd0 = Dns[Typ[i]] * PCL_DST*PCL_DST*PCL_DST;
			double dd1 = v3_tmp_squaredLength + PCL_DST*PCL_DST*0.1666666666666666666666666;

			InertiaTensor[0][0] -= (v3_tmp0*v3_tmp0 - dd1) * dd0;
			InertiaTensor[0][1] -= (v3_tmp0*v3_tmp1) * dd0;
			InertiaTensor[0][2] -= (v3_tmp0*v3_tmp2) * dd0;
			InertiaTensor[1][0] -= (v3_tmp1*v3_tmp0) * dd0;
			InertiaTensor[1][1] -= (v3_tmp1*v3_tmp1 - dd1) * dd0;
			InertiaTensor[1][2] -= (v3_tmp1*v3_tmp2) * dd0;
			InertiaTensor[2][0] -= (v3_tmp2*v3_tmp0) * dd0;
			InertiaTensor[2][1] -= (v3_tmp2*v3_tmp1) * dd0;
			InertiaTensor[2][2] -= (v3_tmp2*v3_tmp2 - dd1) * dd0;
		}
	}

	double buf;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	A00 = 1.0;		A01 = 0.0;		A02 = 0.0;
	A10 = 0.0;		A11 = 1.0;		A12 = 0.0;
	A20 = 0.0;		A21 = 0.0;		A22 = 1.0;
	double B00, B01, B02, B10, B11, B12, B20, B21, B22;
	B00 = InertiaTensor[0][0];	B01 = InertiaTensor[0][1];	B02 = InertiaTensor[0][2];
	B10 = InertiaTensor[1][0];	B11 = InertiaTensor[1][1];	B12 = InertiaTensor[1][2];
	B20 = InertiaTensor[2][0];	B21 = InertiaTensor[2][1];	B22 = InertiaTensor[2][2];

	buf = 1.0 / B00;
	B01 *= buf;	B02 *= buf;
	A00 *= buf;	A01 *= buf;	A02 *= buf;
	B11 -= B01*B10;	B12 -= B02*B10;
	B21 -= B01*B20;	B22 -= B02*B20;
	A10 -= A00*B10;	A11 -= A01*B10;	A12 -= A02*B10;
	A20 -= A00*B20;	A21 -= A01*B20;	A22 -= A02*B20;

	buf = 1.0 / B11;
	B12 *= buf;
	A10 *= buf;	A11 *= buf;	A12 *= buf;
	B02 -= B12*B01;	B22 -= B12*B21;
	A00 -= A10*B01;	A01 -= A11*B01;	A02 -= A12*B01;
	A20 -= A10*B21;	A21 -= A11*B21;	A22 -= A12*B21;

	buf = 1.0 / B22;
	A20 *= buf;	A21 *= buf;	A22 *= buf;
	A00 -= A20*B02;	A01 -= A21*B02;	A02 -= A22*B02;
	A10 -= A20*B12;	A11 -= A21*B12;	A12 -= A22*B12;
	InertiaTensorInv2[0][0] = A00; InertiaTensorInv2[0][1] = A01; InertiaTensorInv2[0][2] = A02;
	InertiaTensorInv2[1][0] = A10; InertiaTensorInv2[1][1] = A11; InertiaTensorInv2[1][2] = A12;
	InertiaTensorInv2[2][0] = A20; InertiaTensorInv2[2][1] = A21; InertiaTensorInv2[2][2] = A22;
}
*/
void MkBkt(void) {
	for (int i = 0; i< nBxyz; i++){ bfst[i] = -1; }
	for (int i = 0; i< nBxyz; i++){ blst[i] = -1; }
	for (int i = 0; i< nP; i++){ nxt[i] = -1; }
	for (int i = 0; i<nP; i++){
		if (Typ[i] == GST)continue;
		int ix = (int)((Pos[i * 3] - MIN_X)*DBinv) + 1;
		int iy = (int)((Pos[i * 3 + 1] - MIN_Y)*DBinv) + 1;
		int iz = (int)((Pos[i * 3 + 2] - MIN_Z)*DBinv) + 1;
		int ib = iz*nBxy + iy*nBx + ix;
		int j = blst[ib];
		blst[ib] = i;
		if (j == -1){ bfst[ib] = i; }
		else{ nxt[j] = i; }
	}
}
double min(double a, double b, double c, double d, double e){
	double num[6] = { a, b, c, d, e };
	double p = num[0];
	for (int i = 1; i < 5; i++)
	{
		if (num[i] <= p)
			p = num[i];
	}
	return p;
}

double fi(double k){
	double p;
	if (k < -0.35 * kt && k >= -0.4 *kt)
	{
		p = double(20 * pt*(22.3482 - 29.726)*(k + 0.35 *kt) + 22.3482); return p;
	}
	else if (k < -0.3 * kt  && k >= -0.35 *kt)
	{
		p = double(20 * pt* (18.4537 - 22.3482)*(k + 0.3 *kt) + 18.4537); return p;
	}
	else if (k < -0.25 * kt  && k >= -0.3 *kt)
	{
		p = double(20 * pt* (15.9436 - 18.4537)*(k + 0.25 *kt) + 15.9436); return p;
	}
	else if (k < -0.2 * kt  && k >= -0.25 *kt)
	{
		p = double(20 * pt* (14.12 - 15.9436)*(k + 0.2 *kt) + 14.12); return p;
	}
	else if (k < -0.15 *kt  && k >= -0.2 * kt)
	{
		p = double(20 * pt* (12.683 - 14.12)*(k + 0.15 *kt) + 12.683); return p;
	}
	else if (k < -0.1 *kt  && k >= -0.15 * kt)
	{
		p = double(20 * pt* (11.487 - 12.683)*(k + 0.10 *kt) + 11.487); return p;
	}
	else if (k < -0.05 && k >= -0.1 *kt)
	{
		p = double(20 * pt*  (10.4489 - 11.487)*(k + 0.05 *kt) + 10.4489); return p;
	}
	else if (k < 0 && k >= -0.05 *kt)
	{
		p = double(20 * pt*  (9.52 - 10.4489)*k + 9.52); return p;
	}
	else if (k < 0.05 * kt  && k >= 0.0)
	{
		p = double(20 * pt*  (8.676 - 9.52)*(k - 0.05 *kt) + 8.676); return p;
	}
	else if (k < 0.1 * kt  && k >= 0.05 *kt)
	{
		p = double(20 * pt*  (7.96 - 8.676)*(k - 0.1 *kt) + 7.96); return p;
	}
	else if (k < 0.15 * kt  && k >= 0.1 *kt)
	{
		p = double(20 * pt*  (7.2971 - 7.96)*(k - 0.15 *kt) + 7.2971); return p;
	}
	else if (k < 0.2 * kt  && k >= 0.15 *kt)
	{
		p = double(20 * pt*  (6.7 - 7.2971)*(k - 0.2 *kt) + 6.7); return p;
	}
	else if (k < 0.25 * kt  && k >= 0.2 *kt)
	{
		p = double(20 * pt*(6.135 - 6.7)*(k - 0.25 *kt) + 6.135); return p;
	}
	else if (k < 0.3 * kt  && k >= 0.25 *kt)
	{
		p = double(20 * pt*(5.6 - 6.135)*(k - 0.3 *kt) + 5.6); return p;
	}
	else if (k < 0.35 *kt  && k >= 0.3 * kt)
	{
		p = double(20 * pt* (5.097 - 5.6)*(k - 0.35 *kt) + 5.097); return p;
	}
	else if (k < 0.4 *kt  && k >= 0.35 * kt)
	{
		p = double(20 * pt* (4.69 - 5.097)*(k - 0.4 *kt) + 4.69); return p;
	}
	else if (k < 0.45 *kt  && k >= 0.4 *kt)
	{
		p = double(20 * pt* (4.308 - 4.69)*(k - 0.45 *kt) + 4.308); return p;
	}
	else if (k < 0.5 *kt  && k >= 0.45 *kt)
	{
		p = double(20 * pt* (3.94 - 4.308)*(k - 0.5 *kt) + 3.94); return p;
	}
	else if (k < 0.55 *kt  && k >= 0.5 * kt)
	{
		p = double(20 * pt*  (3.586 - 3.94)*(k - 0.55 *kt) + 3.586); return p;
	}
	else if (k < 0.6 *kt  && k >= 0.55 * kt)
	{
		p = double(20 * pt*  (3.248 - 3.586)*(k - 0.6 *kt) + 3.248); return p;
	}
	else if (k < 0.65 *kt  && k >= 0.6 * kt)
	{
		p = double(20 * pt*(2.946 - 3.248)*(k - 0.65 *kt) + 2.946); return p;
	}
	else if (k < 0.7 *kt  && k >= 0.65 * kt)
	{
		p = double(20 * pt*(2.65 - 2.946)*(k - 0.7 *kt) + 2.65); return p;
	}
	else if (k < 0.75 *kt  && k >= 0.7 *kt)
	{
		p = double(20 * pt* (2.378 - 2.65)*(k - 0.75 *kt) + 2.378); return p;
	}
	else if (k < 0.8 *kt  && k >= 0.75 *kt)
	{
		p = double(20 * pt* (2.1098 - 2.378)*(k - 0.8 *kt) + 2.1098); return p;
	}
	else if (k < 0.85 *kt  && k >= 0.8 * kt)
	{
		p = double(20 * pt*  (1.8519 - 2.1098)*(k - 0.85 *kt) + 1.8519); return p;
	}
	else if (k < 0.9 *kt  && k >= 0.85 * kt)
	{
		p = double(20 * pt*  (1.60356 - 1.8519)*(k - 0.9 *kt) + 1.60356); return p;
	}
	else if (k < 0.95 *kt  && k >= 0.9 *kt)
	{
		p = double(20 * pt*  (1.3644 - 1.60356)*(k - 0.95 *kt) + 1.3644); return p;
	}
	else if (k < 1.0 *kt  && k >= 0.95 *kt)
	{
		p = double(20 * pt*  (1.134 - 1.3644)*(k - 1.0 *kt) + 1.134); return p;
	}
	else if (k < 1.1 *kt  && k >= 1.0 *kt)
	{
		p = double(10 * pt*  (0.7645 - 1.134)*(k - 1.1 *kt) + 0.7645); return p;
	}
	else if (k < 1.2 *kt  && k >= 1.1 *kt)
	{
		p = double(10 * pt* (0.494 - 0.7645)*(k - 1.2 *kt) + 0.494); return p;
	}
	else if (k < 1.3 * kt  && k >= 1.2 * kt)
	{
		p = double(10 * pt* (0.246 - 0.494)*(k - 1.3 *kt) + 0.246); return p;
	}
	else if (k < 1.4 * kt  && k >= 1.3 *kt)
	{
		p = double(10 * pt* (0.105 - 0.246)*(k - 1.4 *kt) + 0.105); return p;
	}
	else if (k < 1.5 * kt  && k >= 1.4 *kt)
	{
		p = double(10 * pt*  (0.05 - 0.105)*(k - 1.5 *kt) + 0.05); return p;
	}
	else if (k < 1.6 *kt  && k >= 1.5 *kt)
	{
		p = double(10 * pt* (0.0 - 0.05)*(k - 1.6 *kt) + 0.0); return p;
	}
	else if (k >= 1.6 * kt) return 0.0;
	else return 29.72;
}

void VscTrm(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
		if (Typ[i] == FLD || Typ[i] == RIGID0 || Typ[i] == RIGID1 || Typ[i] == RIGID2){
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dist2 = v0*v0+v1*v1+v2*v2;
				if(dist2<r2){
				if(j!=i && Typ[j]!=GST){
					double dist = sqrt(dist2);
					double w =  WEI(dist, r);
					Acc_x +=(Vel[j*3  ]-vec_ix)*w;
					Acc_y +=(Vel[j*3+1]-vec_iy)*w;
					Acc_z +=(Vel[j*3+2]-vec_iz)*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		/*
		if (Pos[i * 3] <= 1.6*PCL_DST || Pos[i * 3 + 1] <= 1.6*PCL_DST || Pos[i * 3 + 2] <= 1.6*PCL_DST || Pos[i * 3 + 1] > (MAX_Y - 4.6*PCL_DST) || Pos[i * 3] > (MAX_X - 4.6*PCL_DST)){
			double PL = double(kt / PCL_DST);
			int nx = int((MAX_X - 2.9*PCL_DST)*PL);
			int ny = int((MAX_Y - 2.9*PCL_DST)*PL);
			int nz = int((MAX_Z - 2.9*PCL_DST)*PL);

			int tx = int(pos_ix *PL + 0.0004);
			int ty = int(pos_iy *PL + 0.0004);
			int tz = int(pos_iz *PL + 0.0004);
			double a = pos_ix - tx *PCL_DST*pt;
			double b = (1.0 + tx)*PCL_DST*pt - pos_ix;
			double c = pos_iz - (tz + 0.0)*PCL_DST*pt;
			double d = (1.0 + tz)*PCL_DST*pt - pos_iz;
			double e = pos_iy - (ty + 0.0)*PCL_DST*pt;
			double f = (1.0 + ty)*PCL_DST*pt - pos_iy;
			int r1 = min(tx, ty, tz + 1, nx - tx, ny - ty);
			int r2 = min(tx + 1, ty, tz + 1, nx - tx - 1, ny - ty);
			int r3 = min(tx, ty, tz, nx - tx, ny - ty);
			int r4 = min(tx + 1, ty, tz, nx - tx - 1, ny - ty);
			int r5 = min(tx, ty + 1, tz + 1, nx - tx, ny - ty - 1);
			int r6 = min(tx + 1, ty + 1, tz + 1, nx - tx - 1, ny - ty - 1);
			int r7 = min(tx, ty + 1, tz, nx - tx, ny - ty - 1);
			int r8 = min(tx + 1, ty + 1, tz, nx - tx - 1, ny - ty - 1);

			double riw = double((b*c*f*r1 + a*c*f*r2 + b*d*f*r3 + a*d*f*r4 + b*c*e*r5 + a*c*e*r6 + b*d*e*r7 + a*d*e*r8) *PL *PL *PL);
			double w = fi(riw);
			Acc_x += (-2 * vec_ix)*w;
			Acc_y += (-2 * vec_iy)*w;
			Acc_z += (-2 * vec_iz)*w;
		}
		*/
		Acc[i*3  ]=Acc_x*A1 + G_X;
		Acc[i*3+1]=Acc_y*A1 + G_Y;
		Acc[i*3+2]=Acc_z*A1 + G_Z;
	}}
}

void UpPcl1(){
	for(int i=0;i<nP;i++){
		if (Typ[i] == FLD || Typ[i] == RIGID0 || Typ[i] == RIGID1 || Typ[i] == RIGID2){
			Vel[i*3  ] +=Acc[i*3  ]*DT;	Vel[i*3+1] +=Acc[i*3+1]*DT;	Vel[i*3+2] +=Acc[i*3+2]*DT;
			Pos[i*3  ] +=Vel[i*3  ]*DT;		Pos[i*3+1] +=Vel[i*3+1]*DT;		Pos[i*3+2] +=Vel[i*3+2]*DT;
			Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
			ChkPcl(i);
		}
	}
}

void ChkCol(){
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if (Typ[i] == FLD || Typ[i] == RIGID0 || Typ[i] == RIGID1 || Typ[i] == RIGID2){
		double mi = Dns[Typ[i]];
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double vec_ix2 = Vel[i*3  ];double vec_iy2 = Vel[i*3+1];double vec_iz2 = Vel[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dist2 = v0*v0+v1*v1+v2*v2;
				if(dist2<rlim2){
				if(j!=i && Typ[j]!=GST){
					double fDT = (vec_ix-Vel[j*3  ])*v0+(vec_iy-Vel[j*3+1])*v1+(vec_iz-Vel[j*3+2])*v2;
					if(fDT > 0.0){
						double mj = Dns[Typ[j]];
						fDT *= COL*mj/(mi+mj)/dist2;
						vec_ix2 -= v0*fDT;		vec_iy2 -= v1*fDT;		vec_iz2 -= v2*fDT;
					}
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		Acc[i*3  ]=vec_ix2;	Acc[i*3+1]=vec_iy2;	Acc[i*3+2]=vec_iz2;
	}}
	for(int i=0;i<nP;i++){
		Vel[i*3  ]=Acc[i*3  ];	Vel[i*3+1]=Acc[i*3+1];	Vel[i*3+2]=Acc[i*3+2];
	}
}

void MkPrs(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] != GST){
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double ni = 0.0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dist2 = v0*v0+v1*v1+v2*v2;
				if(dist2<r2){
				if(j!=i && Typ[j]!=GST){
					double dist = sqrt(dist2);
					double w =  WEI(dist, r);
					ni += w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		/*
		if (Pos[i * 3] <= 1.6*PCL_DST || Pos[i * 3 + 1] <= 1.6*PCL_DST || Pos[i * 3 + 2] <= 1.6*PCL_DST || Pos[i * 3 + 1] >(MAX_Y - 4.6*PCL_DST) || Pos[i * 3] > (MAX_X - 4.6*PCL_DST)){
			double PCL = double(kt / PCL_DST);
			int nx = int((MAX_X - 2.9*PCL_DST)*PCL);
			int ny = int((MAX_Y - 2.9*PCL_DST)*PCL);
			int nz = int((MAX_Z - 2.9*PCL_DST)*PCL);

			int tx = int(pos_ix *PCL + 0.0004);
			int ty = int(pos_iy *PCL + 0.0004);
			int tz = int(pos_iz *PCL + 0.0004);
			double a = pos_ix - tx *PCL_DST*pt;
			double b = (1.0 + tx)*PCL_DST*pt - pos_ix;
			double c = pos_iz - (tz + 0.0)*PCL_DST*pt;
			double d = (1.0 + tz)*PCL_DST*pt - pos_iz;
			double e = pos_iy - (ty + 0.0)*PCL_DST*pt;
			double f = (1.0 + ty)*PCL_DST*pt - pos_iy;
			int r1 = min(tx, ty, tz + 1, nx - tx, ny - ty);
			int r2 = min(tx + 1, ty, tz + 1, nx - tx - 1, ny - ty);
			int r3 = min(tx, ty, tz, nx - tx, ny - ty);
			int r4 = min(tx + 1, ty, tz, nx - tx - 1, ny - ty);
			int r5 = min(tx, ty + 1, tz + 1, nx - tx, ny - ty - 1);
			int r6 = min(tx + 1, ty + 1, tz + 1, nx - tx - 1, ny - ty - 1);
			int r7 = min(tx, ty + 1, tz, nx - tx, ny - ty - 1);
			int r8 = min(tx + 1, ty + 1, tz, nx - tx - 1, ny - ty - 1);
			double riw = double((b*c*f*r1 + a*c*f*r2 + b*d*f*r3 + a*d*f*r4 + b*c*e*r5 + a*c*e*r6 + b*d*e*r7 + a*d*e*r8) *PCL*PCL*PCL);
			double w = fi(riw);
			ni += w;
		}
		*/
		double mi = Dns[Typ[i]];
		double pressure = (ni > n0)*(ni - n0) * A2 * mi;
		Prs[i] = pressure;
	}}
}

void PrsGrdTrm(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
		if (Typ[i] == FLD || Typ[i] == RIGID0 || Typ[i] == RIGID1 || Typ[i] == RIGID2){
		double a[3][3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double Acc_1x = 0.0;
		double Acc_1y = 0.0;
		double Acc_1z = 0.0;
		double  m = 0.0;
			int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		int l = 0;
		for (int jz = iz - 1; jz <= iz + 1; jz++){
			for (int jy = iy - 1; jy <= iy + 1; jy++){
				for (int jx = ix - 1; jx <= ix + 1; jx++){
					int jb = jz*nBxy + jy*nBx + jx;
					int j = bfst[jb];
					if (j == -1) continue;
					for (;;){
						double v0 = Pos[j * 3] - pos_ix;
						double v1 = Pos[j * 3 + 1] - pos_iy;
						double v2 = Pos[j * 3 + 2] - pos_iz;
						double dist2 = v0*v0 + v1*v1 + v2*v2;
						if (dist2<r2){
							if (j != i && Typ[j] != GST){
								double dist = sqrt(dist2);
								double w = WEI(dist, r);
								/*
								m += w;

								a[0][0] += DIM * (w*v0*v0 / dist2);  a[0][1] += DIM * (w*v0*v1 / dist2); a[0][2] += DIM * (w*v0*v2 / dist2);
								a[1][0] += DIM * (w*v0*v1 / dist2);  a[1][1] += DIM * (w*v1*v1 / dist2); a[1][2] += DIM * (w*v1*v2 / dist2);
								a[2][0] += DIM * (w*v0*v2 / dist2);  a[2][1] += DIM * (w*v2*v1 / dist2); a[2][2] += DIM * (w*v2*v2 / dist2);
								l++;
								*/
								w *= (Prs[j] + Prs[i]) / dist2;
								Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
							}
						}
						j = nxt[j];
						if (j == -1) break;
					}
				}
			}
		}
		/*
		double a0; double a1; double a2; double a3; double a4; double a5; double a6; double a7; double a8;
			if (l <= 12){
		a0 = 1.0;  a1 = 0.0; a2 = 0.0;
		a3 = 0.0;  a4 = 1.0; a5 = 0.0;
		a6 = 0.0;  a7 = 0.0; a8 = 1.0;
			}
			else {
				inMatrix(a, 3);
				a0 = a[0][0] * m;  a1 = a[0][1] * m; a2 = a[0][2] * m;
				a3 = a[1][0] * m;  a4 = a[1][1] * m; a5 = a[1][2] * m;
				a6 = a[2][0] * m;  a7 = a[2][1] * m; a8 = a[2][2] * m;
			}
			
		Acc[i * 3] = (Acc_x*a0 + Acc_y*a3 + Acc_z*a6)*mi * A3;
		Acc[i * 3 + 1] = (Acc_x*a1 + Acc_y*a4 + Acc_z*a7)*mi * A3;
		Acc[i * 3 + 2] = (Acc_x*a2 + Acc_y*a5 + Acc_z*a8)*mi * A3;
		*/
								double mi = invDns[Typ[i]];
								Acc[i * 3] = Acc_x*mi * A3;
								Acc[i * 3 + 1] = Acc_y*mi * A3;
								Acc[i * 3 + 2] = Acc_z*mi * A3;
	
	}}
}

void UpPcl2(void){
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if (Typ[i] == FLD || Typ[i] == RIGID0 || Typ[i] == RIGID1 || Typ[i] == RIGID2){
			Vel[i*3  ] +=Acc[i*3  ]*DT;
			Vel[i*3+1] +=Acc[i*3+1]*DT;
			Vel[i*3+2] +=Acc[i*3+2]*DT;
			Pos[i*3  ] +=Acc[i*3  ]*DT*DT;
			Pos[i*3+1] +=Acc[i*3+1]*DT*DT;
			Pos[i*3+2] +=Acc[i*3+2]*DT*DT;
			Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
			ChkPcl(i);
		}
	}
}

void Rigid0(void){
	double CenterofGravity[3];
	double ChangeofCenterofGravity[3];
	double torque[3];
	double CG[3] = { 0.0, 0.0, 0.0 };
	double CCG[3] = { 0.0, 0.0, 0.0 };
	double tq[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID0){
			CG[0] += Pos[i * 3];
			CG[1] += Pos[i * 3 + 1];
			CG[2] += Pos[i * 3 + 2];
			CCG[0] += Pos[i * 3] - PrePos[i * 3];
			CCG[1] += Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			CCG[2] += Pos[i * 3 + 2] - PrePos[i * 3 + 2];
		}
	}

	CenterofGravity[0] = CG[0] / nr0;
	CenterofGravity[1] = CG[1] / nr0;
	CenterofGravity[2] = CG[2] / nr0;

	fopen_s(&fp, "X0.txt", "a+");
	if ((iLP*2) % OPT_FQC == 0){
		fprintf(fp, "%lf\n", CenterofGravity[0]);
	}
	fclose(fp);
	fopen_s(&fp, "z0.txt", "a+");
	if ((iLP * 2) % OPT_FQC == 0){
		fprintf(fp, "%lf\n", CenterofGravity[2]);
	}
	fclose(fp);
	
	ChangeofCenterofGravity[0] = CCG[0] / nr0;
	ChangeofCenterofGravity[1] = CCG[1] / nr0;
	ChangeofCenterofGravity[2] = CCG[2] / nr0;

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID0){

			double ChangeofPosition_0 = Pos[i * 3] - PrePos[i * 3];
			double ChangeofPosition_1 = Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			double ChangeofPosition_2 = Pos[i * 3 + 2] - PrePos[i * 3 + 2];
			double RefVec_ip0 = PrePos[i * 3] - CenterofGravity[0];
			double RefVec_ip1 = PrePos[i * 3 + 1] - CenterofGravity[1];
			double RefVec_ip2 = PrePos[i * 3 + 2] - CenterofGravity[2];
			double tmp = Dns[Typ[i]] * PCL_DST*PCL_DST*PCL_DST;
			tq[0] += tmp * (RefVec_ip1*ChangeofPosition_2 - RefVec_ip2*ChangeofPosition_1);
			tq[1] += tmp * (RefVec_ip2*ChangeofPosition_0 - RefVec_ip0*ChangeofPosition_2);
			tq[2] += tmp * (RefVec_ip0*ChangeofPosition_1 - RefVec_ip1*ChangeofPosition_0);

		}
	}
	torque[0] = tq[0]; torque[1] = tq[1]; torque[2] = tq[2];

	double ps2 = qs0;
	double px2 = qx0;
	double py2 = qy0;
	double pz2 = qz0;
	double RotationMatrix[3][3];
	RotationMatrix[0][0] = 1.0 - 2.0*(qy0*qy0 + qz0*qz0);
	RotationMatrix[0][1] = 2.0*(qx0*qy0 - qs0*qz0);
	RotationMatrix[0][2] = 2.0*(qx0*qz0 + qs0*qy0);
	RotationMatrix[1][0] = 2.0*(qx0*qy0 + qs0*qz0);
	RotationMatrix[1][1] = 1.0 - 2.0*(qx0*qx0 + qz0*qz0);
	RotationMatrix[1][2] = 2.0*(qy0*qz0 - qs0*qx0);
	RotationMatrix[2][0] = 2.0*(qx0*qz0 - qs0*qy0);
	RotationMatrix[2][1] = 2.0*(qy0*qz0 + qs0*qx0);
	RotationMatrix[2][2] = 1.0 - 2.0*(qx0*qx0 + qy0*qy0);

	double tmp1_00 = RotationMatrix[0][0] * InertiaTensorInv0[0][0] + RotationMatrix[0][1] * InertiaTensorInv0[1][0] + RotationMatrix[0][2] * InertiaTensorInv0[2][0];
	double tmp1_01 = RotationMatrix[0][0] * InertiaTensorInv0[0][1] + RotationMatrix[0][1] * InertiaTensorInv0[1][1] + RotationMatrix[0][2] * InertiaTensorInv0[2][1];
	double tmp1_02 = RotationMatrix[0][0] * InertiaTensorInv0[0][2] + RotationMatrix[0][1] * InertiaTensorInv0[1][2] + RotationMatrix[0][2] * InertiaTensorInv0[2][2];
	double tmp1_10 = RotationMatrix[1][0] * InertiaTensorInv0[0][0] + RotationMatrix[1][1] * InertiaTensorInv0[1][0] + RotationMatrix[1][2] * InertiaTensorInv0[2][0];
	double tmp1_11 = RotationMatrix[1][0] * InertiaTensorInv0[0][1] + RotationMatrix[1][1] * InertiaTensorInv0[1][1] + RotationMatrix[1][2] * InertiaTensorInv0[2][1];
	double tmp1_12 = RotationMatrix[1][0] * InertiaTensorInv0[0][2] + RotationMatrix[1][1] * InertiaTensorInv0[1][2] + RotationMatrix[1][2] * InertiaTensorInv0[2][2];
	double tmp1_20 = RotationMatrix[2][0] * InertiaTensorInv0[0][0] + RotationMatrix[2][1] * InertiaTensorInv0[1][0] + RotationMatrix[2][2] * InertiaTensorInv0[2][0];
	double tmp1_21 = RotationMatrix[2][0] * InertiaTensorInv0[0][1] + RotationMatrix[2][1] * InertiaTensorInv0[1][1] + RotationMatrix[2][2] * InertiaTensorInv0[2][1];
	double tmp1_22 = RotationMatrix[2][0] * InertiaTensorInv0[0][2] + RotationMatrix[2][1] * InertiaTensorInv0[1][2] + RotationMatrix[2][2] * InertiaTensorInv0[2][2];

	double RIRT00 = tmp1_00 * RotationMatrix[0][0] + tmp1_01 * RotationMatrix[0][1] + tmp1_02 * RotationMatrix[0][2];
	double RIRT01 = tmp1_00 * RotationMatrix[1][0] + tmp1_01 * RotationMatrix[1][1] + tmp1_02 * RotationMatrix[1][2];
	double RIRT02 = tmp1_00 * RotationMatrix[2][0] + tmp1_01 * RotationMatrix[2][1] + tmp1_02 * RotationMatrix[2][2];
	double RIRT10 = tmp1_10 * RotationMatrix[0][0] + tmp1_11 * RotationMatrix[0][1] + tmp1_12 * RotationMatrix[0][2];
	double RIRT11 = tmp1_10 * RotationMatrix[1][0] + tmp1_11 * RotationMatrix[1][1] + tmp1_12 * RotationMatrix[1][2];
	double RIRT12 = tmp1_10 * RotationMatrix[2][0] + tmp1_11 * RotationMatrix[2][1] + tmp1_12 * RotationMatrix[2][2];
	double RIRT20 = tmp1_20 * RotationMatrix[0][0] + tmp1_21 * RotationMatrix[0][1] + tmp1_22 * RotationMatrix[0][2];
	double RIRT21 = tmp1_20 * RotationMatrix[1][0] + tmp1_21 * RotationMatrix[1][1] + tmp1_22 * RotationMatrix[1][2];
	double RIRT22 = tmp1_20 * RotationMatrix[2][0] + tmp1_21 * RotationMatrix[2][1] + tmp1_22 * RotationMatrix[2][2];

	double v3_tmp0 = torque[0];
	double v3_tmp1 = torque[1];
	double v3_tmp2 = torque[2];
	double torque0 = RIRT00 * v3_tmp0 + RIRT01 * v3_tmp1 + RIRT02 * v3_tmp2;
	double torque1 = RIRT10 * v3_tmp0 + RIRT11 * v3_tmp1 + RIRT12 * v3_tmp2;
	double torque2 = RIRT20 * v3_tmp0 + RIRT21 * v3_tmp1 + RIRT22 * v3_tmp2;
	torque0 *= dt_inv;
	torque1 *= dt_inv;
	torque2 *= dt_inv;
	double torque_length = sqrt(torque0*torque0 + torque1*torque1 + torque2*torque2);
	double torque_length_inv = 1.0 / torque_length;
	//		if(pInHDDM->my_rank == pInHDDM->source){printf("iR:%d  torque_length:%20.10e \n",iR,torque_length);}

	double theta;
	double axis0;
	double axis1;
	double axis2;
	if (torque_length < 0.0001){

		theta = 0.0;
		axis0 = 0.0;
		axis1 = 0.0;
		axis2 = 0.0;
	}
	else{
		theta = torque_length * DT;
		axis0 = torque0*torque_length_inv;
		axis1 = torque1*torque_length_inv;
		axis2 = torque2*torque_length_inv;
	}

	qs0 = cos(theta*0.5);
	qx0 = axis0 * sin(theta*0.5);
	qy0 = axis1 * sin(theta*0.5);
	qz0 = axis2 * sin(theta*0.5);

	double qqs = qs0*ps2 - qx0*px2 - qy0*py2 - qz0*pz2;
	double qqx = qs0*px2 + qx0*ps2 + qy0*pz2 - qz0*py2;
	double qqy = qs0*py2 - qx0*pz2 + qy0*ps2 + qz0*px2;
	double qqz = qs0*pz2 + qx0*py2 - qy0*px2 + qz0*ps2;

	double deltaRot00 = -2.0*(qy0*qy0 + qz0*qz0);
	double deltaRot01 = 2.0*(qx0*qy0 - qs0*qz0);
	double deltaRot02 = 2.0*(qx0*qz0 + qs0*qy0);
	double deltaRot10 = 2.0*(qx0*qy0 + qs0*qz0);
	double deltaRot11 = -2.0*(qx0*qx0 + qz0*qz0);
	double deltaRot12 = 2.0*(qy0*qz0 - qs0*qx0);
	double deltaRot20 = 2.0*(qx0*qz0 - qs0*qy0);
	double deltaRot21 = 2.0*(qy0*qz0 + qs0*qx0);
	double deltaRot22 = -2.0*(qx0*qx0 + qy0*qy0);

	qs0 = qqs; qx0 = qqx; qy0=qqy ; qz0=qqz;

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID0){

			double pre_x = PrePos[i * 3];
			double pre_y = PrePos[i * 3 + 1];
			double pre_z = PrePos[i * 3 + 2];
			double RefVec_ip0 = pre_x - CenterofGravity[0];
			double RefVec_ip1 = pre_y - CenterofGravity[1];
			double RefVec_ip2 = pre_z - CenterofGravity[2];
			double Movement_ip0 = ChangeofCenterofGravity[0] + deltaRot00 * RefVec_ip0 + deltaRot01 * RefVec_ip1 + deltaRot02 * RefVec_ip2;
			double Movement_ip1 = ChangeofCenterofGravity[1] + deltaRot10 * RefVec_ip0 + deltaRot11 * RefVec_ip1 + deltaRot12 * RefVec_ip2;
			double Movement_ip2 = ChangeofCenterofGravity[2] + deltaRot20 * RefVec_ip0 + deltaRot21 * RefVec_ip1 + deltaRot22 * RefVec_ip2;

			double vec2_i0 = Movement_ip0 *dt_inv;
			double vec2_i1 = Movement_ip1 *dt_inv;
			double vec2_i2 = Movement_ip2 *dt_inv;

			Vel[i * 3] = vec2_i0;
			Vel[i * 3 + 1] = vec2_i1;
			Vel[i * 3 + 2] = vec2_i2;

			Pos[i * 3] = pre_x + Movement_ip0;
			Pos[i * 3 + 1] = pre_y + Movement_ip1;
			Pos[i * 3 + 2] = pre_z + Movement_ip2;

			PrePos[i * 3] = Pos[i * 3];
			PrePos[i * 3 + 1] = Pos[i * 3 + 1];
			PrePos[i * 3 + 2] = Pos[i * 3 + 2];

			ChkPcl(i);
		}
	}
}

void Rigid1(void){
	double CenterofGravity[3];
	double ChangeofCenterofGravity[3];
	double torque[3];
	double CG[3] = { 0.0, 0.0, 0.0 };
	double CCG[3] = { 0.0, 0.0, 0.0 };
	double tq[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID1){
			CG[0] += Pos[i * 3];
			CG[1] += Pos[i * 3 + 1];
			CG[2] += Pos[i * 3 + 2];
			CCG[0] += Pos[i * 3] - PrePos[i * 3];
			CCG[1] += Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			CCG[2] += Pos[i * 3 + 2] - PrePos[i * 3 + 2];
		}
	}

	CenterofGravity[0] = CG[0] / nr1;
	CenterofGravity[1] = CG[1] / nr1;
	CenterofGravity[2] = CG[2] / nr1;

	fopen_s(&fp, "X1.txt", "a+");
	if ((iLP * 2) % OPT_FQC == 0){
		fprintf(fp, "%lf\n", CenterofGravity[0]);
	}
	fclose(fp);
	fopen_s(&fp, "z1.txt", "a+");
	if ((iLP * 2) % OPT_FQC == 0){
		fprintf(fp, "%lf\n", CenterofGravity[2]);
	}
	fclose(fp);
	
	ChangeofCenterofGravity[0] = CCG[0] / nr1;
	ChangeofCenterofGravity[1] = CCG[1] / nr1;
	ChangeofCenterofGravity[2] = CCG[2] / nr1;

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID1){
			double ChangeofPosition_0 = Pos[i * 3] - PrePos[i * 3];
			double ChangeofPosition_1 = Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			double ChangeofPosition_2 = Pos[i * 3 + 2] - PrePos[i * 3 + 2];
			double RefVec_ip0 = PrePos[i * 3] - CenterofGravity[0];
			double RefVec_ip1 = PrePos[i * 3 + 1] - CenterofGravity[1];
			double RefVec_ip2 = PrePos[i * 3 + 2] - CenterofGravity[2];
			double tmp = Dns[Typ[i]] * PCL_DST*PCL_DST*PCL_DST;
			tq[0] += tmp * (RefVec_ip1*ChangeofPosition_2 - RefVec_ip2*ChangeofPosition_1);
			tq[1] += tmp * (RefVec_ip2*ChangeofPosition_0 - RefVec_ip0*ChangeofPosition_2);
			tq[2] += tmp * (RefVec_ip0*ChangeofPosition_1 - RefVec_ip1*ChangeofPosition_0);
		}
	}
	torque[0] = tq[0]; torque[1] = tq[1]; torque[2] = tq[2];

	double ps2 = qs1;
	double px2 = qx1;
	double py2 = qy1;
	double pz2 = qz1;
	double RotationMatrix[3][3];
	RotationMatrix[0][0] = 1.0 - 2.0*(qy1*qy1 + qz1*qz1);
	RotationMatrix[0][1] = 2.0*(qx1*qy1 - qs1*qz1);
	RotationMatrix[0][2] = 2.0*(qx1*qz1 + qs1*qy1);
	RotationMatrix[1][0] = 2.0*(qx1*qy1 + qs1*qz1);
	RotationMatrix[1][1] = 1.0 - 2.0*(qx1*qx1 + qz1*qz1);
	RotationMatrix[1][2] = 2.0*(qy1*qz1 - qs1*qx1);
	RotationMatrix[2][0] = 2.0*(qx1*qz1 - qs1*qy1);
	RotationMatrix[2][1] = 2.0*(qy1*qz1 + qs1*qx1);
	RotationMatrix[2][2] = 1.0 - 2.0*(qx1*qx1 + qy1*qy1);

	double tmp1_00 = RotationMatrix[0][0] * InertiaTensorInv1[0][0] + RotationMatrix[0][1] * InertiaTensorInv1[1][0] + RotationMatrix[0][2] * InertiaTensorInv1[2][0];
	double tmp1_01 = RotationMatrix[0][0] * InertiaTensorInv1[0][1] + RotationMatrix[0][1] * InertiaTensorInv1[1][1] + RotationMatrix[0][2] * InertiaTensorInv1[2][1];
	double tmp1_02 = RotationMatrix[0][0] * InertiaTensorInv1[0][2] + RotationMatrix[0][1] * InertiaTensorInv1[1][2] + RotationMatrix[0][2] * InertiaTensorInv1[2][2];
	double tmp1_10 = RotationMatrix[1][0] * InertiaTensorInv1[0][0] + RotationMatrix[1][1] * InertiaTensorInv1[1][0] + RotationMatrix[1][2] * InertiaTensorInv1[2][0];
	double tmp1_11 = RotationMatrix[1][0] * InertiaTensorInv1[0][1] + RotationMatrix[1][1] * InertiaTensorInv1[1][1] + RotationMatrix[1][2] * InertiaTensorInv1[2][1];
	double tmp1_12 = RotationMatrix[1][0] * InertiaTensorInv1[0][2] + RotationMatrix[1][1] * InertiaTensorInv1[1][2] + RotationMatrix[1][2] * InertiaTensorInv1[2][2];
	double tmp1_20 = RotationMatrix[2][0] * InertiaTensorInv1[0][0] + RotationMatrix[2][1] * InertiaTensorInv1[1][0] + RotationMatrix[2][2] * InertiaTensorInv1[2][0];
	double tmp1_21 = RotationMatrix[2][0] * InertiaTensorInv1[0][1] + RotationMatrix[2][1] * InertiaTensorInv1[1][1] + RotationMatrix[2][2] * InertiaTensorInv1[2][1];
	double tmp1_22 = RotationMatrix[2][0] * InertiaTensorInv1[0][2] + RotationMatrix[2][1] * InertiaTensorInv1[1][2] + RotationMatrix[2][2] * InertiaTensorInv1[2][2];

	double RIRT00 = tmp1_00 * RotationMatrix[0][0] + tmp1_01 * RotationMatrix[0][1] + tmp1_02 * RotationMatrix[0][2];
	double RIRT01 = tmp1_00 * RotationMatrix[1][0] + tmp1_01 * RotationMatrix[1][1] + tmp1_02 * RotationMatrix[1][2];
	double RIRT02 = tmp1_00 * RotationMatrix[2][0] + tmp1_01 * RotationMatrix[2][1] + tmp1_02 * RotationMatrix[2][2];
	double RIRT10 = tmp1_10 * RotationMatrix[0][0] + tmp1_11 * RotationMatrix[0][1] + tmp1_12 * RotationMatrix[0][2];
	double RIRT11 = tmp1_10 * RotationMatrix[1][0] + tmp1_11 * RotationMatrix[1][1] + tmp1_12 * RotationMatrix[1][2];
	double RIRT12 = tmp1_10 * RotationMatrix[2][0] + tmp1_11 * RotationMatrix[2][1] + tmp1_12 * RotationMatrix[2][2];
	double RIRT20 = tmp1_20 * RotationMatrix[0][0] + tmp1_21 * RotationMatrix[0][1] + tmp1_22 * RotationMatrix[0][2];
	double RIRT21 = tmp1_20 * RotationMatrix[1][0] + tmp1_21 * RotationMatrix[1][1] + tmp1_22 * RotationMatrix[1][2];
	double RIRT22 = tmp1_20 * RotationMatrix[2][0] + tmp1_21 * RotationMatrix[2][1] + tmp1_22 * RotationMatrix[2][2];

	double v3_tmp0 = torque[0];
	double v3_tmp1 = torque[1];
	double v3_tmp2 = torque[2];
	double torque0 = RIRT00 * v3_tmp0 + RIRT01 * v3_tmp1 + RIRT02 * v3_tmp2;
	double torque1 = RIRT10 * v3_tmp0 + RIRT11 * v3_tmp1 + RIRT12 * v3_tmp2;
	double torque2 = RIRT20 * v3_tmp0 + RIRT21 * v3_tmp1 + RIRT22 * v3_tmp2;
	torque0 *= dt_inv;
	torque1 *= dt_inv;
	torque2 *= dt_inv;
	double torque_length = sqrt(torque0*torque0 + torque1*torque1 + torque2*torque2);
	double torque_length_inv = 1.0 / torque_length;
	//		if(pInHDDM->my_rank == pInHDDM->source){printf("iR:%d  torque_length:%20.10e \n",iR,torque_length);}

	double theta;
	double axis0;
	double axis1;
	double axis2;
	if (torque_length < 0.0001){
		theta = 0.0;
		axis0 = 0.0;
		axis1 = 0.0;
		axis2 = 0.0;
	}
	else{
		theta = torque_length * DT;
		axis0 = torque0*torque_length_inv;
		axis1 = torque1*torque_length_inv;
		axis2 = torque2*torque_length_inv;
	}

	qs1 = cos(theta*0.5);
	qx1 = axis0 * sin(theta*0.5);
	qy1 = axis1 * sin(theta*0.5);
	qz1 = axis2 * sin(theta*0.5);

	double qqs = qs1*ps2 - qx1*px2 - qy1*py2 - qz1*pz2;
	double qqx = qs1*px2 + qx1*ps2 + qy1*pz2 - qz1*py2;
	double qqy = qs1*py2 - qx1*pz2 + qy1*ps2 + qz1*px2;
	double qqz = qs1*pz2 + qx1*py2 - qy1*px2 + qz1*ps2;

	double deltaRot00 = -2.0*(qy1*qy1 + qz1*qz1);
	double deltaRot01 = 2.0*(qx1*qy1 - qs1*qz1);
	double deltaRot02 = 2.0*(qx1*qz1 + qs1*qy1);
	double deltaRot10 = 2.0*(qx1*qy1 + qs1*qz1);
	double deltaRot11 = -2.0*(qx1*qx1 + qz1*qz1);
	double deltaRot12 = 2.0*(qy1*qz1 - qs1*qx1);
	double deltaRot20 = 2.0*(qx1*qz1 - qs1*qy1);
	double deltaRot21 = 2.0*(qy1*qz1 + qs1*qx1);
	double deltaRot22 = -2.0*(qx1*qx1 + qy1*qy1);

	qs1 = qqs; qx1 = qqx; qy1 = qqy; qz1 = qqz;

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID1){

			double pre_x = PrePos[i * 3];
			double pre_y = PrePos[i * 3 + 1];
			double pre_z = PrePos[i * 3 + 2];
			double RefVec_ip0 = pre_x - CenterofGravity[0];
			double RefVec_ip1 = pre_y - CenterofGravity[1];
			double RefVec_ip2 = pre_z - CenterofGravity[2];
			double Movement_ip0 = ChangeofCenterofGravity[0] + deltaRot00 * RefVec_ip0 + deltaRot01 * RefVec_ip1 + deltaRot02 * RefVec_ip2;
			double Movement_ip1 = ChangeofCenterofGravity[1] + deltaRot10 * RefVec_ip0 + deltaRot11 * RefVec_ip1 + deltaRot12 * RefVec_ip2;
			double Movement_ip2 = ChangeofCenterofGravity[2] + deltaRot20 * RefVec_ip0 + deltaRot21 * RefVec_ip1 + deltaRot22 * RefVec_ip2;

			double vec2_i0 = Movement_ip0 *dt_inv;
			double vec2_i1 = Movement_ip1 *dt_inv;
			double vec2_i2 = Movement_ip2 *dt_inv;

			Vel[i * 3] = vec2_i0;
			Vel[i * 3 + 1] = vec2_i1;
			Vel[i * 3 + 2] = vec2_i2;

			Pos[i * 3] = pre_x + Movement_ip0;
			Pos[i * 3 + 1] = pre_y + Movement_ip1;
			Pos[i * 3 + 2] = pre_z + Movement_ip2;

			PrePos[i * 3] = Pos[i * 3];
			PrePos[i * 3 + 1] = Pos[i * 3 + 1];
			PrePos[i * 3 + 2] = Pos[i * 3 + 2];

			ChkPcl(i);
		}
	}
}
void Rigid2(void){
	double CenterofGravity[3];
	double ChangeofCenterofGravity[3];
	double torque[3];
	double CG[3] = { 0.0, 0.0, 0.0 };
	double CCG[3] = { 0.0, 0.0, 0.0 };
	double tq[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID2){
			CG[0] += Pos[i * 3];
			CG[1] += Pos[i * 3 + 1];
			CG[2] += Pos[i * 3 + 2];
			CCG[0] += Pos[i * 3] - PrePos[i * 3];
			CCG[1] += Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			CCG[2] += Pos[i * 3 + 2] - PrePos[i * 3 + 2];
		}
	}

	CenterofGravity[0] = CG[0] / nr2;
	CenterofGravity[1] = CG[1] / nr2;
	CenterofGravity[2] = CG[2] / nr2;

	fopen_s(&fp, "X2.txt", "a+");
	if ((iLP * 2) % OPT_FQC == 0){
		fprintf(fp, "%lf\n", CenterofGravity[0]);
	}
	fclose(fp);
	fopen_s(&fp, "z2.txt", "a+");
	if ((iLP * 2) % OPT_FQC == 0){
		fprintf(fp, "%lf\n", CenterofGravity[2]);
	}
	fclose(fp);

	ChangeofCenterofGravity[0] = CCG[0] / nr2;
	ChangeofCenterofGravity[1] = CCG[1] / nr2;
	ChangeofCenterofGravity[2] = CCG[2] / nr2;

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID2){
			double ChangeofPosition_0 = Pos[i * 3] - PrePos[i * 3];
			double ChangeofPosition_1 = Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			double ChangeofPosition_2 = Pos[i * 3 + 2] - PrePos[i * 3 + 2];
			double RefVec_ip0 = PrePos[i * 3] - CenterofGravity[0];
			double RefVec_ip1 = PrePos[i * 3 + 1] - CenterofGravity[1];
			double RefVec_ip2 = PrePos[i * 3 + 2] - CenterofGravity[2];
			double tmp = Dns[Typ[i]]*PCL_DST*PCL_DST*PCL_DST;
			tq[0] += tmp * (RefVec_ip1*ChangeofPosition_2 - RefVec_ip2*ChangeofPosition_1);
			tq[1] += tmp * (RefVec_ip2*ChangeofPosition_0 - RefVec_ip0*ChangeofPosition_2);
			tq[2] += tmp * (RefVec_ip0*ChangeofPosition_1 - RefVec_ip1*ChangeofPosition_0);

		}
	}

	torque[0] = tq[0]; torque[1] = tq[1]; torque[2] = tq[2];
	double ps2 = qs2;
	double px2 = qx2;
	double py2 = qy2;
	double pz2 = qz2;
	double RotationMatrix[3][3];
	RotationMatrix[0][0] = 1.0 - 2.0*(qy2*qy2 + qz2*qz2);
	RotationMatrix[0][1] = 2.0*(qx2*qy2 - qs2*qz2);
	RotationMatrix[0][2] = 2.0*(qx2*qz2 + qs2*qy2);
	RotationMatrix[1][0] = 2.0*(qx2*qy2 + qs2*qz2);
	RotationMatrix[1][1] = 1.0 - 2.0*(qx2*qx2 + qz2*qz2);
	RotationMatrix[1][2] = 2.0*(qy2*qz2 - qs2*qx2);
	RotationMatrix[2][0] = 2.0*(qx2*qz2 - qs2*qy2);
	RotationMatrix[2][1] = 2.0*(qy2*qz2 + qs2*qx2);
	RotationMatrix[2][2] = 1.0 - 2.0*(qx2*qx2 + qy2*qy2);

	double tmp1_00 = RotationMatrix[0][0] * InertiaTensorInv2[0][0] + RotationMatrix[0][1] * InertiaTensorInv2[1][0] + RotationMatrix[0][2] * InertiaTensorInv2[2][0];
	double tmp1_01 = RotationMatrix[0][0] * InertiaTensorInv2[0][1] + RotationMatrix[0][1] * InertiaTensorInv2[1][1] + RotationMatrix[0][2] * InertiaTensorInv2[2][1];
	double tmp1_02 = RotationMatrix[0][0] * InertiaTensorInv2[0][2] + RotationMatrix[0][1] * InertiaTensorInv2[1][2] + RotationMatrix[0][2] * InertiaTensorInv2[2][2];
	double tmp1_10 = RotationMatrix[1][0] * InertiaTensorInv2[0][0] + RotationMatrix[1][1] * InertiaTensorInv2[1][0] + RotationMatrix[1][2] * InertiaTensorInv2[2][0];
	double tmp1_11 = RotationMatrix[1][0] * InertiaTensorInv2[0][1] + RotationMatrix[1][1] * InertiaTensorInv2[1][1] + RotationMatrix[1][2] * InertiaTensorInv2[2][1];
	double tmp1_12 = RotationMatrix[1][0] * InertiaTensorInv2[0][2] + RotationMatrix[1][1] * InertiaTensorInv2[1][2] + RotationMatrix[1][2] * InertiaTensorInv2[2][2];
	double tmp1_20 = RotationMatrix[2][0] * InertiaTensorInv2[0][0] + RotationMatrix[2][1] * InertiaTensorInv2[1][0] + RotationMatrix[2][2] * InertiaTensorInv2[2][0];
	double tmp1_21 = RotationMatrix[2][0] * InertiaTensorInv2[0][1] + RotationMatrix[2][1] * InertiaTensorInv2[1][1] + RotationMatrix[2][2] * InertiaTensorInv2[2][1];
	double tmp1_22 = RotationMatrix[2][0] * InertiaTensorInv2[0][2] + RotationMatrix[2][1] * InertiaTensorInv2[1][2] + RotationMatrix[2][2] * InertiaTensorInv2[2][2];

	double RIRT00 = tmp1_00 * RotationMatrix[0][0] + tmp1_01 * RotationMatrix[0][1] + tmp1_02 * RotationMatrix[0][2];
	double RIRT01 = tmp1_00 * RotationMatrix[1][0] + tmp1_01 * RotationMatrix[1][1] + tmp1_02 * RotationMatrix[1][2];
	double RIRT02 = tmp1_00 * RotationMatrix[2][0] + tmp1_01 * RotationMatrix[2][1] + tmp1_02 * RotationMatrix[2][2];
	double RIRT10 = tmp1_10 * RotationMatrix[0][0] + tmp1_11 * RotationMatrix[0][1] + tmp1_12 * RotationMatrix[0][2];
	double RIRT11 = tmp1_10 * RotationMatrix[1][0] + tmp1_11 * RotationMatrix[1][1] + tmp1_12 * RotationMatrix[1][2];
	double RIRT12 = tmp1_10 * RotationMatrix[2][0] + tmp1_11 * RotationMatrix[2][1] + tmp1_12 * RotationMatrix[2][2];
	double RIRT20 = tmp1_20 * RotationMatrix[0][0] + tmp1_21 * RotationMatrix[0][1] + tmp1_22 * RotationMatrix[0][2];
	double RIRT21 = tmp1_20 * RotationMatrix[1][0] + tmp1_21 * RotationMatrix[1][1] + tmp1_22 * RotationMatrix[1][2];
	double RIRT22 = tmp1_20 * RotationMatrix[2][0] + tmp1_21 * RotationMatrix[2][1] + tmp1_22 * RotationMatrix[2][2];

	double v3_tmp0 = torque[0];
	double v3_tmp1 = torque[1];
	double v3_tmp2 = torque[2];
	double torque0 = RIRT00 * v3_tmp0 + RIRT01 * v3_tmp1 + RIRT02 * v3_tmp2;
	double torque1 = RIRT10 * v3_tmp0 + RIRT11 * v3_tmp1 + RIRT12 * v3_tmp2;
	double torque2 = RIRT20 * v3_tmp0 + RIRT21 * v3_tmp1 + RIRT22 * v3_tmp2;
	torque0 *= dt_inv;
	torque1 *= dt_inv;
	torque2 *= dt_inv;
	double torque_length = sqrt(torque0*torque0 + torque1*torque1 + torque2*torque2);
	double torque_length_inv = 1.0 / torque_length;
	//		if(pInHDDM->my_rank == pInHDDM->source){printf("iR:%d  torque_length:%20.10e \n",iR,torque_length);}

	double theta;
	double axis0;
	double axis1;
	double axis2;
	if (torque_length < 0.0001){

		theta = 0.0;
		axis0 = 0.0;
		axis1 = 0.0;
		axis2 = 0.0;
	}
	else{
		theta = torque_length * DT;
		axis0 = torque0*torque_length_inv;
		axis1 = torque1*torque_length_inv;
		axis2 = torque2*torque_length_inv;
	}

	qs2 = cos(theta*0.5);
	qx2 = axis0 * sin(theta*0.5);
	qy2 = axis1 * sin(theta*0.5);
	qz2 = axis2 * sin(theta*0.5);

	double qqs = qs2*ps2 - qx2*px2 - qy2*py2 - qz2*pz2;
	double qqx = qs2*px2 + qx2*ps2 + qy2*pz2 - qz2*py2;
	double qqy = qs2*py2 - qx2*pz2 + qy2*ps2 + qz2*px2;
	double qqz = qs2*pz2 + qx2*py2 - qy2*px2 + qz2*ps2;

	double deltaRot00 = -2.0*(qy2*qy2 + qz2*qz2);
	double deltaRot01 = 2.0*(qx2*qy2 - qs2*qz2);
	double deltaRot02 = 2.0*(qx2*qz2 + qs2*qy2);
	double deltaRot10 = 2.0*(qx2*qy2 + qs2*qz2);
	double deltaRot11 = -2.0*(qx2*qx2 + qz2*qz2);
	double deltaRot12 = 2.0*(qy2*qz2 - qs2*qx2);
	double deltaRot20 = 2.0*(qx2*qz2 - qs2*qy2);
	double deltaRot21 = 2.0*(qy2*qz2 + qs2*qx2);
	double deltaRot22 = -2.0*(qx2*qx2 + qy2*qy2);

	qs2 = qqs; qx2 = qqx; qy2 = qqy; qz2 = qqz;

	for (int i = 0; i < nP; i++){
		if (Typ[i] == RIGID2){

			double pre_x = PrePos[i * 3];
			double pre_y = PrePos[i * 3 + 1];
			double pre_z = PrePos[i * 3 + 2];
			double RefVec_ip0 = pre_x - CenterofGravity[0];
			double RefVec_ip1 = pre_y - CenterofGravity[1];
			double RefVec_ip2 = pre_z - CenterofGravity[2];
			double Movement_ip0 = ChangeofCenterofGravity[0] + deltaRot00 * RefVec_ip0 + deltaRot01 * RefVec_ip1 + deltaRot02 * RefVec_ip2;
			double Movement_ip1 = ChangeofCenterofGravity[1] + deltaRot10 * RefVec_ip0 + deltaRot11 * RefVec_ip1 + deltaRot12 * RefVec_ip2;
			double Movement_ip2 = ChangeofCenterofGravity[2] + deltaRot20 * RefVec_ip0 + deltaRot21 * RefVec_ip1 + deltaRot22 * RefVec_ip2;

			double vec2_i0 = Movement_ip0 *dt_inv;
			double vec2_i1 = Movement_ip1 *dt_inv;
			double vec2_i2 = Movement_ip2 *dt_inv;

			Vel[i * 3] = vec2_i0;
			Vel[i * 3 + 1] = vec2_i1;
			Vel[i * 3 + 2] = vec2_i2;

			Pos[i * 3] = pre_x + Movement_ip0;
			Pos[i * 3 + 1] = pre_y + Movement_ip1;
			Pos[i * 3 + 2] = pre_z + Movement_ip2;

			PrePos[i * 3] = Pos[i * 3];
			PrePos[i * 3 + 1] = Pos[i * 3 + 1];
			PrePos[i * 3 + 2] = Pos[i * 3 + 2];

			ChkPcl(i);
		}
	}
}

void wdataG1(void){
	fopen_s(&fp, "G1.txt", "a+");
	if (iLP % OPT_FQC == 0){
		double maxd2 = 0.0;
		for (int k = 0; k < nP; k++){
			if (Typ[k] == FLD &&  Pos[k * 3] >= 2.5 - 5 * PCL_DST  &&Pos[k * 3] <= 2.5 &&Pos[k * 3 + 1] <= 0.13 + 15 * PCL_DST && Pos[k * 3 + 1] >= 0.13 - 15 * PCL_DST){
				if (Prs[k] >= 0.0 && Pos[k * 3 + 2] > maxd2)
					maxd2 = Pos[k * 3 + 2];
				else continue;
			}
			else continue;
		}
		double maxd3 = 0.0;
		for (int k = 0; k < nP; k++){
			if (Typ[k] == FLD &&  Pos[k * 3] >= 2.5 - 5 * PCL_DST  &&Pos[k * 3] <= 2.5 &&Pos[k * 3 + 1] <= 0.13 + 15 * PCL_DST && Pos[k * 3 + 1] >= 0.13 - 15 * PCL_DST&&Pos[k * 3 + 2] <= maxd2 + 4 * PCL_DST){
				if (Prs[k] >= 0.0 && Pos[k * 3 + 2] > maxd3)
					maxd3 = Pos[k * 3 + 2];
				else continue;
			}
			else continue;
		}
		fprintf(fp, "%lf\n", maxd3);
	}
	fclose(fp);

	fopen_s(&fp, "G2.txt", "a+");
	if (iLP % OPT_FQC == 0){
		double maxd2 = 0.0;
		for (int k = 0; k < nP; k++){
			if (Typ[k] == FLD &&  Pos[k * 3] >= 2.95 - 5 * PCL_DST  &&Pos[k * 3] <= 2.95 &&Pos[k * 3 + 1] <= 0.5 - 0.13 + 15 * PCL_DST && Pos[k * 3 + 1] >= 0.5 - 0.13 - 15 * PCL_DST){
				if (Prs[k] >= 0.0 && Pos[k * 3 + 2] > maxd2)
					maxd2 = Pos[k * 3 + 2];
				else continue;
			}
			else continue;
		}
		double maxd3 = 0.0;
		for (int k = 0; k < nP; k++){
			if (Typ[k] == FLD &&  Pos[k * 3] >= 2.95 - 5 * PCL_DST  &&Pos[k * 3] <= 2.95 &&Pos[k * 3 + 1] <= 00.5 - 0.13 + 15 * PCL_DST && Pos[k * 3 + 1] >= 0.5 - 0.13 - 15 * PCL_DST&&Pos[k * 3 + 2] <= maxd2 + 4 * PCL_DST){
				if (Prs[k] >= 0.0 && Pos[k * 3 + 2] > maxd3)
					maxd3 = Pos[k * 3 + 2];
				else continue;
			}
			else continue;
		}
		fprintf(fp, "%lf\n", maxd3);
	}
	fclose(fp);
}

void GALIFT(void){
	if (TIM >= 1.0){
		for (int i = 0; i < nP; i++){
			if (Typ[i] == GATE){
				Pos[i * 3 + 2] += 0.11*DT;
			}
		}
	}
}
void ClcEMPS(void){
	while(1){
		if (iLP%OPT_FQC == 0){
			int p_num=0;
			for(int i=0;i<nP;i++){if(Typ[i] != GST)p_num++;}
			printf("%5d th TIM: %lf / p_num: %d\n", iLP,TIM,p_num);
		}
		if(iLP%OPT_FQC == 0 ){
			WrtDat();
			if(TIM >= FIN_TIM ){break;}
		}

		MkBkt();		
		VscTrm();
		UpPcl1();
		ChkCol();
		MkPrs();
		PrsGrdTrm();
		UpPcl2();
		MkPrs();
		Rigid0();
//		Rigid1();
//		Rigid2();
		wdataG1();
		GALIFT();

		for (int i=0;i<nP;i++){pav[i] += Prs[i];}
		iLP++;
		TIM += DT;
	}
}

int main( int argc, char** argv) {
	printf("start emps.\n");
	RdDat();
	AlcBkt();
	SetPara();
	init_rigid0();
//	init_rigid1();
//	init_rigid2();

	ClcEMPS();

	free(Acc);	free(Pos);	free(Vel);
	free(Prs);	free(pav);	free(Typ);
	free(bfst);	free(blst);	free(nxt);
	printf("end emps.\n");
	return 0;
}
