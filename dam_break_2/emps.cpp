#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <io.h>
#include <time.h>
#include <malloc.h>
#include "../params.h"

#define IN_FILE "../data/init_position.prof"
#define IN_DIR_VTU "../data"
#define OUT_DIR_VTU "../data"

/*
vtuへの出力はファイル分割するべき
*/

/*
疑問点
剛体の質量、動粘度の扱い？(現コードでは異なる質量、同一動粘度)
preaveの扱いが異なる

AlcBkt(void):クーラン数だけバケット長の延長
init_rigid0(void):論文と異なる剛体の扱い
VscTrm():動粘性係数の調和平均
ChkCol():衝突の際の運動量の交換、修正方法
Presmoothing():なぜ必要
PrsGrdTrm():なぜ違うgrad-modelにしたのか
*/

FILE* fp;
char filename[256];
int iLP,iF;
double TIM;
int nP,nr0,nr1,nr2;
double *Acc,*Pos,*Vel,*Prs,*pav;
double *PrePos;
double qs0, qx0, qy0, qz0;
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
double alcSize;

// vtuへの変換で利用
int NumberOfParticle;
int ProfFileNumber;
double *Position;
double *Velocity;
double *Pressure;
double *pressave;
int *ParticleType;

void ChkPcl(int i) {
	if (Pos[i * 3]>MAX_X || Pos[i * 3]<MIN_X ||
		Pos[i * 3 + 1]>MAX_Y || Pos[i * 3 + 1]<MIN_Y ||
		Pos[i * 3 + 2]>MAX_Z || Pos[i * 3 + 2]<MIN_Z)
	{
		Typ[i] = GHOST;
		Prs[i] = Vel[i * 3] = Vel[i * 3 + 1] = Vel[i * 3 + 2] = 0.0;
	}
}

void RdDat(void) {
	fopen_s(&fp, IN_FILE, "r");
	fscanf_s(fp, "%d %d", &nP, &nr0);
	printf("nP: %d\n", nP);
	alcSize = 0;
	Acc = (double*)malloc(sizeof(double)*nP * 3);	alcSize += _msize(Acc);		//粒子の加速度
	Pos = (double*)malloc(sizeof(double)*nP * 3);	alcSize += _msize(Pos);		//粒子の座標
	PrePos = (double*)malloc(sizeof(double)*nP * 3); alcSize += _msize(PrePos);	//粒子の座標
	Vel = (double*)malloc(sizeof(double)*nP * 3);	alcSize += _msize(Vel);		//粒子の速度
	Prs = (double*)malloc(sizeof(double)*nP); alcSize += _msize(Prs);		//粒子の圧力
	pav = (double*)malloc(sizeof(double)*nP); alcSize += _msize(pav);		//時間平均された粒子の圧力
	Typ = (int*)malloc(sizeof(int)*nP);	alcSize += _msize(Typ);				//粒子の種類
	for (int i = 0; i<nP; i++) {
		int a[2];
		double b[8];
		fscanf_s(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7]);
		Typ[i] = a[1];
		Pos[i * 3] = b[0];	Pos[i * 3 + 1] = b[1];	Pos[i * 3 + 2] = b[2];
		Vel[i * 3] = b[3];	Vel[i * 3 + 1] = b[4];	Vel[i * 3 + 2] = b[5];
		Prs[i] = b[6];		pav[i] = b[7];
	}
	fclose(fp);
	for (int i = 0; i<nP; i++) { ChkPcl(i); }
	for (int i = 0; i<nP * 3; i++) { Acc[i] = 0.0; }
	for (int i = 0; i<nP * 3; i++) { PrePos[i] = Pos[i]; }
}

void WrtDat(void) {
	char outout_filename[256];
	sprintf_s(outout_filename, "../data/output%05d.prof", iF);
	fopen_s(&fp, outout_filename, "w");
	fprintf(fp, "%d\n", nP);
	for (int i = 0; i<nP; i++) {
		int a[2];
		double b[8];
		a[0] = i;	a[1] = Typ[i];
		b[0] = Pos[i * 3];	b[1] = Pos[i * 3 + 1];	b[2] = Pos[i * 3 + 2];
		b[3] = Vel[i * 3];	b[4] = Vel[i * 3 + 1];	b[5] = Vel[i * 3 + 2];
		b[6] = Prs[i];		b[7] = pav[i] / OPT_FQC;
		if (Typ[i] == GHOST)continue;
		fprintf(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", a[0], a[1], b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
		pav[i] = 0.0;
	}
	fclose(fp);
	iF++;
}

void AlcBkt(void) {
	r = PARTICLE_DISTANCE*3.1;		//影響半径
	r2 = r*r;
	DB = r*(1.0 + CRT_NUM);	//バケット1辺の長さ
	DB2 = DB*DB;
	DBinv = 1.0 / DB;
	nBx = (int)((MAX_X - MIN_X)*DBinv) + 3;//解析領域内のx方向のバケット数
	nBy = (int)((MAX_Y - MIN_Y)*DBinv) + 3;//解析領域内のy方向のバケット数
	nBz = (int)((MAX_Z - MIN_Z)*DBinv) + 3;//解析領域内のz方向のバケット数
	nBxy = nBx*nBy;
	nBxyz = nBx*nBy*nBz;		//解析領域内のバケット数
	printf("nBx:%d  nBy:%d  nBz:%d  nBxy:%d  nBxyz:%d\n", nBx, nBy, nBz, nBxy, nBxyz);
	bfst = (int*)malloc(sizeof(int) * nBxyz);	alcSize += _msize(bfst);		//バケットに格納された先頭の粒子番号
	blst = (int*)malloc(sizeof(int) * nBxyz);	alcSize += _msize(blst);		//バケットに格納された最後尾の粒子番号
	nxt = (int*)malloc(sizeof(int) * nP);		alcSize += _msize(nxt);		//同じバケット内の次の粒子番号
	printf("allocate %f GB in heap\n", alcSize / (1024 * 1024 * 1024));
}

void SetPara(void) {
	double tn0 = 0.0;
	double tlmd = 0.0;
	for (int ix = -4; ix < 5; ix++) {
		for (int iy = -4; iy < 5; iy++) {
			for (int iz = -4; iz < 5; iz++) {
				double x = PARTICLE_DISTANCE* (double)ix;
				double y = PARTICLE_DISTANCE* (double)iy;
				double z = PARTICLE_DISTANCE* (double)iz;
				double dist2 = x*x + y*y + z*z;
				if (dist2 <= r2) {
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
	A1 = 2.0*KNM_VSC_FRUID*DIM / n0 / lmd;//粘性項の計算に用いる係数
	A2 = SND*SND / n0;				//圧力の計算に用いる係数
	A3 = -DIM / n0;					//圧力勾配項の計算に用いる係数
	Dns[FLUID] = DNS_FLUID;
	Dns[WALL] = DNS_WALL;
	Dns[RIGID0] = DNS_RIGID0;
	invDns[FLUID] = 1.0 / DNS_FLUID;
	invDns[WALL] = 1.0 / DNS_WALL;
	invDns[RIGID0] = 1.0 / DNS_RIGID0;
	rlim = PARTICLE_DISTANCE * DST_LMT_RAT;//これ以上の粒子間の接近を許さない距離
	rlim2 = rlim*rlim;
	COL = 1.0 + COL_RAT;
	iLP = 0;			//反復数
	iF = 0;			//ファイル番号
	TIM = 0.0;		//時刻
}

void init_rigid0(void) {

	double CenterofGravity[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++) {
		if (Typ[i] == RIGID0) {
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

	for (int i = 0; i < nP; i++) {
		if (Typ[i] == RIGID0) {
			double  v3_tmp0 = Pos[3 * i] - CenterofGravity[0];
			double  v3_tmp1 = Pos[3 * i + 1] - CenterofGravity[1];
			double  v3_tmp2 = Pos[3 * i + 2] - CenterofGravity[2];
			double v3_tmp_squaredLength = v3_tmp0*v3_tmp0 + v3_tmp1*v3_tmp1 + v3_tmp2*v3_tmp2;
			double dd0 = Dns[Typ[i]] * PARTICLE_DISTANCE*PARTICLE_DISTANCE*PARTICLE_DISTANCE;
			double dd1 = v3_tmp_squaredLength + PARTICLE_DISTANCE*PARTICLE_DISTANCE*0.1666666666666666666666666;

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
	InertiaTensorInv0[0][0] = A00; InertiaTensorInv0[0][1] = A01; InertiaTensorInv0[0][2] = A02;
	InertiaTensorInv0[1][0] = A10; InertiaTensorInv0[1][1] = A11; InertiaTensorInv0[1][2] = A12;
	InertiaTensorInv0[2][0] = A20; InertiaTensorInv0[2][1] = A21; InertiaTensorInv0[2][2] = A22;
}

void MkBkt(void) {
	for (int i = 0; i< nBxyz; i++) { bfst[i] = -1; }
	for (int i = 0; i< nBxyz; i++) { blst[i] = -1; }
	for (int i = 0; i< nP; i++) { nxt[i] = -1; }
	for (int i = 0; i<nP; i++) {
		if (Typ[i] == GHOST)continue;
		int ix = (int)((Pos[i * 3] - MIN_X)*DBinv) + 1;
		int iy = (int)((Pos[i * 3 + 1] - MIN_Y)*DBinv) + 1;
		int iz = (int)((Pos[i * 3 + 2] - MIN_Z)*DBinv) + 1;
		int ib = iz*nBxy + iy*nBx + ix;
		int j = blst[ib];
		blst[ib] = i;
		if (j == -1) { bfst[ib] = i; }
		else { nxt[j] = i; }
	}
}

void VscTrm() {
#pragma omp parallel for schedule(dynamic,64)
	for (int i = 0; i<nP; i++) {
		if (Typ[i] == FLUID || Typ[i] == RIGID0) {
			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i * 3];	double pos_iy = Pos[i * 3 + 1];	double pos_iz = Pos[i * 3 + 2];
			double vec_ix = Vel[i * 3];	double vec_iy = Vel[i * 3 + 1];	double vec_iz = Vel[i * 3 + 2];
			int ix = (int)((pos_ix - MIN_X)*DBinv) + 1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) + 1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz*nBxy + jy*nBx + jx;
						int j = bfst[jb];
						if (j == -1) continue;
						for (;;) {
							double v0 = Pos[j * 3] - pos_ix;
							double v1 = Pos[j * 3 + 1] - pos_iy;
							double v2 = Pos[j * 3 + 2] - pos_iz;
							double dist2 = v0*v0 + v1*v1 + v2*v2;
							if (dist2<r2) {
								if (j != i && Typ[j] != GHOST) {
									double dist = sqrt(dist2);
									double w = WEI(dist, r);
									Acc_x += (Vel[j * 3] - vec_ix)*w;
									Acc_y += (Vel[j * 3 + 1] - vec_iy)*w;
									Acc_z += (Vel[j * 3 + 2] - vec_iz)*w;
								}
							}
							j = nxt[j];
							if (j == -1) break;
						}
					}
				}
			}
			Acc[i * 3] = Acc_x*A1 + G_X;
			Acc[i * 3 + 1] = Acc_y*A1 + G_Y;
			Acc[i * 3 + 2] = Acc_z*A1 + G_Z;
		}
	}
}

void UpPcl1() {
	for (int i = 0; i<nP; i++) {
		if (Typ[i] == FLUID || Typ[i] == RIGID0) {
			Vel[i * 3] += Acc[i * 3] * DT;	Vel[i * 3 + 1] += Acc[i * 3 + 1] * DT;	Vel[i * 3 + 2] += Acc[i * 3 + 2] * DT;
			Pos[i * 3] += Vel[i * 3] * DT;		Pos[i * 3 + 1] += Vel[i * 3 + 1] * DT;		Pos[i * 3 + 2] += Vel[i * 3 + 2] * DT;
			Acc[i * 3] = Acc[i * 3 + 1] = Acc[i * 3 + 2] = 0.0;
			ChkPcl(i);
		}
	}
}

void ChkCol() {
#pragma omp parallel for
	for (int i = 0; i<nP; i++) {
		if (Typ[i] == FLUID || Typ[i] == RIGID0) {
			double mi = Dns[Typ[i]];
			double pos_ix = Pos[i * 3];	double pos_iy = Pos[i * 3 + 1];	double pos_iz = Pos[i * 3 + 2];
			double vec_ix = Vel[i * 3];	double vec_iy = Vel[i * 3 + 1];	double vec_iz = Vel[i * 3 + 2];
			double vec_ix2 = Vel[i * 3]; double vec_iy2 = Vel[i * 3 + 1]; double vec_iz2 = Vel[i * 3 + 2];
			int ix = (int)((pos_ix - MIN_X)*DBinv) + 1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) + 1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz*nBxy + jy*nBx + jx;
						int j = bfst[jb];
						if (j == -1) continue;
						for (;;) {
							double v0 = Pos[j * 3] - pos_ix;
							double v1 = Pos[j * 3 + 1] - pos_iy;
							double v2 = Pos[j * 3 + 2] - pos_iz;
							double dist2 = v0*v0 + v1*v1 + v2*v2;
							if (dist2<rlim2) {
								if (j != i && Typ[j] != GHOST) {
									double fDT = (vec_ix - Vel[j * 3])*v0 + (vec_iy - Vel[j * 3 + 1])*v1 + (vec_iz - Vel[j * 3 + 2])*v2;
									if (fDT > 0.0) {
										double mj = Dns[Typ[j]];
										fDT *= COL*mj / (mi + mj) / dist2;
										vec_ix2 -= v0*fDT;		vec_iy2 -= v1*fDT;		vec_iz2 -= v2*fDT;
									}
								}
							}
							j = nxt[j];
							if (j == -1) break;
						}
					}
				}
			}
			Acc[i * 3] = vec_ix2;	Acc[i * 3 + 1] = vec_iy2;	Acc[i * 3 + 2] = vec_iz2;
		}
	}
	for (int i = 0; i<nP; i++) {
		Vel[i * 3] = Acc[i * 3];	Vel[i * 3 + 1] = Acc[i * 3 + 1];	Vel[i * 3 + 2] = Acc[i * 3 + 2];
	}
}

void MkPrs() {
#pragma omp parallel for schedule(dynamic,64)
	for (int i = 0; i<nP; i++) {
		if (Typ[i] != GHOST) {
			double pos_ix = Pos[i * 3];	double pos_iy = Pos[i * 3 + 1];	double pos_iz = Pos[i * 3 + 2];
			double ni = 0.0;
			int ix = (int)((pos_ix - MIN_X)*DBinv) + 1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) + 1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz*nBxy + jy*nBx + jx;
						int j = bfst[jb];
						if (j == -1) continue;
						for (;;) {
							double v0 = Pos[j * 3] - pos_ix;
							double v1 = Pos[j * 3 + 1] - pos_iy;
							double v2 = Pos[j * 3 + 2] - pos_iz;
							double dist2 = v0*v0 + v1*v1 + v2*v2;
							if (dist2<r2) {
								if (j != i && Typ[j] != GHOST) {
									double dist = sqrt(dist2);
									double w = WEI(dist, r);
									ni += w;
								}
							}
							j = nxt[j];
							if (j == -1) break;
						}
					}
				}
			}
			double mi = Dns[Typ[i]];
			double pressure = (ni > n0)*(ni - n0) * A2 * mi;
			Prs[i] = pressure;
		}
	}
}

void PrsGrdTrm() {
#pragma omp parallel for schedule(dynamic,64)
	for (int i = 0; i<nP; i++) {
		if (Typ[i] == FLUID || Typ[i] == RIGID0) {
			double a[3][3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i * 3];	double pos_iy = Pos[i * 3 + 1];	double pos_iz = Pos[i * 3 + 2];
			double Acc_1x = 0.0;
			double Acc_1y = 0.0;
			double Acc_1z = 0.0;
			double  m = 0.0;
			int ix = (int)((pos_ix - MIN_X)*DBinv) + 1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) + 1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) + 1;
			int l = 0;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz*nBxy + jy*nBx + jx;
						int j = bfst[jb];
						if (j == -1) continue;
						for (;;) {
							double v0 = Pos[j * 3] - pos_ix;
							double v1 = Pos[j * 3 + 1] - pos_iy;
							double v2 = Pos[j * 3 + 2] - pos_iz;
							double dist2 = v0*v0 + v1*v1 + v2*v2;
							if (dist2<r2) {
								if (j != i && Typ[j] != GHOST) {
									double dist = sqrt(dist2);
									double w = WEI(dist, r);
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
			double mi = invDns[Typ[i]];
			Acc[i * 3] = Acc_x*mi * A3;
			Acc[i * 3 + 1] = Acc_y*mi * A3;
			Acc[i * 3 + 2] = Acc_z*mi * A3;
		}
	}
}

void UpPcl2(void) {
#pragma omp parallel for
	for (int i = 0; i<nP; i++) {
		if (Typ[i] == FLUID || Typ[i] == RIGID0) {
			Vel[i * 3] += Acc[i * 3] * DT;
			Vel[i * 3 + 1] += Acc[i * 3 + 1] * DT;
			Vel[i * 3 + 2] += Acc[i * 3 + 2] * DT;
			Pos[i * 3] += Acc[i * 3] * DT*DT;
			Pos[i * 3 + 1] += Acc[i * 3 + 1] * DT*DT;
			Pos[i * 3 + 2] += Acc[i * 3 + 2] * DT*DT;
			Acc[i * 3] = Acc[i * 3 + 1] = Acc[i * 3 + 2] = 0.0;
			ChkPcl(i);
		}
	}
}

void Rigid0(void) {
	double CenterofGravity[3];
	double ChangeofCenterofGravity[3];
	double torque[3];
	double CG[3] = { 0.0, 0.0, 0.0 };
	double CCG[3] = { 0.0, 0.0, 0.0 };
	double tq[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < nP; i++) {
		if (Typ[i] == RIGID0) {
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

	/*
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
	*/

	ChangeofCenterofGravity[0] = CCG[0] / nr0;
	ChangeofCenterofGravity[1] = CCG[1] / nr0;
	ChangeofCenterofGravity[2] = CCG[2] / nr0;

	for (int i = 0; i < nP; i++) {
		if (Typ[i] == RIGID0) {
			double ChangeofPosition_0 = Pos[i * 3] - PrePos[i * 3];
			double ChangeofPosition_1 = Pos[i * 3 + 1] - PrePos[i * 3 + 1];
			double ChangeofPosition_2 = Pos[i * 3 + 2] - PrePos[i * 3 + 2];
			double RefVec_ip0 = PrePos[i * 3] - CenterofGravity[0];
			double RefVec_ip1 = PrePos[i * 3 + 1] - CenterofGravity[1];
			double RefVec_ip2 = PrePos[i * 3 + 2] - CenterofGravity[2];
			double tmp = Dns[Typ[i]] * PARTICLE_DISTANCE*PARTICLE_DISTANCE*PARTICLE_DISTANCE;
			//? 単位おかしくね？
			tq[0] += tmp * (RefVec_ip1*ChangeofPosition_2 - RefVec_ip2*ChangeofPosition_1);
			tq[1] += tmp * (RefVec_ip2*ChangeofPosition_0 - RefVec_ip0*ChangeofPosition_2);
			tq[2] += tmp * (RefVec_ip0*ChangeofPosition_1 - RefVec_ip1*ChangeofPosition_0);
		}
	}
	torque[0] = tq[0]; torque[1] = tq[1]; torque[2] = tq[2];

	// I(t)^-1を計算
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
	// temp=R*I(init)^-1
	double tmp1_00 = RotationMatrix[0][0] * InertiaTensorInv0[0][0] + RotationMatrix[0][1] * InertiaTensorInv0[1][0] + RotationMatrix[0][2] * InertiaTensorInv0[2][0];
	double tmp1_01 = RotationMatrix[0][0] * InertiaTensorInv0[0][1] + RotationMatrix[0][1] * InertiaTensorInv0[1][1] + RotationMatrix[0][2] * InertiaTensorInv0[2][1];
	double tmp1_02 = RotationMatrix[0][0] * InertiaTensorInv0[0][2] + RotationMatrix[0][1] * InertiaTensorInv0[1][2] + RotationMatrix[0][2] * InertiaTensorInv0[2][2];
	double tmp1_10 = RotationMatrix[1][0] * InertiaTensorInv0[0][0] + RotationMatrix[1][1] * InertiaTensorInv0[1][0] + RotationMatrix[1][2] * InertiaTensorInv0[2][0];
	double tmp1_11 = RotationMatrix[1][0] * InertiaTensorInv0[0][1] + RotationMatrix[1][1] * InertiaTensorInv0[1][1] + RotationMatrix[1][2] * InertiaTensorInv0[2][1];
	double tmp1_12 = RotationMatrix[1][0] * InertiaTensorInv0[0][2] + RotationMatrix[1][1] * InertiaTensorInv0[1][2] + RotationMatrix[1][2] * InertiaTensorInv0[2][2];
	double tmp1_20 = RotationMatrix[2][0] * InertiaTensorInv0[0][0] + RotationMatrix[2][1] * InertiaTensorInv0[1][0] + RotationMatrix[2][2] * InertiaTensorInv0[2][0];
	double tmp1_21 = RotationMatrix[2][0] * InertiaTensorInv0[0][1] + RotationMatrix[2][1] * InertiaTensorInv0[1][1] + RotationMatrix[2][2] * InertiaTensorInv0[2][1];
	double tmp1_22 = RotationMatrix[2][0] * InertiaTensorInv0[0][2] + RotationMatrix[2][1] * InertiaTensorInv0[1][2] + RotationMatrix[2][2] * InertiaTensorInv0[2][2];
	// RIRT=temp*R_T
	double RIRT00 = tmp1_00 * RotationMatrix[0][0] + tmp1_01 * RotationMatrix[0][1] + tmp1_02 * RotationMatrix[0][2];
	double RIRT01 = tmp1_00 * RotationMatrix[1][0] + tmp1_01 * RotationMatrix[1][1] + tmp1_02 * RotationMatrix[1][2];
	double RIRT02 = tmp1_00 * RotationMatrix[2][0] + tmp1_01 * RotationMatrix[2][1] + tmp1_02 * RotationMatrix[2][2];
	double RIRT10 = tmp1_10 * RotationMatrix[0][0] + tmp1_11 * RotationMatrix[0][1] + tmp1_12 * RotationMatrix[0][2];
	double RIRT11 = tmp1_10 * RotationMatrix[1][0] + tmp1_11 * RotationMatrix[1][1] + tmp1_12 * RotationMatrix[1][2];
	double RIRT12 = tmp1_10 * RotationMatrix[2][0] + tmp1_11 * RotationMatrix[2][1] + tmp1_12 * RotationMatrix[2][2];
	double RIRT20 = tmp1_20 * RotationMatrix[0][0] + tmp1_21 * RotationMatrix[0][1] + tmp1_22 * RotationMatrix[0][2];
	double RIRT21 = tmp1_20 * RotationMatrix[1][0] + tmp1_21 * RotationMatrix[1][1] + tmp1_22 * RotationMatrix[1][2];
	double RIRT22 = tmp1_20 * RotationMatrix[2][0] + tmp1_21 * RotationMatrix[2][1] + tmp1_22 * RotationMatrix[2][2];
	// I(t)^-1の計算ここまで

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

	double theta;
	double axis0;
	double axis1;
	double axis2;
	if (torque_length < 0.0001) {
		theta = 0.0;
		axis0 = 0.0;
		axis1 = 0.0;
		axis2 = 0.0;
	}
	else {
		theta = torque_length * DT;
		axis0 = torque0*torque_length_inv;
		axis1 = torque1*torque_length_inv;
		axis2 = torque2*torque_length_inv;
	}

	//
	qs0 = cos(theta*0.5);
	qx0 = axis0 * sin(theta*0.5);
	qy0 = axis1 * sin(theta*0.5);
	qz0 = axis2 * sin(theta*0.5);

	// クオータニオン更新
	double qqs = qs0*ps2 - qx0*px2 - qy0*py2 - qz0*pz2;
	double qqx = qs0*px2 + qx0*ps2 + qy0*pz2 - qz0*py2;
	double qqy = qs0*py2 - qx0*pz2 + qy0*ps2 + qz0*px2;
	double qqz = qs0*pz2 + qx0*py2 - qy0*px2 + qz0*ps2;

	//?微小変化になってなくね？
	double deltaRot00 = -2.0*(qy0*qy0 + qz0*qz0);
	double deltaRot01 = 2.0*(qx0*qy0 - qs0*qz0);
	double deltaRot02 = 2.0*(qx0*qz0 + qs0*qy0);
	double deltaRot10 = 2.0*(qx0*qy0 + qs0*qz0);
	double deltaRot11 = -2.0*(qx0*qx0 + qz0*qz0);
	double deltaRot12 = 2.0*(qy0*qz0 - qs0*qx0);
	double deltaRot20 = 2.0*(qx0*qz0 - qs0*qy0);
	double deltaRot21 = 2.0*(qy0*qz0 + qs0*qx0);
	double deltaRot22 = -2.0*(qx0*qx0 + qy0*qy0);

	/*比較
	RotationMatrix[0][0] = 1.0 - 2.0*(qy0*qy0 + qz0*qz0);
	RotationMatrix[0][1] = 2.0*(qx0*qy0 - qs0*qz0);
	RotationMatrix[0][2] = 2.0*(qx0*qz0 + qs0*qy0);
	RotationMatrix[1][0] = 2.0*(qx0*qy0 + qs0*qz0);
	RotationMatrix[1][1] = 1.0 - 2.0*(qx0*qx0 + qz0*qz0);
	RotationMatrix[1][2] = 2.0*(qy0*qz0 - qs0*qx0);
	RotationMatrix[2][0] = 2.0*(qx0*qz0 - qs0*qy0);
	RotationMatrix[2][1] = 2.0*(qy0*qz0 + qs0*qx0);
	RotationMatrix[2][2] = 1.0 - 2.0*(qx0*qx0 + qy0*qy0);
	*/

	qs0 = qqs; qx0 = qqx; qy0 = qqy; qz0 = qqz;

	for (int i = 0; i < nP; i++) {
		if (Typ[i] == RIGID0) {

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

/*何かの測定。こいつのせいでエラー吐くので注意*/
/*
void wdataG1(void) {
	fopen_s(&fp, "G1.txt", "a+");
	if (iLP % OPT_FQC == 0) {
		double maxd2 = 0.0;
		for (int k = 0; k < nP; k++) {
			if (Typ[k] == FLUID &&  Pos[k * 3] >= 2.5 - 5 * PARTICLE_DISTANCE  &&Pos[k * 3] <= 2.5 &&Pos[k * 3 + 1] <= 0.13 + 15 * PARTICLE_DISTANCE && Pos[k * 3 + 1] >= 0.13 - 15 * PARTICLE_DISTANCE) {
				if (Prs[k] >= 0.0 && Pos[k * 3 + 2] > maxd2)
					maxd2 = Pos[k * 3 + 2];
				else continue;
			}
			else continue;
		}
		double maxd3 = 0.0;
		for (int k = 0; k < nP; k++) {
			if (Typ[k] == FLUID &&  Pos[k * 3] >= 2.5 - 5 * PARTICLE_DISTANCE  &&Pos[k * 3] <= 2.5 &&Pos[k * 3 + 1] <= 0.13 + 15 * PARTICLE_DISTANCE && Pos[k * 3 + 1] >= 0.13 - 15 * PARTICLE_DISTANCE&&Pos[k * 3 + 2] <= maxd2 + 4 * PARTICLE_DISTANCE) {
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
	if (iLP % OPT_FQC == 0) {
		double maxd2 = 0.0;
		for (int k = 0; k < nP; k++) {
			if (Typ[k] == FLUID &&  Pos[k * 3] >= 2.95 - 5 * PARTICLE_DISTANCE  &&Pos[k * 3] <= 2.95 &&Pos[k * 3 + 1] <= 0.5 - 0.13 + 15 * PARTICLE_DISTANCE && Pos[k * 3 + 1] >= 0.5 - 0.13 - 15 * PARTICLE_DISTANCE) {
				if (Prs[k] >= 0.0 && Pos[k * 3 + 2] > maxd2)
					maxd2 = Pos[k * 3 + 2];
				else continue;
			}
			else continue;
		}
		double maxd3 = 0.0;
		for (int k = 0; k < nP; k++) {
			if (Typ[k] == FLUID &&  Pos[k * 3] >= 2.95 - 5 * PARTICLE_DISTANCE  &&Pos[k * 3] <= 2.95 &&Pos[k * 3 + 1] <= 00.5 - 0.13 + 15 * PARTICLE_DISTANCE && Pos[k * 3 + 1] >= 0.5 - 0.13 - 15 * PARTICLE_DISTANCE&&Pos[k * 3 + 2] <= maxd2 + 4 * PARTICLE_DISTANCE) {
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
*/

void ClcEMPS(void) {
	while (1) {
		if (iLP%OPT_FQC == 0) {
			int p_num = 0;
			for (int i = 0; i<nP; i++) { if (Typ[i] != GHOST)p_num++; }
			printf("%5d th TIM: %lf / p_num: %d\n", iLP, TIM, p_num);
		}
		if (iLP%OPT_FQC == 0) {
			WrtDat();
			if (TIM >= FIN_TIM) { break; }
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
		//wdataG1();

		for (int i = 0; i<nP; i++) { pav[i] += Prs[i]; }
		iLP++;
		TIM += DT;
	}
}

void read_data(int iFile) {
	sprintf_s(filename, "%s/output%05d.prof", IN_DIR_VTU, iFile);
	fopen_s(&fp, filename, "r");
	fscanf_s(fp, "%d", &NumberOfParticle);
	printf("NumberOfParticle: %d\n", NumberOfParticle);

	Position = (double*)malloc(sizeof(double)*NumberOfParticle * 3);
	Velocity = (double*)malloc(sizeof(double)*NumberOfParticle * 3);
	Pressure = (double*)malloc(sizeof(double)*NumberOfParticle);
	pressave = (double*)malloc(sizeof(double)*NumberOfParticle);
	ParticleType = (int*)malloc(sizeof(int)*NumberOfParticle);

	for (int i = 0; i<NumberOfParticle; i++) {
		int a[2];
		double b[8];
		fscanf_s(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7]);
		ParticleType[i] = a[1];
		Position[i * 3] = b[0];	Position[i * 3 + 1] = b[1];	Position[i * 3 + 2] = b[2];
		Velocity[i * 3] = b[3];	Velocity[i * 3 + 1] = b[4];	Velocity[i * 3 + 2] = b[5];
		Pressure[i] = b[6];
		pressave[i] = b[7];
	}
	fclose(fp);
}

void mk_vtu(int iFile) {
	sprintf_s(filename, "%s/particle_%05d.vtu", OUT_DIR_VTU, iFile);
	printf("Creating %s ... ", filename);

	FILE *fp;
	fopen_s(&fp, filename, "w");
	fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(fp, "<UnstructuredGrid>\n");
	fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", NumberOfParticle, NumberOfParticle);

	fprintf(fp, "<Points>\n");
	fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++)fprintf(fp, "%lf %lf %lf\n", Position[i * 3], Position[i * 3 + 1], Position[i * 3 + 2]);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</Points>\n");

	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++) { fprintf(fp, "%d\n", ParticleType[i]); }
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++) {
		double val = sqrt(Velocity[i * 3] * Velocity[i * 3] + Velocity[i * 3 + 1] * Velocity[i * 3 + 1] + Velocity[i * 3 + 2] * Velocity[i * 3 + 2]);
		fprintf(fp, "%f\n", (float)val);
	}
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='pressave' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++) { fprintf(fp, "%f\n", (float)pressave[i]); }
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</PointData>\n");

	fprintf(fp, "<Cells>\n");
	fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++)fprintf(fp, "%d\n", i);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++)fprintf(fp, "%d\n", i + 1);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++)fprintf(fp, "1\n");
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</Cells>\n");

	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</UnstructuredGrid>\n");
	fprintf(fp, "</VTKFile>\n");

	fclose(fp);
	printf("done.\n");
}

void setProfFileNum() {

	ProfFileNumber = 0;

	struct _finddata_t c_file;
	long   hFile;

	// 最初のファイルを探索
	sprintf_s(filename, "%s/output*.prof", IN_DIR_VTU);
	if ((hFile = _findfirst(filename, &c_file)) == -1L)
		printf("ファイルは存在しません。\n");
	else {
		ProfFileNumber++;	//printf("%s\n", c_file.name);

							// 残りのファイルを探索
		while (_findnext(hFile, &c_file) == 0) {
			ProfFileNumber++;	//printf("%s\n", c_file.name);
		}
		_findclose(hFile);
	}
	printf("%i", ProfFileNumber);
}

int main(int argc, char** argv) {
	printf("start emps.\n");
	RdDat();
	AlcBkt();
	SetPara();
	init_rigid0();

	ClcEMPS();

	free(Acc);	free(Pos);	free(Vel);
	free(Prs);	free(pav);	free(Typ);
	free(bfst);	free(blst);	free(nxt);
	printf("end emps.\n");

	// vtu
	setProfFileNum();
	for (int iFile = 0; iFile<ProfFileNumber; iFile++) {

		read_data(iFile);
		mk_vtu(iFile);

		free(Position);
		free(Velocity);
		free(Pressure);
		free(pressave);
		free(ParticleType);
	}

	printf(" END. \n");
	//getchar();
	return 0;
}
