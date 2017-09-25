#pragma once
/*****************************************************************

	mk_particle.cppで利用

*****************************************************************/

//流体（と剛体）の範囲
#define mk_MIN_X  0.0
#define mk_MIN_Y  0.0
#define mk_MIN_Z  0.0
#define mk_MAX_X  1.0
#define mk_MAX_Y  0.2
#define mk_MAX_Z  0.6

#define GHOST -1
#define WALL  0
#define FLUID 1
#define RIGID0 2
#define NUM_TYP  4		//粒子の種類数

#define RE  3  //3 4 5		PARTICLE_DISTANCE*(RE+0.1)が影響半径であることに注意

#define WAVE_HEIGHT 0.5
#define WAVE_WIDTH 0.25
// 奥行きはmk_MAX_Y

#define PARTICLE_DISTANCE 0.02 

#define CUBE_LENGTH 0.1
// 各面と平行におく
#define CENTER_CUBE_X 0.25
#define CENTER_CUBE_Y mk_MAX_Y/2
#define CENTER_CUBE_Z CUBE_LENGTH/2


/*****************************************************************
	
	emps.cppで利用

*****************************************************************/

//#define PARTICLE_DISTANCE 0.002					//平均粒子間距離

#define MIN_X  (mk_MIN_X - PARTICLE_DISTANCE*3)	//解析領域のx方向の最小値
#define MIN_Y  (mk_MIN_Y - PARTICLE_DISTANCE*3)	//解析領域のy方向の最小値
#define MIN_Z  (mk_MIN_Z - PARTICLE_DISTANCE*3)	//解析領域のz方向の最小値
#define MAX_X  (mk_MAX_X + PARTICLE_DISTANCE*3)	//解析領域のx方向の最大値
#define MAX_Y  (mk_MAX_Y + PARTICLE_DISTANCE*3)	//解析領域のy方向の最大値
#define MAX_Z  (mk_MAX_Z + PARTICLE_DISTANCE*10)	//解析領域のz方向の最大値

/*
#define GHOST -1
#define FLUID 0
#define WALL  1
#define RIGID0 2
*/

//#define NUM_TYP  4		//粒子の種類数

#define DNS_FLUID 1000		//流体粒子の密度
#define DNS_WALL 1000		//壁粒子の密度
#define DNS_RIGID0 500	//剛体密度	

#define DT 0.0005			//時間刻み幅
#define dt_inv   double (1/DT)	
#define FIN_TIM 3.0		//時間の上限
#define SND 22			//音速
#define OPT_FQC 100		//出力間隔を決める反復数
#define KNM_VSC_FRUID 0.000001	//動粘性係数
//#define KNM_VSC_WALL 0.000001	//動粘性係数
//#define KNM_VSC_RIGID0 0.000001	//動粘性係数
#define DIM 3				//次元数
#define CRT_NUM 0.1		//クーラン条件数
#define COL_RAT 0.2		//接近した粒子の反発率
#define DST_LMT_RAT 0.9	//これ以上の粒子間の接近を許さない距離の係数
#define G_X 0.0			//重力加速度のx成分
#define G_Y 0.0			//重力加速度のy成分
#define G_Z -9.8			//重力加速度のz成分
#define WEI(dist, re) ((re/dist) - 1.0)	//重み関数

