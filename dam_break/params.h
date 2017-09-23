#pragma once
/*****************************************************************

	mk_particle.cpp�ŗ��p

*****************************************************************/

//���́i�ƍ��́j�͈̔�
#define mk_MIN_X  0.0
#define mk_MIN_Y  0.0
#define mk_MIN_Z  0.0
#define mk_MAX_X  3.5
#define mk_MAX_Y  0.5
#define mk_MAX_Z  0.75

#define GHOST -1
#define WALL  0
#define FLUID 1
#define RIGID0 2
#define NUM_TYP  4		//���q�̎�ސ�

#define RE  3  //3 4 5		PARTICLE_DISTANCE*(RE+0.1)���e�����a�ł��邱�Ƃɒ���
#define WAVE_HEIGHT 0.35
#define WAVE_WIDTH 0.5
#define PARTICLE_DISTANCE 0.005 
#define CUBE_LENGTH 0.03


/*****************************************************************
	
	emps.cpp�ŗ��p

*****************************************************************/

//#define PARTICLE_DISTANCE 0.002					//���ϗ��q�ԋ���

#define MIN_X  (mk_MIN_X - PARTICLE_DISTANCE*3)	//��͗̈��x�����̍ŏ��l
#define MIN_Y  (mk_MIN_Y - PARTICLE_DISTANCE*3)	//��͗̈��y�����̍ŏ��l
#define MIN_Z  (mk_MIN_Z - PARTICLE_DISTANCE*3)	//��͗̈��z�����̍ŏ��l
#define MAX_X  (mk_MAX_X + PARTICLE_DISTANCE*3)	//��͗̈��x�����̍ő�l
#define MAX_Y  (mk_MAX_Y + PARTICLE_DISTANCE*3)	//��͗̈��y�����̍ő�l
#define MAX_Z  (mk_MAX_Z + PARTICLE_DISTANCE*30)	//��͗̈��z�����̍ő�l

/*
#define GHOST -1
#define FLUID 0
#define WALL  1
#define RIGID0 2
*/

//#define NUM_TYP  4		//���q�̎�ސ�

#define DNS_FLUID 996.51		//���̗��q�̖��x
#define DNS_WALL 996.51		//�Ǘ��q�̖��x
#define DNS_RIGID0 2120.0	//���̖��x	

#define DT 0.0001			//���ԍ��ݕ�
#define dt_inv   double (1/DT)	
#define FIN_TIM 0.5		//���Ԃ̏��
#define SND 18.0			//����
#define OPT_FQC 200		//�o�͊Ԋu�����߂锽����
#define KNM_VSC_FRUID 0.000001	//���S���W��
#define KNM_VSC_WALL 0.000001	//���S���W��
#define KNM_VSC_RIGID0 0.000001	//���S���W��
#define DIM 3				//������
#define CRT_NUM 0.1		//�N�[����������
#define COL_RAT 0.2		//�ڋ߂������q�̔�����
#define DST_LMT_RAT 0.9	//����ȏ�̗��q�Ԃ̐ڋ߂������Ȃ������̌W��
#define G_X 0.0			//�d�͉����x��x����
#define G_Y 0.0			//�d�͉����x��y����
#define G_Z -9.8			//�d�͉����x��z����
#define WEI(dist, re) ((re/dist) - 1.0)	//�d�݊֐�
