#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../params.h"

#define OUTPUT_FILE "../data/init_position.prof"


// ���b�V�����쐬
// ����������wall�̂��߁~2
int nx = (int)((mk_MAX_X - mk_MIN_X) / PARTICLE_DISTANCE) + RE * 2;	// 0 <= x < nx�Ƃ��ė��p
int ny = (int)((mk_MAX_Y - mk_MIN_Y) / PARTICLE_DISTANCE) + RE * 2;	// 0 <= y < ny�Ƃ��ė��p
int nz = (int)((mk_MAX_Z - mk_MIN_Z) / PARTICLE_DISTANCE) + RE * 2;	// 0 <= z < nz�Ƃ��ė��p

int nxy = nx*ny;
int nxyz = nx*ny*nz;	//�����q��

int NumberOfParticle;	//�L�����q��

int *ParticleType;
float *Position;

int main(int argc, char** argv) {

	printf("start mk_particle\n");
	printf("nx:%d ny:%d nz:%d nxyz:%d\n", nx, ny, nz, nxyz);

	ParticleType = (int*)malloc(sizeof(int)*nxyz);
	Position = (float*)malloc(sizeof(float)*nxyz * 3);

	NumberOfParticle = 0;

	// position��ݒ�
	for (int iz = 0; iz<nz; iz++) {
		for (int iy = 0; iy<ny; iy++) {
			for (int ix = 0; ix<nx; ix++) {
				int ip = iz*nxy + iy*nx + ix;
				ParticleType[ip] = GHOST;		//��U���ׂ�GHOST��
				Position[ip * 3] = mk_MIN_X + PARTICLE_DISTANCE*(ix - RE + 0.5);		//���b�V���̒��Ԃɗ��q�𐶐�
				Position[ip * 3 + 1] = mk_MIN_Y + PARTICLE_DISTANCE*(iy - RE + 0.5);
				Position[ip * 3 + 2] = mk_MIN_Z + PARTICLE_DISTANCE*(iz - RE + 0.5);
			}
		}
	}

	int nr0 = 0;		//����0�̐�

	// particleType��ݒ�
	for (int iz = 0; iz<nz; iz++) {
		for (int iy = 0; iy<ny; iy++) {
			for (int ix = 0; ix<nx; ix++) {
				int ip = iz*nxy + iy*nx + ix;
				double x = Position[ip * 3];
				double y = Position[ip * 3 + 1];
				double z = Position[ip * 3 + 2];

				// z������������goast�̂܂�
				if (ix<RE || ix >= nx - RE || iy<RE || iy >= ny - RE || iz<RE) {
					ParticleType[ip] = WALL;
					NumberOfParticle++;
				}
				else if (x >= mk_MAX_X/2 - CUBE_LENGTH/2 && x < mk_MAX_X / 2 + CUBE_LENGTH/2 && y >= mk_MAX_Y / 2 - CUBE_LENGTH/2 && y < mk_MAX_Y / 2 + CUBE_LENGTH/2 && z >= WAVE_HEIGHT - CUBE_LENGTH/2 && z < WAVE_HEIGHT + CUBE_LENGTH/2) {
					ParticleType[ip] = RIGID0;
					NumberOfParticle++;
					nr0++;
				}
				else if (z <= WAVE_HEIGHT) {
					ParticleType[ip] = FLUID;
					NumberOfParticle++;
				}
			}
		}
	}

	printf("NumberOfParticle:     %d\n", NumberOfParticle);

	// goast�ȊO���o��
	FILE* fp;
	fopen_s(&fp, OUTPUT_FILE, "w");
	fprintf(fp, "%d %d\n", NumberOfParticle, nr0);
	int k = 0;
	for (int iz = 0; iz<nz; iz++) {
		for (int iy = 0; iy<ny; iy++) {
			for (int ix = 0; ix<nx; ix++) {
				int ip = iz*nxy + iy*nx + ix;
				if (ParticleType[ip] == GHOST)continue;
				fprintf(fp, " %d %d %lf %lf %lf 0.0 0.0 0.0 0.0 0.0\n", k, ParticleType[ip], Position[ip * 3], Position[ip * 3 + 1], Position[ip * 3 + 2]);
				k++;
			}
		}
	}
	fclose(fp);

	free(ParticleType);	free(Position);
	printf("end mk_particle\n");

	//getchar();

	return 0;
}

