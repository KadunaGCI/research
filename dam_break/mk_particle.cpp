#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../params.h"

#define OUTPUT_FILE "../data/init_position.prof"

// —Ìˆæ
int nx = (int)((mk_MAX_X - mk_MIN_X) / PARTICLE_DISTANCE) + RE * 2;
int ny = (int)((mk_MAX_Y - mk_MIN_Y) / PARTICLE_DISTANCE) + RE * 2;
int nz = (int)((mk_MAX_Z - mk_MIN_Z) / PARTICLE_DISTANCE) + RE * 2;

int tx = (int)((0.055) / PARTICLE_DISTANCE);
int ty = (int)((0.055) / PARTICLE_DISTANCE);
int tz = (int)((0.055) / PARTICLE_DISTANCE);

int nxy = nx*ny;
int nxyz = nx*ny*nz;
int txyz = tx*ty*tz;
int nxyz1 = nxyz + txyz;

int NumberOfParticle;
int *ParticleType;
float *Position;

int main(int argc, char** argv){

	printf("start mk_particle\n");
	printf("nx:%d ny:%d nz:%d nxyz:%d\n", nx, ny, nz,nxyz1);

	ParticleType = (int*)malloc(sizeof(int)*nxyz1);
	Position = (float*)malloc(sizeof(float)*nxyz1*3);
	int NumberOfParticle = 0;

	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		ParticleType[ip] = GHOST;
		Position[ip*3  ] = mk_MIN_X + PARTICLE_DISTANCE*(ix-RE+0.5);
		Position[ip * 3 + 1] = mk_MIN_Y + PARTICLE_DISTANCE*(iy-RE+0.5);
		Position[ip * 3 + 2] = mk_MIN_Z + PARTICLE_DISTANCE*(iz-RE+0.5);
	}}}

	int nr0 = 0;
	int nr1 = 0;
	int nr2 = 0;

	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		double x = Position[ip * 3];
		double y = Position[ip * 3 + 1];
		double z = Position[ip * 3 + 2];

		if (ix < RE || ix >= nx - RE || iy < RE || iy >= ny - RE || iz < RE){
			ParticleType[ip] = GHOST;
		}
		else if (x > 2.5+5*PARTICLE_DISTANCE && x <= 2.5+0.15-5*PARTICLE_DISTANCE && y >= 0.06+5*PARTICLE_DISTANCE && y <= 0.06+0.15-5*PARTICLE_DISTANCE && z <= 0.75-5*PARTICLE_DISTANCE){
			ParticleType[ip] = GHOST;
		}
		else if (x > 2.95 + 5 * PARTICLE_DISTANCE && x <= 2.95+0.15 - 5 * PARTICLE_DISTANCE && y >= 0.5-0.06-0.15 + 5 * PARTICLE_DISTANCE && y <= 0.5-0.06 - 5 * PARTICLE_DISTANCE && z <= 0.75-5* PARTICLE_DISTANCE){
			ParticleType[ip] = GHOST;
		}
		else if (x > 2.5 && x <= 2.5+0.15 && y >= 0.06 && y <= 0.06+0.15 && z <= 0.75){
			ParticleType[ip] = WALL;
			NumberOfParticle++;
		}
		else if (x > 2.95 && x <= 2.95+0.15 && y >= 0.5-0.06-0.15 && y <= 0.5-0.06 && z <= 0.75){
			ParticleType[ip] = WALL;
			NumberOfParticle++;
		}
		
		else if (z <= WAVE_HEIGHT&&x <= WAVE_WIDTH){
			ParticleType[ip] = FLUID;
			NumberOfParticle++;
		}
	}
	}
	}

	for (int iz = 0; iz < tz; iz++){
		for (int iy = 0; iy < ty; iy++){
			for (int ix = 0; ix < tx; ix++){
				int ip = nxyz + iz*tx*ty + iy*tx + ix;
				Position[ip * 3] = 2.532 - 0.055*0.5 + PARTICLE_DISTANCE*(ix + 0.5);
				Position[ip * 3 + 1] = 0.313 - 0.055*0.5 + PARTICLE_DISTANCE*(iy + 0.5);
				Position[ip * 3 + 2] = PARTICLE_DISTANCE*(iz + 0.5);
				ParticleType[ip] = RIGID0;
				nr0++;
				NumberOfParticle++;
			}
		}
	}

	printf("NumberOfParticle:     %d\n",NumberOfParticle);

	FILE* fp;
	fopen_s(&fp, OUTPUT_FILE, "w");
	fprintf(fp,"%d %d %d %d\n",NumberOfParticle,nr0,nr1,nr2);
	int k=0;
	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		if(ParticleType[ip]==GHOST)continue;
		fprintf(fp," %d %d %lf %lf %lf 0.0 0.0 0.0 0.0 0.0\n",k,ParticleType[ip],Position[ip*3],Position[ip*3+1],Position[ip*3+2]);
		k++;
	}}}
	printf("%d\n", k);
	fprintf(fp, "%s\n","‚Í‚°");
	for (int iz = 0; iz < tz; iz++){
		for (int iy = 0; iy < ty; iy++){
			for (int ix = 0; ix < tx; ix++){
				int ip = nxyz + iz*tx*ty + iy*tx + ix;
				if (ParticleType[ip] == GHOST)continue;
				fprintf(fp, " %d %d %lf %lf %lf 0.0 0.0 0.0 0.0 0.0\n", k, ParticleType[ip], Position[ip * 3], Position[ip * 3 + 1], Position[ip * 3 + 2]);
				k++;
			}
		}
	}
	fclose(fp);

	free(ParticleType);	free(Position);
	printf("end mk_particle\n");
	getchar();
	return 0;
}
