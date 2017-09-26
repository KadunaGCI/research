#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../params.h"

//#define OUTPUT_FILE "../data/init_position.prof"
//#define IN_DATA OUTPUT_FILE
//#define OUT_DATA_VTU "../data/chk_init_position.vtu"
char OUTPUT_FILE[100];
char IN_DATA[100];
char OUT_DATA_VTU[100];

// 領域
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

// vtu
FILE* fp;
char filename[256];
//int NumberOfParticle, nr0;
int nr0;
double *Position_vtu;
double *Velocity_vtu;
double *Pressure_vtu;
double *pressave_vtu;
//int *ParticleType;

/*
初期条件を浮遊物に→密度を調整して浮いているところから
流体外に設置するのであれば、別の計算定を作ればよい
*/

void read_data() {

	fopen_s(&fp, IN_DATA, "r");
	fscanf_s(fp, "%d %d", &NumberOfParticle, &nr0);
	printf("NumberOfParticle: %d\n", NumberOfParticle);

	Position_vtu = (double*)malloc(sizeof(double)*NumberOfParticle * 3);
	Velocity_vtu = (double*)malloc(sizeof(double)*NumberOfParticle * 3);
	Pressure_vtu = (double*)malloc(sizeof(double)*NumberOfParticle);
	pressave_vtu = (double*)malloc(sizeof(double)*NumberOfParticle);
	ParticleType = (int*)malloc(sizeof(int)*NumberOfParticle);

	for (int i = 0; i<NumberOfParticle; i++) {
		int a[2];
		double b[8];
		fscanf_s(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7]);
		ParticleType[i] = a[1];
		Position_vtu[i * 3] = b[0];	Position_vtu[i * 3 + 1] = b[1];	Position_vtu[i * 3 + 2] = b[2];
		Velocity_vtu[i * 3] = b[3];	Velocity_vtu[i * 3 + 1] = b[4];	Velocity_vtu[i * 3 + 2] = b[5];
		Pressure_vtu[i] = b[6];
		pressave_vtu[i] = b[7];
	}
	fclose(fp);
}

void mk_vtu() {

	printf("%s", OUT_DATA_VTU);

	FILE *fp;
	fopen_s(&fp, OUT_DATA_VTU, "w");
	fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(fp, "<UnstructuredGrid>\n");
	fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", NumberOfParticle, NumberOfParticle);

	fprintf(fp, "<Points>\n");
	fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++)fprintf(fp, "%lf %lf %lf\n", Position_vtu[i * 3], Position_vtu[i * 3 + 1], Position_vtu[i * 3 + 2]);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</Points>\n");

	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++) { fprintf(fp, "%d\n", ParticleType[i]); }
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++) {
		double val = sqrt(Velocity_vtu[i * 3] * Velocity_vtu[i * 3] + Velocity_vtu[i * 3 + 1] * Velocity_vtu[i * 3 + 1] + Velocity_vtu[i * 3 + 2] * Velocity_vtu[i * 3 + 2]);
		fprintf(fp, "%f\n", (float)val);
	}
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='pressave' format='ascii'>\n");
	for (int i = 0; i<NumberOfParticle; i++) { fprintf(fp, "%f\n", (float)pressave_vtu[i]); }
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

// argcは引数の総個数(プログラム名も含む)(渡す必要なし)
int main(int argc, char** argv){

	// 引数で受け取る（バッチから）
	char *fileNumber;
	double WAVE_HEIGHT;
	double CENTER_CUBE_X;
	double DNS_RIGID0;	//剛体密度
	
	if (argc == 2) {
		
		char sp[] = ",";
		char *tok;

		tok = strtok(argv[1], sp);
		fileNumber = tok;
		tok = strtok(NULL, sp);
		WAVE_HEIGHT = atof(tok);
		tok = strtok(NULL, sp);
		CENTER_CUBE_X = atof(tok);
		tok = strtok(NULL, sp);
		DNS_RIGID0 = atof(tok);
		sprintf(OUTPUT_FILE, "../data/init_position_%s.prof", fileNumber);
		sprintf(IN_DATA, "../data/init_position_%s.prof", fileNumber);
		sprintf(OUT_DATA_VTU, "../data/chk_init_position_%s.vtu", fileNumber);
	}
	else {
		printf("%s", "error in args");
		getchar();
		return 1;
	}

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

	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		double x = Position[ip * 3];
		double y = Position[ip * 3 + 1];
		double z = Position[ip * 3 + 2];

		if (ix < RE || ix >= nx - RE || iy < RE || iy >= ny - RE || iz < RE){
			ParticleType[ip] = WALL;
			NumberOfParticle++;
		}
		else if (x >= CENTER_CUBE_X - CUBE_LENGTH / 2 && x < CENTER_CUBE_X + CUBE_LENGTH / 2 && y >= CENTER_CUBE_Y - CUBE_LENGTH / 2 && y < CENTER_CUBE_Y + CUBE_LENGTH / 2 && z >= CENTER_CUBE_Z - CUBE_LENGTH / 2 && z < CENTER_CUBE_Z + CUBE_LENGTH / 2) {
			ParticleType[ip] = RIGID0;
			NumberOfParticle++;
			nr0++;
		}
		else if (z <= WAVE_HEIGHT&&x <= WAVE_WIDTH){
			ParticleType[ip] = FLUID;
			NumberOfParticle++;
		}
	}
	}
	}

	printf("NumberOfParticle:     %d\n",NumberOfParticle);
	fopen_s(&fp, OUTPUT_FILE, "w");
	fprintf(fp,"%d %d\n",NumberOfParticle,nr0);
	int k=0;
	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		if(ParticleType[ip]==GHOST)continue;
		fprintf(fp," %d %d %lf %lf %lf 0.0 0.0 0.0 0.0 0.0\n",k,ParticleType[ip],Position[ip*3],Position[ip*3+1],Position[ip*3+2]);
		k++;
	}}}

	fclose(fp);

	free(ParticleType);	free(Position);
	printf("end mk_particle\n");
	
	// vtu
	read_data();
	mk_vtu();

	free(Position_vtu);
	free(Velocity_vtu);
	free(Pressure_vtu);
	free(pressave_vtu);
	free(ParticleType);

	return 0;
}
