#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	int i;

	char *fileNumber;
	double WAVE_HEIGHT;
	double CENTER_CUBE_X;
	double DNS_RIGID0;	//���̖��x

	printf("�����̑��� = %d\n", argc);
	// �����Ŏ󂯎��i�o�b�`�����j
	printf("%d\n", argc);
	if (argc==5) {
		//printf("%s, %s, %s, %s, %s", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
		fileNumber = argv[1];
		WAVE_HEIGHT = atof(argv[2]);
		CENTER_CUBE_X = atof(argv[3]);
		DNS_RIGID0 = atof(argv[4]);	//���̖��x
		printf("fileNumber:%s, WAVE_HEIGHT:%f, CENTER_CUBE_X:%f, DNS_RIGID0:%f", fileNumber, WAVE_HEIGHT, CENTER_CUBE_X, DNS_RIGID0);
	}
	else {
		printf("%s", "error in args");
	}
	strcat("../data/init_position", fileNumber, ".prof");
	getchar();

	return 0;
}