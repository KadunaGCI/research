#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	int i;

	char *fileNumber;
	double WAVE_HEIGHT;
	double CENTER_CUBE_X;
	double DNS_RIGID0;	//剛体密度

	printf("引数の総個数 = %d\n", argc);
	// 引数で受け取る（バッチ処理）
	printf("%d\n", argc);
	if (argc==5) {
		//printf("%s, %s, %s, %s, %s", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
		fileNumber = argv[1];
		WAVE_HEIGHT = atof(argv[2]);
		CENTER_CUBE_X = atof(argv[3]);
		DNS_RIGID0 = atof(argv[4]);	//剛体密度
		printf("fileNumber:%s, WAVE_HEIGHT:%f, CENTER_CUBE_X:%f, DNS_RIGID0:%f", fileNumber, WAVE_HEIGHT, CENTER_CUBE_X, DNS_RIGID0);
	}
	else {
		printf("%s", "error in args");
	}
	strcat("../data/init_position", fileNumber, ".prof");
	getchar();

	return 0;
}