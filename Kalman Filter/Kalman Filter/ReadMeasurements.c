#include "Kalman.h"

Matrix ReadMeasurements(char* filePath)
{
	char data;
	int dataNumber = 0;
	int i;
	FILE* fp = fopen(filePath, "r");

	if (fp == NULL)
	{
		printf("Could not open file!\n");
		return;
	}

	while ((data = fgetc(fp)) != EOF)
	{
		if (data == '\n')
			dataNumber++;
	}

	dataNumber++; // ++1 since no '\n'(return line character) for the last data

	fclose(fp);
	FILE* fp2 = fopen(filePath, "r");
	Matrix measurement = CreateEmptyMatrix(dataNumber, 1);

	char* test = (char *)malloc(100 * sizeof(char)); // allocate memory for a string of 100 characters

	for (i = 0; i < dataNumber; i++)
	{
		fscanf(fp2, "%s", test);
		measurement.content[i][0] = atof(test);
	}

	fclose(fp2);
	return measurement;
}