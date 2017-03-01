#ifndef DICTIONARY_GENERATION_H
#define DICTIONARY_GENERATION_H

namespace dpl{

double **InitializeDictionary( int featureDim, int sampleDim ){

	double **Wd = (double**)malloc(featureDim*sizeof(double*));
	for( unsigned int j=0; j<featureDim; j++ )
		Wd[j] = (double*)malloc(sampleDim*sizeof(double));

	return Wd;
}

double **GenerateRandomDictionary( int featureDim, int sampleDim ){

	double **Wd = InitializeDictionary( featureDim, sampleDim );
	srand((unsigned)time(0));
	unsigned int myseed = (unsigned int) RAND_MAX * rand();

    	for( unsigned int j=0; j<featureDim; j++ ){
        	for( unsigned int i=0; i<sampleDim; i++ )
			Wd[j][i] = 2*(double) (rand_r(&myseed) / (RAND_MAX + 1.0))-1;
	}
	return Wd;
}

double **GenerateRandomPatchDictionary( int featureDim, int sampleDim, int sampleNumber, double **sample ){

	double **Wd = InitializeDictionary( featureDim, sampleDim );
	srand((unsigned)time(0));
	unsigned int myseed = (unsigned int) RAND_MAX * rand();

	for( unsigned int i=0; i<featureDim; i++ ){
        	unsigned int index = rand_r(&myseed)%sampleNumber;
        	for( unsigned int j=0; j<sampleDim; j++ )
            		Wd[i][j]=sample[index][j];
	}
	return Wd;
}

void DictionaryNormalization( int featureDim, int sampleDim, double **Wd ){
	for( unsigned int i=0; i<featureDim; i++ ){
		double sum = 0;
		for( unsigned int j=0; j<sampleDim; j++ )
			sum += Wd[i][j]*Wd[i][j];
		sum = sqrt(sum);

		if( sum!=0 ){
			for( unsigned int j=0; j<sampleDim; j++ )
				Wd[i][j] = Wd[i][j]/sum;
		}
	}
}

double **readDictionary( char *FileName, int featureDim, int sampleDim ) {

	double **Wd = InitializeDictionary( featureDim, sampleDim );
	FILE *fp;
    	fp = fopen( FileName, "rw");
    	if( fp == NULL ){
		printf("could not find dictionary file %s\n", FileName);
        	exit(0);
	}
    	for( unsigned int i=0; i<sampleDim; i++ ){
        	for( unsigned int j=0; j<featureDim; j++){
	        	if( fscanf(fp, "%lf", &Wd[j][i])==0 )
				exit(0);
		}
        }
	fclose(fp);
	return Wd;
}

void saveDictionary( int featureDim, int sampleDim, double **Wd, char *dictionaryName ){

	FILE *fp;
        fp = fopen( dictionaryName, "w");
        if( fp == NULL ){
		printf("could not find dictionary file: %s\n", dictionaryName);
            	exit(0);
	}

	for( unsigned int i=0; i<sampleDim; i++ ){
		for( unsigned int j=0; j<featureDim; j++)
	        	fprintf(fp, "%.50lf ", Wd[j][i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("Save dictionary file in %s\n", dictionaryName);
}


void clearDictionary( int featureDim, double **Wd ){
	for( unsigned int i=0; i<featureDim; i++ )
		free(Wd[i]);
	free(Wd);
}

}

#endif /* Dictionary Generation*/
