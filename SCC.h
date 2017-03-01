#ifndef SPARSE_COORDINATE_CODING_H
#define SPARSE_COORDINATE_CODING_H

namespace dpl{

static unsigned int myseed;

double getAbs( double value ){
	if( value < 0 )
		return -1*value;
	else
		return value;
}

double **FeatureInitialization( int featureDim, int sampleNumber ){

	double **feature = (double**)malloc(sampleNumber*sizeof(double*));
        for( unsigned int i=0; i<sampleNumber; i++ ){
		feature[i] = (double*)malloc(featureDim*sizeof(double));
		for( unsigned int j=0; j<featureDim; j++ )
			feature[i][j] = 0;
	}
	return feature;
}

std::vector<int> *NonZeroIndexInitialization( int sampleNumber ){
	std::vector<int> *nonZeroIndex = new std::vector<int> [sampleNumber];
	return nonZeroIndex;
}

double ShrinkageFunction( double value, double theta ){

	if( value < -theta )
		return value+theta;
	else if( value > theta )
		return value-theta;
	else
		return 0;
}

double *Initialize_A_Copy( int featureDim ){

	double *A_Copy = (double*)malloc(featureDim*sizeof(double));
	for( unsigned int i=0; i<featureDim; i++ )
		A_Copy[i]=0;
	return A_Copy;
}

double *Initialize_A( int featureDim ){

	double *A = (double*)malloc(featureDim*sizeof(double));
	for( unsigned int i=0; i<featureDim; i++ )
		A[i]=0;
	return A;
}

void Initialize_A( double *A, double *A_Copy, int featureDim ){
	for( unsigned int i=0; i<featureDim; i++ ){
		A[i]=A_Copy[i];
		A_Copy[i]=0;
    	}
}

void Update_A( double *A, double *A_Copy, double *feature, std::vector<int> &nonZeroIndex ){
	for( unsigned int i=0; i<nonZeroIndex.size(); i++ ){
		A[nonZeroIndex[i]] += feature[nonZeroIndex[i]]*feature[nonZeroIndex[i]];
		A_Copy[nonZeroIndex[i]] += feature[nonZeroIndex[i]]*feature[nonZeroIndex[i]];
	}
}

double getNonNegativeFeature( double featureElement, double optimalT ){
	if( featureElement+optimalT>=0 )
		return optimalT;
	else
		return -1*featureElement;
}

int *getRandomIndex( int size ){

	std::vector<int> index (size);
	int *data=(int*)malloc(size*sizeof(int));
	for( unsigned int i=0; i<size; i++ )
	        index[i] = i;

	for( unsigned int i=0; i<size; i++ ){
    		int randomIndex = rand_r(&myseed)%index.size();
        	data[i] = index[randomIndex];
        	index.erase(index.begin()+randomIndex);
    	}
    	return data;
}

void UpdateFeature( double **Wd, double *sample, double *residuals, double *feature, std::vector<int> &nonZeroIndex, double lambda, int layers, int featureDim, int sampleDim, bool NonNegative ){

    for( unsigned int i = 0; i<sampleDim; i++ ){
  		residuals[i] = -sample[i];
        for( unsigned int j = 0; j<nonZeroIndex.size(); j++ )
            residuals[i] += Wd[nonZeroIndex[j]][i]*feature[nonZeroIndex[j]];
    }

	nonZeroIndex.resize(0);
	int *randomIndex = getRandomIndex(featureDim );

	for ( unsigned int i = 0; i < featureDim; i++ ){

        	double optimalT;
        	double derivative = 0;

        	for (unsigned int j = 0;j < sampleDim; j++)
                	derivative += (residuals[j]*Wd[randomIndex[i]][j]);

		optimalT = ShrinkageFunction( feature[randomIndex[i]]-derivative, lambda )-feature[randomIndex[i]];

		if( NonNegative )
			optimalT = getNonNegativeFeature( feature[randomIndex[i]], optimalT );

		feature[randomIndex[i]] += optimalT;

        	if ( optimalT!=0 ){
            		for (unsigned int j = 0;j < sampleDim; j++)
                		residuals[j] += optimalT*Wd[randomIndex[i]][j];
        	}

		if( feature[randomIndex[i]]!=0 )
			nonZeroIndex.push_back(randomIndex[i]);

	}

	for ( unsigned int k = 1; k < layers; k++ ){
		for ( unsigned int i = 0; i < nonZeroIndex.size(); i++ ){
        		double optimalT;
        		double derivative = 0;
        		for (unsigned int j = 0;j < sampleDim; j++)
                		derivative += (residuals[j]*Wd[nonZeroIndex[i]][j]);

			optimalT = ShrinkageFunction( feature[nonZeroIndex[i]]-derivative, lambda )-feature[nonZeroIndex[i]];

			if( NonNegative )
				optimalT = getNonNegativeFeature( feature[nonZeroIndex[i]], optimalT );

			feature[nonZeroIndex[i]] += optimalT;

        		if ( optimalT!=0 ){
            			for (unsigned int j = 0;j < sampleDim; j++)
                			residuals[j] += optimalT*Wd[nonZeroIndex[i]][j];
        		}
		}
	}

	nonZeroIndex.resize(0);
	for ( unsigned int i = 0; i < featureDim; i++ ){
		if( feature[i]!=0 )
			nonZeroIndex.push_back(i);
	}
	free(randomIndex);
}


void UpdateWd( double **Wd, double *residuals, double *feature, double *A, std::vector<int> &nonZeroIndex, int sampleDim, bool NonNegative ){

    	for ( unsigned int i = 0; i < sampleDim; i++ ){
        	for ( unsigned int j = 0; j < nonZeroIndex.size(); j++ ){
                if( NonNegative && Wd[nonZeroIndex[j]][i] - feature[nonZeroIndex[j]]*residuals[i]*dpl::learningRate(A,nonZeroIndex[j])<0 )
                    Wd[nonZeroIndex[j]][i] = 0;
                else
                    Wd[nonZeroIndex[j]][i] = Wd[nonZeroIndex[j]][i] - feature[nonZeroIndex[j]]*residuals[i]*dpl::learningRate(A,nonZeroIndex[j]);
            }
    	}
}

void NormalizeWd( double **Wd, std::vector<int> &nonZeroIndex, int sampleDim ){
	for( unsigned int i=0; i<nonZeroIndex.size(); i++ ){
		double sum = 0;
		for( unsigned int j=0; j<sampleDim; j++ )
			sum += Wd[nonZeroIndex[i]][j]*Wd[nonZeroIndex[i]][j];
		sum = sqrt(sum);

		if( sum!=0 ){
			for( unsigned int j=0; j<sampleDim; j++ )
				Wd[nonZeroIndex[i]][j] = Wd[nonZeroIndex[i]][j]/sum;
		}
	}
}

void saveFeature( double **feature, char *FeatureFileName, int featureDim, int sampleNumber ){

	printf("Save Features in %s\n", FeatureFileName);

	FILE *fp;
    fp = fopen( FeatureFileName, "w");
    if( fp == NULL ){
		printf("could not find feature file %s\n", FeatureFileName);
        exit(0);
	}

	for( unsigned int i=0; i<featureDim; i++ ){
		for( unsigned int j=0; j<sampleNumber; j++)
	        	fprintf(fp, "%.15lf ", feature[j][i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void saveNonZeroIndex( std::vector<int> *nonZeroIndex, char *IndexFileName, int featureDim, int sampleNumber ){

	printf("Save nonZero index in %s\n", IndexFileName);

	FILE *fp;
    fp = fopen( IndexFileName, "w");
    if( fp == NULL ){
		printf("could not find index file %s\n", IndexFileName);
        exit(0);
	}

	for( unsigned int i=0; i<sampleNumber; i++ ){
		for( unsigned int j=0; j<nonZeroIndex[i].size(); j++)
	        	fprintf(fp, "%d ", nonZeroIndex[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void clearFeature( int sampleNumber, double **feature ){

	for( unsigned int i=0; i<sampleNumber; i++ )
		free(feature[i]);
	free(feature);
}

double computeLassoResult( double **Wd, double *sample, double *feature, double lambda, int sampleDim, int featureDim ){

	double LassoResult = 0;
	double residuals;
	for( unsigned int i=0; i<sampleDim; i++ ){
		residuals = -sample[i];
		for( unsigned int j=0; j<featureDim; j++ )
			residuals += Wd[j][i]*feature[j];

		LassoResult += residuals*residuals;
	}

	double sum_feature = 0;
	for( unsigned int j=0; j<featureDim; j++ )
		sum_feature += getAbs(feature[j]);

    	return 0.5*LassoResult+lambda*sum_feature;
}


void calculateError(  double **Wd,  double **sample, double **feature, double lambda, int sampleNumber, int sampleDim, int featureDim ) {

	double TotalDecError = 0;
	for( unsigned int t=0; t<sampleNumber; t++ ){
		TotalDecError += computeLassoResult( Wd, sample[t], feature[t], lambda, sampleDim, featureDim);
	}
	TotalDecError /= sampleNumber;
	std::cout<<"Total Decode Error is "<<TotalDecError<<std::endl;
}

void trainDecoder( double **Wd, double **feature, double **sample, double lambda, int layers, int featureDim, int sampleNumber, int sampleDim, int iterationNumber, bool NonNegative ){

	double *residuals = (double*)malloc(sampleDim*sizeof(double));
	double *A = Initialize_A( featureDim );
	double *A_Copy = Initialize_A_Copy( featureDim );
	std::vector<int> *nonZeroIndex = NonZeroIndexInitialization( sampleNumber );

	srand((unsigned)time(0));
	myseed = (unsigned int) RAND_MAX * rand();

	std::cout<<"Train decoder"<<std::endl;
	double ComputionalTime = 0;
	clock_t BeginTime = clock();
    	for( unsigned int it=0; it<iterationNumber; it++ ){

		int index = it%sampleNumber;
		if( index==0 )
			Initialize_A( A, A_Copy, featureDim );

		UpdateFeature( Wd, sample[index], residuals, feature[index], nonZeroIndex[index], lambda, layers, featureDim, sampleDim, NonNegative );

		Update_A( A, A_Copy, feature[index], nonZeroIndex[index] );

		UpdateWd( Wd, residuals, feature[index], A, nonZeroIndex[index], sampleDim, NonNegative );
		NormalizeWd( Wd, nonZeroIndex[index], sampleDim );

		if( it%sampleNumber==sampleNumber-1 ){
            		std::cout<<it+1<<" iterations finished"<<std::endl;
			//calculateError(Wd, sample, feature, lambda, sampleNumber, sampleDim, featureDim);
		}
	}
	clock_t EndTime = clock();
	ComputionalTime = (double)(EndTime-BeginTime)/CLOCKS_PER_SEC;

   	std::cout<<"Finish decoding process:"<<std::endl;
	std::cout<<"Train Decode Time is "<<ComputionalTime<<" seconds."<<std::endl;
	calculateError( Wd, sample, feature, lambda, sampleNumber, sampleDim, featureDim );
	free(A_Copy);
	free(A);
	free(residuals);
	delete [] nonZeroIndex;
}


}

#endif /* Sparse Coordinate Coding */
