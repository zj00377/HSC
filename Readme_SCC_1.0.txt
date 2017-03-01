Stochastic Coordinate Coding 

version 1.0

Authors: Qingyang Li  qingyang.li@asu.edu
	 Binbin Lin   binbinlinzju@gmail.com
	 

Input:

1. Sample
Input sample is a matrix, and its size is MxN.
M is the number of elements in one sample.
N is the total number of samples. 

Output:
1. Sparse code
Ouput sparse code is a matrix, and its size is KxN.
K is the number of elements in one sparse code.
N is the number of sparse codes, equal to the number of samples.

2. Dictionary
Output dictionary is a matirx, and its size is MxK.
M is the number of elements in one sample.
K is the number of bases in the dictionary, i.e., the number of elements in one sparse code.

Variables to be configured:

File name:
1. SampleFileName  ( the name of sample file )

2. FeatureFileName ( the name of generated sparse code )

3. initializedDictionaryName ( the name of initialial dictionary )

4. savedDictionaryName ( the name of learned dictionary )

Variable number:
1. featureDim ( number of elements in each sparse code )

2. sampleDim( number of elements in each sample )

3. layers ( number of layers in traing process, default value is 3 )

4. epochNumber ( number of epochs )

5. lambda ( the value of lambda, default is 0.13 )

6. DictionaryGenerationState ( true/false )
   true: generate an initial dictionary 
   false:  use the initial dictionary provided by the user, the name of initialized dictionary file should be the value of initializedDictionaryName.

7. NonNegative ( true/false )
   true: the sample, dictionary and sparse code are non-negative.
   false: no constraint on the sign of the sample, dictionary and sparse code.



Operating Environment:
Linux (ubuntu/red hat)
g++ 4.6 and above version


How to run the program:

Command to compile the program:  g++ run.cpp -o run -O3
command to run the program:   ./run



