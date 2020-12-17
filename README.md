# Kmeans_OpenMP
K-means algorithm implemented in C, paralelized with OpenMP


How to generate random data to be clustered:

python getinput.py k n > dados.txt

where k is the number of clusters and n is the number of points.
Example for k=100 and n=40000:
python getinput.py 100 40000 > dados.txt

How to compile the kmeans parallel:
gcc -o kmeans_par -fopenmp kmeans_par.c -lm -O3 -Wall

How to execute:
./kmeans_par < dados.txt
