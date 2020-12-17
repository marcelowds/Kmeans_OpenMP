# Kmeans_OpenMP
K-means algorithm implemented with OpenMP


How to generate random data to be clustered:
python getinput.py k n > dados.txt
where k is the number of clusters and n is the number of points.

How to compile the kmeans parallel:
gcc -o kmeans_par -fopenmp kmeans_par.c -lm -O3 -Wall

How to execute:
./kmeans_par < dados.txt
