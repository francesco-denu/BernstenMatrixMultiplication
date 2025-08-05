#!/bin/bash
echo -e "\n64" 
echo -e "\n64" > sequential
for i in {1..300}
do
    #mpirun -np 1 ../bernsten matrix64 matrix64 matrixC64_real >> bernsten_1_core
    ../seq_mat_mult matrix64 matrix64 matrixC64_real>> sequential
done

echo -e "\n128" >> sequential
echo -e "\n128" 
for i in {1..300}
do
    #mpirun -np 1 ../bernsten matrix128 matrix128 matrixC128_real >> bernsten_1_core
    ../seq_mat_mult matrix128 matrix128 matrixC128_real >> sequential

done

echo -e "\n256" >> sequential
echo -e "\n256" 
for i in {1..300}
do
    #mpirun -np 1 ../bernsten matrix256 matrix256 matrixC256_real >> bernsten_1_core
    ../seq_mat_mult matrix256 matrix256 matrixC256_real >> sequential

done

echo -e "\n512" >> sequential
echo -e "\n512" 
for i in {1..300}
do
    #mpirun -np 1 ../bernsten matrix512 matrix512 matrixC512_real >> bernsten_1_core
    ../seq_mat_mult matrix512 matrix512 matrixC512_real >> sequential

done


echo -e "\n1024" >> sequential
echo -e "\n1024" 
for i in {1..50}
do
    #mpirun -np 1 ../bernsten matrix1024 matrix1024 matrixC1024_real >> bernsten_1_core
    ../seq_mat_mult matrix1024 matrix1024 matrixC1024_real >> sequential

done

echo -e "\n2048" >> sequential
echo -e "\n2048" 
for i in {1..10}
do
    #mpirun -np 1 ../bernsten matrix2048 matrix2048 matrixC2048_real >> bernsten_1_core
    ../seq_mat_mult matrix2048 matrix2048 matrixC2048_real >> sequential

done

echo -e "\n4096" >> sequential
echo -e "\n4096" 
for i in {1..4}
do
    #mpirun -np 1 ../bernsten matrix4096 matrix4096 matrixC4096_real >> bernsten_1_core
    ../seq_mat_mult matrix4096 matrix4096 matrixC4096_real >> sequential

done

echo -e "\n6000" >> sequential
echo -e "\n6000" 
for i in {1..3}
do
    #mpirun -np 1 ../bernsten matrix6000 matrix6000 matrixC6000_real >> bernsten_1_core
    ../seq_mat_mult matrix6000 matrix6000 matrixC6000_real >> sequential

done

echo -e "\n8190" >> sequential
echo -e "\n8190" 
for i in {1..2}
do
    #mpirun -np 1 ../bernsten matrix8190 matrix8190 matrixC8190_real >> bernsten_1_core
    ../seq_mat_mult matrix8190 matrix8190 matrixC8190_real >> sequential
done