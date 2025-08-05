#!/bin/bash
echo -e "\n64"
echo -e "\n64" > bernsten_8_cores
for i in {1..300}
do
    mpirun -np 8 ../bernsten matrix64 matrix64 matrixC64 >> bernsten_8_cores
done

echo -e "\n128" >> bernsten_8_cores
echo -e "\n128" 
for i in {1..300}
do
    mpirun -np 8 ../bernsten matrix128 matrix128 matrixC128 >> bernsten_8_cores
done

echo -e "\n256" >> bernsten_8_cores
echo -e "\n256" 
for i in {1..300}
do
    mpirun -np 8 ../bernsten matrix256 matrix256 matrixC256 >> bernsten_8_cores
done

echo -e "\n512" >> bernsten_8_cores
echo -e "\n512" 
for i in {1..300}
do
    mpirun -np 8 ../bernsten matrix512 matrix512 matrixC512 >> bernsten_8_cores
done


echo -e "\n1024" >> bernsten_8_cores
echo -e "\n1024" 
for i in {1..50}
do
    mpirun -np 8 ../bernsten matrix1024 matrix1024 matrixC1024 >> bernsten_8_cores
done

echo -e "\n2048" >> bernsten_8_cores
echo -e "\n2048" 
for i in {1..10}
do
    mpirun -np 8 ../bernsten matrix2048 matrix2048 matrixC2048 >> bernsten_8_cores
done

echo -e "\n4096" >> bernsten_8_cores
echo -e "\n4096" 
for i in {1..4}
do
    mpirun -np 8 ../bernsten matrix4096 matrix4096 matrixC4096 >> bernsten_8_cores
done

echo -e "\n6000" >> bernsten_8_cores
echo -e "\n6000" 
for i in {1..3}
do
    mpirun -np 8 ../bernsten matrix6000 matrix6000 matrixC6000 >> bernsten_8_cores
done

echo -e "\n8190" >> bernsten_8_cores
echo -e "\n8190" 
for i in {1..2}
do
    mpirun -np 8 ../bernsten matrix8190 matrix8190 matrixC8190 >> bernsten_8_cores
done
