#!/bin/bash
echo -e "\n66" 
echo -e "\n66" > cannon_9_cores
for i in {1..300}
do
    mpirun --oversubscribe --np 9 ../cannon matrix66 matrix66 matrixC66 >> cannon_9_cores 
done

echo -e "\n129" >> cannon_9_cores
echo -e "\n129" 
for i in {1..300}
do
    mpirun --oversubscribe --np 9 ../cannon matrix129 matrix129 matrixC129 >> cannon_9_cores
done

echo -e "\n258" >> cannon_9_cores
echo -e "\n258" 
for i in {1..300}
do
    mpirun --oversubscribe --np 9  ../cannon matrix258 matrix258 matrixC258 >> cannon_9_cores
done

echo -e "\n513" >> cannon_9_cores
echo -e "\n513" 
for i in {1..300}
do
    mpirun --oversubscribe --np 9  ../cannon matrix513 matrix513 matrixC513 >> cannon_9_cores
done


echo -e "\n1026" >> cannon_9_cores
echo -e "\n1026" 
for i in {1..50}
do
    mpirun --oversubscribe --np 9 ../cannon matrix1026 matrix1026 matrixC1026 >> cannon_9_cores
done

echo -e "\n2049" >> cannon_9_cores
echo -e "\n2049" 
for i in {1..10}
do
    mpirun --oversubscribe --np 9 ../cannon matrix2049 matrix2049 matrixC2049 >> cannon_9_cores
done

echo -e "\n4098" >> cannon_9_cores
echo -e "\n4098" 
for i in {1..4}
do
    mpirun --oversubscribe --np 9 ../cannon matrix4098 matrix4098 matrixC4098 >> cannon_9_cores
done

echo -e "\n6000" >> cannon_9_cores
echo -e "\n6000" 
for i in {1..3}
do
    mpirun --oversubscribe --np 9 ../cannon matrix6000 matrix6000 matrixC6000 >> cannon_9_cores
done

echo -e "\n8190" >> cannon_9_cores
echo -e "\n8190" 
for i in {1..2}
do
    mpirun --oversubscribe --np 9 ../cannon matrix8190 matrix8190 matrixC8190 >> cannon_9_cores
done
