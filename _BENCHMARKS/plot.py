import matplotlib.pyplot as plt

# Specified matrix sizes as x-coordinates
x_points = ['64x64', '128x128', '256x256', '512x512', '1024x1024', '2048x2048', '4096x4096', '6000x6000']

# Given y-coordinates for the green and blue points
y_points_green = [1.25, 1.64, 3.12, 4.49, 6.53, 6.178, 11.33, 8.61]
y_points_blue = [1.91, 3.21, 4.22, 5.14, 6.68, 6.80, 15.31, 9.46]

plt.figure(figsize=(10, 6))
plt.plot(x_points, y_points_green, 'go-', label='Bernsten')
plt.plot(x_points, y_points_blue, 'bo-', label='Cannon')

# Adding the straight lines y = 1 and y = 8
plt.axhline(y=1, color='r', linestyle='-', label='Minimum Speedup')
plt.axhline(y=8, color='r', linestyle='-', label='Maxinum Speedup')

plt.xlabel('Size of matrices')
plt.ylabel('Speedup')
plt.title('Plot of Bernsten and Cannon Speedups')
plt.legend()
plt.grid(True)

plt.savefig('plot.svg', format='svg')



plt.show()
