import pandas as pd
import matplotlib.pyplot as plt

# Read CSV file
df = pd.read_csv('heat_conduction_results.csv')

# Plotting
for index, row in df.iterrows():
    plt.plot(row[::2], row[1::2], label=f'Time step {index}')

plt.title('1D Unsteady Heat Conduction')
plt.xlabel('x')
plt.ylabel('Temperature')
plt.legend()
plt.show()
