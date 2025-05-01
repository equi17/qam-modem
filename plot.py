import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv('ber.csv')
data = data[data.iloc[:,1] > 0]

plt.figure(figsize=(10, 6))
plt.semilogy(data.iloc[:,0], data.iloc[:,1], 'b-o')
plt.xlabel('Дисперсия шума $\\sigma^2$', fontsize=12)
plt.ylabel('Вероятность ошибки на бит (BER)', fontsize=12)
plt.title('Зависимость BER от дисперсии шума при QAM16', fontsize=14)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.savefig('ber_plot.png', dpi=300, bbox_inches='tight')
print("График сохранён в ber_plot.png")