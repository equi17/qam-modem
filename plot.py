import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import erfc

def qpsk_ber_theoretical(snr_db):
    snr = 10**(snr_db/10)
    return 0.5 * erfc(np.sqrt(snr))

def qam16_ber_theoretical(snr_db):
    snr = 10**(snr_db/10)
    return 0.375 * erfc(np.sqrt(snr/5))

def qam64_ber_theoretical(snr_db):
    snr = 10**(snr_db/10)
    return (7/12) * erfc(np.sqrt(snr/21))

qpsk_data = pd.read_csv('qpsk_ber.csv')
qam16_data = pd.read_csv('qam16_ber.csv')
qam64_data = pd.read_csv('qam64_ber.csv')

qpsk_data = qpsk_data[qpsk_data.iloc[:,1] > 0]
qam16_data = qam16_data[qam16_data.iloc[:,1] > 0]
qam64_data = qam64_data[qam64_data.iloc[:,1] > 0]

plt.figure(figsize=(12, 8))

plt.semilogy(qpsk_data.iloc[:,0], qpsk_data.iloc[:,1], 'b-o', label='QPSK (Simulated)', markersize=5)
plt.semilogy(qam16_data.iloc[:,0], qam16_data.iloc[:,1], 'r-s', label='16-QAM (Simulated)', markersize=5)
plt.semilogy(qam64_data.iloc[:,0], qam64_data.iloc[:,1], 'g-d', label='64-QAM (Simulated)', markersize=5)

snr_range = np.linspace(min(qpsk_data.iloc[:,0].min(), qam16_data.iloc[:,0].min(), qam64_data.iloc[:,0].min()),
                        max(qpsk_data.iloc[:,0].max(), qam16_data.iloc[:,0].max(), qam64_data.iloc[:,0].max()), 100)

plt.semilogy(snr_range, qpsk_ber_theoretical(snr_range), 'b--', label='QPSK (Theoretical)', linewidth=2)
plt.semilogy(snr_range, qam16_ber_theoretical(snr_range), 'r--', label='16-QAM (Theoretical)', linewidth=2)
plt.semilogy(snr_range, qam64_ber_theoretical(snr_range), 'g--', label='64-QAM (Theoretical)', linewidth=2)

plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('Bit Error Rate (BER)', fontsize=12)
plt.title('BER Performance of QPSK, 16-QAM and 64-QAM', fontsize=14)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(fontsize=10, loc='upper right')
plt.ylim(1e-6, 1)
plt.xlim(snr_range[0], snr_range[-1])

plt.savefig('ber_comparison.png', dpi=300, bbox_inches='tight')
print("График сохранён в ber_comparison.png")