"""
eggAnalysis.py - Python version of eggAnalysis.m

Performs EGG (electrogastrography) signal processing:
  - Bandpass filtering (3rd-order Butterworth, 0.03-0.25 Hz)
  - FFT-based dominant frequency (DF) calculation
  - Manual peak corrections
  - Paired-sample t-tests (fasting vs postprandial, by BMI group)

Input:  EGG-database/*.txt  (20 subjects x 2 conditions x 3 channels)
Output: df.csv  (dominant frequencies in cycles-per-minute)

Original Octave code by:
  Nenad B. Popovic (nenad.pop92@gmail.com)
  Nadica Miljkovic (nadica.miljkovic@etf.rs)

References:
  [1] Popovic, N.B., Miljkovic, N. and Popovic M.B., 2019. Simple gastric
      motility assessment method with a single-channel electrogastrogram.
      Biomedical Engineering/Biomedizinsche Technik, 64(2), pp.177-185.
  [2] Popovic, N.B., Miljkovic, N. and Popovic M.B., 2020. Three-channel
      surface EGG dataset recorded during fasting and post-prandial states
      in 20 healthy individuals. Zenodo, doi: 10.5281/zenodo.3730617.
"""

import numpy as np
from scipy import signal
from scipy.stats import ttest_rel
import os

# ==================== Parameters ====================
fs = 2          # Hz, sampling frequency
N = 4096        # number of points for FFT analysis

DATA_DIR = 'EGG-database'   # adjust if needed

# ==================== Bandpass filter ====================
# 3rd-order Butterworth bandpass, [0.03, 0.25] Hz
b, a = signal.butter(3, [0.03, 0.25], btype='bandpass', fs=fs)

# ==================== Dominant Frequency (DF) calculation ====================
df = np.zeros((20, 6))   # columns: [fastCH1, fastCH2, fastCH3, postCH1, postCH2, postCH3]

for ind in range(1, 21):          # subject ID 1..20
    # --- FASTING ---
    file_f = os.path.join(DATA_DIR, f'ID{ind}_fasting.txt')
    dat_f = np.loadtxt(file_f)

    ch1_f = signal.filtfilt(b, a, dat_f[:2400, 0])
    ch2_f = signal.filtfilt(b, a, dat_f[:2400, 1])
    ch3_f = signal.filtfilt(b, a, dat_f[:2400, 2])

    # FFT power spectrum (one-sided)
    p1f = np.abs(np.fft.fft(ch1_f, N))**2; fft1f = p1f[:N//2 + 1]
    p2f = np.abs(np.fft.fft(ch2_f, N))**2; fft2f = p2f[:N//2 + 1]
    p3f = np.abs(np.fft.fft(ch3_f, N))**2; fft3f = p3f[:N//2 + 1]

    # DF = peak bin index / 2048  (+1 because Octave is 1-indexed)
    df[ind - 1, 0] = (np.argmax(fft1f) + 1) / 2048
    df[ind - 1, 1] = (np.argmax(fft2f) + 1) / 2048
    df[ind - 1, 2] = (np.argmax(fft3f) + 1) / 2048

    # --- POSTPRANDIAL ---
    file_p = os.path.join(DATA_DIR, f'ID{ind}_postprandial.txt')
    dat_p = np.loadtxt(file_p)

    ch1_p = signal.filtfilt(b, a, dat_p[:2400, 0])
    ch2_p = signal.filtfilt(b, a, dat_p[:2400, 1])
    ch3_p = signal.filtfilt(b, a, dat_p[:2400, 2])

    p1p = np.abs(np.fft.fft(ch1_p, N))**2; fft1p = p1p[:N//2 + 1]
    p2p = np.abs(np.fft.fft(ch2_p, N))**2; fft2p = p2p[:N//2 + 1]
    p3p = np.abs(np.fft.fft(ch3_p, N))**2; fft3p = p3p[:N//2 + 1]

    df[ind - 1, 3] = (np.argmax(fft1p) + 1) / 2048
    df[ind - 1, 4] = (np.argmax(fft2p) + 1) / 2048
    df[ind - 1, 5] = (np.argmax(fft3p) + 1) / 2048

    print(f'  ID{ind} processed')

# Convert Hz -> cpm (cycles per minute)
df = df * 60

# ==================== Manual corrections (visual inspection) ====================
df[3, 3] = 3.1934    # ID4  postprandial CH1
df[3, 4] = 3.1348    # ID4  postprandial CH2
df[3, 5] = 3.1348    # ID4  postprandial CH3
df[5, 2] = 2.4900     # ID6  fasting     CH3
df[14, 3] = 2.2560   # ID15 postprandial CH1
df[16, 0] = 2.9592   # ID17 fasting     CH1
df[16, 1] = 3.0469   # ID17 fasting     CH2
df[16, 2] = 3.0762   # ID17 fasting     CH3

# ==================== Save DF values ====================
np.savetxt('df.csv', df, delimiter=',')
print('\nDF values saved to df.csv')
print('\nDF matrix (cpm):')
print('  Columns: fastCH1  fastCH2  fastCH3  postCH1  postCH2  postCH3')
for i in range(20):
    print(f'  ID{i+1:2d}:  ' + '  '.join(f'{v:7.4f}' for v in df[i]))

# ==================== Paired-sample t-tests ====================
print('\n' + '='*60)
print('Paired-sample t-tests (fasting vs postprandial)')
print('='*60)

# All subjects ID1-ID20
t1, p1 = ttest_rel(df[:, 0], df[:, 3])
t2, p2 = ttest_rel(df[:, 1], df[:, 4])
t3, p3 = ttest_rel(df[:, 2], df[:, 5])

# Lower BMI subjects: ID1,2,3,9,12,14,16,17,18,19
low_idx = [0, 1, 2, 8, 11, 13, 15, 16, 17, 18]
df_low = df[low_idx, :]
t4, p4 = ttest_rel(df_low[:, 0], df_low[:, 3])
t5, p5 = ttest_rel(df_low[:, 1], df_low[:, 4])
t6, p6 = ttest_rel(df_low[:, 2], df_low[:, 5])

# Higher BMI subjects: ID4,5,6,7,8,10,11,13,15,20
high_idx = [3, 4, 5, 6, 7, 9, 10, 12, 14, 19]
df_high = df[high_idx, :]
t7, p7 = ttest_rel(df_high[:, 0], df_high[:, 3])
t8, p8 = ttest_rel(df_high[:, 1], df_high[:, 4])
t9, p9 = ttest_rel(df_high[:, 2], df_high[:, 5])

p_vals = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
h_vals = [1 if p < 0.05 else 0 for p in p_vals]

print(f'\n--- All subjects (ID1-ID20) ---')
print(f'  CH1: t={t1:8.4f}  p={p1:.6f}  h={h_vals[0]}')
print(f'  CH2: t={t2:8.4f}  p={p2:.6f}  h={h_vals[1]}')
print(f'  CH3: t={t3:8.4f}  p={p3:.6f}  h={h_vals[2]}')

print(f'\n--- Lower BMI subjects ---')
print(f'  CH1: t={t4:8.4f}  p={p4:.6f}  h={h_vals[3]}')
print(f'  CH2: t={t5:8.4f}  p={p5:.6f}  h={h_vals[4]}')
print(f'  CH3: t={t6:8.4f}  p={p6:.6f}  h={h_vals[5]}')

print(f'\n--- Higher BMI subjects ---')
print(f'  CH1: t={t7:8.4f}  p={p7:.6f}  h={h_vals[6]}')
print(f'  CH2: t={t8:8.4f}  p={p8:.6f}  h={h_vals[7]}')
print(f'  CH3: t={t9:8.4f}  p={p9:.6f}  h={h_vals[8]}')

print(f'\np-values: {p_vals}')
print(f'h (reject H0 at alpha=0.05): {h_vals}')
