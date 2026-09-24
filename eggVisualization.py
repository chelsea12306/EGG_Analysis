"""
eggVisualization.py - 可视化模块

绘制以下图表：
  1. EGG 原始信号 vs 滤波后信号（以 ID1 为例）
  2. FFT 功率谱 + 主频标注
  3. 全体受试者空腹 vs 餐后主频柱状图
  4. t 检验结果汇总表

依赖：matplotlib（Google Colab 预装）
"""

import numpy as np
from scipy import signal as sig
import matplotlib.pyplot as plt
import os

# ==================== Parameters ====================
fs = 2
N = 4096
DATA_DIR = 'EGG-database'

b, a = sig.butter(3, [0.03, 0.25], btype='bandpass', fs=fs)

# ==================== 图1: 原始信号 vs 滤波后 ====================
dat = np.loadtxt(os.path.join(DATA_DIR, 'ID1_fasting.txt'))
t = np.arange(2400) / fs  # 时间轴（秒）

raw_ch1 = dat[:2400, 0]
filt_ch1 = sig.filtfilt(b, a, raw_ch1)

fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
axes[0].plot(t, raw_ch1, color='#888888', linewidth=0.8)
axes[0].set_title('ID1 Fasting - Channel 1 (Raw)', fontsize=13)
axes[0].set_ylabel('Amplitude')
axes[0].grid(True, alpha=0.3)

axes[1].plot(t, filt_ch1, color='#2563eb', linewidth=1.2)
axes[1].set_title('ID1 Fasting - Channel 1 (Filtered 0.03-0.25 Hz)', fontsize=13)
axes[1].set_ylabel('Amplitude')
axes[1].set_xlabel('Time (s)')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fig1_signal.png', dpi=150)
plt.show()
print('Figure 1 saved: fig1_signal.png')

# ==================== 图2: FFT 功率谱 ====================
fft_power = np.abs(np.fft.fft(filt_ch1, N))**2
fft_half = fft_power[:N//2 + 1]
freqs = np.arange(len(fft_half)) / 2048 * 60  # 转换为 cpm
peak_idx = np.argmax(fft_half)
peak_cpm = (peak_idx + 1) / 2048 * 60

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(freqs, fft_half, color='#2563eb', linewidth=1.2)
ax.axvline(peak_cpm, color='#dc2626', linestyle='--', linewidth=1.5, label=f'DF = {peak_cpm:.2f} cpm')
ax.set_xlim(0, 10)
ax.set_title('ID1 Fasting CH1 - FFT Power Spectrum', fontsize=13)
ax.set_xlabel('Frequency (cpm)')
ax.set_ylabel('Power')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig2_fft.png', dpi=150)
plt.show()
print(f'Figure 2 saved: fig2_fft.png  (DF = {peak_cpm:.2f} cpm)')

# ==================== 图3: 全体 DF 柱状图 ====================
# 先重新计算 DF（复用 eggAnalysis.py 的逻辑）
df = np.zeros((20, 6))
for ind in range(1, 21):
    dat_f = np.loadtxt(os.path.join(DATA_DIR, f'ID{ind}_fasting.txt'))
    dat_p = np.loadtxt(os.path.join(DATA_DIR, f'ID{ind}_postprandial.txt'))
    for ch in range(3):
        cf = sig.filtfilt(b, a, dat_f[:2400, ch])
        cp = sig.filtfilt(b, a, dat_p[:2400, ch])
        pf = np.abs(np.fft.fft(cf, N))**2; ff = pf[:N//2+1]
        pp = np.abs(np.fft.fft(cp, N))**2; fp = pp[:N//2+1]
        df[ind-1, ch] = (np.argmax(ff) + 1) / 2048
        df[ind-1, ch+3] = (np.argmax(fp) + 1) / 2048
df = df * 60

# 手动修正
df[3, 3] = 3.1934; df[3, 4] = 3.1348; df[3, 5] = 3.1348
df[5, 2] = 2.4900
df[14, 3] = 2.2560
df[16, 0] = 2.9592; df[16, 1] = 3.0469; df[16, 2] = 3.0762

# 绘图：CH1 空腹 vs 餐后
fig, ax = plt.subplots(figsize=(14, 5))
x = np.arange(20)
width = 0.35
ax.bar(x - width/2, df[:, 0], width, label='Fasting', color='#60a5fa')
ax.bar(x + width/2, df[:, 3], width, label='Postprandial', color='#f87171')
ax.set_xticks(x)
ax.set_xticklabels([f'ID{i+1}' for i in range(20)], fontsize=9)
ax.set_ylabel('Dominant Frequency (cpm)')
ax.set_title('Dominant Frequency: Fasting vs Postprandial (Channel 1)', fontsize=13)
ax.legend(fontsize=11)
ax.axhline(3.0, color='#888888', linestyle=':', alpha=0.5, label='Normal (3 cpm)')
ax.grid(True, alpha=0.2, axis='y')
plt.tight_layout()
plt.savefig('fig3_df_bar.png', dpi=150)
plt.show()
print('Figure 3 saved: fig3_df_bar.png')

# ==================== 图4: t 检验汇总表 ====================
from scipy.stats import ttest_rel

# All subjects
t1, p1 = ttest_rel(df[:, 0], df[:, 3])
t2, p2 = ttest_rel(df[:, 1], df[:, 4])
t3, p3 = ttest_rel(df[:, 2], df[:, 5])

low_idx = [0, 1, 2, 8, 11, 13, 15, 16, 17, 18]
df_low = df[low_idx, :]
t4, p4 = ttest_rel(df_low[:, 0], df_low[:, 3])
t5, p5 = ttest_rel(df_low[:, 1], df_low[:, 4])
t6, p6 = ttest_rel(df_low[:, 2], df_low[:, 5])

high_idx = [3, 4, 5, 6, 7, 9, 10, 12, 14, 19]
df_high = df[high_idx, :]
t7, p7 = ttest_rel(df_high[:, 0], df_high[:, 3])
t8, p8 = ttest_rel(df_high[:, 1], df_high[:, 4])
t9, p9 = ttest_rel(df_high[:, 2], df_high[:, 5])

groups = ['All (n=20)', 'All (n=20)', 'All (n=20)',
          'Low BMI (n=10)', 'Low BMI (n=10)', 'Low BMI (n=10)',
          'High BMI (n=10)', 'High BMI (n=10)', 'High BMI (n=10)']
channels = ['CH1', 'CH2', 'CH3'] * 3
t_vals = [t1, t2, t3, t4, t5, t6, t7, t8, t9]
p_vals = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
sig = ['***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.' for p in p_vals]

fig, ax = plt.subplots(figsize=(10, 4))
ax.axis('off')
col_labels = ['Group', 'Channel', 't-value', 'p-value', 'Sig.']
table_data = [[groups[i], channels[i], f'{t_vals[i]:.3f}', f'{p_vals[i]:.6f}', sig[i]] for i in range(9)]
colors = [['#dbeafe' if p < 0.05 else '#f3f4f6' for _ in range(5)] for p in p_vals]
tbl = ax.table(cellText=table_data, colLabels=col_labels, cellColours=colors, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(11)
tbl.scale(1.2, 1.8)
ax.set_title('Paired t-test Results (Fasting vs Postprandial)', fontsize=14, pad=20)
plt.tight_layout()
plt.savefig('fig4_ttest.png', dpi=150, bbox_inches='tight')
plt.show()
print('Figure 4 saved: fig4_ttest.png')
print('\nDone! All figures saved.')
