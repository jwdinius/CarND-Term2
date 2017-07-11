"""
# script to read and plot output from udacity self-driving car nanodegree UKF project
# to generate the needed file, run (from project root level) "./build/UnscentedKF out.txt"

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

csv_file = "out.txt"

df = pd.read_csv(csv_file)

vxe_np = np.array(df["vx_est"])
vye_np = np.array(df["vy_est"])
xe_np  = np.array(df["x_est"])
ye_np  = np.array(df["y_est"])
vxt_np = np.array(df["vx_true"])
vyt_np = np.array(df["vy_true"])
xt_np  = np.array(df["x_true"])
yt_np  = np.array(df["y_true"])
nis_l  = []
nis_r  = []
for i,stype in enumerate(df["sensor"]):
	if (stype == 'L'):
		nis_l.append(df["nis"][i])
	elif (stype == 'R'):
		nis_r.append(df["nis"][i])
nis_l_np = np.array(nis_l)
nis_r_np = np.array(nis_r)
k_np     = np.arange(vxe_np.shape[0])

plt.figure()
plt.subplot(2,2,1)
plt.plot(k_np, xt_np, 'r--', label='truth')
plt.plot(k_np, xe_np, 'b--', label='estimate')
plt.legend(loc='best')
plt.title("X pos vs. timestep")
plt.ylim((-30,30))
plt.ylabel("X pos(m)")

plt.subplot(2,2,2)
plt.plot(k_np, yt_np, 'r--', label='truth')
plt.plot(k_np, ye_np, 'b--', label='estimate')
plt.title("Y pos vs. timestep")
plt.ylim((-30,30))
plt.ylabel("Y pos(m)")

plt.subplot(2,2,3)
plt.plot(k_np, vxt_np, 'r--', label='truth')
plt.plot(k_np, vxe_np, 'b--', label='estimate')
plt.title("X vel vs. timestep")
plt.ylim((-10,10))
plt.ylabel("X vel(m/s)")
plt.xlabel("timestep")

plt.subplot(2,2,4)
plt.plot(k_np, vyt_np, 'r--', label='truth')
plt.plot(k_np, vye_np, 'b--', label='estimate')
plt.title("Y vel vs. timestep")
plt.ylim((-10,10))
plt.ylabel("Y vel(m/s)")
plt.xlabel("timestep")

plt.savefig("state_comp.png")

plt.figure()
plt.subplot(2,1,1)
plt.plot(nis_l_np,'ro')
chi2_2_95 = 5.99
plt.plot([0,nis_l_np.shape[0]],[chi2_2_95,chi2_2_95],'b--')
plt.title("Laser NIS")
plt.ylabel("NIS")
plt.ylim((0,12))

plt.subplot(2,1,2)
plt.plot(nis_r_np,'ro')
chi2_3_95 = 7.81
plt.plot([0,nis_r_np.shape[0]],[chi2_3_95,chi2_3_95],'b--')
plt.title("Radar NIS")
plt.ylabel("NIS")
plt.ylim((0,12))

plt.savefig("nis.png")
