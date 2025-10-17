const Float = Float32

#Neurons
#membrane potential
Cm = Float(240e-3)  #*1e9
τe = Float(1.1e-3)
τi = Float(1.1e-3)
τw = Float(144e-3)
gL = Float(4.19)
Vpeak = Float(0)
Ee = Float(10e-3)
Ei = Float(-75e-3)
El = Float(-70.6e-3)
V_reset_e = Float(-60e-3)
V_reset_i = Float(-60e-3)
ΔT = Float(2.0e-3)

#Firing threshold
μθe = Float(-50.4e-3)
μθi = Float(-50.4e-3)
σθe = Float(0e-3)
σθi = Float(0e-3)
rt_e = 0.45
rt_i = 2.7  #target firing rates
ηIP_e = 0.05e-3
ηIP_i = 0

#adaptation
a_w = Float(4.0)
b_w = Float(80.5e-3)

#random inputs
wn_dt = Float(1e-3)
μ_wn = Float(62.5e-3)
σ_wn = Float(21.4e-3)
ξext = Float(1e-3)
τext = Float(4e-3)

#Networks
Ne = 1200
Ni = 240
p_ee = 0.06
p_ei = 0.2
p_ie = 0.2
p_ii = 0.12
d_ee = Float(1.5e-3)
d_ei = Float(1.5e-3)
d_ie = Float(1.5e-3)
d_ii = Float(1.5e-3)

#STDP
AeLTP = Float(1.6)
AeLTD = Float(-0.32)
Aipre = Float(1.6)
Aipost = Float(1.6)
τLTP = Float(15e-3)
τLTD = Float(30e-3)
τpre = Float(15e-3)
τpost = Float(15e-3)
LTDalpha = 2*(Aipre*τpre + Aipost*τpost)*rt_e

#Initial weights
w0ee = Float(4.2)
w0ei = Float(4.2)
w0ie = Float(13.5)
w0ii = Float(13.0)

#Weight normalization
min_SN_correction = Float(0.01)

tstim = 3.1
stim_period = 0.4
