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
corr_t_e = 100.0
gap_t_e = 5.0

#adaptation
a_w = Float(4.0)
b_w = Float(80.5e-3)

#random inputs
wn_dt = Float(1e-3)
μ_wn = Float(62.5e-3)
σ_wn = Float(21.4e-3)
ξext_e = Float(20e-3)
τext_e = Float(400e-3)
ξext_i = Float(20e-3)
τext_i = Float(300e-3)

#Networks
Ne = 1200
Ni = 240
p_ee = 0.17/3.0
p_ei = 0.57/3.0
p_ie = 0.61/3.0
p_ii = 0.32/3.0
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
μTee = Float(297.0)
μTei = Float(1012.0)
μTie = Float(643.0)
μTii = Float(373.0)
σTee = Float(0.0)
σTei = Float(0.0)
σTie = Float(0.0)
σTii = Float(0.0)

#Weight normalization
min_SN_correction = Float(0.01)
