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

#adaptation
a_w = Float(4.0)
b_w = Float(80.5e-3)

#random inputs
wn_dt = Float(1e-3)
μ_wn = Float(62.5e-3)
σ_wn = Float(21.4e-3)

#stimulation
tstim = [0.4, 0.8, 1.2, 1.6]
Ne = 1000
effw = 1
V0 = -63e-3
