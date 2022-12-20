import osiris
import matplotlib.pyplot as plt
import numpy as np
mydata = osiris.RamsesData(-1,scale="au",verbose=True)
mumh='1.67e-24*'+str(mydata.info['mu_gas'])
mydata.new_field(name="Nh",operation="density/"+mumh,unit="cm$^{-3}$",label="N_H")

Z_15 = mydata.get("passive_scalar_cons_14",only_leafs=True)
Z_30 = mydata.get("passive_scalar_cons_31",only_leafs=True)
Z_100= mydata.get("passive_scalar_cons_41",only_leafs=True)
Z_250= mydata.get("passive_scalar_cons_50",only_leafs=True)

n_i= mydata.get("passive_scalar_cons_51",only_leafs=True)
n_e= mydata.get("passive_scalar_cons_52",only_leafs=True)
nH = mydata.get("Nh",only_leafs=True)

fig, ax1 = plt.subplots()
ax1.loglog(nH,Z_15,marker='.' , linestyle='none', color='orange',label='15 nm' )
ax1.loglog(nH,Z_30,marker='.' , linestyle='none', color='cyan'  ,label='30 nm' )
ax1.loglog(nH,Z_100,marker='.', linestyle='none', color='green' ,label='100 nm')
ax1.loglog(nH,Z_250,marker='.', linestyle='none', color='purple',label='250 nm')

ax2=ax1.twinx()
ax2.loglog(nH,n_i,marker='.',color='k',linestyle='none',label='ions')
ax2.loglog(nH,n_e,marker='.',color='k',linestyle='none',label='electron')

ax1.set_ylim(1e-4,10)
ax2.set_ylim(1e-4,0.1)

ax1.set_xlabel("$n_{\\mathrm{H}}$ (cm$^{-3}$)")
ax1.set_ylabel("-$Z$")
ax2.set_ylabel("$n_{\\mathrm{i}}$,$n_{\\mathrm{e}}$ (cm$^{-3}$)")

ax1.legend()
ax2.legend()

plt.show()
