get_ipython().run_line_magic('run', '-i common.py')

def find_epsilon(reds, mstar, rho_star, exfac=1.0):
   # take measured rho_star, m_star & find the \epsilon that it implies.
   # exfac allows you to change normalization of rho(>M) at fixed M (for EDE)
   from scipy.optimize import fsolve
   my_mf.update(z=reds)
   spl_mrho=US(log(my_mf.m/h0), log(my_mf.rho_gtm*h0**2), s=0)
   def my_func(eps, exfac=exfac):
      return eps-rho_star/fbary/exp(spl_mrho(log(mstar/eps/fbary)))/exfac
   
   output=fsolve(my_func, 1.0, args=(exfac))
   return output

# z=9
xv9=rhostar_data_z9['mstar'][-1].copy()
yv9=rhostar_data_z9['rhostar'][-1].copy()/vol_z9

# z=7.5:
xv7=rhostar_data_z8['mstar'][-1].copy()
yv7=rhostar_data_z8['rhostar'][-1].copy()/vol_z8



eps_z7=find_epsilon(zz8, xv7, yv7)
eps_z9=find_epsilon(zz9, xv9, yv9)
print('implied epsilon at z=7.5: {0:.2f}'.format(eps_z7[0]))
print('implied epsilon at z=9: {0:.2f}'.format(eps_z9[0]))

# considering errors:
yv9_low=rhostar_data_z9['rhostar_lo'][-1].copy()/vol_z9
eps_z9_low=find_epsilon(zz9, xv9, yv9_low)
print('implied epsilon at z=9 (1 sig low): {0:.2f}'.format(eps_z9_low[0]))


#######################################################################
# also compute implied epsilon for EDE at fixed M_star:

my_mstar=xv9

my_mf.update(z=zz9)

my_mf_ede=hmf.MassFunction(Mmin=min_mval_log10, hmf_model='SMT',
                           cosmo_model=planck2020_ede,
                           sigma_8=sig8_ede, n=ns_ede,
                           transfer_model=hmf.density_field.transfer_models.CAMB,
                           transfer_params={'extrapolate_with_eh':True})
my_mf_ede.update(z=zz9)

print('are the mass arrays equal?', numpy.allclose(my_mf.m, my_mf_ede.m))


# spline from m_halo to rho_bary(>f_bary mhalo):
spl_rho_pl=US(log(my_mf.m/h0*fbary), log(fbary*my_mf.rho_gtm*h0**2), s=0)
spl_rho_ede=US(log(my_mf_ede.m/h0_ede*fbary_ede),
               log(fbary_ede*my_mf_ede.rho_gtm*h0_ede**2), s=0)

ede_planck_ratio=exp(spl_rho_ede(log(my_mstar))-spl_rho_pl(log(my_mstar)))

# ######################################################################
# version that just uses nearest neighbor interpolation, for comparison:
# ######################################################################
# ind_planck=my_mf.m.searchsorted(my_mstar*h0/fbary)
# ind_ede=my_mf_ede.m.searchsorted(my_mstar*h0_ede/fbary_ede)
#  
# ede_planck_ratio=((fbary_ede*my_mf_ede.rho_gtm[ind_ede]*h0_ede**2)/
#                   (fbary*my_mf.rho_gtm[ind_planck]*h0**2))
# ######################################################################

print('*'*70)
print('for EDE (where rho(M_halo) is factor of {0:.2f} higher at z=9):'.format(ede_planck_ratio))
# need to update data points to use volume associated with EDE cosmology, not Planck:
eps_ede=find_epsilon(zz9, xv9, yv9/vol_z9_ede*vol_z9, exfac=ede_planck_ratio)[0]
print('implied epsilon at z=9 is {0:.2f}'.format(eps_ede))
