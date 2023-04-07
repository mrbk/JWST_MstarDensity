my_filename=__file__

get_ipython().run_line_magic('run', '-i common.py')

my_cmap_rev=get_cmap('plasma')
my_xlims=r_[20, 5]
my_ylims=r_[1e8, 1e12]

######################################################################
num_z=35
zdat=numpy.linspace(5, 20, num_z)
nuvals=r_[3:7][::-1]
# nu and rho_gt_m map 1-to-1, independent of z:
splnu_rev=US(my_mf.nu**0.5, log(my_mf.rho_gtm*h0**2), s=0)
# get equivalent log(rho_gtm) for each nu:
rho_gtm_nu_log=splnu_rev(nuvals)


######################################################################
# plot z vs M* with lines for fixed cumulative mass density
outdat_nu=numpy.zeros((zdat.size, nuvals.size))
# add data for lower number densities in a different color
for i in range(zdat.size):
    my_mf.update(z=zdat[i])
    rho_gtm_i=my_mf.rho_gtm*h0**2
    ngtm_i=my_mf.ngtm*h0**3
    m_i=my_mf.m/h0
    
    # only do spline to maxind to avoid issues with nonsensically low abundances:
    try:
        maxind=(ngtm_i < 0).nonzero()[0][0]-1
    except IndexError:
        maxind=ngtm_i.size

    if (ngtm_i[maxind-1] == 0):
        maxind=(ngtm_i > 0).nonzero()[0].size

    # spline that goes from log(rho(>M)) to log(M) at this redshift
    myspl_rho=US(log(rho_gtm_i[:maxind][::-1]), log(m_i[:maxind][::-1]), s=0)

    outdat_nu[i]=myspl_rho(rho_gtm_nu_log)
    
# note the data are in terms of m_halo; galaxies are a factor of fbary lower in
# (stellar) mass, so we need to multiply by a factor of fbary to get the
# cumulative number density corresponding to a given galaxy.
my_size=outdat_nu.shape[1]
my_zinds=zdat.searchsorted(r_[18, 15.5, 12.65, 9.6])
for i in range(my_size):
    my_rhodens=exp(splnu_rev(nuvals[i]))*fbary
    my_rhodens_log10=log10(my_rhodens)
    my_exponent=numpy.floor(my_rhodens_log10)
    my_coeff=my_rhodens/10.**my_exponent
    my_ls='-'
    semilogy(zdat, exp(outdat_nu[:,i])*fbary, lw=4,
             color=my_cmap_rev(i/(my_size-0.5)), ls=my_ls,
             label=r'${1:.1f} \times 10^{{{0}}}$'.format(int(my_exponent),
                                                         my_coeff))

    myzind=my_zinds[i]
    if nuvals[i] > 5:
        text(zdat[myzind]-0.5, exp(outdat_nu[myzind,i])*fbary*1.5,
             r'$\nu={0:.0f}$'.format(nuvals[i]),
             color=my_cmap_rev(i/(my_size-0.5)), ha='right')
    else:
        text(zdat[myzind],
             exp(outdat_nu[myzind,i])*fbary*1.2, 
             r'${0:.0f}$'.format(nuvals[i]),
             color=my_cmap_rev(i/(my_size-0.5)), ha='right')


xlabel(r'redshift '+r'$z$')
ylabel(r'$M_{\star}\;\;{\rm or}\;\;\epsilon\,f_{\rm b}\,M_{\rm halo}\;[M_{\odot}]$')


######################################################################
# add data from Labbe et al.:

xvs=mstar_data['z'][:2]
xe_low=xvs-mstar_data['z_low'][:2]
xe_upp=mstar_data['z_hi'][:2]-xvs

yvs=10.**mstar_data['log10_mstar'][:2]
ye_low=yvs-10.**mstar_data['log10_mstar_lo'][:2]
ye_upp=10.**mstar_data['log10_mstar_hi'][:2]-yvs

errorbar(xvs, yvs,
         yerr=c_[ye_low, ye_upp].T,
         xerr=c_[xe_low, xe_upp].T, 
         color=top_color, ecolor=top_color,
         ls='none', capsize=5, elinewidth=2, marker='*', ms=14)


legend(frameon=False, loc=(0.04, 0.56))
text(19.5, 4.2e11, r'$\rho_{\rm b}(>f_{\rm b}^{-1} \epsilon^{-1} \, M_{\star})\;[M_{\odot}\,{\rm Mpc}^{-3}]:$')

# ######################################################################
# add annotation for \epsilon < 1
text(7.5, 1.7*1e9,  r'$\downarrow$', fontsize=30,
     color='k', va='top', ha='center')
text(7.1, 0.9*1.7*1e9, r'$\epsilon < 1$', color='k',
     va='top', ha='left')

xlim(*my_xlims)
ylim(*my_ylims)


plot(19.15, 1.275e10, marker='*', ms=16, color=top_color, linestyle='none', mew=0)
text(18.55, 1.2e10, r'${\rm Labb\'e\;et\;al. [13]}$', color=top_color,
     va='center', ha='left')


savefig(my_filename.split('.')[0] + '.pdf')












######################################################################
######################################################################
# also want to output, at z=9.9 and z=7.6 (the Labbe et al. points), the
# comoving number densities for \epsilon = 0.316 and 0.1:

######################################################################
eff_vals=r_[1, 10**-0.5, 0.1]

# at z=zz9:
my_mf.update(z=zz9)
myspl_zz9=US(log(my_mf.m/h0), log(my_mf.ngtm*h0**3), s=0, k=5)
splnu_output_zz9=US(my_mf.nu**0.5, log(my_mf.m/h0*fbary), s=0)
splnu_output_rev_zz9=US(log(my_mf.m/h0*fbary), my_mf.nu**0.5, s=0)

n_zz9_eff1, n_zz9_eff2, n_zz9_eff3=splnu_output_rev_zz9(log(10.**mstar_data['log10_mstar'][1]/eff_vals))
print('peak height at z={0:.1f} for efficiency={1:.2f}: {2:.2f}'.format(zz9, eff_vals[0], n_zz9_eff1))
print('peak height at z={0:.1f} for efficiency={1:.2f}: {2:.2f}'.format(zz9, eff_vals[1], n_zz9_eff2))
print('peak height at z={0:.1f} for efficiency={1:.2f}: {2:.2f}\n'.format(zz9, eff_vals[2], n_zz9_eff3))

# also want to quote frac of collapsed mass for the peak height implied by this
# galaxy (assuming epsilon=1):
splnu_rhogtm=US(my_mf.nu**0.5,  log(my_mf.rho_gtm*h0**2), s=0, k=5)
frac_zz9=exp(splnu_rhogtm(n_zz9_eff1))/rhom0_noh
print('fraction of collapsed mass above nu={0:.1f} at z={1:.1f}: {2:.3e}\n'.format(n_zz9_eff1, zz9, frac_zz9))
######################################################################
# at z=zz8:
my_mf.update(z=zz8)
myspl_zz8=US(log(my_mf.m/h0), log(my_mf.ngtm*h0**3), s=0, k=5)
splnu_output_zz8=US(my_mf.nu**0.5, log(my_mf.m/h0*fbary), s=0)
splnu_output_rev_zz8=US(log(my_mf.m/h0*fbary), my_mf.nu**0.5, s=0)

n_zz8_eff1, n_zz8_eff2, n_zz8_eff3=splnu_output_rev_zz8(log(10.**mstar_data['log10_mstar'][0]/eff_vals))
print('peak height at z={0:.1f} for efficiency={1:.2f}: {2:.2f}'.format(zz8, eff_vals[0], n_zz8_eff1))
print('peak height at z={0:.1f} for efficiency={1:.2f}: {2:.2f}'.format(zz8, eff_vals[1], n_zz8_eff2))
print('peak height at z={0:.1f} for efficiency={1:.2f}: {2:.2f}\n'.format(zz8, eff_vals[2], n_zz8_eff3))

frac_zz8=exp(splnu_rhogtm(n_zz8_eff1))/rhom0_noh
print('fraction of collapsed mass above nu={0:.1f} at z={1:.1f}: {2:.3e}\n'.format(n_zz8_eff1, zz8, frac_zz8))


# get log10 cumulative number densities of the Labbe galaxies:
print('number densities: ')
print(log10(exp(myspl_zz9(log(10.**mstar_data['log10_mstar'][1]/fbary)))))
print(log10(exp(myspl_zz8(log(10.**mstar_data['log10_mstar'][0]/fbary)))))
print()

# print M_halo for nu corresponding to z=7.5 galaxy, assuming eff=1:
my_mf_nu.update(z=0)
nu_ind=(my_mf_nu.nu**0.5).searchsorted(n_zz8_eff1)+1
print('mass at z=0 corresponding to nu={0:.2f}: {1:.2e}'.format((my_mf_nu.nu**0.5)[nu_ind], my_mf_nu.m[nu_ind]/h0))
