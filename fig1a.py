my_filename=__file__

get_ipython().run_line_magic('run', '-i common.py')

my_xlims=r_[20, 5]
my_ylims=r_[1e8, 1e12]

######################################################################
# plot z vs M* with lines for fixed cumulative number density
numdens=logspace(-8, -2, 7)
num_z=25
zdat=numpy.linspace(5, 20, num_z)
outdat=numpy.zeros((zdat.size, numdens.size))
# add data for lower number densities in a different color
numdens2=logspace(-10, -9, 2)
outdat2=numpy.zeros((zdat.size, 2))
for i in range(zdat.size):
    my_mf.update(z=zdat[i])
    ngtm_i=my_mf.ngtm*h0**3
    m_i=my_mf.m/h0

    # only do spline to maxind to avoid issues with nonsensically low abundances:
    try:
        maxind=(ngtm_i < 0).nonzero()[0][0]-1
    except IndexError:
        maxind=ngtm_i.size

    if (ngtm_i[maxind-1] == 0):
        maxind=(ngtm_i > 0).nonzero()[0].size

    myspl=US(log(ngtm_i[:maxind][::-1]), log(m_i[:maxind][::-1]), s=0)
    outdat[i]=myspl(log(numdens))
    outdat2[i]=myspl(log(numdens2))
    
# note the data are in terms of m_halo; galaxies are at least a factor of f_bary
# lower in (stellar) mass, so we need to multiply by a factor of f_bary to get
# the cumulative number density corresponding to a given galaxy.
my_size=outdat.shape[1]
for i in range(my_size):
    if i % 2 == 0:
        semilogy(zdat, exp(outdat[:,i])*fbary, lw=4, ls='-',
                 color=my_cmap(i/(my_size-0.5)), 
                 label=r'$10^{{{0}}}$'.format(int(log10(numdens[i]))))
    else:
        semilogy(zdat, exp(outdat[:,i])*fbary, lw=4, ls=(0, (5, 1)), 
                 color=my_cmap(i/(my_size-0.5)))

# ######################################################################
# also plot 2 lower values of number densities
semilogy(zdat, exp(outdat2[:,0])*fbary, lw=2.5, color='k', ls='-', alpha=0.2)
semilogy(zdat, exp(outdat2[:,1])*fbary, lw=2, color='k', ls='-', alpha=0.5)


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


# ######################################################################
# add annotation for \epsilon < 1
text(7.5, 0.15*1e10,  r'$\downarrow$', fontsize=30,
     color='k', va='top', ha='center')
text(7.1, 0.9*0.15*1e10, r'$\epsilon < 1$', color='k',
     va='top', ha='left')

xlim(*my_xlims)
ylim(*my_ylims)
legend(frameon=False, loc=(0.04, 0.56))
text(19.5, 4.2e11,
     r'$n(>f_{\rm b}^{-1}\epsilon^{-1}\,M_{\star})\;[{\rm Mpc}^{-3}]:$')

# ######################################################################
# add annotation for "extra" number densities

text(14.6, 0.925*2.995e11, r'$10^{-10}$', va='top', ha='left', color='k',
     alpha=0.35, fontsize=16)
text(14.6, 0.925*1.5e11, r'$10^{-9}$', va='top', ha='left', color='k',
     alpha=0.6, fontsize=16)


# this placement is actually right; it looks a little off on screen but ends up
# right in the PDF.
plot(10.95, 1.95e8, marker='*', ms=16, color=top_color, linestyle='none', mew=0)
text(5.5, 1.95e8/1.0625, r'${\rm Labb\'e\;et\;al. [13]}$', color=top_color,
     va='center', ha='right')

savefig(my_filename.split('.')[0] + '.pdf')










######################################################################
######################################################################
# also want to output, at z=9.1 and z=7.5 (the Labbe et al. points), the
# comoving number densities for \epsilon = 0.316 and 0.1:
######################################################################
eff_vals=r_[1, 10**-0.5, 0.1]

######################################################################
# at z=zz9:
my_mf.update(z=zz9)
myspl_zz9=US(log(my_mf.m/h0), log(my_mf.ngtm*h0**3), s=0, k=5)

n_zz9_eff1, n_zz9_eff2, n_zz9_eff3=exp(myspl_zz9(log(10.**mstar_data['log10_mstar'][1]/eff_vals/fbary)))
print('n/mpc^-3 at z={0:.1f} for efficiency={1:.2f}: {2:.2e}'.format(zz9, eff_vals[0], n_zz9_eff1))
print('n/mpc^-3 at z={0:.1f} for efficiency={1:.2f}: {2:.2e}'.format(zz9, eff_vals[1], n_zz9_eff2))
print('n/mpc^-3 at z={0:.1f} for efficiency={1:.2f}: {2:.2e}\n'.format(zz9, eff_vals[2], n_zz9_eff3))

######################################################################
# at z=zz8:
my_mf.update(z=zz8)
myspl_zz8=US(log(my_mf.m/h0), log(my_mf.ngtm*h0**3), s=0, k=5)

n_zz8_eff1, n_zz8_eff2, n_zz8_eff3=exp(myspl_zz8(log(10.**mstar_data['log10_mstar'][0]/eff_vals/fbary)))
print('n/mpc^-3 at z={0:.1f} for efficiency={1:.2f}: {2:.2e}'.format(zz8, eff_vals[0], n_zz8_eff1))
print('n/mpc^-3 at z={0:.1f} for efficiency={1:.2f}: {2:.2e}'.format(zz8, eff_vals[1], n_zz8_eff2))
print('n/mpc^-3 at z={0:.1f} for efficiency={1:.2f}: {2:.2e}'.format(zz8, eff_vals[2], n_zz8_eff3))


# get log10 cumulative number densities of the Labbe galaxies:
print('number densities: ')
print(log10(exp(myspl_zz9(log(yvs[1]/fbary)))))
print(log10(exp(myspl_zz8(log(yvs[0]/fbary)))))
