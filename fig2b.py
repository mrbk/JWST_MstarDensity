my_filename=__file__
rcParams['legend.labelspacing']=0.5
rcParams['legend.handlelength']=2.

get_ipython().run_line_magic('run', '-i common.py')

# peak of PDF for Labbe et al. z~7.5 most massive galaxy
my_mf.update(z=zz8)
rho_gtm=my_mf.rho_gtm*h0**2


# data from Labbe / van Dokkum
xvs=rhostar_data_z8['mstar'][-2:].copy()
yvs_log=log10(rhostar_data_z8['rhostar'][-2:].copy()/vol_z8)


yvs=10.**yvs_log
ye_low=yvs-rhostar_data_z8['rhostar_lo'][-2:].copy()/vol_z8
ye_upp=rhostar_data_z8['rhostar_hi'][-2:].copy()/vol_z8-yvs

my_xlims=r_[1e8, 1e12]
my_ylims=r_[1e2, 1e8]

######################################################################
figure()
errorbar(xvs, yvs,
         yerr=c_[ye_low, ye_upp].T,
         # xerr=c_[xe_low, xe_upp].T, 
         color='k', ecolor='k',
         ls='none', capsize=5, elinewidth=3, marker='s', ms=8, 
         label=r'${\rm Labb\'e\;et\;al. [13]}$', zorder=100)
loglog(mvals_noh*fbary, fbary*rho_gtm, 'k-', lw=4, label=r'$\epsilon=1.0$')
# if conversion efficiency is 31.6%:
loglog(mvals_noh*fbary*10**-0.5, fbary*rho_gtm*10**-0.5,
       ls='-', color=bo_color, lw=4, label=r'$\epsilon=0.32$')
# if conversion efficiency is 10%
loglog(mvals_noh*fbary*0.1, fbary*rho_gtm*0.1,
       ls='-', color='0.5', lw=4, label=r'$\epsilon=0.1$')

fill_between(mvals_noh*fbary, fbary*rho_gtm, 1e9, facecolor='b',
             alpha=0.25, edgecolors='face', linestyle='None', 
             lw=0)

xlim(*my_xlims)
ylim(*my_ylims)
xlabel(r'$M_{\star} \;\;{\rm or}\;\; \epsilon\,f_{\rm b}\,M_{\rm halo}\;[M_{\odot}]$')
ylabel(r'$\rho_{\star}(>M_{\star}) \;\;{\rm or}\;\; \epsilon \,f_{\rm b}\,\rho_{\rm m}(>M_{\rm halo})\;\,[M_{\odot}\,{\rm Mpc}^{-3}]$')
legend(loc=(0.02, 0.035), frameon=False)
text(6.5e11, 3.5e7, 'more stellar mass than', color='b', ha='right')
text(6.5e11, 1.5e7, 'available baryons', color='b', ha='right')
text(1.25e11, 4e6, r'$z = 7.5$', fontsize=18, ha='center')
savefig(my_filename.split('.')[0] + '.pdf')
