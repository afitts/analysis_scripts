def tvir(mvir,rvir):
    G = 4.3e-6 #in kpc/M_sun (km/s)^2
    kb = 1.3806e-26 #in km^2*g/s^2 
    mu = 0.59
    mp = 1.6726219e-24
    tv = G*mvir/rvir*(0.5*mu*mp)/kb
    return tv
hdf['props']
rvir = hdf['props']['rvir']/(hdf['props']['redshift']+1)*1e3
rvir
rvir = float(rvir)
rvir
mvir = sum(hdf['particles/dm']['mass'].as_matrix())+sum(hdf['particles/gas']['mass'].as_matrix())
mvir
tvir(mvir,rvir)
