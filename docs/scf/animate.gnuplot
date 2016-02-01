#plot sin(x+i*dt) w l lt 1 lw 1.5 title sprintf("t=%i",i)
#set logscale x
ind=3
indd=11
j=i+1
#plot "./G_Imp_Re_".i u 1:ind w l lw 2 t 'Impurity GF', "./G_Lat_Re_".i u 1:ind w l t 'Lattice before', "./G_Lat_Re_".j u 1:ind w l t 'Lattice after'
plot "./G_Imp_".i u 1:ind w l lw 2 t 'Impurity GF iter='.i, "./G_Lat_".i u 1:ind w l lw 2 t 'Lattice GF iter='.i
#,"./G_Imp_Fit_".i u 1:indd w l lw 2 t 'Impurity GF', "./G_Lat_Fit_".i u 1:indd w l t 'Lattice before', "./G_Lat_Fit_".j u 1:indd w l t 'Lattice after'
#plot 'fort.'.i w l lt 1
i=i+1
pause -1
#pause mouse
#pause 5
#if (i < n) reread
reread
