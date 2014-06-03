from math import sqrt

#Units of angstrom
BB_BB_diam = {'NH CH':'1.46', 'CH CO':'1.51', 'CO NH':'1.33'}
BB_BB_diam.update({'CH NH':BB_BB_diam['NH CH'], 'CO CH':BB_BB_diam['CH CO'], 'NH CO':BB_BB_diam['CO NH']})

BB_pseudo = {'NH CO':'2.45', 'CH NH':'2.41', 'CO CH': '2.45', 'CH CH': '3.80'}
BB_pseudo.update( {'CO NH':'2.45', 'NH CH':'2.41', 'CH CO': '2.45'} )
BB_SC_diam = {}
BB_SC_well_diam = {}
BB_SC_pseudo_diam = {}
BB_SC_pseudo_well_diam = {}

#My own parameter set, see PRIME20_parametrisation dir
BB_SC_ideal_pseudobondlens = {}
BB_SC_ideal_bondlens = {
    'CH A':'1.57', 'A CH':'1.57',
    'CH R':'4.23', 'R CH':'4.23',
    'CH N':'2.53', 'N CH':'2.53',
    'CH D':'2.53', 'D CH':'2.53',
    'CH C':'2.37', 'C CH':'2.37',
    'CH Q':'3.11', 'Q CH':'3.11',
    'CH E':'3.24', 'E CH':'3.24',
    'CH H':'3.16', 'H CH':'3.16',
    'CH I':'2.36', 'I CH':'2.36',
    'CH L':'2.64', 'L CH':'2.64',
    'CH K':'3.61', 'K CH':'3.61',
    'CH M':'3.14', 'M CH':'3.14',
    'CH F':'3.41', 'F CH':'3.41',
    'CH P':'1.91', 'P CH':'1.91',
    'CH S':'1.96', 'S CH':'1.96',
    'CH T':'1.97', 'T CH':'1.97',
    'CH W':'3.93', 'W CH':'3.93',
    'CH Y':'3.86', 'Y CH':'3.86',
    'CH V':'1.99', 'V CH':'1.99'
}

#numbers from http://0-onlinelibrary.wiley.com.pugwash.lib.warwick.ac.uk/doi/10.1002/prot.1100/full
#(Alpha helix formation)

##rationale: keep the ratio intact
#for SC in list('ACDEFHIKLMNPQRSTVWY'):
    #BB_SC_ideal_pseudobondlens['NH ' + SC] = str(2.44*(float(BB_SC_ideal_bondlens['CH ' + SC])/1.531))
    #BB_SC_ideal_pseudobondlens[SC + ' NH'] = str(2.44*(float(BB_SC_ideal_bondlens['CH ' + SC])/1.531))

    #BB_SC_ideal_pseudobondlens['CO ' + SC] = str(2.49*(float(BB_SC_ideal_bondlens['CH ' + SC])/1.531))
    #BB_SC_ideal_pseudobondlens[SC + ' CO'] = str(2.49*(float(BB_SC_ideal_bondlens['CH ' + SC])/1.531))

#Rationale, keep angles CO-CH-SC and NH-CH-SC intact
for SC in list('ACDEFHIKLMNPQRSTVWY'):

    CH_SC_dist = float(BB_SC_ideal_bondlens['CH ' + SC])

    BB_SC_ideal_pseudobondlens['NH ' + SC] = sqrt(2.1316 + CH_SC_dist**2 + (CH_SC_dist*1.05) )
    BB_SC_ideal_pseudobondlens[SC + 'NH '] = BB_SC_ideal_pseudobondlens['NH ' + SC]

    BB_SC_ideal_pseudobondlens['CO ' + SC] = sqrt(2.2801 + CH_SC_dist**2 + (CH_SC_dist*1.12) )
    BB_SC_ideal_pseudobondlens[SC + 'CO '] = BB_SC_ideal_pseudobondlens['CO ' + SC]



#################################
#  create diams and well diams  #
#################################
for key in BB_SC_ideal_bondlens:
    BB_SC_diam[key] = str(float(BB_SC_ideal_bondlens[key])*0.97625)
    BB_SC_well_diam[key] = str(float(BB_SC_ideal_bondlens[key])*1.02375)

for key in BB_SC_ideal_pseudobondlens:
    BB_SC_pseudo_diam[key] = str(float(BB_SC_ideal_pseudobondlens[key])*0.97625)
    BB_SC_pseudo_well_diam[key] = str(float(BB_SC_ideal_pseudobondlens[key])*1.02375)

#For printing out LAMMPS parameters:
if __name__ == "__main__":
    def plus1(i):
        while True:
            i+=1
            yield i

    go=plus1(7)

    for SC in list('ACDEFHIKLMNPQRSTVWY'):
        for BB in ['NH','CH','CO']:
            try:
                len=BB_SC_ideal_bondlens[BB+' '+SC]
            except:
                len=BB_SC_ideal_pseudobondlens[BB+' '+SC]
            print "bond_coeff   {3:2}   90.0         {0:.4}   #{1}-{2}".format(len,SC,BB, go.next())

