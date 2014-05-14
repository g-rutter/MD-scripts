#Filled out in the order of table 2 in http://onlinelibrary.wiley.com/doi/10.1002/prot.22817/full
SC_SC_well_depth = {'L N': '0.015', 'L M': '-0.2', 'L L': '-0.2', 'L K': '0.015', 'L I': '-0.2', 'L H': '0.015', 'L F': '-0.203', 'L E': '0.015', 'L D': '0.015', 'N Q': '-0.080', 'N N': '-0.080', 'N H': '-0.080', 'L Y': '-0.203', 'L W': '-0.203', 'L V': '-0.2', 'L T': '0.015', 'L S': '0.015', 'L R': '0.015', 'L Q': '0.015', 'L P': '0.015', 'P S': '0.074', 'P R': '0.074', 'W N': '-0.086', 'P P': '0.074', 'W H': '-0.086', 'W K': '0.015', 'W D': '-0.086', 'W E': '-0.086', 'P E': '0.074', 'P D': '0.074', 'W T': '-0.086', 'W W': '-0.205', 'W P': '0.015', 'W Q': '-0.086', 'W R': '0.015', 'W S': '-0.086', 'Q H': '-0.080', 'K R': '0.073', 'K S': '-0.086', 'F T': '0.015', 'D E': '0.253', 'D D': '0.253', 'F Q': '0.015', 'F P': '0.015', 'F S': '0.015', 'F R': '0.015', 'D N': '-0.086', 'D K': '-0.136', 'D H': '-0.086', 'F E': '0.015', 'F D': '0.015', 'D T': '-0.086', 'D S': '-0.086', 'D R': '-0.136', 'D Q': '-0.086', 'F N': '0.015', 'F H': '0.015', 'F K': '0.015', 'H H': '-0.080', 'T Q': '-0.086', 'C D': '-0.116', 'H Q': '-0.080', 'R K': '0.073', 'T T': '-0.086', 'T S': '-0.086', 'Y D': '-0.086', 'Y E': '-0.086', 'Y K': '-0.086', 'Y H': '-0.086', 'Y N': '-0.086', 'P Q': '0.074', 'C C': '-0.139', 'Y R': '-0.086', 'Y S': '-0.086', 'Y P': '0.015', 'Y Q': '-0.086', 'Y W': '-0.201', 'Y T': '-0.086', 'T N': '-0.086', 'Y Y': '-0.201', 'T H': '-0.086', 'P T': '0.074', 'R S': '-0.086', 'R R': '0.073', 'R T': '-0.086', 'M W': '-0.210', 'M T': '-0.116', 'Q Q': '-0.080', 'M R': '-0.116', 'M S': '-0.116', 'M P': '0.015', 'M Q': '-0.116', 'K K': '0.073', 'F W': '-0.205', 'M Y': '-0.210', 'M F': '-0.203', 'M D': '0.015', 'M E': '0.015', 'K T': '-0.086', 'M N': '-0.116', 'M M': '-0.2', 'M K': '-0.116', 'M H': '-0.116', 'V E': '0.015', 'V D': '0.015', 'V F': '-0.203', 'V M': '-0.2', 'V L': '-0.2', 'V N': '0.015', 'V I': '-0.2', 'V H': '0.015', 'V K': '0.015', 'V T': '0.015', 'V W': '-0.203', 'V V': '-0.2', 'V Q': '0.015', 'V P': '0.015', 'V S': '0.015', 'V R': '0.015', 'Q N': '-0.080', 'V Y': '-0.203', 'C E': '-0.116', 'C Y': '-0.116', 'F Y': '-0.205', 'P K': '0.074', 'P H': '0.074', 'P N': '0.074', 'F F': '-0.205', 'E N': '-0.086', 'A A': '-0.084', 'E K': '-0.136', 'E H': '-0.086', 'A F': '-0.148', 'E D': '0.253', 'E E': '0.253', 'C V': '-0.139', 'E T': '-0.086', 'C W': '-0.116', 'E R': '-0.136', 'E S': '-0.086', 'E Q': '-0.086', 'I R': '0.015', 'I S': '0.015', 'I P': '0.015', 'I Q': '0.015', 'I V': '-0.2', 'I W': '-0.203', 'I T': '0.015', 'I Y': '-0.203', 'A Y': '-0.148', 'I F': '-0.203', 'I D': '0.015', 'I E': '0.015', 'I K': '0.015', 'I H': '0.015', 'I I': '-0.2', 'I N': '0.015', 'I L': '-0.2', 'I M': '-0.2', 'S H': '-0.086', 'S N': '-0.086', 'C A': '-0.139', 'C F': '-0.139', 'S Q': '-0.086', 'S S': '-0.086', 'S T': '-0.086', 'H N': '-0.080', 'A K': '0.074', 'A H': '0.074', 'A I': '-0.148', 'A N': '0.074', 'A L': '-0.148', 'A M': '-0.148', 'C P': '0.015', 'C Q': '-0.116', 'C R': '-0.116', 'C S': '-0.116', 'C T': '-0.116', 'A D': '0.074', 'A E': '0.074', 'C H': '-0.116', 'C I': '-0.139', 'C K': '-0.116', 'C L': '-0.139', 'C M': '-0.139', 'C N': '-0.116', 'A R': '0.074', 'A S': '0.074', 'A P': '0.074', 'A Q': '0.074', 'A V': '-0.148', 'A W': '-0.148', 'A T': '0.074'}

for aa1 in ['K','R']:
    for aa2 in ['N','H','Q']:
        SC_SC_well_depth[aa1 + ' ' + aa2] = '-0.086'

#Supplemental information table 2 and 3 from main PRIME20 paper
SC_SC_diam = { 'A A':'2.7', 'A V':'2.7', 'A P':'2.9', 'A T':'2.6', 'A S':'2.3', 'A N':'2.8', 'A D':'2.6', 'A R':'3.0', 'A K':'3.3', 'A E':'2.9', 'A Q':'3.0', 'A L':'2.7', 'A I':'2.9', 'A F':'2.4', 'A Y':'2.7', 'A W':'2.7', 'A M':'2.9', 'A C':'2.8', 'A H':'3.1', 'V V':'3.3', 'V P':'3.3', 'V T':'2.8', 'V S':'2.8', 'V N':'3.1', 'V D':'3.0', 'V R':'3.1', 'V K':'3.1', 'V E':'3.1', 'V Q':'3.3', 'V L':'3.0', 'V I':'3.3', 'V F':'3.2', 'V Y':'3.0', 'V W':'2.9', 'V M':'3.0', 'V C':'2.9', 'V H':'3.1', 'P P':'3.1', 'P T':'2.6', 'P S':'3.2', 'P N':'3.3', 'P D':'3.2', 'P R':'3.0', 'P K':'3.6', 'P E':'3.5', 'P Q':'3.6', 'P L':'3.5', 'P I':'3.5', 'P F':'3.1', 'P Y':'3.3', 'P W':'3.4', 'P M':'3.7', 'P C':'3.0', 'P H':'3.7', 'T T':'2.9', 'T S':'2.9', 'T N':'3.1', 'T D':'3.1', 'T R':'3.2', 'T K':'3.1', 'T E':'3.1', 'T Q':'3.3', 'T L':'3.2', 'T I':'3.0', 'T F':'2.8', 'T Y':'3.2', 'T W':'3.3', 'T M':'3.6', 'T C':'2.7', 'T H':'2.9', 'S S':'2.5', 'S N':'3.0', 'S D':'2.8', 'S R':'3.0', 'S K':'3.0', 'S E':'2.9', 'S Q':'2.7', 'S L':'3.0', 'S I':'2.6', 'S F':'2.9', 'S Y':'2.9', 'S W':'2.7', 'S M':'3.2', 'S C':'2.8', 'S H':'2.6', 'N N':'3.3', 'N D':'3.2', 'N R':'2.9', 'N K':'3.2', 'N E':'3.1', 'N Q':'3.5', 'N L':'3.4', 'N I':'2.8', 'N F':'2.7', 'N Y':'3.3', 'N W':'2.8', 'N M':'3.5', 'N C':'3.1', 'N H':'3.4', 'D D':'3.4', 'D R':'3.0', 'D K':'3.0', 'D E':'2.9', 'D Q':'2.8', 'D L':'3.0', 'D I':'3.4', 'D F':'3.1', 'D Y':'2.8', 'D W':'3.2', 'D M':'3.6', 'D C':'3.2', 'D H':'2.8', 'R R':'3.2', 'R K':'3.9', 'R E':'3.1', 'R Q':'3.6', 'R L':'3.4', 'R I':'3.6', 'R F':'3.3', 'R Y':'3.1', 'R W':'3.0', 'R M':'3.7', 'R C':'3.3', 'R H':'3.5', 'K K':'3.5', 'K E':'3.4', 'K Q':'3.4', 'K L':'3.5', 'K I':'2.9', 'K F':'3.5', 'K Y':'3.5', 'K W':'3.5', 'K M':'3.7', 'K C':'2.7', 'K H':'3.4', 'E E':'3.2', 'E Q':'2.9', 'E L':'3.3', 'E I':'3.2', 'E F':'3.3', 'E Y':'3.3', 'E W':'3.5', 'E M':'3.3', 'E C':'2.7', 'E H':'3.3', 'Q Q':'3.6', 'Q L':'3.5', 'Q I':'3.1', 'Q F':'3.3', 'Q Y':'3.4', 'Q W':'3.4', 'Q M':'3.4', 'Q C':'3.1', 'Q H':'3.3', 'L L':'3.4', 'L I':'3.4', 'L F':'3.4', 'L Y':'3.2', 'L W':'3.4', 'L M':'3.6', 'L C':'3.4', 'L H':'3.2', 'I I':'3.3', 'I F':'3.4', 'I Y':'3.0', 'I W':'3.2', 'I M':'3.6', 'I C':'3.3', 'I H':'3.1', 'F F':'3.3', 'F Y':'3.2', 'F W':'3.4', 'F M':'3.2', 'F C':'3.2', 'F H':'2.9', 'Y Y':'3.0', 'Y W':'3.2', 'Y M':'3.2', 'Y C':'2.9', 'Y H':'3.1', 'W W':'3.7', 'W M':'3.2', 'W C':'3.3', 'W H':'3.2', 'M M':'3.7', 'M C':'3.4', 'M H':'3.6', 'C C':'2.1', 'C H':'2.8', 'H H':'3.4' }
SC_SC_well_diam = { 'A A':'5.4','A V':'6.1','A P':'6.2','A T':'6.2','A S':'5.9','A N':'5.6','A D':'5.6','A R':'6.1','A K':'6.0','A E':'5.9','A Q':'5.8','A L':'5.6','A I':'5.7','A F':'5.9','A Y':'5.7','A W':'5.5','A M':'5.8','A C':'5.9','A H':'5.5', 'V V':'6.3','V P':'6.3','V T':'6.4','V S':'6.2','V N':'6.3','V D':'6.3','V R':'6.8','V K':'6.6','V E':'6.5','V Q':'6.5','V L':'6.2','V I':'6.4','V F':'6.5','V Y':'6.5','V W':'6.6','V M':'6.4','V C':'6.0','V H':'6.2', 'P P':'6.5','P T':'6.6','P S':'6.1','P N':'6.2','P D':'6.3','P R':'6.8','P K':'6.7','P E':'6.4','P Q':'6.5','P L':'6.3','P I':'6.4','P F':'6.5','P Y':'6.4','P W':'6.3','P M':'6.2','P C':'6.0','P H':'6.3', 'T T':'6.5','T S':'6.0','T N':'6.3','T D':'6.2','T R':'6.8','T K':'6.5','T E':'6.4','T Q':'6.4','T L':'6.2','T I':'6.4','T F':'6.6','T Y':'6.4','T W':'6.5','T M':'6.4','T C':'6.1','T H':'6.3', 'S S':'6.4','S N':'6.2','S D':'6.1','S R':'6.3','S K':'6.1','S E':'6.0','S Q':'6.0','S L':'6.3','S I':'6.4','S F':'6.2','S Y':'6.5','S W':'6.3','S M':'6.4','S C':'6.3','S H':'6.3', 'N N':'6.3','N D':'6.5','N R':'6.6','N K':'6.5','N E':'6.4','N Q':'6.4','N L':'6.4','N I':'6.6','N F':'6.5','N Y':'6.7','N W':'6.9','N M':'6.4','N C':'6.2','N H':'6.5', 'D D':'6.5','D R':'6.5','D K':'6.3','D E':'6.6','D Q':'6.3','D L':'6.5','D I':'6.5','D F':'6.7','D Y':'6.9','D W':'6.9','D M':'6.7','D C':'6.2','D H':'6.6', 'R R':'7.2','R K':'6.8','R E':'6.6','R Q':'6.9','R L':'6.8','R I':'6.7','R F':'6.9','R Y':'7.0','R W':'6.9','R M':'6.6','R C':'6.3','R H':'6.9', 'K K':'6.9','K E':'6.4','K Q':'6.7','K L':'6.5','K I':'6.7','K F':'6.9','K Y':'6.7','K W':'6.5','K M':'6.4','K C':'6.4','K H':'6.6', 'E E':'6.7','E Q':'6.6','E L':'6.4','E I':'6.6','E F':'6.8','E Y':'6.8','E W':'6.9','E M':'6.4','E C':'6.1','E H':'6.4', 'Q Q':'6.6','Q L':'6.3','Q I':'6.6','Q F':'6.6','Q Y':'6.7','Q W':'6.7','Q M':'6.4','Q C':'6.1','Q H':'6.6', 'L L':'6.4','L I':'6.5','L F':'6.6','L Y':'6.7','L W':'6.9','L M':'6.5','L C':'6.1','L H':'6.5', 'I I':'6.6','I F':'6.6','I Y':'6.8','I W':'6.8','I M':'6.7','I C':'6.4','I H':'6.6', 'F F':'6.8','F Y':'6.8','F W':'7.0','F M':'6.5','F C':'6.4','F H':'6.5', 'Y Y':'7.0','Y W':'7.0','Y M':'6.6','Y C':'6.5','Y H':'6.9', 'W W':'7.4','W M':'7.0','W C':'6.4','W H':'7.1', 'M M':'6.7','M C':'6.3','M H':'6.5', 'C C':'6.2','C H':'6.2', 'H H':'6.7' }

#Make sure that the parameter will be found no matter which way round it is queried:
for dictionary in [SC_SC_diam,SC_SC_well_diam,SC_SC_well_depth]:
    for key,value in dictionary.items():
        dictionary[key[::-1]] = value

BB_BB_diam = {'NH NH':'3.3','NH CH':'3.5','NH CO':'3.65','CH CH':'3.7','CH CO':'3.85','CO CO':'4.0'}
BB_BB_diam.update( {'NH NH':'3.3','CH NH':'3.5','CO NH':'3.65','CH CH':'3.7','CO CH':'3.85'} )

close_diam = {}
BB_SC_diam = {}
BB_SC_well_diam = {}
BB_SC_well_depth = {}

for bb in ['NH', 'CH', 'CO']:
    for res in list('ACDEFHIKLMNPQRSTVWY'):
        pair                    = bb + ' ' + res
        pair2                   = res + ' ' + bb

        #Geometric mixing rules based on 2010 Hall paper for SC and 2001 paper parameters for BB.
        BB_SC_diam[pair]  = str (0.5*( float(BB_BB_diam[bb + ' ' + bb]) + float(SC_SC_diam[res + ' ' + res]) ))
        BB_SC_diam[pair2] = BB_SC_diam[pair]

        if bb == 'NH' and ['M','Y','C','E','D','S','T','Q','N','H'].count(res) == 1: #Residues that can Hbond with NH
            BB_SC_well_depth[pair] = '-0.15'
            BB_SC_well_depth[pair2] = '-0.15'

            BB_SC_well_diam[pair]  = str ( 2.1 + 0.5*float(SC_SC_well_diam[res + ' ' + res]) )
            BB_SC_well_diam[pair2] = BB_SC_diam[pair]

        elif bb == 'CO' and ['Y','W','C','K','R','S','T','Q','N','H'].count(res) == 1: #Residues that can Hbond with CO
            BB_SC_well_depth[pair] = '-0.15'
            BB_SC_well_depth[pair2] = '-0.15'

            BB_SC_well_diam[pair]  = str ( 2.1 + 0.5*float(SC_SC_well_diam[res + ' ' + res]) )
            BB_SC_well_diam[pair2] = BB_SC_diam[pair]

for diam_dict in [BB_BB_diam, BB_SC_diam, SC_SC_diam]:
    for pair, value in diam_dict.items():
        close_diam[pair] = str( 0.75*float(value) )
