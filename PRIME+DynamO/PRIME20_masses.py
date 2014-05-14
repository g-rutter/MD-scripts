#From: http://www.cell.com/biophysj/fulltext/S0006-3495%2811%2901018-6
mass = {'CH':'0.866', 'NH':'0.999', 'CO':'1.863', 'K':'4.865', 'L':'3.799', 'V':'2.866', 'F':'6.061', 'A':'1.000', 'E':'4.793'}

#From: http://onlinelibrary.wiley.com/doi/10.1002/pro.2141/full
mass.update({'Q':'4.795', 'I':'3.799', 'Y':'7.126'})

#Calculated by asking Wolfram Alpha the molecular weight, taking it in AMU, subtracting 74 (which is the weight in AMU of the backbone atoms) and dividing the result by 15.
mass.update({'R':'6.666', 'N':'3.866', 'D':'3.933', 'C':'3.133', 'H':'5.400', 'L':'3.800', 'M':'5.000',
    'P':'2.733', 'S':'2.066', 'T':'3.000', 'W':'8.666'})

#The units of these are currently 1.0 = 15 AMU.
