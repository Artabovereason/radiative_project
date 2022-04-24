print('Processus names :')
print('— Radio')
print('— Exoplanet detection')
print('— Obstacle\n')
print('Checks names :')
print('— 3D plot')
print('— Planck sampling\n')
print('Enter the name of either a processus or a check you want to execute.')
name = input('What do you want to execute ?\n')

if   name == 'Radio':
    print('')
    print('This processus takes lot\'s of photons to show something.')
    print('The figure shown in the report used near 100000 photons,')
    print('which took about 1 hour to run.')
    print('We suggest you to run it with 20000 photons.')
    number_photons= int(input('How much photons do you want in your simulation ? (20000 takes 300s)\n'))#10000
    exec(open('radio.py').read())

elif name == 'Exoplanet detection':
    print('')
    print('This processus takes lot\'s of photons to show something.')
    print('It\'ll run for about 20 minutes')
    exec(open('exoplanet.py').read())

elif name == 'Obstacle':
    print('')
    exec(open('obstacle.py').read())
    print('If you want to reproduce the Figure (7) of the report,')
    print('You just have to modify the first parameter of the line 239 of obstacle.py')
    print('to different values and you\'ll be able to plot the Figure (7)')


elif name == '3D plot':
    exec(open('3Dplot.py').read())

elif name == 'Planck sampling': 
    print('')
    print('This processus takes about 100s to run.')
    exec(open('plancksampling.py').read())

else:
    print('The processus you wanted to execute is not in the list.')
