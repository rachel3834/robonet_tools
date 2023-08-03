from pyLIMA.magnification import magnification_VBB


params = [0.5029366596327286, 0.2522489899965984, -368.5663343290413, -1606.0751474552462, 0.0005305607885400907]

params2 = [0.5029366596327286, 0.2522489899965984, -293.9209277989914, -1280.1654233196082, 0.0005305607885400907]

mags = []

for i in range(10000):

    if i%2==0:
    
        par = params
        
    else:
    
        par = params2


    mago = magnification_VBB.magnification_USBL([par[0]],par[1],[par[2]],[par[3]],par[4])
    mags.append(mago[0])

print(mago)