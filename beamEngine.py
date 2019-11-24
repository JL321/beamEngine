import numpy as np

#Calculates failure load of baldwin loads using a PI beam configuration

def evalBeam(w_init, h2_init, h3_init, printFlag = False):
    
    t = 1.27 #Thickness
    w = w_init #Mini rectangle width
    h2 = h2_init #Length of webbing leading up to top mini rectangle - ADD TO SMALL RECTANGLE HEIGHT FOR TOTAL WEB HEIGHT
    h3 = h3_init #Mini rectangle height
    b = 100 #Top width - < 100
    h1 = 1.27 #Height of top component
    aVal = 100 #A value for shear buckling
    E = 4000
    availableExt = b - 2*w - 2*t
    sideL = availableExt/5 #Manual computation
    middleL = availableExt/5*3
    maxComp = 6
    maxShear = 4
    maxL =[maxComp, maxShear]
    
    #print("Verified- middle and side L: {}, {}".format(middleL, sideL))
    
    total_h= h1+h2+h3
    
    #Note restriction: Total height of bridge - h1+h2+h3 - should be smaller than 100
    
    I_1 = t*(h2+h3)**3/12
    I_2 = w*h3**3/12
    I_3 = h1**3*b/12
    I = np.array([I_1*2, I_2*2, I_3])
    
    A_1 = (h2+h3)*t
    A_2 = w*h3
    A_3 = h1*b 
    A = np.array([2*A_1, 2*A_2, A_3])
    y_1 = (h2+h3)/2
    y_2 = h2+h3/2
    y_3 = h2+h3+h1/2
    Y = np.array([y_1, y_2, y_3])
    
    y_cent = np.sum(A*Y)/np.sum(A)
    I = np.sum(I+A*(Y-y_cent)**2)
    #print("I value: {}".format(I))
    if y_cent > h2:
        raise Exception('Centroidal axis is above h2!')
    
    Q_cent = 2*(y_cent-y_cent/2)*(y_cent*t)
    local_y = (A_2*y_2+A_1*y_1)/(A_2+A_1)
    
    Q_glue_bot = 2*(A_1+A_2)*(y_cent - local_y) #Second term represents the rectangular sections towards the top of
    Q_glue_top = (A_3)*(total_h-h1/2-y_cent) #the pi beam, while the first represents the rectangular component 
    
    Q_glue = np.maximum(Q_glue_bot, Q_glue_top)
    #print("Bot and Top: {}, {}".format(Q_glue_bot, Q_glue_top))
    #print("Q Glue and Q Cent: {} {}".format(Q_glue, Q_cent))
    
    maxFlexTNeg = I*30/(.19*(total_h-y_cent)*10**3)
    maxFlexTPos = I*30/(.166*(y_cent)*10**3) #Account for m to mm conversion
    
    #All Flexural Failure Cases - What we really care about for the PI beam however is the case of positive moments
    maxShearGlue = 2*I*(2*t+2*w)/(Q_glue)
    
    #Shear Failure cases for glue and centroidal
    
    pbc = np.pi**2*E/(12*(1-.2**2)) #Plate buckling constant
    plateBuckFixed = 4*pbc*(h1/middleL)**2
    plateBuckFree = 0.425*pbc*(h1/sideL)**2
    plateBuckWeb = 6*pbc*(t/((h2+h3)/2))**2    
    plateBuckShearWeb = 5*pbc*((t/(h2+h3))**2 + (t/aVal)**2)

    newComp = np.min(np.array([plateBuckFixed, plateBuckFree, plateBuckWeb, maxComp]))
    if newComp != maxComp and printFlag:
        print("New Comp: {}".format(newComp))
    newShear = np.min(np.array([plateBuckShearWeb, maxShear]))
    if newShear != maxShear and printFlag:
        print("New Shear: {}".format(newShear))

    maxFlexCNeg = I*maxComp/(.19*(y_cent)*10**3)
    maxFlexCPos = I*maxComp/(.166*(total_h-y_cent)*10**3)
    if printFlag:
        print("Base Flexural C (+M): {}".format(maxFlexCPos))
        print("Base Flexural C (-M): {}".format(maxFlexCNeg))
    maxShearCent = newShear*I*(2*t)/(Q_cent) #Max shear is 1P in the denominator
    maxFlexCNeg = I*newComp/(.19*(y_cent)*10**3)
    maxFlexCPos = I*newComp/(.166*(total_h-y_cent)*10**3)
    if printFlag:
        print("Max Flexurals: Negative Moment")
        print("T, C: {}, {}".format(maxFlexTNeg, maxFlexCNeg))
        print("Max Flexurals: Positive Moment")
        print("T, C: {}, {}".format(maxFlexTPos, maxFlexCPos))
        print("Shear Failures (Glue,Cent) {}, {}".format(maxShearGlue, maxShearCent))
        print("Plate Bucklings")
        print("Plate Fixed and Free (Stress): {}, {}".format(plateBuckFixed, plateBuckFree))
        print("Plate Web and Shear (Stress): {}, {}".format(plateBuckWeb, plateBuckShearWeb))
    
    failureLoads = np.array([maxShearCent, maxFlexCPos, maxFlexTPos])
    idx = np.argmin(failureLoads)
    if idx == 0:
        cause = 'Shear'
    elif idx == 1:
        cause = 'Compression'
    elif idx == 2:
        cause = 'Tension'
    
    return failureLoads[idx], cause, (w, h2, h3)
    
def main():
    
    bestLoad = 0
    
    for w in range(5, 30):
        for h2 in range(110, 131):
            for h3 in [1.27, 2.54]:
                fLoad, fCause, params = evalBeam(w, h2, h3)
                if fLoad > bestLoad:
                    bestLoad = fLoad
                    assCause = fCause
                    bestParam = params
                    
    print("Best Case Failure Load: {}".format(fLoad))
    print("Cause of said load: {}".format(assCause))
    print("Set of best params: {}".format(bestParam))
    
if __name__ == '__main__':
    
    main()
