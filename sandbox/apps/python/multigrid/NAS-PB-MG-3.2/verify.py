def setVerification(dataDict, isPeriodic):
    verifyDict = {}
    verifyDict['verified'] = False
    verifyDict['epsilon'] = 1.0e-8

    verifyValue = 0.0

    Class = dataDict['probClass']

    if Class == 'S':
        if isPeriodic == True:
            verifyValue = 0.530770700573e-04
        else:
            verifyValue = 0.115028215110e-01
    elif Class == 'W':
        if isPeriodic == True:
            verifyValue = 0.646732937534e-05
        else:
            verifyValue = 0.731581264884e-01
    elif Class == 'A':
        if isPeriodic == True:
            verifyValue = 0.243336530907e-05
        else:
            verifyValue = 0.979991065870e-01
    elif Class == 'B':
        if isPeriodic == True:
            verifyValue = 0.180056440136e-05
        else:
            verifyValue = 0.184378110108e+19
    elif Class == 'C':
        if isPeriodic == True:
            verifyValue = 0.570673228574e-06
        else:
            verifyValue = 0.474623829181e+24
    elif Class == 'D':
        if isPeriodic == True:
            verifyValue = 0.158327506043e-09
        else:
            verifyValue = 0.183242406095e+84
    #fi

    verifyDict['verifyValue'] = verifyValue

    dataDict['verifyDict'] = verifyDict

    return

def verifyNorm(dataDict):
    verifyDict = dataDict['verifyDict']
    rnm2 = dataDict['rnm2']

    verifyValue = verifyDict['verifyValue']
    epsilon = verifyDict['epsilon']

    err = abs(rnm2 - verifyValue) / verifyValue

    if err <= epsilon:
        print "[verify]: VERIFICATION SUCCESSFUL"
        print "[verify]: L2 norm            =",rnm2
        print "[verify]: Error              =",err
    else:
        print "[verify]: VERIFICATION FAILED"
        print "[verify]: L2 norm            =",rnm2
        print "[verify]: Correct L2 norm    =",verifyValue
    #fi

    verifyDict['verified'] = True

    return
