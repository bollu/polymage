def set_verification(app_data, is_periodic):
    print_line()
    print("[verify]: Setting the verification values ...")

    verify_data = {}
    verify_data['verified'] = False
    verify_data['epsilon'] = 1.0e-8

    verify_value = 0.0

    Class = app_data['prob_class']

    if Class == 'S':
        if is_periodic == True:
            verify_value = 0.530770700573e-04
        else:
            verify_value = 0.115028215110e-01
    elif Class == 'W':
        if is_periodic == True:
            verify_value = 0.646732937534e-05
        else:
            verify_value = 0.731581264884e-01
    elif Class == 'A':
        if is_periodic == True:
            verify_value = 0.243336530907e-05
        else:
            verify_value = 0.979991065870e-01
    elif Class == 'B':
        if is_periodic == True:
            verify_value = 0.180056440136e-05
        else:
            verify_value = 0.184378110108e+19
    elif Class == 'C':
        if is_periodic == True:
            verify_value = 0.570673228574e-06
        else:
            verify_value = 0.474623829181e+24
    elif Class == 'D':
        if is_periodic == True:
            verify_value = 0.158327506043e-09
        else:
            verify_value = 0.183242406095e+84
    #fi

    verify_data['verify_value'] = verify_value

    app_data['verify_data'] = verify_data

    return

def verify_norm(app_data):
    print_line()
    print("[verify]: Verifying the results ...")

    verify_data = app_data['verify_data']
    rnm2 = app_data['rnm2']

    verify_value = verify_data['verify_value']
    epsilon = verify_data['epsilon']

    err = abs(rnm2 - verify_value) / verify_value

    if err <= epsilon:
        print("[verify]: VERIFICATION SUCCESSFUL")
        print("[verify]: L2 norm            =", rnm2)
        print("[verify]: Error              =", err)
    else:
        print("[verify]: VERIFICATION FAILED")
        print("[verify]: L2 norm            =", rnm2)
        print("[verify]: Correct L2 norm    =", verify_value)
    #fi

    verify_data['verified'] = True

    return
